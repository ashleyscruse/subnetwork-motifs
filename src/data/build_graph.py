"""
Build the consolidated yeast regulatory genomics graph artifact.

This script joins the five processed CSV outputs in data/ into a single
PyTorch Geometric Data object stored at data/processed/yeast_graph.pt. It is
the orchestrator that implements the schema specified in section 3 of
docs/data.md. All node-feature computation lives in features.py; this module
owns the I/O and the join logic.

Usage:
    python -m src.data.build_graph
    python -m src.data.build_graph --node-table   # also write sidecar CSV

Sources (paths fixed by repo layout):
    data/harbison/edges.csv             regulatory edges
    data/wapinski/paralogs.csv          paralog pairs + duplication-fate labels
    data/families/families.csv          gene -> family string ID
    data/families/sequence_similarity.csv  paralog pair percent identity
    data/expression/expression.csv      gene x condition expression matrix

Output:
    data/processed/yeast_graph.pt       PyG Data object (single artifact)
    data/processed/node_table.csv       optional human-readable sidecar
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch_geometric.data import Data

from src.data.features import (
    compute_degree_features,
    compute_family_ids,
    compute_paralog_similarities,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
HARBISON_PATH = REPO_ROOT / "data" / "harbison" / "edges.csv"
PARALOGS_PATH = REPO_ROOT / "data" / "wapinski" / "paralogs.csv"
FAMILIES_PATH = REPO_ROOT / "data" / "families" / "families.csv"
SIMILARITY_PATH = REPO_ROOT / "data" / "families" / "sequence_similarity.csv"
EXPRESSION_PATH = REPO_ROOT / "data" / "expression" / "expression.csv"

DEFAULT_OUTPUT = REPO_ROOT / "data" / "processed" / "yeast_graph.pt"
DEFAULT_NODE_TABLE = REPO_ROOT / "data" / "processed" / "node_table.csv"

# CV fold count and random seed are fixed at build time so the artifact's
# splits are reproducible across machines and re-runs of the pipeline.
NUM_FOLDS = 5
RANDOM_SEED = 0


def load_sources():
    """Read the five processed CSVs into pandas DataFrames.

    No transformation is applied here; downstream functions expect the raw
    output schemas documented in docs/data.md section 2.

    Returns
    -------
    edges, paralogs, families, similarity, expression : pandas.DataFrame
    """
    edges = pd.read_csv(HARBISON_PATH)
    paralogs = pd.read_csv(PARALOGS_PATH)
    families = pd.read_csv(FAMILIES_PATH)
    similarity = pd.read_csv(SIMILARITY_PATH)
    expression = pd.read_csv(EXPRESSION_PATH)
    return edges, paralogs, families, similarity, expression


def build_node_order(edges, paralogs, expression):
    """Compute the canonical node order as the sorted union of genes seen.

    A gene becomes a node if it appears in any of: a regulatory edge (as TF
    or target), a paralog pair, or the expression matrix. Sorting yields a
    deterministic order, so the same input CSVs always produce the same
    node indices.

    Returns
    -------
    list of str
        Gene systematic names ordered by node index.
    """
    genes = (
        set(edges["tf"]) | set(edges["target"])
        | set(paralogs["gene1"]) | set(paralogs["gene2"])
        | set(expression["gene"])
    )
    return sorted(genes)


def build_edge_index(edges, gene_to_idx):
    """Translate (tf, target) string pairs to integer node-index pairs.

    Returns
    -------
    edge_index : numpy.ndarray, shape (2, num_edges), dtype int64
        PyG convention: row 0 = source, row 1 = destination.
    edge_attr : numpy.ndarray, shape (num_edges, 2), dtype float64
        Columns are [min_pvalue, n_conditions_bound] from the Harbison
        output. Preserved on each edge so confidence-weighted variants of
        downstream code can use them; ordinary code can ignore.
    """
    src = edges["tf"].map(gene_to_idx).to_numpy(dtype=np.int64)
    dst = edges["target"].map(gene_to_idx).to_numpy(dtype=np.int64)
    edge_index = np.stack([src, dst])
    edge_attr = edges[["min_pvalue", "n_conditions_bound"]].to_numpy(
        dtype=np.float64
    )
    return edge_index, edge_attr


def build_paralog_edges_and_labels(paralogs, gene_to_idx):
    """Translate paralog pairs to integer indices; extract per-pair labels.

    Pair-level storage (rather than per-gene projection) preserves every
    paralog relationship individually. See docs/data.md section 3.3 for
    the rationale.

    Returns
    -------
    paralog_edges : numpy.ndarray, shape (2, num_pairs), dtype int64
        Same PyG edge format as edge_index. Indexes into the same node set
        as the regulatory edges.
    paralog_labels : numpy.ndarray, shape (num_pairs,), dtype int64
        regulation_retained per pair: 1 = retained, 0 = lost / rewired.
    """
    g1 = paralogs["gene1"].map(gene_to_idx).to_numpy(dtype=np.int64)
    g2 = paralogs["gene2"].map(gene_to_idx).to_numpy(dtype=np.int64)
    paralog_edges = np.stack([g1, g2])
    paralog_labels = paralogs["regulation_retained"].to_numpy(dtype=np.int64)
    return paralog_edges, paralog_labels


def build_expression_matrix(genes, expression):
    """Align the expression matrix to the canonical gene order.

    Genes absent from `expression` get an all-NaN row from `reindex`; those
    NaNs are then imputed to 0.0 (the no-change baseline in log2-ratio
    space) per the missing-value rule in docs/data.md section 3.2. This
    treats "no expression measurement" as "no observed change," which is
    the most conservative imputation in this representation.

    Returns
    -------
    numpy.ndarray, shape (len(genes), num_conditions), dtype float64
    """
    matrix = (
        expression.set_index("gene").reindex(genes).to_numpy(dtype=np.float64)
    )
    return np.nan_to_num(matrix, nan=0.0)


def assign_family_folds(
    family_ids, paralog_edges, num_folds=NUM_FOLDS, seed=RANDOM_SEED
):
    """Assign each gene to a CV fold via its family.

    Only families containing at least one paralog pair (i.e. families that
    contribute supervisable pairs) are assigned to folds 0..num_folds-1.
    Genes whose families have no paralog pairs (singletons and any other
    unlabeled families) receive fold = -1, marking them as not used in
    cross-validation. Family-level assignment prevents leakage: every gene
    in a family lands in the same fold, so a paralog pair never has one
    member in train and the other in validation.

    Parameters
    ----------
    family_ids : numpy.ndarray, shape (num_nodes,), int dtype
        Output of compute_family_ids.
    paralog_edges : numpy.ndarray, shape (2, num_pairs), int dtype
        Output of build_paralog_edges_and_labels.
    num_folds : int
    seed : int
        Fixed at build time so fold assignment is reproducible.

    Returns
    -------
    numpy.ndarray, shape (num_nodes,), dtype int64
    """
    rng = np.random.default_rng(seed)

    # Genes that participate in at least one paralog pair.
    pair_genes = np.unique(paralog_edges)
    # Families those genes belong to. By construction (families are
    # connected components of the paralog graph), both members of any pair
    # share a family, so every pair lands in a single fold.
    families_with_pairs = np.unique(family_ids[pair_genes])

    # Random shuffle, then round-robin assignment to folds.
    shuffled = rng.permutation(families_with_pairs)
    family_to_fold = {
        int(fam): int(i % num_folds) for i, fam in enumerate(shuffled)
    }

    # Propagate from family fold to per-gene fold.
    family_fold = np.full(len(family_ids), -1, dtype=np.int64)
    for i, fam in enumerate(family_ids):
        family_fold[i] = family_to_fold.get(int(fam), -1)
    return family_fold


def build_graph():
    """Assemble the consolidated PyTorch Geometric Data object.

    See docs/data.md section 3 for the artifact schema. Steps:
      1. Load all five source CSVs.
      2. Define the canonical node order (sorted union of genes).
      3. Translate regulatory edges and paralog pairs to integer indices.
      4. Compute per-node attributes via features.py.
      5. Pack continuous features into a single x tensor for PyG
         convention; family_id stays separate as a categorical attribute.
      6. Compute family-level CV folds (reproducible from RANDOM_SEED).
      7. Wrap everything in a Data object.
    """
    edges, paralogs, families, similarity, expression = load_sources()

    genes = build_node_order(edges, paralogs, expression)
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    num_nodes = len(genes)

    edge_index, edge_attr = build_edge_index(edges, gene_to_idx)
    paralog_edges, paralog_labels = build_paralog_edges_and_labels(
        paralogs, gene_to_idx
    )

    family_id = compute_family_ids(genes, families)
    in_degree, out_degree = compute_degree_features(edge_index, num_nodes)
    expression_matrix = build_expression_matrix(genes, expression)
    mean_paralog_identity = compute_paralog_similarities(genes, similarity)

    # Pack continuous per-node features into x in a fixed, documented order.
    # data.feature_names records the column meaning so downstream code does
    # not have to hard-code positions.
    x = np.concatenate(
        [
            in_degree[:, None].astype(np.float64),
            out_degree[:, None].astype(np.float64),
            expression_matrix,
            mean_paralog_identity[:, None],
        ],
        axis=1,
    )
    feature_names = (
        ["in_degree", "out_degree"]
        + [f"expression_{i}" for i in range(expression_matrix.shape[1])]
        + ["mean_paralog_identity"]
    )

    family_fold = assign_family_folds(family_id, paralog_edges)

    data = Data(
        x=torch.from_numpy(np.array(x)).float(),
        edge_index=torch.from_numpy(np.array(edge_index)).long(),
        edge_attr=torch.from_numpy(np.array(edge_attr)).float(),
        paralog_edges=torch.from_numpy(np.array(paralog_edges)).long(),
        paralog_labels=torch.from_numpy(np.array(paralog_labels)).long(),
        family_id=torch.from_numpy(np.array(family_id)).long(),
        family_fold=torch.from_numpy(np.array(family_fold)).long(),
    )
    # Non-tensor metadata. Pickled with the rest of the Data object.
    data.gene_names = genes
    data.feature_names = feature_names
    return data


def build_node_table(data):
    """Build a human-readable per-gene sidecar table for spot-checking.

    Includes the scalar attributes only (degrees, family ID, identity,
    fold). The 173 expression columns are omitted to keep the table
    inspectable; the full matrix is already in data/expression/.
    """
    return pd.DataFrame(
        {
            "node_idx": np.arange(data.num_nodes),
            "gene": data.gene_names,
            "family_id": data.family_id.numpy(),
            "in_degree": data.x[:, 0].numpy().astype(np.int64),
            "out_degree": data.x[:, 1].numpy().astype(np.int64),
            "mean_paralog_identity": data.x[:, -1].numpy(),
            "family_fold": data.family_fold.numpy(),
        }
    )


def report(data):
    """Print build-time stats; useful for verifying sanity numbers."""
    print(f"Nodes:                  {data.num_nodes:,}")
    print(f"Regulatory edges:       {data.edge_index.shape[1]:,}")
    print(f"Paralog pairs:          {data.paralog_edges.shape[1]:,}")
    print(f"  retained (label=1):   {(data.paralog_labels == 1).sum().item():,}")
    print(f"  lost     (label=0):   {(data.paralog_labels == 0).sum().item():,}")
    print(f"Continuous feature dim: {data.x.shape[1]:,}")
    print(f"Unique family IDs:      {len(torch.unique(data.family_id)):,}")
    in_cv = (data.family_fold >= 0).sum().item()
    print(f"Genes assigned to CV:   {in_cv:,}  ({in_cv / data.num_nodes:.1%})")


def main():
    parser = argparse.ArgumentParser(
        description="Build the consolidated yeast regulatory graph artifact."
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--node-table",
        action="store_true",
        help="Also write a human-readable sidecar CSV.",
    )
    parser.add_argument(
        "--node-table-path", type=Path, default=DEFAULT_NODE_TABLE
    )
    args = parser.parse_args()

    data = build_graph()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    torch.save(data, args.output)
    print(f"Wrote {args.output}")
    report(data)

    if args.node_table:
        table = build_node_table(data)
        args.node_table_path.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(args.node_table_path, index=False)
        print(f"Wrote {args.node_table_path}")


if __name__ == "__main__":
    main()
