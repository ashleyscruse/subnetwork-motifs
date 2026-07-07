"""
Per-node feature computation for the consolidated yeast regulatory genomics
dataset.

This module contains pure functions that transform already-loaded source data
into per-node attribute arrays aligned to a fixed node order. Loading the
source CSVs and orchestrating the calls is the job of build_graph.py.

Every function takes a `genes` argument (or an aligned input) defining the
canonical node order 0..N-1, and returns a numpy array of length N (or shape
[N, k]). This contract keeps each function independently testable with
synthetic inputs and decouples feature logic from I/O.

Documented behaviors here are the implementation of the per-node attribute
schema specified in docs/data.md section 3.
"""

import numpy as np
import pandas as pd


def compute_family_ids(genes, families_df):
    """Map each gene to an integer family ID aligned to the input order.

    Wapinski-derived families from data/families/families.csv (with string
    IDs like "FAM0001") are encoded as integers 0..M-1 in sorted-string
    order. Genes not present in families_df are singletons and receive
    unique integer IDs starting at M, in the order they appear in `genes`.
    This guarantees:
      (1) no collision between multi-gene families and singletons,
      (2) every gene gets a well-defined family ID,
      (3) the encoding is deterministic for a fixed (genes, families_df).

    Parameters
    ----------
    genes : sequence of str
        Ordered gene systematic names. Position i in the output corresponds
        to genes[i].
    families_df : pandas.DataFrame
        Loaded data/families/families.csv with at least the columns 'gene'
        and 'family_id'.

    Returns
    -------
    numpy.ndarray of shape (len(genes),), dtype int64
        Integer family IDs aligned to the `genes` order.
    """
    # Sort the multi-gene Wapinski family IDs and assign 0..M-1 so the
    # encoding is reproducible across runs of the build pipeline.
    multi_family_strings = sorted(families_df["family_id"].unique())
    string_to_int = {s: i for i, s in enumerate(multi_family_strings)}
    num_multi_families = len(multi_family_strings)

    # Per-gene lookup table for genes that appear in any Wapinski pair.
    gene_to_family_string = dict(
        zip(families_df["gene"], families_df["family_id"])
    )

    # Singleton IDs start after the last multi-gene-family ID, so the two
    # ranges are disjoint and downstream code can distinguish "in a Wapinski
    # family" (id < num_multi_families) from "singleton" (id >=).
    family_ids = np.empty(len(genes), dtype=np.int64)
    next_singleton_id = num_multi_families
    for i, gene in enumerate(genes):
        if gene in gene_to_family_string:
            family_ids[i] = string_to_int[gene_to_family_string[gene]]
        else:
            family_ids[i] = next_singleton_id
            next_singleton_id += 1
    return family_ids


def compute_degree_features(edge_index, num_nodes):
    """Count incoming and outgoing edges for each node.

    Given an edge_index in PyG convention (row 0 = source, row 1 =
    destination), returns the in-degree (count of edges where the node is
    the destination) and out-degree (count of edges where the node is the
    source) per node.

    Parameters
    ----------
    edge_index : numpy.ndarray of shape (2, num_edges), integer dtype
        Edge list. edge_index[0, k] is the source node index of edge k;
        edge_index[1, k] is the destination node index.
    num_nodes : int
        Total node count. Sets the length of the output arrays so nodes
        with no edges still appear with degree 0.

    Returns
    -------
    in_degree : numpy.ndarray of shape (num_nodes,), dtype int64
    out_degree : numpy.ndarray of shape (num_nodes,), dtype int64
    """
    # np.bincount with minlength ensures every node index 0..num_nodes-1
    # has an entry, including nodes that never appear in edge_index.
    in_degree = np.bincount(edge_index[1], minlength=num_nodes).astype(np.int64)
    out_degree = np.bincount(edge_index[0], minlength=num_nodes).astype(np.int64)
    return in_degree, out_degree


def normalize_expression(expression_matrix):
    """Z-score normalize each condition column of an expression matrix.

    For each column: subtract the column mean and divide by the column
    standard deviation. NaN values are ignored when computing column
    statistics (nanmean / nanstd), so missing measurements do not bias the
    standardization. NaN positions in the input are preserved as NaN in
    the output; imputation policy is the caller's responsibility, kept
    separate to match the convention in docs/data.md section 3.2.

    Constant columns (zero standard deviation) are centered without
    scaling, avoiding division by zero.

    Parameters
    ----------
    expression_matrix : numpy.ndarray of shape (num_genes, num_conditions)

    Returns
    -------
    numpy.ndarray of shape (num_genes, num_conditions), dtype float64
    """
    matrix = np.asarray(expression_matrix, dtype=np.float64)
    means = np.nanmean(matrix, axis=0, keepdims=True)
    stds = np.nanstd(matrix, axis=0, keepdims=True)
    # Replace zero stds with 1 so constant columns return as
    # (column - mean) / 1 = centered column, with no warning or NaN.
    safe_stds = np.where(stds == 0, 1.0, stds)
    return (matrix - means) / safe_stds


def compute_paralog_similarities(genes, similarity_df):
    """Compute mean percent identity to paralogs, per gene.

    For each gene, averages percent_identity across all paralog pairs in
    which the gene appears. The similarity table records each pair once as
    (gene1, gene2, percent_identity); both members of the pair receive the
    pair's contribution, so we symmetrize before grouping.

    Genes that do not appear in any paralog pair (singletons), and genes
    whose pairs were excluded during alignment (e.g. ORFs missing from the
    SGD reference proteome, listed in docs/data.md section 2.5), receive
    0.0. This matches the missing-value rule for `mean_paralog_identity`
    in docs/data.md section 3.2.

    Parameters
    ----------
    genes : sequence of str
        Ordered gene systematic names defining the output position.
    similarity_df : pandas.DataFrame
        Loaded data/families/sequence_similarity.csv with at least the
        columns 'gene1', 'gene2', 'percent_identity'.

    Returns
    -------
    numpy.ndarray of shape (len(genes),), dtype float64
    """
    # Symmetrize: each pair (a, b, sim) contributes once to gene a's
    # average and once to gene b's average.
    a = similarity_df[["gene1", "percent_identity"]].rename(
        columns={"gene1": "gene"}
    )
    b = similarity_df[["gene2", "percent_identity"]].rename(
        columns={"gene2": "gene"}
    )
    flat = pd.concat([a, b], ignore_index=True)

    # Mean identity per gene across all its paralog pairs.
    per_gene_mean = flat.groupby("gene")["percent_identity"].mean()

    # Align to the requested gene order; genes absent from the table get
    # 0.0 per the missing-value rule documented in docs/data.md.
    return per_gene_mean.reindex(genes, fill_value=0.0).to_numpy()
