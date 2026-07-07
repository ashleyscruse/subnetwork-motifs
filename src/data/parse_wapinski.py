"""
Parse the Wapinski 2007 PARALOGS sheet into duplication labels and families.

Input:
    data/wapinski/raw/MOESM281_table_S5_paralog_pairs.xls
    (also accepts the unrenamed MOESM281.xls)

    The "PARALOGS" sheet is laid out with two header rows:
    Row 0 group titles: PARALOGS | DUPLICATION POINT | GENE SET MIGRATION (4) |
                        PROTEIN INTERACTION NETWORK (3 genetic + 3 biochemical)
    Row 1 column names: GENE1, GENE2, (blank), GO-component, GO-process,
                        GO-function, Transcription Module, FRACTION GENETIC PI
                        SHARED NEIGHBORHOOD, CONS. SIGNIFICANCE, DIV.
                        SIGNIFICANCE, FRACTION BIOCHEMICAL PI SHARED NEIGHBOR-
                        HOOD, CONS. SIGNIFICANCE, DIV. SIGNIFICANCE
    Body: 1 row per S. cerevisiae paralog pair (~801 rows total).

Outputs:
    data/wapinski/paralogs.csv
        Per-pair table with all 13 cleaned columns, one row per paralog pair.
        Adds a derived 'regulation_retained' column: 1 if the two paralogs
        share the same transcription module (col 'transcription_module' == 0),
        0 if they diverged (col == 1). This is the binary label the GAT
        learns to predict (per plan.md Phase 2/3).

    data/families/families.csv
        Per-gene family assignment derived from connected components of the
        paralog graph. Two genes are in the same family if they appear in the
        same paralog pair, or are connected through a chain of pairs. Genes
        not in any pair get singleton family IDs.
        Columns: gene, family_id, family_size

Usage:
    python -m src.data.parse_wapinski
"""

import argparse
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = REPO_ROOT / "data" / "wapinski" / "raw"
DEFAULT_OUTPUT_PARALOGS = REPO_ROOT / "data" / "wapinski" / "paralogs.csv"
DEFAULT_OUTPUT_FAMILIES = REPO_ROOT / "data" / "families" / "families.csv"

PARALOG_COLUMNS = [
    "gene1",
    "gene2",
    "duplication_point",
    "migration_go_component",
    "migration_go_process",
    "migration_go_function",
    "transcription_module",
    "frac_genetic_shared",
    "genetic_cons_significance",
    "genetic_div_significance",
    "frac_biochemical_shared",
    "biochemical_cons_significance",
    "biochemical_div_significance",
]


def find_input(raw_dir: Path) -> Path:
    candidates = [
        raw_dir / "MOESM281_table_S5_paralog_pairs.xls",
        raw_dir / "MOESM281.xls",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise SystemExit(
        f"Wapinski paralog table not found. Looked for:\n  "
        + "\n  ".join(str(c) for c in candidates)
        + "\nRun download_wapinski.py first."
    )


def load_paralogs(path: Path) -> pd.DataFrame:
    df = pd.read_excel(
        path,
        sheet_name="PARALOGS",
        engine="xlrd",
        header=None,
        skiprows=2,
        names=PARALOG_COLUMNS,
    )
    df = df.dropna(subset=["gene1", "gene2"]).reset_index(drop=True)
    df["gene1"] = df["gene1"].astype(str).str.strip()
    df["gene2"] = df["gene2"].astype(str).str.strip()
    df["duplication_point"] = df["duplication_point"].astype(str).str.strip()
    for c in PARALOG_COLUMNS[3:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["regulation_retained"] = 1 - df["transcription_module"]
    return df


def derive_families(paralogs: pd.DataFrame) -> pd.DataFrame:
    """Connected-components family assignment over the paralog edge set."""
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        parent.setdefault(a, a)
        parent.setdefault(b, b)
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for a, b in zip(paralogs["gene1"], paralogs["gene2"]):
        union(a, b)

    rows = [(g, find(g)) for g in parent]
    fam = pd.DataFrame(rows, columns=["gene", "family_root"])
    family_ids = (
        fam["family_root"]
        .drop_duplicates()
        .sort_values()
        .reset_index(drop=True)
        .reset_index()
        .rename(columns={"index": "family_id", "family_root": "family_root"})
    )
    family_ids["family_id"] = family_ids["family_id"].apply(lambda i: f"FAM{i:04d}")
    fam = fam.merge(family_ids, on="family_root", how="left")
    sizes = fam.groupby("family_id").size().rename("family_size").reset_index()
    fam = fam.merge(sizes, on="family_id", how="left")
    return fam[["gene", "family_id", "family_size"]].sort_values("gene").reset_index(drop=True)


def report(paralogs: pd.DataFrame, families: pd.DataFrame) -> None:
    n_pairs = len(paralogs)
    n_wgd = (paralogs["duplication_point"] == "WGD").sum()
    n_with_label = paralogs["transcription_module"].notna().sum()
    n_retained = (paralogs["regulation_retained"] == 1).sum()
    n_lost = (paralogs["regulation_retained"] == 0).sum()
    n_genes = families["gene"].nunique()
    n_families = families["family_id"].nunique()
    multi = families["family_size"] > 1
    n_multi_families = families.loc[multi, "family_id"].nunique()
    n_genes_in_multi = multi.sum()

    print(f"Paralog pairs:          {n_pairs:,}")
    print(f"  WGD pairs:            {n_wgd:,}")
    print(f"  with labeled module:  {n_with_label:,}")
    print(f"    regulation kept:    {n_retained:,}")
    print(f"    regulation lost:    {n_lost:,}")
    print(f"Genes in any pair:      {n_genes:,}")
    print(f"Families derived:       {n_families:,}")
    print(f"  multi-gene families:  {n_multi_families:,}")
    print(f"  genes in multi fams:  {n_genes_in_multi:,}")


def main() -> None:
    p = argparse.ArgumentParser(description="Parse Wapinski 2007 paralog table")
    p.add_argument("--input", type=Path, default=None)
    p.add_argument("--output-paralogs", type=Path, default=DEFAULT_OUTPUT_PARALOGS)
    p.add_argument("--output-families", type=Path, default=DEFAULT_OUTPUT_FAMILIES)
    args = p.parse_args()

    src = args.input if args.input else find_input(RAW_DIR)
    paralogs = load_paralogs(src)
    families = derive_families(paralogs)

    args.output_paralogs.parent.mkdir(parents=True, exist_ok=True)
    args.output_families.parent.mkdir(parents=True, exist_ok=True)
    paralogs.to_csv(args.output_paralogs, index=False)
    families.to_csv(args.output_families, index=False)

    report(paralogs, families)
    print(f"Wrote {args.output_paralogs}")
    print(f"Wrote {args.output_families}")


if __name__ == "__main__":
    main()
