"""
Parse the Harbison 2004 binding p-value text files into a TF-target edge list.

Inputs:
    data/harbison/raw/pvalbygene_forpaper_ypd_txt_withheader.txt
    data/harbison/raw/pvalbygene_forpaper_otherconditions_txt_withheader.txt

File format (both):
    Row 1: experiment IDs (numeric, ignored).
    Row 2: '\\t\\t\\t' then TF_condition labels (ABF1_YPD, ARG81_SM, ...).
    Body: ORF (systematic, e.g. YAL001C) \\t common name \\t description
          \\t p-value per TF_condition column. Missing values: 'NaN'.

Output:
    data/harbison/edges.csv
    Columns: tf, target, min_pvalue, n_conditions_bound
    One row per (TF, target) pair where any condition crossed the threshold.
    For TFs profiled under multiple conditions, min_pvalue = best p-value
    seen, n_conditions_bound = number of conditions where p <= threshold.

Threshold:
    Default p <= 0.001, the canonical Harbison cutoff.

Usage:
    python -m src.data.parse_harbison
    python -m src.data.parse_harbison --threshold 0.005
"""

import argparse
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = REPO_ROOT / "data" / "harbison" / "raw"
DEFAULT_INPUTS = (
    RAW_DIR / "pvalbygene_forpaper_ypd_txt_withheader.txt",
    RAW_DIR / "pvalbygene_forpaper_otherconditions_txt_withheader.txt",
)
DEFAULT_OUTPUT = REPO_ROOT / "data" / "harbison" / "edges.csv"
DEFAULT_THRESHOLD = 0.001

ID_COLUMNS = ["target", "common_name", "description"]


def load_file(path: Path) -> pd.DataFrame:
    """Load one Harbison p-value file in long format.

    Returns columns: target (ORF), tf, condition, pvalue.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        skiprows=1,             # drop the row of numeric experiment IDs
        header=0,               # second row is now the header
        names=None,
        na_values=["NaN", "nan", "NA", ""],
        low_memory=False,
    )
    df.columns = [*ID_COLUMNS, *df.columns[len(ID_COLUMNS):]]
    long = df.melt(
        id_vars=ID_COLUMNS,
        var_name="tf_condition",
        value_name="pvalue",
    ).dropna(subset=["pvalue"])
    split = long["tf_condition"].str.rsplit("_", n=1, expand=True)
    long["tf"] = split[0].str.strip()
    long["condition"] = split[1].str.strip()
    return long[["target", "tf", "condition", "pvalue"]]


def matrix_to_edges(long: pd.DataFrame, threshold: float) -> pd.DataFrame:
    bound = long[long["pvalue"] <= threshold]
    edges = (
        bound.groupby(["tf", "target"], as_index=False)
        .agg(min_pvalue=("pvalue", "min"), n_conditions_bound=("condition", "nunique"))
        .sort_values(["tf", "min_pvalue"])
        .reset_index(drop=True)
    )
    return edges


def report(long: pd.DataFrame, edges: pd.DataFrame, threshold: float) -> None:
    n_orfs = long["target"].nunique()
    n_tf_conditions = long.groupby(["tf", "condition"]).ngroups
    n_tfs_total = long["tf"].nunique()
    n_tfs_bound = edges["tf"].nunique()
    n_targets_bound = edges["target"].nunique()
    print(f"ORFs assayed:        {n_orfs:,}")
    print(f"TF_condition pairs:  {n_tf_conditions:,}")
    print(f"Unique TFs assayed:  {n_tfs_total:,}  (expected ~203)")
    print(f"Threshold:           p <= {threshold}")
    print(f"Edges:               {len(edges):,}")
    print(f"  TFs with edges:    {n_tfs_bound:,}")
    print(f"  Targets reached:   {n_targets_bound:,}  (expected ~6,000)")


def main() -> None:
    p = argparse.ArgumentParser(description="Parse Harbison 2004 p-value files to edge list")
    p.add_argument("--inputs", nargs="+", type=Path, default=list(DEFAULT_INPUTS))
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    args = p.parse_args()

    missing = [str(p) for p in args.inputs if not p.exists()]
    if missing:
        raise SystemExit(
            "Input(s) not found:\n  " + "\n  ".join(missing)
            + "\nRun download_harbison.py first."
        )

    long = pd.concat([load_file(p) for p in args.inputs], ignore_index=True)
    edges = matrix_to_edges(long, args.threshold)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    edges.to_csv(args.output, index=False)
    report(long, edges, args.threshold)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
