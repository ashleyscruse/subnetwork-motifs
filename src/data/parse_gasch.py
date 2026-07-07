"""
Parse the Gasch 2000 expression matrix into a clean ORF x condition CSV.

Input:
    data/expression/raw/gasch_complete_dataset.txt
    Tab-separated, single header row. Columns:
      [0] UID    - ORF systematic name (e.g. YAL001C)
      [1] NAME   - free-form description (gene name + ontology + SGDID)
      [2] GWEIGHT - gene weight (always 1)
      [3..] log2 expression ratios under labeled conditions

Output:
    data/expression/expression.csv
    Wide: rows = ORFs, columns = condition labels, values = log2 ratios.
    First column header: 'gene' (renamed from 'UID' for consistency with
    the rest of the pipeline). Missing values stay as empty cells.

Optional:
    --long  also emit data/expression/expression_long.csv (gene,
            condition, value), skipping NaN entries.

Usage:
    python -m src.data.parse_gasch
    python -m src.data.parse_gasch --long
"""

import argparse
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = REPO_ROOT / "data" / "expression" / "raw" / "gasch_complete_dataset.txt"
DEFAULT_OUTPUT = REPO_ROOT / "data" / "expression" / "expression.csv"
DEFAULT_OUTPUT_LONG = REPO_ROOT / "data" / "expression" / "expression_long.csv"

DROP_COLUMNS = ["NAME", "GWEIGHT"]


def load_and_clean(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = df.drop(columns=[c for c in DROP_COLUMNS if c in df.columns])
    df = df.rename(columns={"UID": "gene"})
    df["gene"] = df["gene"].astype(str).str.strip()
    df = df[df["gene"].str.len() > 0].drop_duplicates(subset=["gene"]).reset_index(drop=True)
    return df


def report(df: pd.DataFrame) -> None:
    n_genes = len(df)
    n_conditions = df.shape[1] - 1
    cell_total = n_genes * n_conditions
    cell_present = df.iloc[:, 1:].notna().sum().sum()
    coverage = cell_present / cell_total if cell_total else 0.0
    print(f"Genes:              {n_genes:,}")
    print(f"Conditions:         {n_conditions:,}")
    print(f"Cells present:      {cell_present:,} of {cell_total:,} ({coverage:.1%})")


def main() -> None:
    p = argparse.ArgumentParser(description="Parse Gasch 2000 expression matrix")
    p.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--long", action="store_true", help="Also write a long-format CSV")
    p.add_argument("--output-long", type=Path, default=DEFAULT_OUTPUT_LONG)
    args = p.parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input not found: {args.input}\nRun download_gasch.py first.")

    df = load_and_clean(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"Wrote {args.output}")

    if args.long:
        long = (
            df.melt(id_vars=["gene"], var_name="condition", value_name="log2_ratio")
              .dropna(subset=["log2_ratio"])
        )
        long.to_csv(args.output_long, index=False)
        print(f"Wrote {args.output_long} ({len(long):,} rows)")

    report(df)


if __name__ == "__main__":
    main()
