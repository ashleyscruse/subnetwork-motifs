"""
Compute pairwise sequence identity for each Wapinski paralog pair.

Inputs:
    data/wapinski/paralogs.csv          (gene1, gene2, ...)
    data/families/raw/orf_trans.fasta   (SGD S288C protein sequences,
                                         FASTA headers keyed by systematic
                                         ORF name)

Method:
    For each (gene1, gene2) pair, look up both protein sequences from the
    FASTA, run a global pairwise alignment with BLOSUM62 + standard gap
    penalties (open=-10, extend=-0.5), and compute percent identity over
    the alignment length. Pairs where one or both sequences are missing
    from the FASTA are reported but excluded from the output.

Output:
    data/families/sequence_similarity.csv
    Columns: gene1, gene2, length1, length2, aln_length, n_matches,
             percent_identity

Usage:
    python -m src.data.compute_paralog_similarity
"""

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PARALOGS = REPO_ROOT / "data" / "wapinski" / "paralogs.csv"
DEFAULT_FASTA = REPO_ROOT / "data" / "families" / "raw" / "orf_trans.fasta"
DEFAULT_OUTPUT = REPO_ROOT / "data" / "families" / "sequence_similarity.csv"


def load_proteome(fasta: Path) -> dict[str, str]:
    """SGD orf_trans.fasta headers look like '>YAL001C TFC3 SGDID:S0000001 ...'.
    Key by the first whitespace token (systematic name).
    """
    records: dict[str, str] = {}
    for rec in SeqIO.parse(str(fasta), "fasta"):
        gene = rec.id.split()[0]
        seq = str(rec.seq).rstrip("*")
        records[gene] = seq
    return records


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner


def percent_identity(seq1: str, seq2: str, aligner: PairwiseAligner) -> tuple[int, int, float]:
    aln = aligner.align(seq1, seq2)[0]
    a, b = str(aln[0]), str(aln[1])
    aln_len = len(a)
    matches = sum(1 for x, y in zip(a, b) if x == y and x != "-")
    pct = 100.0 * matches / aln_len if aln_len else 0.0
    return aln_len, matches, pct


def main() -> None:
    p = argparse.ArgumentParser(description="Compute identity per Wapinski paralog pair")
    p.add_argument("--paralogs", type=Path, default=DEFAULT_PARALOGS)
    p.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--progress-every", type=int, default=100)
    args = p.parse_args()

    if not args.paralogs.exists():
        raise SystemExit(f"Missing {args.paralogs}; run parse_wapinski.py first.")
    if not args.fasta.exists():
        raise SystemExit(f"Missing {args.fasta}; run download_yeast_proteome.py first.")

    pairs = pd.read_csv(args.paralogs, usecols=["gene1", "gene2"])
    proteome = load_proteome(args.fasta)
    print(f"Loaded {len(proteome):,} protein sequences from FASTA")

    aligner = make_aligner()
    rows = []
    missing = []
    for i, (g1, g2) in enumerate(zip(pairs["gene1"], pairs["gene2"]), start=1):
        s1 = proteome.get(g1)
        s2 = proteome.get(g2)
        if s1 is None or s2 is None:
            missing.append((g1, g2, s1 is None, s2 is None))
            continue
        aln_len, matches, pct = percent_identity(s1, s2, aligner)
        rows.append({
            "gene1": g1,
            "gene2": g2,
            "length1": len(s1),
            "length2": len(s2),
            "aln_length": aln_len,
            "n_matches": matches,
            "percent_identity": round(pct, 3),
        })
        if i % args.progress_every == 0:
            print(f"  aligned {i}/{len(pairs)} pairs")

    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)

    print(f"\nPairs aligned:     {len(out):,}")
    if len(out):
        print(f"Mean % identity:   {out['percent_identity'].mean():.1f}")
        print(f"Median % identity: {out['percent_identity'].median():.1f}")
    if missing:
        print(f"Pairs skipped (missing in FASTA): {len(missing)}")
        for g1, g2, m1, m2 in missing[:10]:
            print(f"  {g1} (missing: {m1})  {g2} (missing: {m2})")
        if len(missing) > 10:
            print(f"  ...and {len(missing) - 10} more")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
