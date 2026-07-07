"""
Download the S. cerevisiae S288C reference protein FASTA from SGD.

Source:
    Saccharomyces Genome Database (https://www.yeastgenome.org/), file
    sequence/S288C_reference/orf_protein/orf_trans.fasta.gz on
    sgd-archive.yeastgenome.org.

Why this file:
    FASTA headers are keyed by systematic ORF names (e.g. >YAL001C), which
    matches every other dataset in this project (Harbison, Wapinski). No
    cross-database ID mapping needed.

Usage:
    python -m src.data.download_yeast_proteome
"""

import argparse
import gzip
import shutil
from pathlib import Path
from urllib.error import URLError
from urllib.request import Request, urlopen

DEFAULT_URL = (
    "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/"
    "orf_protein/orf_trans.fasta.gz"
)

REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = REPO_ROOT / "data" / "families" / "raw"
DEFAULT_GZ = RAW_DIR / "orf_trans.fasta.gz"
DEFAULT_FASTA = RAW_DIR / "orf_trans.fasta"

USER_AGENT = "subnetwork-motifs/0.1 (academic research)"


def download(url: str, out_gz: Path) -> None:
    out_gz.parent.mkdir(parents=True, exist_ok=True)
    print(f"GET {url}")
    req = Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(req, timeout=120) as resp:
            out_gz.write_bytes(resp.read())
    except (URLError, TimeoutError) as e:
        raise SystemExit(f"Download failed: {e}")
    print(f"Wrote {out_gz.stat().st_size:,} bytes to {out_gz}")


def gunzip(gz_path: Path, out_path: Path) -> None:
    with gzip.open(gz_path, "rb") as fin, out_path.open("wb") as fout:
        shutil.copyfileobj(fin, fout)
    print(f"Extracted {out_path} ({out_path.stat().st_size:,} bytes)")


def main() -> None:
    p = argparse.ArgumentParser(description="Download SGD yeast protein FASTA")
    p.add_argument("--url", default=DEFAULT_URL)
    p.add_argument("--gz", type=Path, default=DEFAULT_GZ)
    p.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    p.add_argument("--keep-gz", action="store_true")
    p.add_argument("--force", action="store_true")
    args = p.parse_args()

    if args.fasta.exists() and not args.force:
        print(f"Already present: {args.fasta} (use --force to re-download)")
        return

    download(args.url, args.gz)
    gunzip(args.gz, args.fasta)
    if not args.keep_gz:
        args.gz.unlink()


if __name__ == "__main__":
    main()
