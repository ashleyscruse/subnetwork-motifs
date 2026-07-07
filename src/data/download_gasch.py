"""
Download the Gasch et al. 2000 yeast environmental stress expression matrix.

Source paper:
    Gasch et al. 2000, "Genomic expression programs in the response of yeast
    cells to environmental changes." Mol Biol Cell 11:4241-4257.
    DOI: 10.1091/mbc.11.12.4241

Hosting note:
    The original companion site at genome-www.stanford.edu/yeast_stress/ is
    offline. The complete dataset file is preserved in the Internet Archive
    at the Wayback URL below (snapshot 2021-10-19).

What gets downloaded:
    complete_dataset.txt (~6 MB), tab-separated:
      Header row: UID, NAME, GWEIGHT, then 173 condition labels
      Body: one row per ORF (~6,150 rows). UID is systematic name.
            Values are log2 expression ratios; some cells empty (missing).

Usage:
    python -m src.data.download_gasch
"""

import argparse
from pathlib import Path
from urllib.error import URLError
from urllib.request import Request, urlopen

DEFAULT_URL = (
    "http://web.archive.org/web/20211019175531if_/"
    "http://genome-www.stanford.edu/yeast_stress/data/rawdata/complete_dataset.txt"
)

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO_ROOT / "data" / "expression" / "raw" / "gasch_complete_dataset.txt"

USER_AGENT = "subnetwork-motifs/0.1 (academic research)"


def download(url: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"GET {url}")
    req = Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(req, timeout=120) as resp:
            out_path.write_bytes(resp.read())
    except (URLError, TimeoutError) as e:
        raise SystemExit(f"Download failed: {e}")
    print(f"Wrote {out_path.stat().st_size:,} bytes to {out_path}")


def main() -> None:
    p = argparse.ArgumentParser(description="Download Gasch 2000 stress-response data")
    p.add_argument("--url", default=DEFAULT_URL)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--force", action="store_true")
    args = p.parse_args()
    if args.out.exists() and not args.force:
        print(f"Already present: {args.out} (use --force to re-download)")
        return
    download(args.url, args.out)


if __name__ == "__main__":
    main()
