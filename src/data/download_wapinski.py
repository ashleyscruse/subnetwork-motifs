"""
Download the Wapinski et al. 2007 supplementary tables from Nature.

Source paper:
    Wapinski, Pfeffer, Friedman & Regev 2007,
    "Natural history and evolutionary principles of gene duplication in fungi"
    Nature 449:54-61. DOI: 10.1038/nature06107

Hosting note:
    The original full orthogroup dataset was hosted at the Broad Institute
    (broadinstitute.org/regev/orthogroups/), which is now offline and the
    bulk text dump is not preserved in the Wayback Machine. The Nature
    supplementary tables remain available; for our GAT we need MOESM281
    (the paralog pair classifications), MOESM278 (transcription modules,
    optional features), and MOESM276 (supplementary methods, for reference).

What gets downloaded:
    MOESM281.xls -- PARALOGS sheet: ~801 yeast paralog pairs with binary
                    indicators for GO/transcription-module migration and
                    interaction-network divergence. THIS IS THE KEY FILE.
    MOESM278.xls -- Transcription modules per ORF (optional GAT features).
    MOESM276.pdf -- Supplementary methods (reference only).

Usage:
    python -m src.data.download_wapinski
    python -m src.data.download_wapinski --all   # download all 7 supplements
"""

import argparse
from pathlib import Path
from urllib.error import URLError
from urllib.request import Request, urlopen

DOI_PATH = "art%3A10.1038%2Fnature06107/MediaObjects"
SPRINGER_BASE = f"https://static-content.springer.com/esm/{DOI_PATH}"

# (MOESM number, filename, extension, default include?)
SUPPLEMENTS = [
    (276, "supplementary_methods", "pdf", True),
    (277, "table_S1_GO_volatility", "xls", False),
    (278, "table_S2_transcription_modules", "xls", True),
    (279, "table_S3_module_volatility", "xls", False),
    (280, "table_S4_geneset_coherence", "xls", False),
    (281, "table_S5_paralog_pairs", "xls", True),
    (282, "table_S6_GEO_accessions", "pdf", False),
]

REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = REPO_ROOT / "data" / "wapinski" / "raw"

USER_AGENT = "subnetwork-motifs/0.1 (academic research)"


def fetch(moesm: int, ext: str, out_path: Path) -> None:
    url = f"{SPRINGER_BASE}/41586_2007_BFnature06107_MOESM{moesm}_ESM.{ext}"
    print(f"GET {url}")
    req = Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(req, timeout=120) as resp:
            out_path.write_bytes(resp.read())
    except (URLError, TimeoutError) as e:
        raise SystemExit(f"Download failed for MOESM{moesm}: {e}")
    print(f"  -> {out_path} ({out_path.stat().st_size:,} bytes)")


def main() -> None:
    p = argparse.ArgumentParser(description="Download Wapinski 2007 supplements")
    p.add_argument("--all", action="store_true", help="Download all 7 supplements")
    p.add_argument("--dest", type=Path, default=RAW_DIR)
    p.add_argument("--force", action="store_true", help="Re-download existing files")
    args = p.parse_args()

    args.dest.mkdir(parents=True, exist_ok=True)
    targets = [s for s in SUPPLEMENTS if args.all or s[3]]

    for moesm, name, ext, _ in targets:
        out = args.dest / f"MOESM{moesm}_{name}.{ext}"
        if out.exists() and not args.force:
            print(f"Skip (exists): {out.name}")
            continue
        fetch(moesm, ext, out)


if __name__ == "__main__":
    main()
