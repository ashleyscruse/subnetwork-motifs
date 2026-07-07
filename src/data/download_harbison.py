"""
Download the Harbison et al. 2004 yeast ChIP-chip binding p-value text files.

Source paper:
    Harbison et al. 2004, "Transcriptional regulatory code of a eukaryotic genome"
    Nature 431:99-104. DOI: 10.1038/nature02800

Hosting history:
    The data was originally distributed from the Young Lab at MIT
    (http://younglab.wi.mit.edu/regulatory_code/, mirrored at
    http://jura.wi.mit.edu/young_public/regulatory_code/). Both URLs are
    now dead. The text-format zip is preserved in the Internet Archive at
    the Wayback URL below, captured 2015-09-16.

What gets downloaded:
    files_for_paper_txt.zip (~92 MB), which contains 22 files. We extract
    the two we need:
        pvalbygene_forpaper_ypd_txt_withheader.txt   (TFs profiled in YPD)
        pvalbygene_forpaper_otherconditions_txt_withheader.txt
                                                     (TFs profiled in stress)
    These together cover the ~203 TFs from the paper across all conditions.

File format (both files):
    Tab-separated, two header rows + body.
    Header row 1: experiment IDs (numeric, ignored).
    Header row 2: '\\t\\t\\t' then TF_condition labels (ABF1_YPD, ARG81_SM, ...)
    Body rows: ORF (systematic, e.g. YAL001C) \\t common name \\t description
               \\t p-value per TF_condition column. Missing values: 'NaN'.

Usage:
    python -m src.data.download_harbison
    python -m src.data.download_harbison --keep-zip  # keep the 92 MB archive
"""

import argparse
import zipfile
from pathlib import Path
from urllib.error import URLError
from urllib.request import Request, urlopen

DEFAULT_URL = (
    "http://web.archive.org/web/20150916025747if_/"
    "http://jura.wi.mit.edu/young_public/regulatory_code/files_for_paper_txt.zip"
)

EXTRACT_FILES = (
    "pvalbygene_forpaper_ypd_txt_withheader.txt",
    "pvalbygene_forpaper_otherconditions_txt_withheader.txt",
)

REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = REPO_ROOT / "data" / "harbison" / "raw"
ZIP_PATH = RAW_DIR / "files_for_paper_txt.zip"

USER_AGENT = "subnetwork-motifs/0.1 (academic research; +https://github.com/)"

MANUAL_INSTRUCTIONS = """
Could not download automatically. To proceed manually:
  1. Locate files_for_paper_txt.zip from any of:
     - The Internet Archive Wayback Machine (search the original Young Lab
       URL: http://jura.wi.mit.edu/young_public/regulatory_code/)
     - YEASTRACT (http://www.yeastract.com/) for derivative edge lists
     - SGD or BioGRID for downstream curated versions
  2. Save the zip to: {zip}
  3. Run this script again (it will skip download if the zip is already there).
"""


def download(url: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"GET {url}")
    req = Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(req, timeout=120) as resp:
            out_path.write_bytes(resp.read())
    except (URLError, TimeoutError) as e:
        raise SystemExit(f"Download failed: {e}\n{MANUAL_INSTRUCTIONS.format(zip=out_path)}")
    print(f"Wrote {out_path.stat().st_size:,} bytes to {out_path}")


def extract(zip_path: Path, files: tuple[str, ...], dest_dir: Path) -> None:
    with zipfile.ZipFile(zip_path) as zf:
        members = set(zf.namelist())
        for fname in files:
            if fname not in members:
                raise SystemExit(f"Expected {fname!r} not found in {zip_path.name}")
            zf.extract(fname, dest_dir)
            out = dest_dir / fname
            print(f"Extracted {out} ({out.stat().st_size:,} bytes)")


def main() -> None:
    p = argparse.ArgumentParser(description="Download Harbison 2004 ChIP-chip data")
    p.add_argument("--url", default=DEFAULT_URL)
    p.add_argument("--zip", type=Path, default=ZIP_PATH)
    p.add_argument("--dest", type=Path, default=RAW_DIR)
    p.add_argument("--keep-zip", action="store_true", help="Keep the zip after extraction")
    p.add_argument("--force", action="store_true", help="Re-download even if zip exists")
    args = p.parse_args()

    if args.zip.exists() and not args.force:
        print(f"Zip already present: {args.zip} (use --force to re-download)")
    else:
        download(args.url, args.zip)

    extract(args.zip, EXTRACT_FILES, args.dest)

    if not args.keep_zip:
        args.zip.unlink()
        print(f"Removed {args.zip} (pass --keep-zip to retain)")


if __name__ == "__main__":
    main()
