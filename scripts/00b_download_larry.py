"""Download LARRY in vitro haematopoiesis data (Weinreb et al., Science 2020; ANALYSIS_PLAN.md
amendment A4.6) from GEO GSE140802, sample GSM4185642 (stateFate_inVitro), and record SHA-256.

Files go to data/raw/larry/ (git-ignored). On the first run the SHA-256 of each file is
written to results/24_larry/download_sha256.csv; later runs verify against that table and stop
on a mismatch.

Usage: python scripts/00b_download_larry.py
"""
import csv
import hashlib
import shutil
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4185nnn/GSM4185642/suppl/"
NAMES = ["cell_barcodes.txt.gz", "clone_matrix.mtx.gz", "gene_names.txt.gz", "library_names.txt.gz",
         "metadata.txt.gz", "normed_counts.mtx.gz"]
DEST = ROOT / "data/raw/larry"
TABLE = ROOT / "results/24_larry/download_sha256.csv"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    known = {}
    if TABLE.exists():
        known = {r["file"]: r["sha256"] for r in csv.DictReader(open(TABLE))}
    rows, bad = [], 0
    for n in NAMES:
        fn = f"GSM4185642_stateFate_inVitro_{n}"
        dest = DEST / fn
        if not dest.exists():
            print("downloading", fn, flush=True)
            tmp = dest.with_suffix(dest.suffix + ".part")
            with urllib.request.urlopen(BASE + fn) as r, open(tmp, "wb") as f:
                shutil.copyfileobj(r, f, length=1 << 20)
            tmp.replace(dest)
        h = sha256(dest)
        if fn in known and known[fn] != h:
            print("MISMATCH", fn, h, "expected", known[fn])
            bad += 1
        rows.append(dict(file=fn, url=BASE + fn, bytes=dest.stat().st_size, sha256=h))
        print("OK" if fn not in known or known[fn] == h else "BAD", fn, h)
    if not TABLE.exists():
        TABLE.parent.mkdir(parents=True, exist_ok=True)
        with open(TABLE, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["file", "url", "bytes", "sha256"])
            w.writeheader()
            w.writerows(rows)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
