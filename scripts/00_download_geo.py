"""Download the public SISBAR data from GEO and verify it against the recorded SHA-256 values.

Files (You et al., Cell Stem Cell 2023; GEO GSE207921 = replicate 1, GSE221592 = replicates
2-3 plus the clone table):
  data/raw/GSE207921/GSE207921_RAW.tar
  data/raw/GSE221592/GSE221592_RAW.tar
  data/raw/GSE221592/GSE221592_Virusbar_processed_data_combined.xlsx
The expected SHA-256 values are the ones recorded in results/*/provenance.txt by the scripts
that used these files (01_clone_counts.py, 02_library_map.py, 03_build_objects.py). A file
that already exists with the right checksum is not downloaded again. The script stops with an
error if a checksum does not match (GEO occasionally repackages RAW tars; if so, compare the
per-library matrices instead and record the new value).

The .tar files are not extracted here; 02_library_map.py and 03_build_objects.py read them
directly. About 3 GB in total.

Usage: python scripts/00_download_geo.py [--check-only]
"""
import argparse
import hashlib
import shutil
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series"
FILES = [
    ("GSE207921", "GSE207921_RAW.tar",
     "8409ac36991415e2f3ff803829423fb6b3ffb6c8d5e7d7e0d3028f995f67a554"),
    ("GSE221592", "GSE221592_RAW.tar",
     "f74f7eb925f7f990433987fbea78896ce203faf85119ac18f5f91c2ef1f2b2fa"),
    ("GSE221592", "GSE221592_Virusbar_processed_data_combined.xlsx",
     "ec64e23ff7ccf2c0e3d245bdb14e16444a948b70d52918892554e56a3961eb2b"),
]


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def url_for(gse, name):
    return f"{GEO}/{gse[:-3]}nnn/{gse}/suppl/{name}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-only", action="store_true", help="only verify files already present")
    a = ap.parse_args()
    bad = 0
    for gse, name, expected in FILES:
        dest = ROOT / "data/raw" / gse / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists() or sha256(dest) != expected:
            if a.check_only:
                print(f"MISSING/CHANGED {dest.relative_to(ROOT)}")
                bad += 1
                continue
            url = url_for(gse, name)
            print(f"downloading {url}", flush=True)
            tmp = dest.with_suffix(dest.suffix + ".part")
            with urllib.request.urlopen(url) as r, open(tmp, "wb") as f:
                shutil.copyfileobj(r, f, length=1 << 20)
            tmp.replace(dest)
        got = sha256(dest)
        ok = got == expected
        bad += not ok
        print(f"{'OK ' if ok else 'BAD'} {dest.relative_to(ROOT)} sha256 {got}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
