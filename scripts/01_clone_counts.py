"""Sanity check: reproduce published cross-stage clone counts (You et al. 2023).

Targets (CLAUDE.md): I->II 768 clones / 2,762 cells; II->III 1,213 / 5,856;
III->IV 559 / 2,816.

Input: GSE221592_Virusbar_processed_data_combined.xlsx (one sheet per stage x Exp).
Viral barcode suffix '-R-E' = replicate R, experiment E.

Clone definitions tried:
  A  each unique viral barcode (per replicate) is a clone
  B  connected components of barcodes that co-occur in any cell (union-find)
  C  exact barcode *set* of a cell is the clone identity
Cell counts = cells (both stages) belonging to cross-stage clones.
"""
import hashlib
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "data/raw/GSE221592/GSE221592_Virusbar_processed_data_combined.xlsx"
OUT = ROOT / "results/01_clone_counts"
OUT.mkdir(parents=True, exist_ok=True)

TARGETS = {1: (768, 2762), 2: (1213, 5856), 3: (559, 2816)}
PAIRS = {1: ("I", "II"), 2: ("II", "III"), 3: ("III", "IV")}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def long_table(sheets):
    rows = []
    for name, df in sheets.items():
        stage, exp = name.replace("Clone_stage_", "").split("_Exp#")
        d = df.assign(stage=stage, exp=int(exp))
        d["cell"] = stage + "|" + d["Cell.Barcode"]
        d["bc"] = d["Virus.Barcode"].str.split(";")
        rows.append(d.explode("bc")[["cell", "stage", "exp", "bc", "Cluster.Id"]])
    t = pd.concat(rows, ignore_index=True)
    t["rep"] = t["bc"].str.extract(r"-(\d+)-\d+$")[0].astype(int)
    return t


class UF:
    def __init__(self):
        self.p = {}

    def find(self, a):
        self.p.setdefault(a, a)
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        self.p[self.find(a)] = self.find(b)


def clone_ids(t, how, max_cells=None):
    """Return DataFrame(cell, stage, clone) for one experiment's long table."""
    t = t.copy()
    if max_cells is not None:  # drop promiscuous barcodes (ambient / contamination)
        n = t.groupby("bc")["cell"].nunique()
        t = t[t["bc"].isin(n[n <= max_cells].index)]
    if how == "A":
        return t.rename(columns={"bc": "clone"})[["cell", "stage", "clone"]]
    if how == "B":
        uf = UF()
        for _, g in t.groupby("cell")["bc"]:
            b = g.tolist()
            uf.find(b[0])
            for x in b[1:]:
                uf.union(b[0], x)
        t["clone"] = t["bc"].map(uf.find)
        return t[["cell", "stage", "clone"]].drop_duplicates()
    if how == "C":
        s = t.groupby(["cell", "stage"])["bc"].apply(lambda b: ";".join(sorted(set(b))))
        return s.rename("clone").reset_index()
    raise ValueError(how)


def cross_stage(c, early, late):
    st = c.groupby("clone")["stage"].agg(set)
    shared = st[st.apply(lambda s: {early, late} <= s)].index
    cells = c[c["clone"].isin(shared)]["cell"].nunique()
    return len(shared), cells


def main():
    sheets = pd.read_excel(XLSX, sheet_name=None, dtype=str)
    t = long_table(sheets)
    res = []
    for exp, (e, l) in PAIRS.items():
        te = t[t["exp"] == exp]
        for how in "ABC":
            for cap in (None, 50, 20):
                for rep in ("all", 1, 2, 3):
                    tr = te if rep == "all" else te[te["rep"] == rep]
                    n_cl, n_cells = cross_stage(clone_ids(tr, how, cap), e, l)
                    res.append(dict(exp=exp, transition=f"{e}->{l}", definition=how,
                                    max_cells_per_bc=cap or "none", replicate=rep,
                                    clones=n_cl, cells=n_cells,
                                    target_clones=TARGETS[exp][0],
                                    target_cells=TARGETS[exp][1]))
    r = pd.DataFrame(res)
    r.to_csv(OUT / "clone_counts.csv", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                            capture_output=True, text=True).stdout.strip()
    (OUT / "provenance.txt").write_text(
        f"input {XLSX.name} sha256 {sha256(XLSX)}\ngit_commit {commit}\n"
        f"python {sys.version.split()[0]} pandas {pd.__version__}\n")
    pd.set_option("display.width", 200)
    print(r.to_string(index=False))


if __name__ == "__main__":
    main()
