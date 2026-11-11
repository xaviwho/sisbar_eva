"""Map clone-table cells to their GEO 10x library and check the Exp3 (III->IV) gap.

Task 1 (CLAUDE.md "Next data tasks"): map each Cell.Barcode `_N` suffix, and for
stage IV the replicate from the viral barcode, to its GSM library, by barcode overlap.
Task 2: one time-boxed pass on the Exp3 discrepancy (555/2,692 vs published
559/2,816). No filter tuning: the clone definition is fixed (exact barcode set).

Inputs: clone-table xlsx; *_barcodes.tsv.gz extracted from both GSE RAW tars into
data/interim/barcodes/.
"""
import gzip
import hashlib
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "data/raw/GSE221592/GSE221592_Virusbar_processed_data_combined.xlsx"
TARS = [ROOT / "data/raw/GSE207921/GSE207921_RAW.tar",
        ROOT / "data/raw/GSE221592/GSE221592_RAW.tar"]
BCDIR = ROOT / "data/interim/barcodes"
OUT = ROOT / "results/02_library_map"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_libraries():
    rx = re.compile(r"(GSM\d+)_Stage_(I{1,3}V?|IV)_(?:Rep(\d)_?)?(\d)?_?Exp(\d)_barcodes")
    libs, bcs = [], {}
    for f in sorted(BCDIR.glob("*_barcodes.tsv.gz")):
        m = rx.match(f.name)
        gsm, stage, rep, sub, exp = m.groups()
        lib = f.name.replace("_barcodes.tsv.gz", "")
        with gzip.open(f, "rt") as fh:
            core = {l.strip().split("-")[0] for l in fh if l.strip()}
        bcs[lib] = core
        libs.append(dict(library=lib, gsm=gsm, stage=stage, rep=int(rep or 1),
                         sublib=int(sub) if sub else None, exp=int(exp),
                         n_cells_matrix=len(core)))
    return pd.DataFrame(libs), bcs


def load_cells():
    sheets = pd.read_excel(XLSX, sheet_name=None, dtype=str)
    rows = []
    for name, df in sheets.items():
        stage, exp = name.replace("Clone_stage_", "").split("_Exp#")
        d = df.copy()
        d["stage"], d["exp"] = stage, int(exp)
        d["core"] = d["Cell.Barcode"].str.extract(r"^([ACGTN]+)")[0]
        d["suffix"] = d["Cell.Barcode"].str.extract(r"(-1(?:_\d+)?)$")[0].fillna("none")
        reps = d["Virus.Barcode"].str.findall(r"-(\d)-\d+(?:;|$)").apply(set)
        assert (reps.apply(len) == 1).all(), "cell with barcodes from >1 replicate"
        d["rep"] = reps.apply(lambda s: int(next(iter(s))))
        d["clone"] = d["Virus.Barcode"].str.split(";").apply(lambda b: ";".join(sorted(set(b))))
        rows.append(d)
    return pd.concat(rows, ignore_index=True)


def main():
    libs, bcs = load_libraries()
    cells = load_cells()

    # ---- Task 1: suffix (+replicate) group -> library by barcode overlap
    maps = []
    for (stage, exp, suf, rep), g in cells.groupby(["stage", "exp", "suffix", "rep"]):
        cand = libs[(libs.stage == stage) & (libs.exp == exp)]
        core = set(g["core"])
        for _, L in cand.iterrows():
            maps.append(dict(stage=stage, exp=exp, suffix=suf, rep_from_virus=rep,
                             n_cells=len(core), library=L.library, gsm=L.gsm,
                             lib_rep=L.rep, frac_found=len(core & bcs[L.library]) / len(core)))
    m = pd.DataFrame(maps)
    best = m.loc[m.groupby(["stage", "exp", "suffix", "rep_from_virus"])["frac_found"].idxmax()]
    best = best.assign(rep_consistent=best.rep_from_virus == best.lib_rep)
    m.to_csv(OUT / "suffix_library_overlap_all.csv", index=False)
    best.to_csv(OUT / "suffix_library_map.csv", index=False)

    # ---- per-cell library assignment (stage IV: library chosen by viral-barcode replicate)
    key = best.set_index(["stage", "exp", "suffix", "rep_from_virus"])["library"]
    cells["library"] = [key[(s, e, x, r)] for s, e, x, r in
                        zip(cells.stage, cells.exp, cells.suffix, cells.rep)]
    cells["in_matrix"] = [c in bcs[l] for c, l in zip(cells.core, cells.library)]
    iv_libs = libs[libs.stage == "IV"].library.tolist()
    cells["n_stageIV_libs_with_bc"] = [
        sum(c in bcs[l] for l in iv_libs) if s == "IV" else pd.NA
        for c, s in zip(cells.core, cells.stage)]
    cells[["stage", "exp", "Cell.Barcode", "core", "suffix", "rep", "library",
           "in_matrix", "Cluster.Id", "clone"]].to_csv(OUT / "cell_library_map.csv.gz", index=False)

    # ---- Task 2: Exp3 discrepancy, one pass
    e3 = cells[cells.exp == 3]
    st = e3.groupby("clone")["stage"].agg(set)
    shared = st[st.apply(lambda s: {"III", "IV"} <= s)].index
    x = e3[e3.clone.isin(shared)]
    iv = cells[cells.stage == "IV"]
    gap = dict(
        clones_all=len(shared), cells_all=len(x),
        clones_in_matrix=x[x.in_matrix].groupby("clone")["stage"].agg(set)
        .apply(lambda s: {"III", "IV"} <= s).sum(),
        cells_in_matrix=int(x.in_matrix.sum()),
        stageIV_cells=len(iv), stageIV_not_in_own_matrix=int((~iv.in_matrix).sum()),
        stageIV_bc_in_gt1_IV_library=int((iv.n_stageIV_libs_with_bc > 1).sum()),
        stageIV_bc_in_0_IV_library=int((iv.n_stageIV_libs_with_bc == 0).sum()),
        published_clones=559, published_cells=2816)
    pd.Series(gap).to_csv(OUT / "exp3_gap_check.csv", header=["value"])

    # ---- coverage summary
    cov = cells.groupby(["stage", "exp", "rep", "library"]).agg(
        n_table=("core", "size"), n_in_matrix=("in_matrix", "sum")).reset_index()
    cov = cov.merge(libs[["library", "gsm", "n_cells_matrix"]], on="library")
    cov.to_csv(OUT / "library_coverage.csv", index=False)

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                            capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT,
                           capture_output=True, text=True).stdout.strip()
    prov = [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
            f"python {sys.version.split()[0]} pandas {pd.__version__}",
            f"input {XLSX.name} sha256 {sha256(XLSX)}"]
    prov += [f"input {t.name} sha256 {sha256(t)}" for t in TARS]
    (OUT / "provenance.txt").write_text("\n".join(prov) + "\n")

    pd.set_option("display.width", 220)
    print("total cells in matrices:", libs.n_cells_matrix.sum())
    print(best[["stage", "exp", "suffix", "rep_from_virus", "n_cells", "gsm",
                "library", "frac_found", "rep_consistent"]].to_string(index=False))
    print(cov.to_string(index=False))
    print(pd.Series(gap).to_string())


if __name__ == "__main__":
    main()
