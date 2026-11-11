"""Apply the chosen QC rule (UPDATES.md Q2, option b) and score doublets.

Rule, fixed before running (author decision 2026-11-11):
  keep = n_genes >= 1,000  AND  pct_mt <= 15  AND  not a predicted doublet.
  No UMI cap (the paper's 30,000 cap removes 40% of stage-IV clone cells).

Doublets: Scrublet (scanpy.pp.scrublet) per 10x library, on raw counts of cells that
pass the gene and mt filters, random_state=0. Expected doublet rate per library is set
from the 10x Chromium guideline of ~0.8% per 1,000 cells loaded-and-recovered, using
the library's barcode count: rate = 0.008 * n_barcodes / 1000. Scrublet's automatic
threshold is used. If Scrublet cannot set one (no bimodal score distribution), the
pre-specified fallback flags the top `rate` fraction of cells by doublet score; such
libraries are listed in the output.

Inputs : data/processed/libs/*.h5ad, results/03_objects/cell_obs.csv.gz
Outputs: results/04_qc/cell_qc.csv.gz      cell_id, doublet_score, predicted_doublet,
                                            threshold_source, keep
         results/04_qc/qc_by_library.csv
         results/04_qc/clone_counts.csv     cross-stage clones under the final rule
         results/04_qc/summary.csv, provenance.txt
"""
import hashlib
import subprocess
import sys
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
LIBS = ROOT / "data/processed/libs"
OBS = ROOT / "results/03_objects/cell_obs.csv.gz"
OUT = ROOT / "results/04_qc"
MIN_GENES, MAX_MT, RATE_PER_1000 = 1000, 15.0, 0.008
PUBLISHED = {"I->II": (768, 2762), "II->III": (1213, 5856), "III->IV": (559, 2816)}
PAIRS = {"I->II": ("I", "II", 1), "II->III": ("II", "III", 2), "III->IV": ("III", "IV", 3)}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def clone_counts(obs, mask, label):
    rows = []
    for name, (s0, s1, e) in PAIRS.items():
        d = obs[mask & obs.in_clone_table & (obs.exp == e) & obs.stage.isin([s0, s1])]
        st = d.groupby("clone")["stage"].agg(set)
        shared = st[st.apply(lambda x: {s0, s1} <= x)].index
        rows.append(dict(transition=name, rule=label, clones=len(shared),
                         cells=int(d.clone.isin(shared).sum()),
                         published_clones=PUBLISHED[name][0], published_cells=PUBLISHED[name][1]))
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    obs = pd.read_csv(OBS, low_memory=False)
    obs["pre_filter"] = (obs.n_genes_by_counts >= MIN_GENES) & (obs.pct_counts_mt <= MAX_MT)

    recs, libsum = [], []
    for lib in sorted(obs.library.unique()):
        a = ad.read_h5ad(LIBS / f"{lib}.h5ad")
        n_bc = a.n_obs
        a = a[obs.set_index("cell_id").loc[a.obs_names, "pre_filter"].to_numpy()].copy()
        rate = RATE_PER_1000 * n_bc / 1000
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sc.pp.scrublet(a, expected_doublet_rate=rate, random_state=0, verbose=False)
        thr = a.uns["scrublet"].get("threshold")
        src = "scrublet_auto"
        if thr is None or not np.isfinite(thr):
            thr = float(np.quantile(a.obs.doublet_score, 1 - rate))
            a.obs["predicted_doublet"] = a.obs.doublet_score > thr
            src = "fallback_top_rate"
        r = a.obs[["doublet_score", "predicted_doublet"]].copy()
        r["threshold_source"] = src
        recs.append(r.reset_index(names="cell_id"))
        libsum.append(dict(library=lib, n_barcodes=n_bc, n_pre_filter=a.n_obs,
                           expected_rate=round(rate, 4), threshold=round(float(thr), 4),
                           threshold_source=src, n_doublets=int(a.obs.predicted_doublet.sum()),
                           frac_doublets=round(float(a.obs.predicted_doublet.mean()), 4),
                           warnings="; ".join(sorted({str(x.message)[:80] for x in w}))))
        print(lib, libsum[-1]["n_doublets"], src, flush=True)
        del a

    q = pd.concat(recs, ignore_index=True)
    obs = obs.merge(q, on="cell_id", how="left")
    obs["predicted_doublet"] = obs.predicted_doublet.fillna(False).astype(bool)
    obs["keep"] = obs.pre_filter & ~obs.predicted_doublet
    obs[["cell_id", "doublet_score", "predicted_doublet", "threshold_source", "keep"]].to_csv(
        OUT / "cell_qc.csv.gz", index=False)

    ls = pd.DataFrame(libsum)
    k = obs.groupby("library").agg(n_keep=("keep", "sum"),
                                   n_clone_table=("in_clone_table", "sum"),
                                   n_clone_table_keep=("keep", lambda s: int((s & obs.loc[s.index, "in_clone_table"]).sum())))
    ls = ls.merge(k.reset_index(), on="library")
    ls.to_csv(OUT / "qc_by_library.csv", index=False)

    rows = clone_counts(obs, pd.Series(True, index=obs.index), "all_clone_table_cells")
    rows += clone_counts(obs, obs.pre_filter, "genes_mt_only")
    rows += clone_counts(obs, obs.keep, "final_genes_mt_doublets")
    cc = pd.DataFrame(rows)
    cc.to_csv(OUT / "clone_counts.csv", index=False)

    ct = obs.in_clone_table
    summary = pd.Series(dict(
        rule=f"genes>={MIN_GENES}; mt%<={MAX_MT}; Scrublet doublets removed; no UMI cap",
        barcodes=len(obs), pass_genes_mt=int(obs.pre_filter.sum()),
        doublets_removed=int((obs.pre_filter & obs.predicted_doublet).sum()),
        cells_kept=int(obs.keep.sum()), published_hq_cells=173564,
        clone_table_cells=int(ct.sum()), clone_table_kept=int((ct & obs.keep).sum()),
        clone_table_removed_genes_mt=int((ct & ~obs.pre_filter).sum()),
        clone_table_removed_doublet=int((ct & obs.pre_filter & obs.predicted_doublet).sum()),
        libraries_using_fallback=int((ls.threshold_source != "scrublet_auto").sum())))
    summary.to_csv(OUT / "summary.csv", header=["value"])

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    prov = [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
            f"python {sys.version.split()[0]} scanpy {sc.__version__} anndata {ad.__version__} "
            f"pandas {pd.__version__} numpy {np.__version__}",
            "environment: Claude cloud workspace (Linux)",
            f"input {OBS.relative_to(ROOT)} sha256 {sha256(OBS)}"]
    (OUT / "provenance.txt").write_text("\n".join(prov) + "\n")

    pd.set_option("display.width", 220)
    print(summary.to_string())
    print(ls.drop(columns="warnings").to_string(index=False))
    print(cc.to_string(index=False))


if __name__ == "__main__":
    main()
