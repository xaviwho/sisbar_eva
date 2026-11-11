"""Final QC cell set and doublet scores (UPDATES.md Q2 option b, Q3 option i).

Rules, fixed before running (author decisions 2026-11-11):
  PRIMARY cell set:  keep = n_genes >= 1,000 AND pct_mt <= 15.  No UMI cap, no doublet
                     removal (Scrublet is unstable on these libraries; see Q3).
  SENSITIVITY flag:  doublet_flag = cell is in its library's top `rate` fraction by
                     Scrublet doublet score, rate = 0.008 * n_barcodes / 1000 (10x
                     guideline ~0.8% per 1,000 cells). RQ1-RQ3 are re-run with these
                     cells excluded as a pre-specified sensitivity analysis.

Scrublet settings (fixed): scanpy.pp.scrublet on raw counts of the library's primary
cells, expected_doublet_rate = rate, mean_center=False (sparse PCA; required to fit in
memory), random_state=0. Scrublet's own threshold/predicted_doublet are recorded but
NOT used. Each library runs in its own process to bound memory.

Inputs : data/processed/libs/*.h5ad, results/03_objects/cell_obs.csv.gz
Outputs: results/04_qc/cell_qc.csv.gz   cell_id, keep, doublet_score, doublet_flag,
                                         scrublet_predicted (recorded only)
         results/04_qc/qc_by_library.csv, clone_counts.csv, summary.csv, provenance.txt
"""
import hashlib
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LIBS = ROOT / "data/processed/libs"
OBS = ROOT / "results/03_objects/cell_obs.csv.gz"
TMP = ROOT / "data/interim/scrublet"
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


def primary(obs):
    return (obs.n_genes_by_counts >= MIN_GENES) & (obs.pct_counts_mt <= MAX_MT)


def score_library(lib):
    """Child process: Scrublet on one library's primary cells."""
    import warnings
    import anndata as ad
    import scanpy as sc
    warnings.filterwarnings("ignore")
    obs = pd.read_csv(OBS, usecols=["cell_id", "library", "n_genes_by_counts", "pct_counts_mt"])
    obs = obs[obs.library == lib].set_index("cell_id")
    a = ad.read_h5ad(LIBS / f"{lib}.h5ad")
    n_bc = a.n_obs
    a = a[primary(obs.loc[a.obs_names]).to_numpy()].copy()
    rate = RATE_PER_1000 * n_bc / 1000
    sc.pp.scrublet(a, expected_doublet_rate=rate, mean_center=False, random_state=0, verbose=False)
    s = a.obs.doublet_score
    cut = np.quantile(s, 1 - rate)
    r = pd.DataFrame({"cell_id": a.obs_names, "doublet_score": s.to_numpy(),
                      "doublet_flag": (s > cut).to_numpy(),
                      "scrublet_predicted": a.obs.predicted_doublet.to_numpy()})
    r.to_csv(TMP / f"{lib}.csv.gz", index=False)
    meta = dict(library=lib, n_barcodes=n_bc, n_primary=a.n_obs, expected_rate=round(rate, 4),
                flag_cutoff=round(float(cut), 4), n_flagged=int(r.doublet_flag.sum()),
                scrublet_threshold=round(float(a.uns["scrublet"]["threshold"]), 4),
                scrublet_n_predicted=int(r.scrublet_predicted.sum()))
    pd.Series(meta).to_csv(TMP / f"{lib}.meta.csv", header=["value"])


def clone_counts(obs, mask, label):
    rows = []
    for name, (s0, s1, e) in PAIRS.items():
        d = obs[mask & obs.in_clone_table & (obs.exp == e) & obs.stage.isin([s0, s1])]
        st = d.groupby("clone")["stage"].agg(set)
        shared = st[st.apply(lambda x: {s0, s1} <= x)].index
        rows.append(dict(transition=name, cell_set=label, clones=len(shared),
                         cells=int(d.clone.isin(shared).sum()),
                         published_clones=PUBLISHED[name][0], published_cells=PUBLISHED[name][1]))
    return rows


def main():
    TMP.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    obs = pd.read_csv(OBS, low_memory=False)
    for lib in sorted(obs.library.unique()):
        if not (TMP / f"{lib}.meta.csv").exists():
            subprocess.run([sys.executable, "-W", "ignore", __file__, "--lib", lib], check=True)
        print("scored", lib, flush=True)

    q = pd.concat([pd.read_csv(TMP / f"{l}.csv.gz") for l in sorted(obs.library.unique())])
    meta = pd.DataFrame([pd.read_csv(TMP / f"{l}.meta.csv", index_col=0)["value"]
                         for l in sorted(obs.library.unique())])
    obs = obs.merge(q, on="cell_id", how="left")
    obs["keep"] = primary(obs)
    obs["doublet_flag"] = obs.doublet_flag.astype("boolean").fillna(False).astype(bool)
    obs["scrublet_predicted"] = obs.scrublet_predicted.astype("boolean").fillna(False).astype(bool)
    assert obs.loc[obs.keep, "doublet_score"].notna().all()
    obs[["cell_id", "keep", "doublet_score", "doublet_flag", "scrublet_predicted"]].to_csv(
        OUT / "cell_qc.csv.gz", index=False)

    ct = obs.in_clone_table
    k = obs.groupby("library").agg(
        n_keep=("keep", "sum"), n_clone_table=("in_clone_table", "sum"),
        n_clone_table_keep=("keep", lambda s: int((s & ct.loc[s.index]).sum())),
        n_clone_table_flagged=("doublet_flag", lambda s: int((s & ct.loc[s.index]).sum())),
        median_umi_flagged=("total_counts", lambda s: s[obs.loc[s.index, "doublet_flag"]].median()),
        median_umi_keep=("total_counts", lambda s: s[obs.loc[s.index, "keep"]].median()))
    meta = meta.merge(k.reset_index(), on="library")
    meta.to_csv(OUT / "qc_by_library.csv", index=False)

    rows = clone_counts(obs, pd.Series(True, index=obs.index), "all_clone_table_cells")
    rows += clone_counts(obs, obs.keep, "primary")
    rows += clone_counts(obs, obs.keep & ~obs.doublet_flag, "sensitivity_doublets_excluded")
    cc = pd.DataFrame(rows)
    cc.to_csv(OUT / "clone_counts.csv", index=False)

    summary = pd.Series(dict(
        primary_rule=f"genes>={MIN_GENES}; mt%<={MAX_MT}; no UMI cap; no doublet removal",
        sensitivity_rule="exclude top 0.8%-per-1,000-barcodes fraction by Scrublet score per library",
        barcodes=len(obs), primary_cells=int(obs.keep.sum()),
        sensitivity_cells=int((obs.keep & ~obs.doublet_flag).sum()),
        published_hq_cells=173564, clone_table_cells=int(ct.sum()),
        clone_table_primary=int((ct & obs.keep).sum()),
        clone_table_sensitivity=int((ct & obs.keep & ~obs.doublet_flag).sum()),
        scrublet_auto_predicted_total=int(obs.scrublet_predicted.sum()),
        doublet_score_vs_umi_spearman=round(float(obs.loc[obs.keep, ["doublet_score", "total_counts"]]
                                                  .corr(method="spearman").iloc[0, 1]), 3)))
    summary.to_csv(OUT / "summary.csv", header=["value"])

    import anndata as ad
    import scanpy as sc
    import scipy
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    prov = [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
            f"python {sys.version.split()[0]} scanpy {sc.__version__} anndata {ad.__version__} "
            f"pandas {pd.__version__} numpy {np.__version__} scipy {scipy.__version__}",
            "environment: Claude cloud workspace (Linux)",
            f"input {OBS.relative_to(ROOT)} sha256 {sha256(OBS)}"]
    (OUT / "provenance.txt").write_text("\n".join(prov) + "\n")

    pd.set_option("display.width", 220)
    print(summary.to_string())
    print(meta.to_string(index=False))
    print(cc.to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--lib":
        score_library(sys.argv[2])
    else:
        main()
