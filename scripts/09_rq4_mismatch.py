"""RQ4 — consequence of pairing errors for perturbation models (ANALYSIS_PLAN.md §4).

Simulation, per replicate, on primary II->III later-stage cells (stage III, Exp 2):
  * Expression: CP10k + log1p; 2,000 HVGs (seurat_v3 on raw counts, batch = library) and
    a 30-PC PCA of the scaled HVGs (random_state 0), both fitted on these cells.
  * Effect genes: the top 200 HVGs by seurat_v3 rank.
  * True per-cell effect e_i = 0.5 + 0.25 * z_i, z_i = standardised PC1 score
    (mean shift 0.5 log units, SD across cells 0.25), added to every effect gene:
    treated_i = control_i + e_i on those genes ("treated" = same cells + effect).
  * Pairing with mismatch m: each treated cell is paired with its own control, except a
    random fraction m that is paired with a control drawn uniformly from the other cells.
    Estimated effect = mean over effect genes of (treated_i - control_partner).
  * m in {0, 0.1, ..., 1.0}, plus the mismatch rates measured in RQ1 (OT balanced,
    eps = 0.05, primary set, pooled, per transition): primary 1 - own-clone mass,
    secondary 1 - top-1 recall, each also at its 95% CI bounds.
  * 100 simulation seeds per (replicate, m). Outputs: per-cell effect RMSE, Pearson r with
    the true effect, and population mean-effect bias. Uncertainty: 2.5/97.5 percentiles
    over the 100 seeds (the plan's "bootstrap over clones" has no clone unit here; the
    RQ1 rate's clone-bootstrap CI is propagated by evaluating at its bounds).

Usage: python scripts/09_rq4_mismatch.py grid | summarise
"""
import argparse
import hashlib
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LIBS = ROOT / "data/processed/libs"
OBS = ROOT / "results/03_objects/cell_obs.csv.gz"
QC = ROOT / "results/04_qc/cell_qc.csv.gz"
RQ1 = ROOT / "results/06_rq1/summary.csv"
CACHE = ROOT / "data/interim/rq4"
OUT = ROOT / "results/09_rq4"
N_HVG, N_EFF, N_SIM = 2000, 200, 100
GRID = [round(x, 1) for x in np.arange(0, 1.01, 0.1)]


def prepare(rep):
    import anndata as ad
    import scanpy as sc
    ck = CACHE / f"r{rep}.npz"
    if ck.exists():
        z = np.load(ck)
        return z["L"], z["e"]
    obs = pd.read_csv(OBS, low_memory=False).merge(pd.read_csv(QC)[["cell_id", "keep"]], on="cell_id")
    d = obs[obs.keep & (obs.exp == 2) & (obs.stage == "III") & (obs.rep == rep)]
    parts = [ad.read_h5ad(LIBS / f"{lib}.h5ad")[g.cell_id.to_numpy()].copy()
             for lib, g in d.groupby("library")]
    a = ad.concat(parts)
    a.obs["library"] = np.concatenate([g.library.to_numpy() for _, g in d.groupby("library")])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=N_HVG, batch_key="library")
        rank = a.var.highly_variable_rank.fillna(np.inf).to_numpy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        eff = np.argsort(rank)[:N_EFF]
        L = a.X[:, eff].toarray().astype(np.float32)
        h = a[:, a.var.highly_variable].copy()
        sc.pp.scale(h, max_value=10)
        sc.pp.pca(h, n_comps=30, random_state=0, svd_solver="arpack")
    pc1 = h.obsm["X_pca"][:, 0]
    e = (0.5 + 0.25 * (pc1 - pc1.mean()) / pc1.std()).astype(np.float32)
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez(ck, L=L, e=e)
    return L, e


def simulate(L, e, m, seed):
    rng = np.random.default_rng(seed)
    n = len(e)
    partner = np.arange(n)
    mis = rng.random(n) < m
    r = rng.integers(0, n - 1, size=mis.sum())
    r = r + (r >= np.where(mis)[0])          # uniform over the other cells
    partner[mis] = r
    est = e + (L - L[partner]).mean(1)       # mean over genes of treated_i - control_partner
    return dict(rmse=float(np.sqrt(np.mean((est - e) ** 2))),
                r=float(np.corrcoef(est, e)[0, 1]),
                bias=float(est.mean() - e.mean()), n_cells=n, frac_mismatched=float(mis.mean()))


def measured_rates():
    s = pd.read_csv(RQ1)
    s = s[(s.cellset == "primary") & (~s.uncapped.astype(bool)) & (s.direction == "forward")
          & (s.rep == "pooled") & (s.coupling == "OT_bal_eps0.05")]
    out = []
    for _, r in s.iterrows():
        for kind, col in [("1-own_clone_mass", "m"), ("1-top1_recall", "top1")]:
            for bound, v in [("estimate", r[col]), ("ci_lo", r[col + "_hi"]), ("ci_hi", r[col + "_lo"])]:
                out.append(dict(transition=r.transition, rate_kind=kind, bound=bound, m=1 - v))
    return pd.DataFrame(out)


def run(points):
    rows = []
    for rep in (1, 2, 3):
        L, e = prepare(rep)
        for _, p in points.iterrows():
            for seed in range(N_SIM):
                rows.append(dict(rep=rep, **p.to_dict(), seed=seed, **simulate(L, e, p.m, seed)))
    return pd.DataFrame(rows)


def summarise_df(df, keys):
    g = df.groupby(keys + ["rep"])
    out = g[["rmse", "r", "bias"]].mean()
    for c in ["rmse", "r", "bias"]:
        out[c + "_lo"] = g[c].quantile(0.025)
        out[c + "_hi"] = g[c].quantile(0.975)
    return out.reset_index()


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def provenance(extra):
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    lines = [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
             "plan: ANALYSIS_PLAN.md approved at 65cd70d",
             f"python {sys.version.split()[0]} numpy {np.__version__} pandas {pd.__version__}",
             f"input {OBS.relative_to(ROOT)} sha256 {sha256(OBS)}",
             f"input {QC.relative_to(ROOT)} sha256 {sha256(QC)}"] + extra
    (OUT / "provenance.txt").write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["grid", "summarise"])
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if a.cmd == "grid":
        df = run(pd.DataFrame({"m": GRID, "point": "grid"}))
        df.to_csv(OUT / "grid_per_seed.csv.gz", index=False)
        summarise_df(df, ["m"]).to_csv(OUT / "grid_summary.csv", index=False)
        provenance([])
    else:
        pts = measured_rates()
        pts["point"] = "rq1_measured"
        df = run(pts)
        df.to_csv(OUT / "measured_per_seed.csv.gz", index=False)
        summarise_df(df, ["transition", "rate_kind", "bound", "m"]).to_csv(OUT / "measured_summary.csv", index=False)
        provenance([f"input {RQ1.relative_to(ROOT)} sha256 {sha256(RQ1)}"])


if __name__ == "__main__":
    main()
