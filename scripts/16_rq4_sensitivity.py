"""A3.6 — RQ4 effect-model sensitivity (ANALYSIS_PLAN.md amendment A3.6; post-review, exploratory).

As scripts/09_rq4_mismatch.py (same cells, HVGs, PCA, mismatch procedure, m grid, RQ1-measured
rates with CI bounds, 100 seeds, outputs), except for the effect model:
  (a) sparse      e_i = 0.5 + 0.25 z1_i, added to the top 20 of the 200 effect genes only
                  (the estimate averages over those 20 genes);
  (b) nonlinear   e_i = 0.5 + 0.25 (z1_i^2 - 1) / sqrt(2);
  (c) multi_axis  e_i = 0.5 + 0.25 (z1_i + z2_i + z3_i) / sqrt(3);
  (d) reversed    e_i = -0.5 + 0.25 z1_i.
z_k = standardised PC k score (population SD, as in 09). Effects (b)-(d) use all 200 genes.
The 09 baseline ("linear") is recomputed alongside as a reproduction check: its per-seed
values must equal results/09_rq4 exactly, or the script stops.

Usage: python scripts/16_rq4_sensitivity.py run | summarise
"""
import importlib.util
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rq4", ROOT / "scripts/09_rq4_mismatch.py")
rq4 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rq4)

CACHE = ROOT / "data/interim/rq4_sens"
OUT = ROOT / "results/16_rq4_sens"
MODELS = ["linear", "sparse", "nonlinear", "multi_axis", "reversed"]


def prepare(rep):
    """09's preparation, also keeping PCs 1-3. Checked against the 09 cache."""
    import anndata as ad
    import scanpy as sc
    ck = CACHE / f"r{rep}.npz"
    if ck.exists():
        z = np.load(ck)
        return z["L"], z["Z"]
    obs = pd.read_csv(rq4.OBS, low_memory=False).merge(pd.read_csv(rq4.QC)[["cell_id", "keep"]], on="cell_id")
    d = obs[obs.keep & (obs.exp == 2) & (obs.stage == "III") & (obs.rep == rep)]
    parts = [ad.read_h5ad(rq4.LIBS / f"{lib}.h5ad")[g.cell_id.to_numpy()].copy()
             for lib, g in d.groupby("library")]
    a = ad.concat(parts)
    a.obs["library"] = np.concatenate([g.library.to_numpy() for _, g in d.groupby("library")])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=rq4.N_HVG, batch_key="library")
        rank = a.var.highly_variable_rank.fillna(np.inf).to_numpy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        eff = np.argsort(rank)[:rq4.N_EFF]          # ordered by seurat_v3 rank
        L = a.X[:, eff].toarray().astype(np.float32)
        h = a[:, a.var.highly_variable].copy()
        sc.pp.scale(h, max_value=10)
        sc.pp.pca(h, n_comps=30, random_state=0, svd_solver="arpack")
    P = h.obsm["X_pca"][:, :3]
    Z = ((P - P.mean(0)) / P.std(0)).astype(np.float64)
    old = np.load(rq4.CACHE / f"r{rep}.npz")
    assert np.allclose(old["L"], L) and np.allclose(old["e"], (0.5 + 0.25 * Z[:, 0]).astype(np.float32), atol=1e-6), rep
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez(ck, L=L, Z=Z)
    return L, Z


def effect(model, Z):
    z1, z2, z3 = Z[:, 0], Z[:, 1], Z[:, 2]
    if model in ("linear", "sparse"):
        e = 0.5 + 0.25 * z1
    elif model == "nonlinear":
        e = 0.5 + 0.25 * (z1 ** 2 - 1) / np.sqrt(2)
    elif model == "multi_axis":
        e = 0.5 + 0.25 * (z1 + z2 + z3) / np.sqrt(3)
    else:
        e = -0.5 + 0.25 * z1
    return e.astype(np.float32)


def run():
    pts = pd.concat([pd.DataFrame({"m": rq4.GRID, "point": "grid"}),
                     rq4.measured_rates().assign(point="rq1_measured")], ignore_index=True)
    ref_g = pd.read_csv(rq4.OUT / "grid_per_seed.csv.gz")
    ref_m = pd.read_csv(rq4.OUT / "measured_per_seed.csv.gz")
    rows = []
    for rep in (1, 2, 3):
        L, Z = prepare(rep)
        for model in MODELS:
            e = effect(model, Z)
            Lm = L[:, :20] if model == "sparse" else L
            for _, p in pts.iterrows():
                for seed in range(rq4.N_SIM):
                    rows.append(dict(rep=rep, effect_model=model, **p.to_dict(), seed=seed,
                                     **rq4.simulate(Lm, e, p.m, seed)))
        print("rep", rep, "done", flush=True)
    df = pd.DataFrame(rows)
    # reproduction check: linear model == 09 outputs
    lin = df[df.effect_model == "linear"]
    for ref, sel in ((ref_g, lin.point == "grid"), (ref_m, lin.point == "rq1_measured")):
        a = lin[sel].reset_index(drop=True)
        assert len(a) == len(ref) and np.allclose(a.rmse, ref.rmse) and np.allclose(a.r, ref.r) \
            and np.allclose(a.bias, ref.bias, atol=1e-7), "linear effect does not reproduce 09"
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "per_seed.csv.gz", index=False)


def summarise():
    df = pd.read_csv(OUT / "per_seed.csv.gz")
    for c in ("transition", "rate_kind", "bound"):
        df[c] = df[c].fillna("")
    out = rq4.summarise_df(df, ["effect_model", "point", "transition", "rate_kind", "bound", "m"])
    out.to_csv(OUT / "summary.csv", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    lines = [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
             "plan: ANALYSIS_PLAN.md amendment A3.6 (2026-11-17, post-review, exploratory)",
             f"python {sys.version.split()[0]} numpy {np.__version__} pandas {pd.__version__}",
             f"input {rq4.OBS.relative_to(ROOT)} sha256 {rq4.sha256(rq4.OBS)}",
             f"input {rq4.QC.relative_to(ROOT)} sha256 {rq4.sha256(rq4.QC)}",
             f"input {rq4.RQ1.relative_to(ROOT)} sha256 {rq4.sha256(rq4.RQ1)}",
             "check: 'linear' effect reproduces results/09_rq4 per-seed values exactly"]
    (OUT / "provenance.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    {"run": run, "summarise": summarise}[sys.argv[1]]()
