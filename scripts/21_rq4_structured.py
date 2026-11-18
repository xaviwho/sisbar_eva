"""A4.1(b) — RQ4 with structured partners (ANALYSIS_PLAN.md amendment A4.1(b)).

Cells, expression and effects as §4 (scripts/09_rq4_mismatch.py): primary II->III later-stage
cells (stage III, Exp 2), per replicate; log1p(CP10k); the 2,000 HVGs and 200 effect genes of
§4 (seurat_v3 on raw counts, batch = library, fitted on the untreated cells; effect genes =
top 200 by rank, the first 200 HVG columns here); linear effect e_i = 0.5 + 0.25 z1_i
(z1 = standardised PC1 of §4) on the 200 effect genes; sparse effect (A3.6) = the same e_i on
the top 20 effect genes only.

Per replicate, seed s (0-99) and effect:
  * cells split at random into a control half (first floor(n/2) of a permutation) and a
    treated half, rng = np.random.default_rng([s, rep]); treated profiles = own profile + e_i
    on the effect genes (the oracle partner of a treated cell is its own pre-effect profile);
  * space: the §4 HVG set (kept fixed; HVGs are not reselected on the treated profiles), scaled
    (scanpy scale, max_value 10) and 30 PCs (arpack, random_state 0) fitted on the pooled
    halves; cost = squared Euclidean treated x control / its median;
  * partners: oracle; balanced entropic OT eps = 0.05 (POT sinkhorn_log, solved at eps 0.1 then
    0.05 warm-started, as RQ1), (i) a partner sampled from the row and (ii) the barycentric
    partner sum_j w_ij L_j; kNN (k = 10) in the same space, sampled and averaged; a uniform
    random control. Sampling draws use rng = np.random.default_rng([s, rep, 7]) in the order
    OT, kNN, random;
  * estimate_i = mean over the effect genes of (treated_i - partner); reported against e_i:
    per-cell RMSE, Pearson r, bias (mean estimate - mean e); OT row statistics: effective
    number of partners exp(H(w_i.)) and weight on the 10 nearest controls.
Summaries per replicate x effect x partner: mean and 2.5/97.5 percentiles over seeds.
Check: the cached expression reproduces the §4 effect genes and effect (09 and A3.6 caches).

Usage: python scripts/21_rq4_structured.py run --rep R --seeds S... | summarise
"""
import argparse
import importlib.util
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


rq4 = _load("rq4", "09_rq4_mismatch.py")
rq1 = _load("rq1", "06_rq1_coupling.py")
CACHE = ROOT / "data/interim/rq4_struct"
OUT = ROOT / "results/21_rq4_structured"
rq1.CACHE = CACHE / "pot"
N_SEEDS = 100


def prepare(rep):
    import anndata as ad
    import scanpy as sc
    ck = CACHE / f"r{rep}.npz"
    if ck.exists():
        z = np.load(ck)
        return z["H"], z["e"]
    obs = pd.read_csv(rq4.OBS, low_memory=False).merge(pd.read_csv(rq4.QC)[["cell_id", "keep"]], on="cell_id")
    d = obs[obs.keep & (obs.exp == 2) & (obs.stage == "III") & (obs.rep == rep)]
    parts = [ad.read_h5ad(rq4.LIBS / f"{lib}.h5ad")[g.cell_id.to_numpy()].copy() for lib, g in d.groupby("library")]
    a = ad.concat(parts)
    a.obs["library"] = np.concatenate([g.library.to_numpy() for _, g in d.groupby("library")])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=rq4.N_HVG, batch_key="library")
        rank = a.var.highly_variable_rank.fillna(np.inf).to_numpy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    order = np.argsort(rank)[:rq4.N_HVG]                  # HVGs ordered by rank
    H = a.X[:, order].toarray().astype(np.float32)
    old = np.load(rq4.CACHE / f"r{rep}.npz")
    Z = np.load(ROOT / "data/interim/rq4_sens" / f"r{rep}.npz")["Z"]
    assert np.allclose(H[:, :rq4.N_EFF], old["L"]), rep
    e = (0.5 + 0.25 * Z[:, 0]).astype(np.float32)
    assert np.allclose(e, old["e"], atol=1e-6), rep
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez(ck, H=H, e=e)
    return H, e


def pcs(Hp):
    import anndata as ad
    import scanpy as sc
    a = ad.AnnData(Hp.copy())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.scale(a, max_value=10)
        sc.pp.pca(a, n_comps=30, random_state=0, svd_solver="arpack")
    return a.obsm["X_pca"].astype(np.float32)


def stats(est, e):
    return dict(rmse=float(np.sqrt(np.mean((est - e) ** 2))), r=float(np.corrcoef(est, e)[0, 1]),
                bias=float(est.mean() - e.mean()))


def run(rep, seeds):
    from scipy.spatial.distance import cdist
    from sklearn.neighbors import NearestNeighbors
    H, e_all = prepare(rep)
    n = len(H)
    rq1.CACHE.mkdir(parents=True, exist_ok=True)
    for s in seeds:
        for effect, G in (("linear", rq4.N_EFF), ("sparse", 20)):
            ck = CACHE / f"r{rep}_s{s}_{effect}.csv"
            if ck.exists():
                continue
            t0 = time.time()
            perm = np.random.default_rng([s, rep]).permutation(n)
            C, T = perm[: n // 2], perm[n // 2:]
            e = e_all[T]
            Hp = np.vstack([H[C], H[T]])
            Hp[len(C):, :G] += e[:, None]
            P = pcs(Hp)
            PC_, PT = P[: len(C)], P[len(C):]
            del Hp
            M = np.empty((len(T), len(C)), dtype=np.float32)
            for k in range(0, len(T), 1000):
                M[k:k + 1000] = cdist(PT[k:k + 1000], PC_, "sqeuclidean")
            M /= np.float32(np.median(M))
            tag = f"r{rep}_s{s}_{effect}"
            sol, warm = None, None
            for eps in (0.1, 0.05):
                if sol is not None:
                    warm = (sol["log_u"] * np.float32(0.1 / eps), sol["log_v"] * np.float32(0.1 / eps))
                sol = rq1.solve_balanced(M, eps, warm, tag, float("inf"))
            conv, niter = bool(sol["converged"]), int(sol["niter"])
            W = rq1.plan_from_potentials(M, 0.05, sol["log_u"], sol["log_v"])
            del M
            W /= W.sum(1, keepdims=True)
            Tg = H[T, :G] + e[:, None]
            Cg = H[C, :G]
            mT = Tg.mean(1)
            nb = NearestNeighbors(n_neighbors=10).fit(PC_).kneighbors(PT, return_distance=False)
            rng = np.random.default_rng([s, rep, 7])
            cw = np.cumsum(W, axis=1)
            u = rng.random(len(T)) * cw[:, -1]
            j_ot = np.minimum((cw < u[:, None]).sum(1), len(C) - 1)
            del cw
            j_knn = nb[np.arange(len(T)), rng.integers(0, 10, size=len(T))]
            j_rand = rng.integers(0, len(C), size=len(T))
            res = {"oracle": mT - H[T, :G].mean(1),
                   "OT_sampled": mT - Cg[j_ot].mean(1),
                   "OT_barycentric": mT - (W @ Cg).mean(1),
                   "kNN_sampled": mT - Cg[j_knn].mean(1),
                   "kNN_averaged": mT - Cg[nb].mean(axis=(1, 2)),
                   "random": mT - Cg[j_rand].mean(1)}
            nz = np.where(W > 0, W, 1.0)
            neff = np.exp(-(W * np.log(nz)).sum(1))
            wnn = W[np.arange(len(T))[:, None], nb].sum(1)
            del W, nz
            rows = [dict(rep=rep, seed=s, effect=effect, partner=k, n_treated=len(T), n_control=len(C),
                         ot_converged=conv, ot_niter=niter, ot_neff_mean=float(neff.mean()),
                         ot_wnn10_mean=float(wnn.mean()), **stats(v.astype(np.float64), e.astype(np.float64)))
                    for k, v in res.items()]
            assert rows[0]["rmse"] < 1e-5          # oracle
            pd.DataFrame(rows).to_csv(ck, index=False)
            for f in rq1.CACHE.glob(f"{tag}_*.npz"):
                f.unlink()                          # potentials not kept (disk)
            print(f"rep {rep} seed {s} {effect} {time.time() - t0:.0f}s", flush=True)


def summarise():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("r*_s*_*.csv"))])
    df.to_csv(OUT / "per_seed.csv.gz", index=False)
    g = df.groupby(["rep", "effect", "partner"])
    out = g[["rmse", "r", "bias", "ot_neff_mean", "ot_wnn10_mean"]].mean()
    for c in ("rmse", "r", "bias"):
        out[c + "_lo"], out[c + "_hi"] = g[c].quantile(0.025), g[c].quantile(0.975)
    out["n_seeds"] = g.size()
    out["all_ot_converged"] = g.ot_converged.all()
    out.reset_index().to_csv(OUT / "summary.csv", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    import ot
    import scanpy
    (OUT / "provenance.txt").write_text("\n".join([
        f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
        "plan: ANALYSIS_PLAN.md amendment A4.1(b) (2026-11-18, second review round, exploratory)",
        f"python {sys.version.split()[0]} scanpy {scanpy.__version__} POT {ot.__version__} numpy {np.__version__}",
        f"input {rq4.OBS.relative_to(ROOT)} sha256 {rq4.sha256(rq4.OBS)}",
        f"input {rq4.QC.relative_to(ROOT)} sha256 {rq4.sha256(rq4.QC)}",
        "check: effect genes and effect equal the §4 (09) and A3.6 caches; oracle error 0"]) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "summarise"])
    ap.add_argument("--rep", type=int, nargs="*", default=[1, 2, 3])
    ap.add_argument("--seeds", type=int, nargs="*", default=list(range(N_SEEDS)))
    a = ap.parse_args()
    if a.cmd == "run":
        for r in a.rep:
            run(r, a.seeds)
    else:
        summarise()


if __name__ == "__main__":
    main()
