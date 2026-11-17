"""A3.2 — inductive RQ2, OT arm (ANALYSIS_PLAN.md amendment A3.2; post-review, exploratory).

Primary cell set, all transitions. Split and sampled early/later cells per seed are those
of 07_rq2 (replayed and saved by scripts/12_rq2_pc_input.py, data/interim/rq2_draws/).

Per seed and replicate, the OT pool is refitted without validation/test information:
  * pool = the RQ1 capped pool (scripts/06_rq1_coupling.py job_cells, background cap 8,000
    per stage, seed 0) minus every cell (both stages) of the seed's validation and test
    clones; i.e. cells of training clones plus background cells;
  * HVGs (seurat_v3, batch = library), normalisation, scaling, 30-PC PCA as RQ1 embed(),
    fitted on the pool; cost = squared Euclidean / its median over the pool;
  * balanced entropic OT solved as in RQ1 (POT sinkhorn_log, eps 0.1 then 0.05 warm-started,
    stopThr 1e-6); potentials cached in data/interim/rq2_inductive/.
Training early cells' OT rows come from the fitted plan. Validation early cells are
projected into the fitted space (pool HVGs, pool scaling mean/SD, clip 10, pool PCA
centring and loadings) and their rows are w_j ∝ exp(log_v_j - C_ij/eps), with C scaled by
the pool's cost median. Targets are drawn from these rows over the pool's later cells (all
of which are training-clone or background cells), rng = np.random.default_rng([seed, 2000]),
training clones then validation clones in table order. Test = the held-out true pairs of 07.
State, min-max scaling and models (ConstMedian, ConstMean, Ridge, CV-tuned RF, Transformer)
as in 07; comparisons vs constant median and vs RF as in 07 (paired over seeds, bootstrap
95% CI, Wilcoxon, Holm within transition, d_z, delta = 0.01). The transductive OT-arm MAE of
07 is carried alongside for context (paired difference inductive - transductive, same seeds).

Seed fallback (A3.2): if seed 0 shows that 20 seeds cannot finish within compute limits,
seeds 0-9 are used. The decision is written to data/interim/rq2_inductive/seed_decision.txt
from seed-0 runtimes only, before any result is summarised.

Usage: python scripts/17_rq2_inductive.py run --transition T --seeds S... | summarise
"""
import argparse
import gc
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


rq1 = _load("rq1", "06_rq1_coupling.py")
rq2 = _load("rq2", "07_rq2_forecast.py")
CELLSET = "primary"
DRAWS = ROOT / "data/interim/rq2_draws"
CACHE = ROOT / "data/interim/rq2_inductive"
OUT = ROOT / "results/17_rq2_inductive"
rq1.CACHE = CACHE            # solve_balanced() caches potentials here


def embed_fit(d_pool, extra_ids):
    """RQ1 embed() on the pool, plus projection of extra cells into the fitted space."""
    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp
    obs = pd.read_csv(rq1.OBS, usecols=["cell_id", "library"]).set_index("cell_id")
    allc = pd.Index(np.concatenate([d_pool.cell_id.to_numpy(), extra_ids]))
    libs = obs.loc[allc, "library"]
    parts = []
    for lib, ids in libs.groupby(libs):
        # backed read: only the needed rows are loaded (memory only; same values)
        a = ad.read_h5ad(rq1.LIBS / f"{lib}.h5ad", backed="r")
        parts.append(a[ids.index.to_numpy()].to_memory())
        a.file.close()
    a_all = ad.concat(parts)[allc.to_numpy()].copy()
    a_all.obs["library"] = libs.to_numpy()
    a = a_all[d_pool.cell_id.to_numpy()].copy()
    a.layers["counts"] = a.X.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=rq1.N_HVG,
                                    batch_key="library", layer="counts")
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        hv = a.var.highly_variable.to_numpy()
        a = a[:, hv].copy()
        sc.pp.scale(a, max_value=10)
        Xs = np.asarray(a.X)
        sc.pp.pca(a, n_comps=rq1.N_PC, random_state=rq1.SEED, svd_solver="arpack")
    Xpool = a.obsm["X_pca"].astype(np.float32)
    mu, sd = a.var["mean"].to_numpy(), a.var["std"].to_numpy()
    centre, PCs = Xs.mean(0), a.varm["PCs"]
    assert np.allclose((Xs - centre) @ PCs, Xpool, atol=1e-3), "projection does not reproduce pool PCs"
    e = a_all[extra_ids].copy()
    M = e.X.astype(np.float32)
    tot = np.asarray(M.sum(1)).ravel()
    L = np.log1p((sp.diags(1e4 / tot) @ M)[:, hv].toarray())
    Z = np.clip((L - mu) / np.where(sd > 0, sd, 1), None, 10)
    return Xpool, ((Z - centre) @ PCs).astype(np.float32)


def run(transition, seeds):
    import torch
    from scipy.spatial.distance import cdist
    torch.set_num_threads(1)
    CACHE.mkdir(parents=True, exist_ok=True)
    obs = rq1.load_obs()
    sc_ = pd.read_csv(rq2.SCORES).set_index("cell_id")
    state = lambda ids: sc_.loc[list(ids), ["Pg", "D"]].to_numpy()
    s0, s1, _ = rq1.TRANSITIONS[transition]
    tt = transition.replace("->", "-")
    for seed in seeds:
        ck = CACHE / f"{CELLSET}_{tt}_s{seed}.csv"
        if ck.exists():
            continue
        t0 = time.time()
        d = pd.read_csv(DRAWS / f"{CELLSET}_{tt}_s{seed}.csv")
        tr, va, te = (d[d.part == p] for p in ("train", "val", "test"))
        rng = np.random.default_rng([seed, 2000])
        tgt_tr, tgt_va, timing = {}, {}, []
        for rep in (1, 2, 3):
            t1 = time.time()
            dp, _, _ = rq1.job_cells(obs, CELLSET, transition, rep, cap=True)
            excl = set(d[(d.rep == rep) & (d.part != "train")].clone)
            dp = dp[~dp.clone.isin(excl)].reset_index(drop=True)
            va_r = va[va.rep == rep]
            Xp, Xva = embed_fit(dp, va_r.early_id.to_numpy())
            gc.collect()        # memory only: release the AnnData objects before OT
            st = dp.stage.to_numpy()
            e_ids, l_ids = dp.cell_id.to_numpy()[st == s0], dp.cell_id.to_numpy()[st == s1]
            X0, X1 = Xp[st == s0], Xp[st == s1]
            # chunked float32 (same values as cdist -> float32; as in 07), memory only
            M = np.empty((len(X0), len(X1)), dtype=np.float32)
            for k in range(0, len(X0), 1000):
                M[k:k + 1000] = cdist(X0[k:k + 1000], X1, "sqeuclidean")
            med = np.float32(np.median(M))
            M /= med
            tag = f"{CELLSET}_{tt}_s{seed}_r{rep}"
            warm, sol = None, None
            for eps in (0.1, 0.05):
                if sol is not None:
                    warm = (sol["log_u"] * np.float32(0.1 / eps), sol["log_v"] * np.float32(0.1 / eps))
                sol = rq1.solve_balanced(M, eps, warm, tag, float("inf"))
            assert bool(sol["converged"]), tag
            lu, lv, niter = sol["log_u"], sol["log_v"], int(sol["niter"])
            epos = pd.Index(e_ids)
            for c, eid in zip(tr[tr.rep == rep].clone, tr[tr.rep == rep].early_id):
                i = epos.get_loc(eid)
                lw = lu[i] - M[i] / np.float32(0.05) + lv
                w = np.exp(lw - lw.max()).astype(np.float64)
                tgt_tr[(rep, c)] = l_ids[rng.choice(len(w), p=w / w.sum())]
            del M, sol
            gc.collect()
            Cva = cdist(Xva, X1, "sqeuclidean").astype(np.float32) / med
            for k, c in enumerate(va_r.clone):
                lw = lv - Cva[k] / np.float32(0.05)
                w = np.exp(lw - lw.max()).astype(np.float64)
                tgt_va[(rep, c)] = l_ids[rng.choice(len(w), p=w / w.sum())]
            timing.append(dict(rep=rep, n_early=len(e_ids), n_later=len(l_ids),
                               niter=niter, seconds=time.time() - t1))
        Btr = np.vstack([state(tr.early_id), state(tr.target_true)])
        lo, hi = Btr.min(0), Btr.max(0)
        scl = lambda A: np.clip((A - lo) / np.where(hi > lo, hi - lo, 1), 0, 1)
        ytr = scl(state([tgt_tr[(r, c)] for r, c in zip(tr.rep, tr.clone)]))
        yva = scl(state([tgt_va[(r, c)] for r, c in zip(va.rep, va.clone)]))
        Xte, yte = scl(state(te.early_id)), scl(state(te.target_true))
        preds, rfp = rq2.fit_models(scl(state(tr.early_id)), ytr, scl(state(va.early_id)), yva, seed)
        res = []
        for m, f in preds.items():
            err = np.abs(f(Xte) - yte)
            res.append(dict(cellset=CELLSET, transition=transition, seed=seed, coupling="OT_inductive",
                            model=m, mae=err.mean(), mae_Pg=err[:, 0].mean(), mae_D=err[:, 1].mean(),
                            n_train=len(tr), n_val=len(va), n_test=len(te),
                            rf_params=str(rfp) if m == "RandomForest" else ""))
        pd.DataFrame(res).to_csv(ck, index=False)
        pd.DataFrame(timing).assign(seed=seed, transition=transition, total_seconds=time.time() - t0).to_csv(
            CACHE / f"timing_{tt}_s{seed}.csv", index=False)
        print(transition, "seed", seed, f"{time.time() - t0:.0f}s", flush=True)


def summarise():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob(f"{CELLSET}_*_s*.csv"))])
    ref = pd.read_csv(ROOT / "results/07_rq2/per_seed.csv")
    ref = ref[(ref.cellset == CELLSET) & (ref.coupling == "OT")]
    df.to_csv(OUT / "per_seed.csv", index=False)
    rows = []
    for tr, g in df.groupby("transition"):
        w = g.pivot(index="seed", columns="model", values="mae")
        wt = ref[ref.transition == tr].pivot(index="seed", columns="model", values="mae").loc[w.index]
        fam = []
        for m in rq2.MODELS:
            if m != "ConstMedian":
                fam.append(dict(model=m, reference="ConstMedian", **rq2.paired_stats(w[m], w["ConstMedian"])))
            if m != "RandomForest":
                fam.append(dict(model=m, reference="RandomForest", **rq2.paired_stats(w[m], w["RandomForest"])))
        adj = rq2.holm([f["p"] for f in fam])
        for f, pa in zip(fam, adj):
            f.update(transition=tr, p_holm=pa, n_seeds=len(w), mae_model=w[f["model"]].mean(),
                     mae_reference=w[f["reference"]].mean(), const_median_mae=w["ConstMedian"].mean(),
                     delta_pct_of_const=100 * f["delta"] / w["ConstMedian"].mean())
            rows.append(f)
        for m in rq2.MODELS:      # context: inductive - transductive, same seeds (no test family)
            st = rq2.paired_stats(w[m], wt[m])
            rows.append(dict(model=m, reference="transductive_OT_same_model", transition=tr,
                             n_seeds=len(w), mae_model=w[m].mean(), mae_reference=wt[m].mean(),
                             const_median_mae=w["ConstMedian"].mean(), p_holm=np.nan,
                             delta_pct_of_const=100 * st["delta"] / w["ConstMedian"].mean(), **st))
    pd.DataFrame(rows).to_csv(OUT / "comparisons.csv", index=False)
    pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("timing_*.csv"))]).to_csv(OUT / "timing.csv", index=False)
    dec = CACHE / "seed_decision.txt"
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    import ot
    import scanpy
    (OUT / "provenance.txt").write_text("\n".join([
        f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
        "plan: ANALYSIS_PLAN.md amendment A3.2 (2026-11-17, post-review, exploratory)",
        f"python {sys.version.split()[0]} scanpy {scanpy.__version__} POT {ot.__version__} "
        f"numpy {np.__version__} pandas {pd.__version__}",
        f"input {rq1.OBS.relative_to(ROOT)} sha256 {rq1.sha256(rq1.OBS)}",
        f"input {rq1.QC.relative_to(ROOT)} sha256 {rq1.sha256(rq1.QC)}",
        f"input {rq2.SCORES.relative_to(ROOT)} sha256 {rq1.sha256(rq2.SCORES)}",
        "input splits/sampled cells: data/interim/rq2_draws/primary_*_s*.csv (07_rq2 replay)",
        "seed decision: " + (dec.read_text().strip() if dec.exists() else "not recorded")]) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "summarise"])
    ap.add_argument("--transition", nargs="*", default=list(rq1.TRANSITIONS))
    ap.add_argument("--seeds", nargs="*", type=int, default=list(range(20)))
    a = ap.parse_args()
    if a.cmd == "run":
        for t in a.transition:
            run(t, a.seeds)
    else:
        summarise()


if __name__ == "__main__":
    main()
