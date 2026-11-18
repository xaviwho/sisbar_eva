"""A4.5 — within-dataset generality of RQ1 (ANALYSIS_PLAN.md amendment A4.5).

Primary cell set, capped RQ1 pools (same cells as RQ1), per transition x replicate.
  * Spaces: 15 PCs of 2,000 HVGs ("hvg2000_pc15"), 50 PCs of 2,000 HVGs ("hvg2000_pc50"),
    30 PCs of 4,000 HVGs ("hvg4000_pc30"); each built exactly as RQ1 embed() (seurat_v3 on raw
    counts, batch = library; CP10k, log1p; scale max 10; arpack PCA, random_state 0) with only
    the HVG number / PC number changed. In each: balanced OT eps 0.05 (POT sinkhorn_log,
    eps 0.1 then 0.05 warm-started, as RQ1) and kNN (k = 10).
  * Growth-informed OT in the RQ1 space (cached RQ1 embedding): row marginal
    a_i ∝ exp(lambda z_i), z_i = standardised (over the pool's early cells) mean log1p(CP10k)
    expression of the Tirosh et al. S-phase and G2/M genes (scanpy cell-cycle tutorial list,
    data/raw/cellcycle/regev_lab_cell_cycle_genes.txt, genes present in the data), lambda in
    {0.5, 1}; column marginal uniform; balanced (POT sinkhorn_log, eps 0.1 then 0.05
    warm-started, same stopping rule as RQ1) and unbalanced tau = 0.9 (RQ1 log-domain solver
    with log a in place of the uniform row marginal), eps = 0.05.
Metrics (all as RQ1, foreground early cells, forward): own-clone mass m vs chance m0, AUROC,
top-k; fate agreement (TVD per early cluster with >= 20 clones) using the RQ1 later-cell labels
(k = 15 transfer in the RQ1 space, so only the coupling differs); A4.1(a) partner error D and
R, computed in the RQ1 30-PC space for every coupling (one common metric space), with the
independent coupling and the sister floor as in scripts/20_partner_error.py. Clone-weighted
summaries, pooled and per replicate, 1,000-resample clone bootstrap (np.random.default_rng(0)
per row).
Check: kNN in the RQ1 space reproduces RQ1's kNN own-clone mass per cell.

Usage: python scripts/23_rq1_generality.py run --transition T --rep R | summarise
"""
import argparse
import importlib.util
import json
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rq1", ROOT / "scripts/06_rq1_coupling.py")
rq1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rq1)
RQ1_CACHE = rq1.CACHE

CELLSET = "primary"
CACHE = ROOT / "data/interim/rq1_generality"
OUT = ROOT / "results/23_rq1_generality"
CC_GENES = ROOT / "data/raw/cellcycle/regev_lab_cell_cycle_genes.txt"
SPACES = {"hvg2000_pc15": (2000, 15), "hvg2000_pc50": (2000, 50), "hvg4000_pc30": (4000, 30)}
LAMBDAS = (0.5, 1.0)
N_BOOT = 1000


def load_counts(d):
    import anndata as ad
    parts = []
    for lib, ids in d.groupby("library").cell_id:
        a = ad.read_h5ad(rq1.LIBS / f"{lib}.h5ad", backed="r")
        parts.append(a[ids.to_numpy()].to_memory())
        a.file.close()
    a = ad.concat(parts)[d.cell_id.to_numpy()].copy()
    a.obs["library"] = d.library.to_numpy()
    return a


def embed(a0, n_hvg, n_pc):
    import scanpy as sc
    a = a0.copy()
    a.layers["counts"] = a.X.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=n_hvg, batch_key="library", layer="counts")
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        a = a[:, a.var.highly_variable].copy()
        sc.pp.scale(a, max_value=10)
        sc.pp.pca(a, n_comps=n_pc, random_state=rq1.SEED, svd_solver="arpack")
    return a.obsm["X_pca"].astype(np.float32)


def cc_score(a0, early_mask):
    import scanpy as sc
    # var_names are Ensembl IDs; the list holds symbols -> match on var['symbol']
    wanted = set(CC_GENES.read_text().split())
    cols = np.where(a0.var["symbol"].isin(wanted).to_numpy())[0]
    assert len(cols) > 0
    a = a0[early_mask].copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    x = np.asarray(a.X[:, cols].mean(1)).ravel()
    return (x - x.mean()) / x.std(), len(cols)


def cost(X0, X1):
    from scipy.spatial.distance import cdist
    M = np.empty((len(X0), len(X1)), dtype=np.float32)
    for k in range(0, len(X0), 1000):
        M[k:k + 1000] = cdist(X0[k:k + 1000], X1, "sqeuclidean")
    M /= np.float32(np.median(M))
    return M


def solve_bal(M, eps, a, warm):
    import ot
    n1 = M.shape[1]
    b = np.full(n1, 1 / n1, dtype=np.float32)
    lu, lv = warm if warm is not None else (np.zeros(M.shape[0], np.float32), np.zeros(n1, np.float32))
    it, err = 0, np.inf
    while it < rq1.MAX_ITER:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, log = ot.sinkhorn(a, b, M, eps, method="sinkhorn_log", numItermax=rq1.CHUNK,
                                 stopThr=rq1.STOP, warmstart=(lu, lv), log=True, warn=False)
        lu, lv = log["log_u"].astype(np.float32), log["log_v"].astype(np.float32)
        it += int(log["niter"]) + 1
        err = float(log["err"][-1]) if len(log["err"]) else np.inf
        if err < rq1.STOP:
            break
    return lu, lv, it, err < rq1.STOP


def solve_bal_path(M, a):
    lu, lv, it1, _ = solve_bal(M, 0.1, a, None)
    lu, lv, it2, conv = solve_bal(M, 0.05, a, (lu * np.float32(2), lv * np.float32(2)))
    return lu, lv, it1 + it2, conv


def solve_unb(M, eps, a):
    from scipy.special import logsumexp
    n0, n1 = M.shape
    la, lb = np.log(a).astype(np.float32), np.float32(-np.log(n1))
    reg_m = eps * rq1.TAU / (1 - rq1.TAU)
    phi = np.float32(reg_m / (reg_m + eps))
    logK = (M / np.float32(-eps)).astype(np.float32)
    lu, lv = np.zeros(n0, np.float32), np.zeros(n1, np.float32)
    it, err = 0, np.inf
    while it < rq1.MAX_ITER:
        lu_new = phi * (la - logsumexp(logK + lv[None, :], axis=1)).astype(np.float32)
        lv_new = phi * (lb - logsumexp(logK + lu_new[:, None], axis=0)).astype(np.float32)
        err = float(max((np.abs(lu_new - lu) / (1 + np.abs(lu))).max(),
                        (np.abs(lv_new - lv) / (1 + np.abs(lv))).max()))
        lu, lv, it = lu_new, lv_new, it + 1
        if err < rq1.STOP:
            break
    del logK
    return lu, lv, it, err < rq1.STOP


def partner_D(W, efg, ce, cl, X1r):
    sq1 = (X1r ** 2).sum(1)
    out, cache = [], {}
    for r, i in enumerate(efg):
        c = ce[i]
        if c not in cache:
            S = X1r[cl == c]
            mu = S.mean(0)
            V = float(((S - mu) ** 2).sum(1).mean())
            cache[c] = (sq1 - 2 * X1r @ mu + mu @ mu, V, 2 * len(S) / (len(S) - 1) * V if len(S) >= 2 else np.nan)
        dm, V, dsis = cache[c]
        w = np.asarray(W[r], dtype=np.float64)
        out.append((float(w @ dm / w.sum() + V), dsis))
    return out


def run_job(obs, transition, rep):
    from sklearn.neighbors import NearestNeighbors
    tag = f"{CELLSET}_{transition.replace('->', '-')}_r{rep}"
    ck = CACHE / f"{tag}_cells.csv.gz"
    if ck.exists():
        return
    t0 = time.time()
    d, s0, s1 = rq1.job_cells(obs, CELLSET, transition, rep, cap=True)
    z = np.load(RQ1_CACHE / f"{tag}_emb.npz", allow_pickle=True)
    assert (z["ids"] == d.cell_id.to_numpy()).all()
    Xr = z["X"].astype(np.float64)
    st = d.stage.to_numpy()
    i0, i1 = np.where(st == s0)[0], np.where(st == s1)[0]
    e, l = d.iloc[i0].reset_index(drop=True), d.iloc[i1].reset_index(drop=True)
    ce = e.clone.fillna("").astype(str).to_numpy()
    cl = np.where(l.fg.to_numpy(), l.clone.fillna("").astype(str).to_numpy(), "__bg__")
    efg = np.where(e.fg.to_numpy())[0]
    late_labels = np.load(RQ1_CACHE / f"{tag}_labels.npy", allow_pickle=True)
    labels = sorted(set(late_labels))
    X1r = Xr[i1]
    a0 = load_counts(d)
    rows, fate, solves = [], [], []

    def score(W, name, space):
        W = np.asarray(W, dtype=np.float64)
        f = pd.DataFrame(rq1.row_metrics(W, efg, ce, cl, name))
        Dd = partner_D(W, efg, ce, cl, X1r)
        f["D"], f["D_sis"] = [x[0] for x in Dd], [x[1] for x in Dd]
        f["space"] = space
        rows.append(f)
        fa = pd.DataFrame(rq1.fate_vectors(W, efg, ce, cl, late_labels, labels, name))
        fa["space"] = space
        fate.append(fa)

    def knn_W(X0, X1):
        nb = NearestNeighbors(n_neighbors=rq1.KNN_K).fit(X1).kneighbors(X0[efg], return_distance=False)
        W = np.zeros((len(efg), len(X1)))
        W[np.repeat(np.arange(len(efg)), rq1.KNN_K), nb.ravel()] = 1 / rq1.KNN_K
        return W

    # RQ1 space: independent (for R) and the kNN reproduction check
    score(np.full((len(efg), len(i1)), 1.0 / len(i1)), "independent", "rq1")
    Wk = knn_W(Xr[i0], X1r)
    mk = np.array([Wk[r][cl == ce[i]].sum() for r, i in enumerate(efg)])
    ref = pd.read_csv(rq1.OUT / "cells_forward.csv.gz")
    ref = ref[(ref.cellset == CELLSET) & (~ref.uncapped.astype(bool)) & (ref.transition == transition)
              & (ref.rep == rep) & (ref.coupling == "kNN")].set_index("cell_id").m
    assert np.allclose(pd.Series(mk, index=e.cell_id.to_numpy()[efg]).loc[ref.index], ref, atol=1e-9), tag
    del Wk
    for space, (nh, npc) in SPACES.items():
        X = embed(a0, nh, npc)
        X0, X1 = X[i0], X[i1]
        score(knn_W(X0, X1), "kNN", space)
        M = cost(X0, X1)
        a = np.full(len(i0), 1 / len(i0), dtype=np.float32)
        lu, lv, it, conv = solve_bal_path(M, a)
        solves.append(dict(space=space, coupling="OT_bal_eps0.05", niter=it, converged=conv))
        score(rq1.plan_from_potentials(M[efg], 0.05, lu[efg], lv), "OT_bal_eps0.05", space)
        del M
        print(tag, space, f"{time.time() - t0:.0f}s", flush=True)
    z_cc, n_cc = cc_score(a0, st == s0)
    del a0
    M = cost(Xr[i0].astype(np.float32), X1r.astype(np.float32))
    for lam in LAMBDAS:
        a = np.exp(lam * z_cc)
        a = (a / a.sum()).astype(np.float32)
        lu, lv, it, conv = solve_bal_path(M, a)
        solves.append(dict(space="rq1", coupling=f"OT_growth_bal_lam{lam}", niter=it, converged=conv))
        score(rq1.plan_from_potentials(M[efg], 0.05, lu[efg], lv), f"OT_growth_bal_lam{lam}", "rq1")
        lu, lv, it, conv = solve_unb(M, 0.05, a)
        solves.append(dict(space="rq1", coupling=f"OT_growth_unb_lam{lam}", niter=it, converged=conv))
        score(rq1.plan_from_potentials(M[efg], 0.05, lu[efg], lv), f"OT_growth_unb_lam{lam}", "rq1")
        print(tag, "growth", lam, f"{time.time() - t0:.0f}s", flush=True)
    del M
    f = pd.concat(rows, ignore_index=True)
    f["clone"] = ce[f.cell]
    f["early_cluster"] = e.cluster.to_numpy()[f.cell]
    f["cell_id"] = e.cell_id.to_numpy()[f.cell]
    fa = pd.concat(fate, ignore_index=True)
    fa["clone"] = ce[fa.cell]
    fa["early_cluster"] = e.cluster.to_numpy()[fa.cell]
    fa["labels"] = json.dumps(labels)
    for t in (f, fa):
        t["transition"], t["rep"] = transition, rep
    CACHE.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(solves).assign(transition=transition, rep=rep, n_cc_genes=n_cc).to_csv(CACHE / f"{tag}_solves.csv", index=False)
    fa.drop(columns="cell").to_csv(CACHE / f"{tag}_fate.csv.gz", index=False)
    f.drop(columns="cell").to_csv(ck, index=False)
    print(tag, "done", f"{time.time() - t0:.0f}s", flush=True)


def summarise():
    OUT.mkdir(parents=True, exist_ok=True)
    f = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_cells.csv.gz"))])
    f["clone"] = f.rep.astype(str) + "|" + f.clone.astype(str)
    f["enrichment"] = f.m / f.m0
    f.to_csv(OUT / "cells_forward.csv.gz", index=False)
    cols = ["m", "m0", "auroc", "top1", "top5", "top10", "enrichment", "D"]
    ind = f[f.coupling == "independent"][["clone", "cell_id", "D"]].rename(columns={"D": "D_ind"})
    rows = []
    for (tr, sp, cp), g in f.groupby(["transition", "space", "coupling"]):
        g = g.merge(ind, on=["clone", "cell_id"])
        for rep, gg in [("pooled", g)] + list(g.groupby("rep")):
            c = gg.groupby("clone")[cols + ["D_sis", "D_ind"]].mean()
            rng = np.random.default_rng(0)
            idx = rng.integers(0, len(c), size=(N_BOOT, len(c)))
            arr = c[cols].to_numpy()
            boots = np.stack([arr[i].mean(0) for i in idx])
            r = dict(transition=tr, space=sp, coupling=cp, rep=rep, n_clones=len(c))
            for k, col in enumerate(cols):
                r[col], r[col + "_lo"], r[col + "_hi"] = c[col].mean(), *np.percentile(boots[:, k], [2.5, 97.5])
            c2 = c[c.D_sis.notna()]
            R = lambda ii: (c2.D.to_numpy()[ii].mean() - c2.D_sis.to_numpy()[ii].mean()) / \
                (c2.D_ind.to_numpy()[ii].mean() - c2.D_sis.to_numpy()[ii].mean())
            idx2 = rng.integers(0, len(c2), size=(N_BOOT, len(c2)))
            rb = np.array([R(i) for i in idx2])
            r["R"], r["R_lo"], r["R_hi"], r["n_clones_floor"] = R(np.arange(len(c2))), *np.percentile(rb, [2.5, 97.5]), len(c2)
            rows.append(r)
    pd.DataFrame(rows).to_csv(OUT / "summary.csv", index=False)
    fa = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_fate.csv.gz"))])
    fa["clone"] = fa.rep.astype(str) + "|" + fa.clone.astype(str)
    frows = []
    for tr, g in fa.groupby("transition"):
        labs = sorted({x for s in g.labels for x in json.loads(s)})
        vec = lambda row, col: [dict(zip(json.loads(row.labels), json.loads(getattr(row, col).replace("nan", "NaN")))).get(x, 0.0) for x in labs]
        g = g.copy()
        g["P"] = [vec(r, "pred") for r in g.itertuples()]
        g["Q"] = [vec(r, "clonal") for r in g.itertuples()]
        nclones = g[(g.coupling == "independent")].groupby("early_cluster").clone.nunique()
        for clu in nclones[nclones >= rq1.FATE_MIN_CLONES].index:
            for (sp, cp), gc in g[g.early_cluster == clu].groupby(["space", "coupling"]):
                P, Q = np.array(gc.P.tolist()), np.array(gc.Q.tolist())
                ids = gc.clone.to_numpy()
                u = np.unique(ids)
                Pc = np.stack([P[ids == c].mean(0) for c in u])
                Qc = np.stack([Q[ids == c].mean(0) for c in u])
                rng = np.random.default_rng(0)
                idx = rng.integers(0, len(u), size=(N_BOOT, len(u)))
                bt = [0.5 * np.abs(Pc[i].mean(0) - Qc[i].mean(0)).sum() for i in idx]
                frows.append(dict(transition=tr, early_cluster=clu, space=sp, coupling=cp, n_clones=len(u),
                                  tvd=0.5 * np.abs(Pc.mean(0) - Qc.mean(0)).sum(),
                                  tvd_lo=np.percentile(bt, 2.5), tvd_hi=np.percentile(bt, 97.5)))
    pd.DataFrame(frows).to_csv(OUT / "fate_agreement.csv", index=False)
    pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_solves.csv"))]).to_csv(OUT / "solves.csv", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    import ot
    import scanpy
    (OUT / "provenance.txt").write_text("\n".join([
        f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
        "plan: ANALYSIS_PLAN.md amendment A4.5 (2026-11-18, second review round, exploratory)",
        f"python {sys.version.split()[0]} scanpy {scanpy.__version__} POT {ot.__version__} numpy {np.__version__}",
        f"input {rq1.OBS.relative_to(ROOT)} sha256 {rq1.sha256(rq1.OBS)}",
        f"input {rq1.QC.relative_to(ROOT)} sha256 {rq1.sha256(rq1.QC)}",
        f"input {CC_GENES.relative_to(ROOT)} sha256 {rq1.sha256(CC_GENES)} (from the scanpy tutorial zip "
        "https://www.dropbox.com/s/3dby3bjsaf5arrw/cell_cycle_vignette_files.zip, sha256 "
        "6557fe2a2761d1550d41edf61b26c6d5ec9b42b4933d56c3161164a595f145ab)",
        "check: kNN in the RQ1 space reproduces RQ1 kNN own-clone mass per cell"]) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "summarise"])
    ap.add_argument("--transition", nargs="*", default=list(rq1.TRANSITIONS))
    ap.add_argument("--rep", type=int, nargs="*", default=[1, 2, 3])
    a = ap.parse_args()
    if a.cmd == "run":
        obs = rq1.load_obs()
        for t in a.transition:
            for r in a.rep:
                run_job(obs, t, r)
    else:
        summarise()


if __name__ == "__main__":
    main()
