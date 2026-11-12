"""RQ1 — coupling recovery (ANALYSIS_PLAN.md §0–1, approved at 65cd70d).

Implements the plan as written. Implementation details the plan leaves open, fixed here
before any run:
  * Sinkhorn stopping: POT `sinkhorn_log` (balanced) with stopThr = 1e-6 on the marginal
    error, run in chunks of CHUNK iterations with warm start (identical to one long run),
    max MAX_ITER iterations; the ε grid is solved from large to small ε, each warm-started
    from the previous ε's potentials rescaled by ε_old/ε_new (changes speed, not the
    solution). Iterations and convergence are recorded per solve.
  * Unbalanced: KL-relaxed entropic OT (same objective as POT `sinkhorn_unbalanced`,
    reg_type 'kl'), reg_m = ε·τ/(1−τ), τ = 0.9, solved by log-domain unbalanced Sinkhorn:
    log_u ← φ(log a − LSE_j(−M/ε + log_v)), log_v ← φ(log b − LSE_i(−M/ε + log_u)),
    φ = reg_m/(reg_m+ε) = τ. Stop when the max relative change |Δ log_u|/(1+|log_u|)
    (and for log_v) < 1e-6; max MAX_ITER. (Revised from an absolute 1e-6, which is below
    float32 resolution for these potentials; solves that met the absolute criterion also
    meet this one, so earlier converged results are unchanged.)
    (Revision before any unbalanced result existed: POT's `sinkhorn_stabilized` did not
    finish ε = 0.1 on the smallest job in 20 min; the objective is unchanged.)
  * HVGs: seurat_v3 on raw counts needs ≥2 batches with variance; batch_key = library.
  * Fate labels of unlabelled later cells: kNN (k = 15) majority vote from clone-table
    cells of the same stage and replicate; ties broken by the smallest summed distance.
  * Metrics are clone-weighted: cell metrics → per-clone mean → mean over clones.

Jobs = cell set (primary, sensitivity) × transition × replicate, plus the pre-specified
uncapped I→II replicate-2 robustness job (balanced ε = 0.05 and kNN only). Every
expensive step is cached in data/interim/rq1/, so the script can be re-run to resume.

Usage:
  python scripts/06_rq1_coupling.py run --jobs primary:III->IV:1 [...] [--budget-min 9]
  python scripts/06_rq1_coupling.py summarise
"""
import argparse
import hashlib
import json
import subprocess
import sys
import time
import warnings
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LIBS = ROOT / "data/processed/libs"
OBS = ROOT / "results/03_objects/cell_obs.csv.gz"
QC = ROOT / "results/04_qc/cell_qc.csv.gz"
CACHE = ROOT / "data/interim/rq1"
OUT = ROOT / "results/06_rq1"

TRANSITIONS = {"I->II": ("I", "II", 1), "II->III": ("II", "III", 2), "III->IV": ("III", "IV", 3)}
EPS = [0.1, 0.05, 0.01, 0.005]          # solved large -> small (warm start)
TAU = 0.9
CAP = 8000
N_HVG, N_PC, KNN_K, LABEL_K = 2000, 30, 10, 15
TOPK = (1, 5, 10)
STOP, CHUNK, MAX_ITER = 1e-6, 100, 20000
N_BOOT, FATE_MIN_CLONES = 1000, 20
SEED = 0


# ----------------------------------------------------------------------------- data
def load_obs():
    obs = pd.read_csv(OBS, low_memory=False)
    qc = pd.read_csv(QC)
    obs = obs.merge(qc[["cell_id", "keep", "doublet_flag"]], on="cell_id")
    obs["sensitivity"] = obs.keep & ~obs.doublet_flag
    obs["primary"] = obs.keep
    return obs


def job_cells(obs, cellset, transition, rep, cap=True):
    s0, s1, e = TRANSITIONS[transition]
    d = obs[obs[cellset] & (obs.exp == e) & (obs.rep == rep) & obs.stage.isin([s0, s1])].copy()
    st = d[d.in_clone_table].groupby("clone").stage.agg(set)
    cross = set(st[st.apply(lambda x: {s0, s1} <= x)].index)
    d["fg"] = d.clone.isin(cross)
    keep = []
    rng = np.random.default_rng(SEED)
    for s in (s0, s1):
        x = d[d.stage == s]
        bg = x[~x.fg].cell_id.to_numpy()
        if cap and len(bg) > CAP:
            bg = rng.choice(bg, CAP, replace=False)
        keep.append(x[x.fg].cell_id.to_numpy())
        keep.append(np.sort(bg))
    d = d[d.cell_id.isin(np.concatenate(keep))]
    return d.reset_index(drop=True), s0, s1


def embed(d):
    """Normalise, HVG (seurat_v3, batch=library), scale, 30-PC PCA on both stages pooled."""
    import anndata as ad
    import scanpy as sc
    parts = []
    for lib, ids in d.groupby("library").cell_id:
        a = ad.read_h5ad(LIBS / f"{lib}.h5ad")
        parts.append(a[ids.to_numpy()].copy())
    a = ad.concat(parts)
    a = a[d.cell_id.to_numpy()].copy()
    a.obs["library"] = d.library.to_numpy()
    a.layers["counts"] = a.X.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=N_HVG,
                                    batch_key="library", layer="counts")
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        a = a[:, a.var.highly_variable].copy()
        sc.pp.scale(a, max_value=10)
        sc.pp.pca(a, n_comps=N_PC, random_state=SEED, svd_solver="arpack")
    return a.obsm["X_pca"].astype(np.float32), list(a.var_names)


def transfer_labels(X, d, stage):
    from sklearn.neighbors import NearestNeighbors
    idx = np.where(d.stage.to_numpy() == stage)[0]
    lab = d.cluster.to_numpy()[idx]
    known = ~pd.isna(lab)
    out = lab.astype(object).copy()
    if (~known).any():
        nn = NearestNeighbors(n_neighbors=LABEL_K).fit(X[idx][known])
        dist, nb = nn.kneighbors(X[idx][~known])
        ref = lab[known]
        res = []
        for dr, nr in zip(dist, nb):
            t = pd.DataFrame({"l": ref[nr], "d": dr}).groupby("l").agg(n=("d", "size"), s=("d", "sum"))
            res.append(t.sort_values(["n", "s"], ascending=[False, True]).index[0])
        out[~known] = res
    return out


# ----------------------------------------------------------------------------- OT
def solve_balanced(M, eps, warm, tag, budget_end):
    import ot
    n0, n1 = M.shape
    a = np.full(n0, 1 / n0, dtype=np.float32)
    b = np.full(n1, 1 / n1, dtype=np.float32)
    ck = CACHE / f"{tag}_bal_{eps}.npz"
    state = dict(np.load(ck)) if ck.exists() else None
    if state is not None and bool(state["done"]):
        return state
    if state is not None:
        lu, lv, it = state["log_u"], state["log_v"], int(state["niter"])
    elif warm is not None:
        lu, lv, it = warm[0], warm[1], 0
    else:
        lu, lv, it = np.zeros(n0, np.float32), np.zeros(n1, np.float32), 0
    err = np.inf
    while it < MAX_ITER:
        if time.time() > budget_end:
            np.savez(ck, log_u=lu, log_v=lv, niter=it, err=err, done=False)
            return None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, log = ot.sinkhorn(a, b, M, eps, method="sinkhorn_log", numItermax=CHUNK,
                                 stopThr=STOP, warmstart=(lu, lv), log=True, warn=False)
        lu, lv = log["log_u"].astype(np.float32), log["log_v"].astype(np.float32)
        it += int(log["niter"]) + 1
        err = float(log["err"][-1]) if len(log["err"]) else np.inf
        if err < STOP:
            break
    state = dict(log_u=lu, log_v=lv, niter=it, err=err, done=True, converged=err < STOP)
    np.savez(ck, **state)
    return state


def plan_from_potentials(M, eps, lu, lv):
    # P_ij = exp(log_u_i - M_ij/eps + log_v_j); returned float32
    P = np.empty_like(M)
    np.divide(M, -eps, out=P)
    P += lu[:, None]
    P += lv[None, :]
    np.exp(P, out=P)
    return P


def solve_unbalanced(M, eps, tag, budget_end):
    from scipy.special import logsumexp
    ck = CACHE / f"{tag}_unb_{eps}.npz"
    if ck.exists():
        z = np.load(ck)
        if bool(z["done"]):
            info = dict(niter=int(z["niter"]), err=float(z["err"]), converged=bool(z["converged"]))
            return plan_from_potentials(M, eps, z["log_u"], z["log_v"]), info
        lu, lv, it = z["log_u"], z["log_v"], int(z["niter"])
    else:
        lu, lv, it = np.zeros(M.shape[0], np.float32), np.zeros(M.shape[1], np.float32), 0
    n0, n1 = M.shape
    la, lb = np.float32(-np.log(n0)), np.float32(-np.log(n1))
    reg_m = eps * TAU / (1 - TAU)
    phi = np.float32(reg_m / (reg_m + eps))
    logK = (M / np.float32(-eps)).astype(np.float32)
    err = np.inf
    while it < MAX_ITER:
        if it % CHUNK == 0 and time.time() > budget_end:
            np.savez(ck, log_u=lu, log_v=lv, niter=it, err=err, done=False, converged=False)
            return None, None
        lu_new = phi * (la - logsumexp(logK + lv[None, :], axis=1)).astype(np.float32)
        lv_new = phi * (lb - logsumexp(logK + lu_new[:, None], axis=0)).astype(np.float32)
        # relative change: an absolute 1e-6 is below float32 resolution for |log u| >~ 8
        err = float(max((np.abs(lu_new - lu) / (1 + np.abs(lu))).max(),
                        (np.abs(lv_new - lv) / (1 + np.abs(lv))).max()))
        lu, lv, it = lu_new, lv_new, it + 1
        if err < STOP:
            break
    del logK
    info = dict(niter=it, err=err, converged=err < STOP)
    np.savez(ck, log_u=lu, log_v=lv, niter=it, err=err, done=True, converged=err < STOP)
    return plan_from_potentials(M, eps, lu, lv), info


# ----------------------------------------------------------------------------- metrics
def mw_auroc(scores, pos):
    """Mann-Whitney AUROC with ties as 1/2."""
    from scipy.stats import rankdata
    n1, n0 = pos.sum(), (~pos).sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(scores)
    return (r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def topk_chance(n, s, k):
    return 1 - comb(n - s, k) / comb(n, k)


def row_metrics(W, early_idx, clone_early, clone_late, coupling):
    """Forward metrics for early cells (rows of W, already restricted to fg early cells)."""
    rows = []
    n = W.shape[1]
    for r, i in enumerate(early_idx):
        w = np.asarray(W[r], dtype=np.float64).ravel()
        tot = w.sum()
        pos = clone_late == clone_early[i]
        s = int(pos.sum())
        rec = dict(cell=i, coupling=coupling, n_later=n, s=s, m0=s / n)
        if tot <= 0 or not np.isfinite(tot):
            rec.update(m=np.nan, auroc=np.nan, **{f"top{k}": np.nan for k in TOPK})
        else:
            w = w / tot
            rec["m"] = float(w[pos].sum())
            rec["auroc"] = float(mw_auroc(w, pos))
            order = np.argsort(-w, kind="stable")[:max(TOPK)]
            for k in TOPK:
                rec[f"top{k}"] = float(pos[order[:k]].any())
        rows.append(rec)
    return rows


def knn_metrics(X0, X1, early_idx, clone_early, clone_late):
    from scipy.stats import rankdata  # noqa: F401
    from sklearn.neighbors import NearestNeighbors
    nb = NearestNeighbors(n_neighbors=KNN_K).fit(X1).kneighbors(X0[early_idx], return_distance=False)
    rows = []
    n = X1.shape[0]
    for r, i in enumerate(early_idx):
        pos = clone_late == clone_early[i]
        s = int(pos.sum())
        q = int(pos[nb[r]].sum())
        w = np.zeros(n)
        w[nb[r]] = 1 / KNN_K
        rec = dict(cell=i, coupling="kNN", n_later=n, s=s, m0=s / n, m=q / KNN_K,
                   auroc=float(mw_auroc(w, pos)))
        for k in TOPK:   # hypergeometric expectation within 10 equally weighted neighbours
            rec[f"top{k}"] = 1 - comb(KNN_K - q, k) / comb(KNN_K, k)
        rows.append(rec)
    return rows, nb


def independent_metrics(early_idx, clone_early, clone_late):
    n = len(clone_late)
    rows = []
    for i in early_idx:
        s = int((clone_late == clone_early[i]).sum())
        rec = dict(cell=i, coupling="independent", n_later=n, s=s, m0=s / n, m=s / n, auroc=0.5)
        for k in TOPK:
            rec[f"top{k}"] = topk_chance(n, s, k)
        rows.append(rec)
    return rows


def backward_metrics(W, late_idx, clone_early_all, clone_late, coupling):
    """Column-normalised: later fg cell -> mass on own-clone early cells."""
    rows = []
    col = np.asarray(W[:, late_idx], dtype=np.float64)
    n0 = W.shape[0]
    for c, j in enumerate(late_idx):
        w = col[:, c]
        tot = w.sum()
        pos = clone_early_all == clone_late[j]
        s = int(pos.sum())
        rec = dict(cell=j, coupling=coupling, n_earlier=n0, s=s, m0=s / n0)
        if tot <= 0 or not np.isfinite(tot):
            rec.update(m=np.nan, auroc=np.nan)
        else:
            w = w / tot
            rec.update(m=float(w[pos].sum()), auroc=float(mw_auroc(w, pos)))
        rows.append(rec)
    return rows


def fate_vectors(W, early_idx, clone_early, clone_late, late_labels, labels, coupling):
    """Per early fg cell: predicted later-label distribution vs its clone's sister labels."""
    L = np.zeros((W.shape[1], len(labels)))
    L[np.arange(W.shape[1]), pd.Index(labels).get_indexer(late_labels)] = 1
    rows = []
    for r, i in enumerate(early_idx):
        w = np.asarray(W[r], dtype=np.float64).ravel()
        tot = w.sum()
        pred = (w / tot) @ L if tot > 0 and np.isfinite(tot) else np.full(len(labels), np.nan)
        sis = L[clone_late == clone_early[i]].mean(0)
        rows.append(dict(cell=i, coupling=coupling, pred=pred.tolist(), clonal=sis.tolist()))
    return rows


# ----------------------------------------------------------------------------- job
def run_job(obs, cellset, transition, rep, budget_end, uncapped=False):
    tag = f"{cellset}_{transition.replace('->', '-')}_r{rep}" + ("_uncapped" if uncapped else "")
    done_flag = CACHE / f"{tag}.done"
    if done_flag.exists():
        return True
    d, s0, s1 = job_cells(obs, cellset, transition, rep, cap=not uncapped)
    emb = CACHE / f"{tag}_emb.npz"
    if emb.exists():
        z = np.load(emb, allow_pickle=True)
        X, ids = z["X"], z["ids"]
        assert (ids == d.cell_id.to_numpy()).all()
    else:
        X, hvg = embed(d)
        np.savez(emb, X=X, ids=d.cell_id.to_numpy())
    st = d.stage.to_numpy()
    i0, i1 = np.where(st == s0)[0], np.where(st == s1)[0]
    X0, X1 = X[i0], X[i1]
    e = d.iloc[i0].reset_index(drop=True)
    l = d.iloc[i1].reset_index(drop=True)
    clone_early = e.clone.fillna("").astype(str).to_numpy()
    clone_late = l.clone.fillna("").astype(str).to_numpy()
    clone_late = np.where(l.fg.to_numpy(), clone_late, "__bg__")       # background never positive
    clone_early_all = np.where(e.fg.to_numpy(), clone_early, "__bg__")
    early_fg = np.where(e.fg.to_numpy())[0]
    late_fg = np.where(l.fg.to_numpy())[0]

    lab_ck = CACHE / f"{tag}_labels.npy"
    if lab_ck.exists():
        late_labels = np.load(lab_ck, allow_pickle=True)
    else:
        late_labels = transfer_labels(X, d, s1)
        np.save(lab_ck, late_labels)
    labels = sorted(set(late_labels))

    from scipy.spatial.distance import cdist
    M = cdist(X0, X1, "sqeuclidean").astype(np.float32)
    M /= np.float32(np.median(M))

    fwd, bwd, fate, solves = [], [], [], []
    rows, _ = knn_metrics(X0, X1, early_fg, clone_early, clone_late)
    fwd += rows
    fwd += independent_metrics(early_fg, clone_early, clone_late)
    # kNN as a matrix for backward/fate metrics
    from sklearn.neighbors import NearestNeighbors
    nb = NearestNeighbors(n_neighbors=KNN_K).fit(X1).kneighbors(X0, return_distance=False)
    Wk = np.zeros(M.shape, dtype=np.float32)
    Wk[np.repeat(np.arange(len(X0)), KNN_K), nb.ravel()] = 1 / KNN_K
    bwd += backward_metrics(Wk, late_fg, clone_early_all, clone_late, "kNN")
    fate += fate_vectors(Wk[early_fg], early_fg, clone_early, clone_late, late_labels, labels, "kNN")
    del Wk
    Wi = np.full((len(early_fg), len(X1)), 1 / len(X1), dtype=np.float32)
    fate += fate_vectors(Wi, early_fg, clone_early, clone_late, late_labels, labels, "independent")
    del Wi

    eps_list = [0.05] if uncapped else EPS
    warm = None
    prev_eps = None
    for eps in (EPS if not uncapped else [0.1, 0.05]):
        if prev_eps is not None and warm is not None:
            warm = (warm[0] * prev_eps / eps, warm[1] * prev_eps / eps)
        stt = solve_balanced(M, eps, warm, tag, budget_end)
        if stt is None:
            return False
        warm, prev_eps = (stt["log_u"], stt["log_v"]), eps
        if eps not in eps_list:
            continue
        name = f"OT_bal_eps{eps}"
        P = plan_from_potentials(M, eps, stt["log_u"], stt["log_v"])
        fwd += row_metrics(P[early_fg], early_fg, clone_early, clone_late, name)
        bwd += backward_metrics(P, late_fg, clone_early_all, clone_late, name)
        fate += fate_vectors(P[early_fg], early_fg, clone_early, clone_late, late_labels, labels, name)
        solves.append(dict(coupling=name, niter=int(stt["niter"]), err=float(stt["err"]),
                           converged=bool(stt["converged"])))
        del P
    if not uncapped:
        for eps in sorted(EPS, reverse=True):
            P, info = solve_unbalanced(M, eps, tag, budget_end)
            if P is None:
                return False
            name = f"OT_unb_eps{eps}"
            fwd += row_metrics(P[early_fg], early_fg, clone_early, clone_late, name)
            bwd += backward_metrics(P, late_fg, clone_early_all, clone_late, name)
            fate += fate_vectors(P[early_fg], early_fg, clone_early, clone_late, late_labels, labels, name)
            solves.append(dict(coupling=name, **info))
            del P

    meta = dict(cellset=cellset, transition=transition, rep=rep, uncapped=uncapped)
    f = pd.DataFrame(fwd)
    f["clone"] = clone_early[f.cell]
    f["early_cluster"] = e.cluster.to_numpy()[f.cell]
    f["cell_id"] = e.cell_id.to_numpy()[f.cell]
    b = pd.DataFrame(bwd)
    b["clone"] = clone_late[b.cell]
    b["cell_id"] = l.cell_id.to_numpy()[b.cell]
    fa = pd.DataFrame(fate)
    fa["clone"] = clone_early[fa.cell]
    fa["early_cluster"] = e.cluster.to_numpy()[fa.cell]
    fa["labels"] = json.dumps(labels)
    for t in (f, b, fa):
        for k, v in meta.items():
            t[k] = v
    f.drop(columns="cell").to_csv(CACHE / f"{tag}_forward.csv.gz", index=False)
    b.drop(columns="cell").to_csv(CACHE / f"{tag}_backward.csv.gz", index=False)
    fa.drop(columns="cell").to_csv(CACHE / f"{tag}_fate.csv.gz", index=False)
    s = pd.DataFrame(solves)
    for k, v in meta.items():
        s[k] = v
    s["n_early"], s["n_later"] = len(X0), len(X1)
    s["n_early_fg"], s["n_later_fg"] = len(early_fg), len(late_fg)
    s["n_clones"] = len(set(clone_early[early_fg]))
    s.to_csv(CACHE / f"{tag}_solves.csv", index=False)
    done_flag.write_text("ok")
    return True


# ----------------------------------------------------------------------------- summary
METRICS = ["m", "m0", "auroc"] + [f"top{k}" for k in TOPK]


def clone_weighted(df, cols, rng, n_boot):
    g = df.groupby("clone")[cols].mean()
    est = g.mean()
    if n_boot == 0:
        return est, None, None
    idx = rng.integers(0, len(g), size=(n_boot, len(g)))
    arr = g.to_numpy()
    boots = np.stack([arr[i].mean(0) for i in idx])
    return est, np.percentile(boots, 2.5, axis=0), np.percentile(boots, 97.5, axis=0)


def summarise():
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    fwd = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_forward.csv.gz"))])
    bwd = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_backward.csv.gz"))])
    fate = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_fate.csv.gz"))])
    solves = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_solves.csv"))])
    for t in (fwd, bwd, fate):
        t["clone"] = t.rep.astype(str) + "|" + t.clone.astype(str)   # clones are replicate-specific
    fwd["enrichment"] = fwd.m / fwd.m0
    fwd["diff"] = fwd.m - fwd.m0
    cols = METRICS + ["enrichment", "diff"]
    rows = []
    keys = ["cellset", "uncapped", "transition", "coupling"]
    for k, g in fwd.groupby(keys):
        for rep, gg in [("pooled", g)] + list(g.groupby("rep")):
            est, lo, hi = clone_weighted(gg, cols, rng, N_BOOT)
            r = dict(zip(keys, k), rep=rep, direction="forward", n_clones=gg.clone.nunique(),
                     n_cells=len(gg))
            for c in cols:
                r[c], r[c + "_lo"], r[c + "_hi"] = est[c], lo[cols.index(c)], hi[cols.index(c)]
            rows.append(r)
    bwd["enrichment"] = bwd.m / bwd.m0
    bwd["diff"] = bwd.m - bwd.m0
    bcols = ["m", "m0", "auroc", "enrichment", "diff"]
    for k, g in bwd.groupby(keys):
        for rep, gg in [("pooled", g)] + list(g.groupby("rep")):
            est, lo, hi = clone_weighted(gg.dropna(subset=["m"]), bcols, rng, N_BOOT)
            r = dict(zip(keys, k), rep=rep, direction="backward", n_clones=gg.clone.nunique(),
                     n_cells=len(gg), n_cells_zero_mass=int(gg.m.isna().sum()))
            for c in bcols:
                r[c], r[c + "_lo"], r[c + "_hi"] = est[c], lo[bcols.index(c)], hi[bcols.index(c)]
            rows.append(r)
    pd.DataFrame(rows).to_csv(OUT / "summary.csv", index=False)

    # fate-level agreement: TVD per early cluster with >= FATE_MIN_CLONES cross-stage clones
    frows = []
    for k, g in fate.groupby(["cellset", "uncapped", "transition"]):
        labs = sorted({x for s in g.labels for x in json.loads(s)})
        def vec(row, col):
            v = dict(zip(json.loads(row.labels), json.loads(row[col]) if isinstance(row[col], str) else row[col]))
            return [v.get(x, 0.0) for x in labs]
        g = g.copy()
        g["P"] = [vec(r, "pred") for r in g.itertuples()]
        g["Q"] = [vec(r, "clonal") for r in g.itertuples()]
        nclones = g[g.coupling == "kNN"].groupby("early_cluster").clone.nunique()
        for cl in nclones[nclones >= FATE_MIN_CLONES].index:
            for cp, gc in g[g.early_cluster == cl].groupby("coupling"):
                P = np.array(gc.P.tolist()); Q = np.array(gc.Q.tolist())
                cl_ids = gc.clone.to_numpy(); u = np.unique(cl_ids)
                Pc = np.stack([P[cl_ids == c].mean(0) for c in u]); Qc = np.stack([Q[cl_ids == c].mean(0) for c in u])
                tvd = 0.5 * np.abs(np.nanmean(Pc, 0) - Qc.mean(0)).sum()
                idx = rng.integers(0, len(u), size=(N_BOOT, len(u)))
                bt = [0.5 * np.abs(np.nanmean(Pc[i], 0) - Qc[i].mean(0)).sum() for i in idx]
                frows.append(dict(cellset=k[0], uncapped=k[1], transition=k[2], early_cluster=cl,
                                  coupling=cp, n_clones=len(u), tvd=tvd,
                                  tvd_lo=np.percentile(bt, 2.5), tvd_hi=np.percentile(bt, 97.5),
                                  clonal_dist=json.dumps(dict(zip(labs, np.round(Qc.mean(0), 4)))),
                                  predicted_dist=json.dumps(dict(zip(labs, np.round(np.nanmean(Pc, 0), 4))))))
    pd.DataFrame(frows).to_csv(OUT / "fate_agreement.csv", index=False)
    solves.to_csv(OUT / "solves.csv", index=False)
    fwd.to_csv(OUT / "cells_forward.csv.gz", index=False)
    bwd.to_csv(OUT / "cells_backward.csv.gz", index=False)
    write_provenance()


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_provenance():
    import ot
    import scanpy as sc
    import scipy
    import sklearn
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    lines = [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
             "plan: ANALYSIS_PLAN.md approved at 65cd70d",
             f"python {sys.version.split()[0]} scanpy {sc.__version__} POT {ot.__version__} "
             f"numpy {np.__version__} scipy {scipy.__version__} sklearn {sklearn.__version__} pandas {pd.__version__}",
             "environment: Claude cloud workspace (Linux, 2 CPU, 7 GB)",
             f"input {OBS.relative_to(ROOT)} sha256 {sha256(OBS)}",
             f"input {QC.relative_to(ROOT)} sha256 {sha256(QC)}"]
    (OUT / "provenance.txt").write_text("\n".join(lines) + "\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("cmd", choices=["run", "summarise"])
    p.add_argument("--jobs", nargs="*", default=[])
    p.add_argument("--budget-min", type=float, default=9.0)
    a = p.parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)
    if a.cmd == "summarise":
        summarise()
        return
    obs = load_obs()
    end = time.time() + a.budget_min * 60
    for j in a.jobs:
        cs, tr, rp = j.split(":")
        unc = rp.endswith("u")
        ok = run_job(obs, cs, tr, int(rp.rstrip("u")), end, uncapped=unc)
        print(j, "done" if ok else "paused", flush=True)
        if not ok:
            sys.exit(3)


if __name__ == "__main__":
    main()
