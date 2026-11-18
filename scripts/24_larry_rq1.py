"""A4.6 — RQ1 on a second lineage dataset: LARRY in vitro haematopoiesis (Weinreb et al. 2020;
ANALYSIS_PLAN.md amendment A4.6). Data: GEO GSM4185642 (scripts/00b_download_larry.py).

Clones and pools (implementation choices the amendment leaves open, fixed here before any run):
  * a clone = a column of the clone matrix (every cell has at most one clone in these data);
    transitions day 2 -> 4 and day 4 -> 6; cross-stage clones have cells at both time points;
  * pool seed r in {0, 1, 2}: cross-stage clones are taken in the order of
    np.random.default_rng(r).permutation; a clone is added if, after adding it, neither time
    point exceeds 6,000 foreground cells (otherwise it is skipped and the next one tried);
  * background = cells of the two time points whose clone is not cross-stage for this
    transition (unbarcoded cells included); background is sampled (same rng, after the clone
    pass) so that each time point holds 8,000 cells in total (foreground + background), or all
    available cells if fewer;
  * the three pools are clone subsamples of one experiment, not biological replicates.
Expression: the deposited counts are total-count normalised; log1p of them; 2,000 HVGs with
scanpy flavor "seurat" (batch = Library); scale (max 10); 30 PCs (arpack, random_state 0);
fitted per pool on all its cells. Cost = squared Euclidean / median.
Couplings and metrics as §1 (balanced and unbalanced tau = 0.9 over eps {0.1, 0.05, 0.01, 0.005},
solved as RQ1; kNN k = 10; independent) plus the A4.4 baselines (kNN k = 5, 20, 50; Gaussian
kernel eps 0.01, 0.05, 0.1; MNN k = 10) and metrics (AP, tie-aware top-k, effective partners,
weight on 10 nearest) and the A4.1(a) partner error, all computed exactly as in
scripts/19_rq1_baselines.py and scripts/20_partner_error.py. Fate agreement: TVD per early
cluster (authors' 'Cell type annotation' at the early time point, >= 20 cross-stage clones)
with the authors' annotation for every later cell (no label transfer).

Usage: python scripts/24_larry_rq1.py build | run --transition T --pool R | summarise
"""
import argparse
import gzip
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


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


rq1 = _load("rq1", "06_rq1_coupling.py")
b19 = _load("b19", "19_rq1_baselines.py")
RAW = ROOT / "data/raw/larry"
PFX = "GSM4185642_stateFate_inVitro_"
CACHE = ROOT / "data/interim/larry"
OUT = ROOT / "results/24_larry"
rq1.CACHE = CACHE / "pot"
TRANSITIONS = {"d2->d4": (2.0, 4.0), "d4->d6": (4.0, 6.0)}
POOLS = (0, 1, 2)
FG_CAP, TOTAL = 6000, 8000
N_BOOT = 1000


def meta():
    import scipy.io
    m = pd.read_csv(RAW / f"{PFX}metadata.txt.gz", sep="\t")
    C = scipy.io.mmread(gzip.open(RAW / f"{PFX}clone_matrix.mtx.gz")).tocsr()
    k = np.asarray((C > 0).sum(1)).ravel()
    assert k.max() <= 1
    clone = np.full(C.shape[0], -1)
    r, c = C.nonzero()
    clone[r] = c
    m["clone"] = clone
    m["row"] = np.arange(len(m))
    return m


def pools(m):
    out = {}
    for tr, (t0, t1) in TRANSITIONS.items():
        x = m[m["Time point"].isin([t0, t1])]
        bc = x[x.clone >= 0].groupby("clone")["Time point"].agg(set)
        cross = np.array(sorted(bc[bc.apply(lambda s: {t0, t1} <= s)].index))
        cnt = x[x.clone.isin(cross)].groupby(["clone", "Time point"]).size().unstack(fill_value=0)
        for r in POOLS:
            rng = np.random.default_rng(r)
            n0 = n1 = 0
            take = []
            for c in cross[rng.permutation(len(cross))]:
                a0, a1 = cnt.loc[c, t0], cnt.loc[c, t1]
                if n0 + a0 <= FG_CAP and n1 + a1 <= FG_CAP:
                    take.append(c)
                    n0, n1 = n0 + a0, n1 + a1
            fg = x[x.clone.isin(take)]
            bgpool = x[~x.clone.isin(cross)]
            rows = [fg.row.to_numpy()]
            for t in (t0, t1):
                b = bgpool[bgpool["Time point"] == t].row.to_numpy()
                need = TOTAL - int((fg["Time point"] == t).sum())
                rows.append(np.sort(rng.choice(b, min(need, len(b)), replace=False)) if need > 0 else b[:0])
            out[(tr, r)] = dict(rows=np.sort(np.concatenate(rows)), fg_clones=np.array(take))
    return out


def build():
    """Stream the normalised counts once, keeping only cells used by some pool."""
    CACHE.mkdir(parents=True, exist_ok=True)
    m = meta()
    P = pools(m)
    need = np.unique(np.concatenate([p["rows"] for p in P.values()]))
    keep = np.zeros(len(m), bool)
    keep[need] = True
    local = np.full(len(m), -1, dtype=np.int64)
    local[need] = np.arange(len(need))          # row index within the needed cells
    rr, cc, vv = [], [], []
    with gzip.open(RAW / f"{PFX}normed_counts.mtx.gz", "rt") as f:
        f.readline(); f.readline()
        n_cells, n_genes, nnz = map(int, f.readline().split())
        for ch in pd.read_csv(f, sep=" ", header=None, names=["r", "c", "v"], chunksize=5_000_000,
                              dtype={"r": np.int32, "c": np.int32, "v": np.float32}):
            k = keep[ch.r.to_numpy() - 1]
            rr.append(local[ch.r.to_numpy()[k] - 1].astype(np.int32)); cc.append(ch.c.to_numpy()[k] - 1)
            vv.append(ch.v.to_numpy()[k])
            print("chunk", sum(len(x) for x in rr), flush=True)
    import scipy.sparse as sp
    # entries are sorted by gene (column): build CSC directly, then CSR once
    cc = np.concatenate(cc)
    indptr = np.concatenate([[0], np.cumsum(np.bincount(cc, minlength=n_genes))])
    del cc
    X = sp.csc_matrix((np.concatenate(vv), np.concatenate(rr), indptr), shape=(len(need), n_genes))
    del rr, vv
    X = X.tocsr()
    for (tr, r), p in P.items():
        tt = tr.replace("->", "-")
        sp.save_npz(CACHE / f"counts_{tt}_r{r}.npz", X[np.searchsorted(need, p["rows"])])
        np.savez(CACHE / f"pool_{tt}_r{r}.npz", rows=p["rows"], fg_clones=p["fg_clones"])
    m.to_csv(CACHE / "meta.csv.gz", index=False)
    print("built", len(need), "cells; nnz", nnz, flush=True)


def embed(Xc, libs):
    import anndata as ad
    import scanpy as sc
    a = ad.AnnData(Xc.astype(np.float32))
    a.obs["Library"] = libs
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.log1p(a)
        sc.pp.highly_variable_genes(a, flavor="seurat", n_top_genes=2000, batch_key="Library")
        a = a[:, a.var.highly_variable].copy()
        sc.pp.scale(a, max_value=10)
        sc.pp.pca(a, n_comps=30, random_state=0, svd_solver="arpack")
    return a.obsm["X_pca"].astype(np.float32)


def run_job(transition, r):
    import scipy.sparse as sp
    from scipy.spatial.distance import cdist
    from scipy.special import logsumexp
    from sklearn.neighbors import NearestNeighbors
    tt = transition.replace("->", "-")
    tag = f"{tt}_r{r}"
    ck = CACHE / f"{tag}_cells.csv.gz"
    if ck.exists():
        return
    t_start = time.time()
    rq1.CACHE.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(CACHE / "meta.csv.gz")
    Xc = sp.load_npz(CACHE / f"counts_{tt}_r{r}.npz")
    p = np.load(CACHE / f"pool_{tt}_r{r}.npz")
    rows, fgc = p["rows"], set(p["fg_clones"].tolist())
    d = m.iloc[rows].reset_index(drop=True)
    d["fg"] = d.clone.isin(fgc)
    t0, t1 = TRANSITIONS[transition]
    emb = CACHE / f"{tag}_emb.npy"
    if emb.exists():
        X = np.load(emb)
    else:
        X = embed(Xc, d.Library.to_numpy())
        np.save(emb, X)
    del Xc
    st = d["Time point"].to_numpy()
    i0, i1 = np.where(st == t0)[0], np.where(st == t1)[0]
    X0, X1 = X[i0], X[i1]
    e, l = d.iloc[i0].reset_index(drop=True), d.iloc[i1].reset_index(drop=True)
    ce = e.clone.astype(str).to_numpy()
    cl = np.where(l.fg.to_numpy(), l.clone.astype(str).to_numpy(), "__bg__")
    efg = np.where(e.fg.to_numpy())[0]
    late_labels = l["Cell type annotation"].astype(str).to_numpy()
    labels = sorted(set(late_labels))
    M = cdist(X0, X1, "sqeuclidean").astype(np.float32)
    M /= np.float32(np.median(M))
    Mf = M[efg].astype(np.float64)
    nn_l = NearestNeighbors(n_neighbors=50).fit(X1).kneighbors(X0[efg], return_distance=False)
    nn10 = nn_l[:, :10]
    X1d = X1.astype(np.float64)
    sq1 = (X1d ** 2).sum(1)
    sis = {}
    for c in np.unique(ce[efg]):
        S = X1d[cl == c]
        mu = S.mean(0)
        V = float(((S - mu) ** 2).sum(1).mean())
        sis[c] = (mu, V, 2 * len(S) / (len(S) - 1) * V if len(S) >= 2 else np.nan)
    recs, fate, solves = [], [], []

    def score(W, name):
        W = np.asarray(W, dtype=np.float64)
        f = pd.DataFrame(b19.row_stats(W, efg, ce, cl, nn10, name))
        Dl, Ds = [], []
        for rr_, i in enumerate(efg):
            mu, V, dsis = sis[ce[i]]
            w = W[rr_]
            tot = w.sum()
            Dl.append(float(w @ (sq1 - 2 * X1d @ mu + mu @ mu) / tot + V) if tot > 0 else np.nan)
            Ds.append(dsis)
        f["D"], f["D_sis"] = Dl, Ds
        recs.append(f)
        ok = W.sum(1) > 0
        if ok.all():
            fate.extend(rq1.fate_vectors(W, efg, ce, cl, late_labels, labels, name))

    for k in (5, 10, 20, 50):
        W = np.zeros((len(efg), len(X1)))
        W[np.repeat(np.arange(len(efg)), k), nn_l[:, :k].ravel()] = 1.0 / k
        score(W, f"kNN{k}")
    nn_e = NearestNeighbors(n_neighbors=10).fit(X0).kneighbors(X1, return_distance=False)
    W = np.zeros((len(efg), len(X1)))
    for rr_, i in enumerate(efg):
        js = [j for j in nn_l[rr_, :10] if i in nn_e[j]]
        if js:
            W[rr_, js] = 1.0 / len(js)
    score(W, "MNN10")
    for eps in (0.01, 0.05, 0.1):
        lw = -Mf / eps
        score(np.exp(lw - logsumexp(lw, axis=1, keepdims=True)), f"Gauss_eps{eps}")
    score(np.full((len(efg), len(X1)), 1.0 / len(X1)), "independent")
    del Mf
    warm, prev = None, None
    for eps in rq1.EPS:
        if prev is not None:
            warm = (warm[0] * prev / eps, warm[1] * prev / eps)
        s = rq1.solve_balanced(M, eps, warm, tag, float("inf"))
        warm, prev = (s["log_u"], s["log_v"]), eps
        solves.append(dict(coupling=f"OT_bal_eps{eps}", niter=int(s["niter"]), converged=bool(s["converged"])))
        score(rq1.plan_from_potentials(M[efg], eps, s["log_u"][efg], s["log_v"]), f"OT_bal_eps{eps}")
    for eps in sorted(rq1.EPS, reverse=True):
        P_, info = rq1.solve_unbalanced(M, eps, tag, float("inf"))
        solves.append(dict(coupling=f"OT_unb_eps{eps}", **info))
        score(P_[efg], f"OT_unb_eps{eps}")
        del P_
    del M
    f = pd.concat(recs, ignore_index=True)
    f["clone"] = ce[f.cell]
    f["early_cluster"] = e["Cell type annotation"].astype(str).to_numpy()[f.cell]
    f["transition"], f["pool"] = transition, r
    fa = pd.DataFrame(fate)
    fa["clone"] = ce[fa.cell]
    fa["early_cluster"] = e["Cell type annotation"].astype(str).to_numpy()[fa.cell]
    fa["labels"] = json.dumps(labels)
    fa["transition"], fa["pool"] = transition, r
    fa.drop(columns="cell").to_csv(CACHE / f"{tag}_fate.csv.gz", index=False)
    pd.DataFrame(solves).assign(transition=transition, pool=r, n_early=len(i0), n_later=len(i1),
                                n_early_fg=len(efg), n_later_fg=int(l.fg.sum()),
                                n_clones=len(set(ce[efg]))).to_csv(CACHE / f"{tag}_solves.csv", index=False)
    f.drop(columns="cell").to_csv(ck, index=False)
    print(tag, "done", f"{time.time() - t_start:.0f}s", flush=True)


def summarise():
    OUT.mkdir(parents=True, exist_ok=True)
    f = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_cells.csv.gz"))])
    f["clone"] = f.pool.astype(str) + "|" + f.clone.astype(str)
    f["enrichment"] = f.m / f.m0
    f.to_csv(OUT / "cells_forward.csv.gz", index=False)
    cols = b19.METRICS + ["enrichment", "D"]
    rows = []
    for (tr, cp), g in f.groupby(["transition", "coupling"]):
        gi = f[(f.transition == tr) & (f.coupling == "independent")]
        for pool, gg in [("all_pools", g)] + list(g.groupby("pool")):
            v = gg[gg.valid.astype(bool)]
            c = v.groupby("clone")[cols + ["D_sis"]].mean()
            ci = (gi if pool == "all_pools" else gi[gi.pool == pool]).groupby("clone").D.mean()
            rng = np.random.default_rng(0)
            idx = rng.integers(0, len(c), size=(N_BOOT, len(c)))
            arr = c[cols].to_numpy()
            boots = np.stack([np.nanmean(arr[i], 0) for i in idx])
            r = dict(transition=tr, coupling=cp, pool=pool, n_clones=len(c), coverage_cells=len(v) / len(gg))
            for k, col in enumerate(cols):
                r[col], r[col + "_lo"], r[col + "_hi"] = c[col].mean(), *np.percentile(boots[:, k], [2.5, 97.5])
            c2 = c[c.D_sis.notna()].join(ci.rename("D_ind"), how="inner")
            R = lambda ii: (c2.D.to_numpy()[ii].mean() - c2.D_sis.to_numpy()[ii].mean()) / \
                (c2.D_ind.to_numpy()[ii].mean() - c2.D_sis.to_numpy()[ii].mean())
            idx2 = rng.integers(0, len(c2), size=(N_BOOT, len(c2)))
            rb = np.array([R(i) for i in idx2])
            r["R"], r["R_lo"], r["R_hi"], r["n_clones_floor"] = R(np.arange(len(c2))), *np.percentile(rb, [2.5, 97.5]), len(c2)
            rows.append(r)
    pd.DataFrame(rows).to_csv(OUT / "summary.csv", index=False)
    fa = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_fate.csv.gz"))])
    fa["clone"] = fa.pool.astype(str) + "|" + fa.clone.astype(str)
    frows = []
    for tr, g in fa.groupby("transition"):
        labs = sorted({x for s in g.labels for x in json.loads(s)})
        vec = lambda row, col: [dict(zip(json.loads(row.labels), json.loads(getattr(row, col).replace("nan", "NaN")))).get(x, 0.0) for x in labs]
        g = g.copy()
        g["P"] = [vec(r, "pred") for r in g.itertuples()]
        g["Q"] = [vec(r, "clonal") for r in g.itertuples()]
        nclones = g[g.coupling == "independent"].groupby("early_cluster").clone.nunique()
        for clu in nclones[nclones >= rq1.FATE_MIN_CLONES].index:
            for cp, gc in g[g.early_cluster == clu].groupby("coupling"):
                P, Q = np.array(gc.P.tolist()), np.array(gc.Q.tolist())
                ids = gc.clone.to_numpy()
                u = np.unique(ids)
                Pc = np.stack([P[ids == c].mean(0) for c in u])
                Qc = np.stack([Q[ids == c].mean(0) for c in u])
                rng = np.random.default_rng(0)
                idx = rng.integers(0, len(u), size=(N_BOOT, len(u)))
                bt = [0.5 * np.abs(Pc[i].mean(0) - Qc[i].mean(0)).sum() for i in idx]
                frows.append(dict(transition=tr, early_cluster=clu, coupling=cp, n_clones=len(u),
                                  tvd=0.5 * np.abs(Pc.mean(0) - Qc.mean(0)).sum(),
                                  tvd_lo=np.percentile(bt, 2.5), tvd_hi=np.percentile(bt, 97.5),
                                  clonal_dist=json.dumps(dict(zip(labs, np.round(Qc.mean(0), 4))))))
    pd.DataFrame(frows).to_csv(OUT / "fate_agreement.csv", index=False)
    pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_solves.csv"))]).to_csv(OUT / "solves.csv", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    import ot
    import scanpy
    sha = pd.read_csv(OUT / "download_sha256.csv")
    (OUT / "provenance.txt").write_text("\n".join(
        [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
         "plan: ANALYSIS_PLAN.md amendment A4.6 (2026-11-18, second review round, exploratory)",
         f"python {sys.version.split()[0]} scanpy {scanpy.__version__} POT {ot.__version__} numpy {np.__version__}"]
        + [f"input data/raw/larry/{r.file} sha256 {r.sha256}" for r in sha.itertuples()]) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build", "run", "summarise"])
    ap.add_argument("--transition", nargs="*", default=list(TRANSITIONS))
    ap.add_argument("--pool", type=int, nargs="*", default=list(POOLS))
    a = ap.parse_args()
    if a.cmd == "build":
        build()
    elif a.cmd == "run":
        for t in a.transition:
            for r in a.pool:
                run_job(t, r)
    else:
        summarise()


if __name__ == "__main__":
    main()
