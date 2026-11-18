"""A4.4 — proximity baselines and ranking metrics for RQ1 (ANALYSIS_PLAN.md amendment A4.4).

Primary cell set, capped RQ1 pools, RQ1 30-PC embedding and cost (squared Euclidean / median),
per transition x replicate. Rows = foreground early cells (forward direction), weights
row-normalised over all later pool cells.

Couplings:
  * kNN, k in {5, 10, 20, 50}: uniform over the k nearest later cells;
  * Gaussian kernel, w_ij ∝ exp(-C_ij / eps), eps in {0.01, 0.05, 0.1} (no column constraint);
  * mutual nearest neighbours, k = 10: later cell j is a partner of early cell i if j is among
    i's 10 nearest later cells and i is among j's 10 nearest early cells (all early pool
    cells); uniform over partners; rows without a partner are counted (coverage) and excluded
    from that coupling's metrics;
  * the existing RQ1 couplings, rebuilt from the cached potentials: balanced and unbalanced
    OT for eps in {0.1, 0.05, 0.01, 0.005}; independent (uniform).
Metrics per row: own-clone mass m and chance m0; Mann-Whitney AUROC (ties 1/2); top-k recall,
k in {1, 5, 10}, as the expectation under uniform tie-breaking (no random draw): later cells
are grouped by equal weight in descending order and P(no sister in the first k) is the product
of hypergeometric terms over the groups the first k positions touch (this equals the §1 rule:
argsort order for untied OT weights, the hypergeometric expectation for kNN, and the analytic
chance for the independent coupling); average precision with the sisters as positives
(sklearn average_precision_score), chance s/n; effective number of partners exp(H(w_i.)), and
the weight on the 10 nearest later cells (A4.1(a) row statistics).
Summaries: clone-weighted (cell -> clone mean -> mean over clones), pooled and per replicate,
95% percentile bootstrap over clones (1,000; np.random.default_rng(0) per row).
Reproduction check: m, m0, AUROC and top-k of the RQ1 couplings equal the RQ1 per-cell
values (results/06_rq1/cells_forward.csv.gz) to 1e-6, or the script stops.

Usage: python scripts/19_rq1_baselines.py run | summarise
"""
import argparse
import importlib.util
import subprocess
import sys
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rq1", ROOT / "scripts/06_rq1_coupling.py")
rq1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rq1)

CELLSET = "primary"
CACHE = ROOT / "data/interim/rq1_baselines"
OUT = ROOT / "results/19_rq1_baselines"
KS = (5, 10, 20, 50)
GAUSS_EPS = (0.01, 0.05, 0.1)
MNN_K = 10
N_BOOT = 1000
METRICS = ["m", "m0", "auroc", "top1", "top5", "top10", "ap", "ap_chance", "n_eff", "w_nn10"]


def topk_expect(w, pos, ks=rq1.TOPK):
    """Expected top-k hit under uniform tie-breaking."""
    vals, inv = np.unique(-w, return_inverse=True)          # groups in descending weight
    gsize = np.bincount(inv)
    gpos = np.bincount(inv, weights=pos.astype(float)).astype(int)
    out = {}
    for k in ks:
        left, p_none = k, 1.0
        for g, p in zip(gsize, gpos):
            if left <= 0:
                break
            r = min(left, g)
            p_none *= comb(g - p, r) / comb(g, r)
            left -= r
            if p_none == 0:
                break
        out[f"top{k}"] = 1 - p_none
    return out


def row_stats(W, efg, ce, cl, nn10, coupling):
    from sklearn.metrics import average_precision_score
    rows = []
    n = W.shape[1]
    for r, i in enumerate(efg):
        w = np.asarray(W[r], dtype=np.float64)
        pos = cl == ce[i]
        s = int(pos.sum())
        rec = dict(cell=i, coupling=coupling, n_later=n, s=s, m0=s / n, ap_chance=s / n)
        tot = w.sum()
        if not np.isfinite(tot) or tot <= 0:
            rec.update(valid=False)
            rows.append(rec)
            continue
        w = w / tot
        nz = w[w > 0]
        rec.update(valid=True, n_pos_w=int((w > 0).sum()), m=float(w[pos].sum()), auroc=float(rq1.mw_auroc(w, pos)),
                   ap=float(average_precision_score(pos, w)),
                   n_eff=float(np.exp(-(nz * np.log(nz)).sum())), w_nn10=float(w[nn10[r]].sum()),
                   **topk_expect(w, pos))
        rows.append(rec)
    return rows


def run_job(obs, transition, rep):
    from scipy.spatial.distance import cdist
    from scipy.special import logsumexp
    from sklearn.neighbors import NearestNeighbors
    tag = f"{CELLSET}_{transition.replace('->', '-')}_r{rep}"
    ck = CACHE / f"{tag}_cells.csv.gz"
    if ck.exists():
        return
    d, s0, s1 = rq1.job_cells(obs, CELLSET, transition, rep, cap=True)
    z = np.load(rq1.CACHE / f"{tag}_emb.npz", allow_pickle=True)
    assert (z["ids"] == d.cell_id.to_numpy()).all()
    X, st = z["X"], d.stage.to_numpy()
    i0, i1 = np.where(st == s0)[0], np.where(st == s1)[0]
    X0, X1 = X[i0], X[i1]
    e, l = d.iloc[i0].reset_index(drop=True), d.iloc[i1].reset_index(drop=True)
    ce = e.clone.fillna("").astype(str).to_numpy()
    cl = np.where(l.fg.to_numpy(), l.clone.fillna("").astype(str).to_numpy(), "__bg__")
    efg = np.where(e.fg.to_numpy())[0]
    M = cdist(X0, X1, "sqeuclidean").astype(np.float32)
    M /= np.float32(np.median(M))
    Mf = M[efg].astype(np.float64)

    nn_l = NearestNeighbors(n_neighbors=max(KS)).fit(X1).kneighbors(X0[efg], return_distance=False)
    nn10 = nn_l[:, :10]
    rows = []
    for k in KS:
        W = np.zeros((len(efg), len(X1)))
        W[np.repeat(np.arange(len(efg)), k), nn_l[:, :k].ravel()] = 1.0 / k
        rows += row_stats(W, efg, ce, cl, nn10, f"kNN{k}")
    # mutual nearest neighbours (k = 10): later -> early neighbours among all early pool cells
    nn_e = NearestNeighbors(n_neighbors=MNN_K).fit(X0).kneighbors(X1, return_distance=False)
    W = np.zeros((len(efg), len(X1)))
    for r, i in enumerate(efg):
        js = [j for j in nn_l[r, :MNN_K] if i in nn_e[j]]
        if js:
            W[r, js] = 1.0 / len(js)
    rows += row_stats(W, efg, ce, cl, nn10, "MNN10")
    for eps in GAUSS_EPS:
        lw = -Mf / eps
        W = np.exp(lw - logsumexp(lw, axis=1, keepdims=True))
        rows += row_stats(W, efg, ce, cl, nn10, f"Gauss_eps{eps}")
    rows += row_stats(np.full((len(efg), len(X1)), 1.0 / len(X1)), efg, ce, cl, nn10, "independent")
    for kind in ("bal", "unb"):
        for eps in rq1.EPS:
            p = np.load(rq1.CACHE / f"{tag}_{kind}_{eps}.npz")
            W = rq1.plan_from_potentials(M[efg], eps, p["log_u"][efg], p["log_v"]).astype(np.float64)
            rows += row_stats(W, efg, ce, cl, nn10, f"OT_{kind}_eps{eps}")
    del M, Mf
    f = pd.DataFrame(rows)
    f["clone"] = ce[f.cell]
    f["cell_id"] = e.cell_id.to_numpy()[f.cell]
    f["transition"], f["rep"] = transition, rep

    # reproduction check against RQ1 per-cell values
    ref = pd.read_csv(rq1.OUT / "cells_forward.csv.gz")
    ref = ref[(ref.cellset == CELLSET) & (~ref.uncapped.astype(bool)) & (ref.transition == transition) & (ref.rep == rep)]
    ref = ref.assign(coupling=ref.coupling.replace({"kNN": "kNN10"}))
    chk = f.merge(ref, on=["cell_id", "coupling"], suffixes=("", "_ref"))
    n_ref = len(ref)
    assert len(chk) == n_ref, (tag, len(chk), n_ref)
    for c in ("m", "m0", "auroc"):
        assert np.allclose(chk[c], chk[c + "_ref"], atol=1e-6, equal_nan=True), (tag, c)
    # top-k: RQ1 broke exact ties by index (stable argsort); here ties use the uniform
    # expectation. They can differ only for OT rows whose top 10 includes exactly tied
    # weights (underflow to 0 at small eps); those rows are exempt and counted.
    tied = chk.coupling.str.startswith("OT_") & (chk.n_pos_w < 10)
    for c in ("top1", "top5", "top10"):
        ok = np.isclose(chk[c], chk[c + "_ref"], atol=1e-6, equal_nan=True)
        assert (ok | tied).all(), (tag, c)
    f["topk_tie_exempt"] = f.cell_id.isin(chk[tied & ~np.isclose(chk.top10, chk.top10_ref, atol=1e-6)].cell_id) & f.coupling.str.startswith("OT_")
    CACHE.mkdir(parents=True, exist_ok=True)
    f.drop(columns="cell").to_csv(ck, index=False)
    print(tag, "done; RQ1 reproduced on", n_ref, "rows", flush=True)


def summarise():
    OUT.mkdir(parents=True, exist_ok=True)
    f = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_cells.csv.gz"))])
    f["clone"] = f.rep.astype(str) + "|" + f.clone.astype(str)
    f["enrichment"] = f.m / f.m0          # per cell, as RQ1
    f.to_csv(OUT / "cells_forward.csv.gz", index=False)
    rows = []
    for (tr, cp), g in f.groupby(["transition", "coupling"]):
        for rep, gg in [("pooled", g)] + list(g.groupby("rep")):
            v = gg[gg.valid.astype(bool)]
            c = v.groupby("clone")[METRICS + ["enrichment"]].mean()
            r = dict(transition=tr, coupling=cp, rep=rep, n_clones=c.shape[0], n_cells=len(v),
                     coverage_cells=len(v) / len(gg), n_clones_all=gg.clone.nunique())
            rng = np.random.default_rng(0)
            idx = rng.integers(0, len(c), size=(N_BOOT, len(c)))
            arr = c.to_numpy()
            boots = np.stack([np.nanmean(arr[i], 0) for i in idx])
            for k, col in enumerate(METRICS + ["enrichment"]):
                r[col] = c[col].mean()
                r[col + "_lo"], r[col + "_hi"] = np.percentile(boots[:, k], [2.5, 97.5])
            rows.append(r)
    pd.DataFrame(rows).to_csv(OUT / "summary.csv", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    import sklearn
    (OUT / "provenance.txt").write_text("\n".join([
        f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
        "plan: ANALYSIS_PLAN.md amendment A4.4 (2026-11-18, second review round, exploratory)",
        f"python {sys.version.split()[0]} sklearn {sklearn.__version__} numpy {np.__version__} pandas {pd.__version__}",
        f"input {rq1.OBS.relative_to(ROOT)} sha256 {rq1.sha256(rq1.OBS)}",
        f"input {rq1.QC.relative_to(ROOT)} sha256 {rq1.sha256(rq1.QC)}",
        "input RQ1 embeddings and potentials: data/interim/rq1/primary_*_{emb,bal_*,unb_*}.npz",
        "check: RQ1 couplings reproduce results/06_rq1/cells_forward.csv.gz per cell (m, m0, AUROC; top-k "
        "except OT rows with exact ties in the top 10, counted in cells_forward 'topk_tie_exempt': "
        f"{int(f.get('topk_tie_exempt', pd.Series(dtype=bool)).fillna(False).astype(bool).sum())} rows)",
        "note: 'enrichment' = clone-weighted mean of per-cell m/m0, as RQ1"]) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "summarise"])
    a = ap.parse_args()
    if a.cmd == "run":
        obs = rq1.load_obs()
        for tr in rq1.TRANSITIONS:
            for rep in (1, 2, 3):
                run_job(obs, tr, rep)
    else:
        summarise()


if __name__ == "__main__":
    main()
