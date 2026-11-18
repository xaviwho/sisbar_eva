"""A4.1(a) — lineage-grounded partner error (ANALYSIS_PLAN.md amendment A4.1(a)).

Primary cell set, capped RQ1 pools, RQ1 30-PC embedding, per transition x replicate.
For each foreground early cell i (all have >= 1 later sister), with S = its later sisters
and w the row-normalised coupling over all later pool cells:
    D(i) = sum_j w_ij * mean_{s in S} ||x_j - x_s||^2        (raw squared Euclidean, 30-PC space)
computed exactly as sum_j w_ij ||x_j - mu_S||^2 + V_S, V_S = mean_{s in S} ||x_s - mu_S||^2.
Couplings: OT balanced eps 0.05, OT unbalanced (tau 0.9) eps 0.05 (from the RQ1 potentials),
kNN (k = 10), independent (uniform).
Floor (clones with >= 2 later sisters): D_sis = mean over pairs of distinct later sisters of
||x_s - x_s'||^2 = 2 |S| / (|S| - 1) * V_S.
Fixed here before any run (the amendment leaves the estimator open): clone-weighted means
(cell -> clone mean -> mean over clones), and the normalised partner error as a ratio of
clone-weighted means over clones with >= 2 later sisters,
    R = (mean D_coupling - mean D_sis) / (mean D_indep - mean D_sis)
(0 = as close as a sister, 1 = as far as a random later cell); a per-clone ratio is unstable
when its denominator is small. 95% percentile bootstrap over clones (1,000 resamples; the
same resampled clones for every quantity of a row; np.random.default_rng(0) per row), pooled
and per replicate. Row statistics: effective number of partners exp(H(w_i.)) and the weight on
the 10 nearest later cells.
Reproduction check: own-clone mass of each coupling equals the RQ1 per-cell values to 1e-6.

Usage: python scripts/20_partner_error.py
"""
import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rq1", ROOT / "scripts/06_rq1_coupling.py")
rq1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rq1)

CELLSET = "primary"
OUT = ROOT / "results/20_partner_error"
N_BOOT = 1000
COUPLINGS = ["OT_bal_eps0.05", "OT_unb_eps0.05", "kNN10", "independent"]


def job(obs, transition, rep):
    from scipy.spatial.distance import cdist
    from sklearn.neighbors import NearestNeighbors
    tag = f"{CELLSET}_{transition.replace('->', '-')}_r{rep}"
    d, s0, s1 = rq1.job_cells(obs, CELLSET, transition, rep, cap=True)
    z = np.load(rq1.CACHE / f"{tag}_emb.npz", allow_pickle=True)
    assert (z["ids"] == d.cell_id.to_numpy()).all()
    X, st = z["X"].astype(np.float64), d.stage.to_numpy()
    i0, i1 = np.where(st == s0)[0], np.where(st == s1)[0]
    X0, X1 = X[i0], X[i1]
    e, l = d.iloc[i0].reset_index(drop=True), d.iloc[i1].reset_index(drop=True)
    ce = e.clone.fillna("").astype(str).to_numpy()
    cl = np.where(l.fg.to_numpy(), l.clone.fillna("").astype(str).to_numpy(), "__bg__")
    efg = np.where(e.fg.to_numpy())[0]
    M = cdist(X0[efg].astype(np.float32), X1.astype(np.float32), "sqeuclidean").astype(np.float32)
    med = np.float32(np.median(cdist(X0.astype(np.float32), X1.astype(np.float32), "sqeuclidean").astype(np.float32)))
    M /= med
    nb = NearestNeighbors(n_neighbors=10).fit(X1).kneighbors(X0[efg], return_distance=False)
    W = {}
    for kind in ("bal", "unb"):
        p = np.load(rq1.CACHE / f"{tag}_{kind}_0.05.npz")
        W[f"OT_{kind}_eps0.05"] = rq1.plan_from_potentials(M, 0.05, p["log_u"][efg], p["log_v"]).astype(np.float64)
    Wk = np.zeros((len(efg), len(X1)))
    Wk[np.repeat(np.arange(len(efg)), 10), nb.ravel()] = 0.1
    W["kNN10"] = Wk
    W["independent"] = np.full((len(efg), len(X1)), 1.0 / len(X1))
    del M
    # sister statistics per clone
    stats = {}
    for c in np.unique(ce[efg]):
        S = X1[cl == c]
        mu = S.mean(0)
        V = float(((S - mu) ** 2).sum(1).mean())
        stats[c] = (mu, V, len(S), 2 * len(S) / (len(S) - 1) * V if len(S) >= 2 else np.nan)
    sq1 = (X1 ** 2).sum(1)
    rows = []
    for r, i in enumerate(efg):
        c = ce[i]
        mu, V, ns, dsis = stats[c]
        dist_mu = sq1 - 2 * X1 @ mu + mu @ mu            # ||x_j - mu_S||^2
        pos = cl == c
        rec = dict(cell_id=e.cell_id.to_numpy()[i], clone=c, n_sisters=ns, D_sis=dsis,
                   early_cluster=e.cluster.to_numpy()[i])
        for cp, Wc in W.items():
            w = Wc[r] / Wc[r].sum()
            nz = w[w > 0]
            rec[f"D_{cp}"] = float(w @ dist_mu + V)
            rec[f"m_{cp}"] = float(w[pos].sum())
            rec[f"neff_{cp}"] = float(np.exp(-(nz * np.log(nz)).sum()))
            rec[f"wnn10_{cp}"] = float(w[nb[r]].sum())
        rows.append(rec)
    f = pd.DataFrame(rows).assign(transition=transition, rep=rep)
    ref = pd.read_csv(rq1.OUT / "cells_forward.csv.gz")
    ref = ref[(ref.cellset == CELLSET) & (~ref.uncapped.astype(bool)) & (ref.transition == transition) & (ref.rep == rep)]
    for cp in COUPLINGS:
        rr = ref[ref.coupling == cp.replace("kNN10", "kNN")].set_index("cell_id").m
        assert np.allclose(f.set_index("cell_id")[f"m_{cp}"].loc[rr.index], rr, atol=1e-6), (tag, cp)
    return f


def summarise_rows(f):
    out = []
    for tr, g in f.groupby("transition"):
        for rep, gg in [("pooled", g)] + list(g.groupby("rep")):
            cols = [c for c in gg.columns if c.startswith(("D_", "neff_", "wnn10_"))]
            c = gg.groupby("clone")[cols].mean()
            c2 = c[c.D_sis.notna()]
            rng = np.random.default_rng(0)
            idx = rng.integers(0, len(c), size=(N_BOOT, len(c)))
            idx2 = rng.integers(0, len(c2), size=(N_BOOT, len(c2)))
            for cp in COUPLINGS:
                r = dict(transition=tr, rep=rep, coupling=cp, n_clones=len(c), n_clones_floor=len(c2))
                for q in ("D", "neff", "wnn10"):
                    v = c[f"{q}_{cp}"].to_numpy()
                    b = np.array([v[i].mean() for i in idx])
                    r[q], r[q + "_lo"], r[q + "_hi"] = v.mean(), *np.percentile(b, [2.5, 97.5])
                a, ind, s = (c2[f"D_{cp}"].to_numpy(), c2["D_independent"].to_numpy(), c2["D_sis"].to_numpy())
                R = lambda ii: (a[ii].mean() - s[ii].mean()) / (ind[ii].mean() - s[ii].mean())
                rb = np.array([R(i) for i in idx2])
                r["D_floor_clones"], r["D_sis"] = a.mean(), s.mean()
                r["R"], r["R_lo"], r["R_hi"] = R(np.arange(len(a))), *np.percentile(rb, [2.5, 97.5])
                out.append(r)
    return pd.DataFrame(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    obs = rq1.load_obs()
    f = []
    for tr in rq1.TRANSITIONS:
        for rep in (1, 2, 3):
            f.append(job(obs, tr, rep))
            print(tr, rep, "done", flush=True)
    f = pd.concat(f, ignore_index=True)
    f["clone"] = f.rep.astype(str) + "|" + f.clone.astype(str)
    f.to_csv(OUT / "cells.csv.gz", index=False)
    summarise_rows(f).to_csv(OUT / "summary.csv", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    (OUT / "provenance.txt").write_text("\n".join([
        f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
        "plan: ANALYSIS_PLAN.md amendment A4.1(a) (2026-11-18, second review round, exploratory)",
        f"python {sys.version.split()[0]} numpy {np.__version__} pandas {pd.__version__}",
        f"input {rq1.OBS.relative_to(ROOT)} sha256 {rq1.sha256(rq1.OBS)}",
        f"input {rq1.QC.relative_to(ROOT)} sha256 {rq1.sha256(rq1.QC)}",
        "input RQ1 embeddings and potentials: data/interim/rq1/primary_*_{emb,bal_0.05,unb_0.05}.npz",
        "check: own-clone mass per cell equals results/06_rq1/cells_forward.csv.gz (1e-6)"]) + "\n")


if __name__ == "__main__":
    main()
