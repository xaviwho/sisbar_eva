"""Validate the block-wise log-domain Sinkhorn against POT before using it (amendment A1).

On a capped job where POT fits in memory (primary, I->II, replicate 2), solve balanced
OT at eps = 0.1 then 0.05 (warm start, as in RQ1) with the blocked solver and compare with
POT's cached eps = 0.05 solution: plan entries, and the RQ1 forward metrics of the
foreground early cells. Pass criterion, fixed before running: mean absolute difference in
own-clone mass < 1e-5 and in AUROC < 1e-3, and max relative plan difference < 1e-2.

Output: results/06_rq1/blocked_solver_validation.csv
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rq1", ROOT / "scripts/06_rq1_coupling.py")
rq1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rq1)


def main():
    obs = rq1.load_obs()
    d, s0, s1 = rq1.job_cells(obs, "primary", "I->II", 2, cap=True)
    tag = "primary_I-II_r2"
    z = np.load(rq1.CACHE / f"{tag}_emb.npz", allow_pickle=True)
    assert (z["ids"] == d.cell_id.to_numpy()).all()
    X, st = z["X"], d.stage.to_numpy()
    X0, X1 = X[st == s0], X[st == s1]
    e, l = d[st == s0].reset_index(drop=True), d[st == s1].reset_index(drop=True)
    from scipy.spatial.distance import cdist
    M = cdist(X0, X1, "sqeuclidean").astype(np.float32)
    M /= np.float32(np.median(M))

    u, v = np.zeros(len(X0), np.float32), np.zeros(len(X1), np.float32)
    u, v, it1, _ = rq1.sinkhorn_log_blocked(M, 0.1, u, v, rq1.STOP, rq1.MAX_ITER)
    u, v = u * np.float32(0.1 / 0.05), v * np.float32(0.1 / 0.05)
    u, v, it2, err = rq1.sinkhorn_log_blocked(M, 0.05, u, v, rq1.STOP, rq1.MAX_ITER)
    pot = np.load(rq1.CACHE / f"{tag}_bal_0.05.npz")
    Pb = rq1.plan_from_potentials(M, 0.05, u, v)
    Pp = rq1.plan_from_potentials(M, 0.05, pot["log_u"], pot["log_v"])
    rel = float(np.abs(Pb - Pp).max() / Pp.max())

    ce = e.clone.fillna("").astype(str).to_numpy()
    cl = np.where(l.fg.to_numpy(), l.clone.fillna("").astype(str).to_numpy(), "__bg__")
    fg = np.where(e.fg.to_numpy())[0]
    mb = pd.DataFrame(rq1.row_metrics(Pb[fg], fg, ce, cl, "blocked"))
    mp = pd.DataFrame(rq1.row_metrics(Pp[fg], fg, ce, cl, "pot"))
    dm = float((mb.m - mp.m).abs().mean())
    da = float((mb.auroc - mp.auroc).abs().mean())
    dt = float((mb.top1 - mp.top1).abs().mean())
    ok = dm < 1e-5 and da < 1e-3 and rel < 1e-2
    out = pd.DataFrame([dict(job=tag, eps=0.05, iters_blocked=it1 + it2 + 2, err_blocked=err,
                             iters_pot=int(pot["niter"]), max_rel_plan_diff=rel,
                             mean_abs_diff_m=dm, mean_abs_diff_auroc=da, mean_abs_diff_top1=dt,
                             mean_m_blocked=mb.m.mean(), mean_m_pot=mp.m.mean(),
                             mean_auroc_blocked=mb.auroc.mean(), mean_auroc_pot=mp.auroc.mean(),
                             passed=ok)])
    out.to_csv(ROOT / "results/06_rq1/blocked_solver_validation.csv", index=False)
    print(out.T.to_string())


if __name__ == "__main__":
    main()
