"""A3.5 — fate agreement on measured labels (ANALYSIS_PLAN.md amendment A3.5; post-review).

Primary cell set, capped RQ1 pools, couplings OT balanced eps = 0.05 (rows rebuilt from the
saved RQ1 potentials), kNN (k = 10) and independent. Per transition x replicate:
  (i)   share of later-pool cells whose label is transferred (no authors' Cluster.Id), and,
        per early cluster, the share of each coupling's transported mass that lands on
        transferred-label cells;
  (ii)  TVD with the predicted distribution renormalised over later cells that carry the
        authors' labels only (clonal distribution from the sisters' authors' labels);
  (iii) TVD with label-transfer k in {5, 15, 30} (RQ1 used k = 15).
Clusters, aggregation and the >= 20-clone rule as in scripts/06_rq1_coupling.py summarise():
per early cluster, clone means then the mean over clones; TVD = 0.5 * L1; 95% percentile
bootstrap over clones (1,000). Bootstrap draws use np.random.default_rng(0) per
(transition, cluster, coupling, variant) row, so CIs differ from RQ1's sequential stream;
point estimates are comparable.
Reproduction check: variant k = 15 must reproduce the RQ1 point TVDs
(results/06_rq1/fate_agreement.csv, primary, capped) to 1e-6, or the script stops.

Usage: python scripts/15_fate_labels.py
"""
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rq1", ROOT / "scripts/06_rq1_coupling.py")
rq1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rq1)

OUT = ROOT / "results/15_fate_labels"
CELLSET = "primary"
KS = (5, 15, 30)
COUPLINGS = ("OT_bal_eps0.05", "kNN", "independent")
N_BOOT = 1000


def job(obs, transition, rep):
    from scipy.spatial.distance import cdist
    from sklearn.neighbors import NearestNeighbors
    tag = f"{CELLSET}_{transition.replace('->', '-')}_r{rep}"
    d, s0, s1 = rq1.job_cells(obs, CELLSET, transition, rep, cap=True)
    z = np.load(rq1.CACHE / f"{tag}_emb.npz", allow_pickle=True)
    assert (z["ids"] == d.cell_id.to_numpy()).all()
    X = z["X"]
    st = d.stage.to_numpy()
    i0, i1 = np.where(st == s0)[0], np.where(st == s1)[0]
    X0, X1 = X[i0], X[i1]
    e, l = d.iloc[i0].reset_index(drop=True), d.iloc[i1].reset_index(drop=True)
    ce = e.clone.fillna("").astype(str).to_numpy()
    cl = np.where(l.fg.to_numpy(), l.clone.fillna("").astype(str).to_numpy(), "__bg__")
    efg = np.where(e.fg.to_numpy())[0]

    measured = l.cluster.to_numpy().astype(object)
    has_meas = ~pd.isna(measured)
    labs = {}
    for k in KS:
        rq1.LABEL_K = k
        labs[k] = rq1.transfer_labels(X, d, s1)
    rq1.LABEL_K = 15
    assert (labs[15] == np.load(rq1.CACHE / f"{tag}_labels.npy", allow_pickle=True)).all(), tag

    # rows of each coupling for the foreground early cells (row-normalised)
    M = cdist(X0, X1, "sqeuclidean").astype(np.float32)
    M /= np.float32(np.median(M))
    pot = np.load(rq1.CACHE / f"{tag}_bal_0.05.npz")
    P = rq1.plan_from_potentials(M[efg], 0.05, pot["log_u"][efg], pot["log_v"]).astype(np.float64)
    del M
    nb = NearestNeighbors(n_neighbors=rq1.KNN_K).fit(X1).kneighbors(X0[efg], return_distance=False)
    Wk = np.zeros((len(efg), len(X1)))
    Wk[np.repeat(np.arange(len(efg)), rq1.KNN_K), nb.ravel()] = 1 / rq1.KNN_K
    W = {"OT_bal_eps0.05": P, "kNN": Wk, "independent": np.full((len(efg), len(X1)), 1 / len(X1))}

    pool = dict(transition=transition, rep=rep, n_later=len(l), n_later_measured=int(has_meas.sum()),
                n_later_fg=int(l.fg.sum()), n_later_fg_measured=int((has_meas & l.fg.to_numpy()).sum()))
    rows = []
    for cp, Wc in W.items():
        Wc = Wc / Wc.sum(1, keepdims=True)
        for r, i in enumerate(efg):
            w = Wc[r]
            sis = cl == ce[i]
            rec = dict(transition=transition, rep=rep, clone=f"{rep}|{ce[i]}",
                       early_cluster=e.cluster.to_numpy()[i], coupling=cp,
                       mass_on_transferred=float(w[~has_meas].sum()))
            for k in KS:
                rec[f"pred_k{k}"] = pd.Series(w).groupby(labs[k]).sum().to_dict()
                rec[f"clonal_k{k}"] = pd.Series(labs[k][sis]).value_counts(normalize=True).to_dict()
            wm = w * has_meas
            rec["pred_measured"] = pd.Series(wm[has_meas] / wm.sum()).groupby(measured[has_meas]).sum().to_dict() \
                if wm.sum() > 0 else {}
            sm = sis & has_meas
            rec["clonal_measured"] = pd.Series(measured[sm]).value_counts(normalize=True).to_dict()
            rec["n_sisters"], rec["n_sisters_measured"] = int(sis.sum()), int(sm.sum())
            rows.append(rec)
    return pool, rows


def tvd_rows(cells):
    out = []
    for tr, g in cells.groupby("transition"):
        nclones = g[g.coupling == "kNN"].groupby("early_cluster").clone.nunique()
        for clu in nclones[nclones >= rq1.FATE_MIN_CLONES].index:
            for cp, gc in g[g.early_cluster == clu].groupby("coupling"):
                for var in [f"k{k}" for k in KS] + ["measured"]:
                    pr, cr = gc[f"pred_{var}"].tolist(), gc[f"clonal_{var}"].tolist()
                    labs = sorted({x for dct in pr + cr for x in dct})
                    Pm = np.array([[dct.get(x, 0.0) for x in labs] for dct in pr])
                    Qm = np.array([[dct.get(x, 0.0) for x in labs] for dct in cr])
                    if var == "measured":     # cells with no measured sister / no measured mass
                        ok = (Qm.sum(1) > 0) & (Pm.sum(1) > 0)
                    else:
                        ok = np.ones(len(Pm), bool)
                    ids = gc.clone.to_numpy()[ok]
                    u = np.unique(ids)
                    Pc = np.stack([Pm[ok][ids == c].mean(0) for c in u])
                    Qc = np.stack([Qm[ok][ids == c].mean(0) for c in u])
                    tvd = 0.5 * np.abs(Pc.mean(0) - Qc.mean(0)).sum()
                    rng = np.random.default_rng(0)
                    idx = rng.integers(0, len(u), size=(N_BOOT, len(u)))
                    bt = [0.5 * np.abs(Pc[i].mean(0) - Qc[i].mean(0)).sum() for i in idx]
                    mt = gc.groupby("clone").mass_on_transferred.mean()
                    out.append(dict(transition=tr, early_cluster=clu, coupling=cp, labels=var,
                                    n_clones=len(u), tvd=tvd, tvd_lo=np.percentile(bt, 2.5),
                                    tvd_hi=np.percentile(bt, 97.5),
                                    mass_on_transferred=mt.mean()))
    return pd.DataFrame(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    obs = rq1.load_obs()
    pools, cells = [], []
    for tr in rq1.TRANSITIONS:
        for rep in (1, 2, 3):
            p, rows = job(obs, tr, rep)
            pools.append(p)
            cells += rows
            print(tr, rep, "done", flush=True)
    pools = pd.DataFrame(pools)
    pooled = pools.groupby("transition")[[c for c in pools.columns if c.startswith("n_")]].sum().reset_index()
    pools = pd.concat([pools, pooled.assign(rep="pooled")], ignore_index=True)
    pools["share_later_transferred"] = 1 - pools.n_later_measured / pools.n_later
    pools["share_fg_later_transferred"] = 1 - pools.n_later_fg_measured / pools.n_later_fg
    pools.to_csv(OUT / "pool_label_shares.csv", index=False)
    cells = pd.DataFrame(cells)
    t = tvd_rows(cells)
    # reproduction check: k = 15 vs RQ1
    ref = pd.read_csv(rq1.OUT / "fate_agreement.csv")
    ref = ref[(ref.cellset == CELLSET) & (~ref.uncapped.astype(bool)) & ref.coupling.isin(COUPLINGS)]
    chk = t[t.labels == "k15"].merge(ref, on=["transition", "early_cluster", "coupling"], suffixes=("", "_rq1"))
    assert len(chk) == len(ref) and np.allclose(chk.tvd, chk.tvd_rq1, atol=1e-6), "k=15 does not reproduce RQ1"
    t.to_csv(OUT / "summary.csv", index=False)
    cells[["transition", "rep", "clone", "early_cluster", "coupling", "mass_on_transferred",
           "n_sisters", "n_sisters_measured"]].to_csv(OUT / "per_cell.csv.gz", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    (OUT / "provenance.txt").write_text("\n".join([
        f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
        "plan: ANALYSIS_PLAN.md amendment A3.5 (2026-11-17, post-review, exploratory)",
        f"python {sys.version.split()[0]} numpy {np.__version__} pandas {pd.__version__}",
        f"input {rq1.OBS.relative_to(ROOT)} sha256 {rq1.sha256(rq1.OBS)}",
        f"input {rq1.QC.relative_to(ROOT)} sha256 {rq1.sha256(rq1.QC)}",
        "input RQ1 embeddings, labels and potentials: data/interim/rq1/primary_*_{emb.npz,labels.npy,bal_0.05.npz}",
        f"check: k=15 reproduces results/06_rq1/fate_agreement.csv point TVDs ({len(chk)} rows)"]) + "\n")


if __name__ == "__main__":
    main()
