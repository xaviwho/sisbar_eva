"""A4.3 — fate labels and per-cluster uncertainty (ANALYSIS_PLAN.md amendment A4.3).

Primary cell set, capped RQ1 pools (per transition x replicate), RQ1 30-PC embedding.

(i) Label-transfer accuracy, leave-one-clone-out. For every later-stage pool cell with an
    authors' label (Cluster.Id), the label is predicted by the RQ1 rule (k = 15 vote among
    labelled cells of the same stage and replicate; ties -> smallest summed distance) with
    every cell of its own clone removed from the reference set. Accuracy, balanced accuracy
    and per-cluster recall by transition x replicate and pooled per transition.
(ii) Paired TVD differences per early cluster: TVD(OT) - TVD(independent) and
    TVD(OT) - TVD(kNN), OT = balanced eps 0.05, kNN k = 10, with all labels (k = 15 transfer,
    as RQ1) and measured labels only (A3.5 rule). TVD as RQ1: per-cell predicted and clonal
    distributions -> clone means -> mean over clones; 0.5 * L1. Paired 95% percentile
    bootstrap over clones (1,000 resamples; the same resampled clones for both couplings;
    np.random.default_rng(0) per row). Pooled (clusters with >= 20 clones) and per replicate
    (clusters with >= 10 clones in that replicate). For a measured-label pair, clones are
    those valid for both couplings (a cell is valid if it has measured sisters and nonzero
    measured mass), so its point TVDs can differ slightly from the unpaired A3.5 values;
    both TVDs are reported with each difference.
(iii) Counts of cross-stage clones and foreground cells (early, later) per transition x
    replicate x early cluster.
Per-cell fate vectors are those of scripts/15_fate_labels.py job() (reused). Reproduction
check: the unpaired pooled TVDs recomputed here equal results/15_fate_labels/summary.csv
(k15 and measured) to 1e-9, or the script stops.

Usage: python scripts/18_fate_uncertainty.py
"""
import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


fl = _load("fl", "15_fate_labels.py")
rq1 = fl.rq1
OUT = ROOT / "results/18_fate_unc"
CELLSET = "primary"
N_BOOT = 1000
PAIRS = [("OT_bal_eps0.05", "independent"), ("OT_bal_eps0.05", "kNN")]


# ------------------------------------------------------------------ (i) label transfer
def loco_labels(X, d, stage):
    from sklearn.neighbors import NearestNeighbors
    idx = np.where(d.stage.to_numpy() == stage)[0]
    lab = d.cluster.to_numpy()[idx]
    known = ~pd.isna(lab)
    Xk, lk = X[idx][known], lab[known]
    ck = d.clone.to_numpy()[idx][known].astype(str)
    maxc = pd.Series(ck).value_counts().max()
    k = min(rq1.LABEL_K + maxc, len(Xk))
    dist, nb = NearestNeighbors(n_neighbors=k).fit(Xk).kneighbors(Xk)
    pred = []
    for q in range(len(Xk)):
        keep = ck[nb[q]] != ck[q]                 # drops the cell itself and its whole clone
        nr, dr = nb[q][keep][:rq1.LABEL_K], dist[q][keep][:rq1.LABEL_K]
        assert len(nr) == rq1.LABEL_K
        t = pd.DataFrame({"l": lk[nr], "d": dr}).groupby("l").agg(n=("d", "size"), s=("d", "sum"))
        pred.append(t.sort_values(["n", "s"], ascending=[False, True]).index[0])
    return pd.DataFrame({"cell_id": d.cell_id.to_numpy()[idx][known], "clone": ck,
                         "true": lk, "pred": pred, "fg": d.fg.to_numpy()[idx][known]})


def acc_rows(df, keys):
    from sklearn.metrics import balanced_accuracy_score
    out, rec = [], []
    for k, g in df.groupby(keys):
        k = k if isinstance(k, tuple) else (k,)
        out.append(dict(zip(keys, k), n_cells=len(g), n_classes=g.true.nunique(),
                        accuracy=(g.true == g.pred).mean(),
                        balanced_accuracy=balanced_accuracy_score(g.true, g.pred),
                        accuracy_fg=(g[g.fg].true == g[g.fg].pred).mean() if g.fg.any() else np.nan))
        for c, gc in g.groupby("true"):
            rec.append(dict(zip(keys, k), cluster=c, n_cells=len(gc), recall=(gc.pred == c).mean()))
    return pd.DataFrame(out), pd.DataFrame(rec)


# ------------------------------------------------------------------ (ii) paired TVD
def clone_matrices(gc, var, labs):
    pr, cr = gc[f"pred_{var}"].tolist(), gc[f"clonal_{var}"].tolist()
    Pm = np.array([[x.get(l, 0.0) for l in labs] for x in pr])
    Qm = np.array([[x.get(l, 0.0) for l in labs] for x in cr])
    ok = ((Qm.sum(1) > 0) & (Pm.sum(1) > 0)) if var == "measured" else np.ones(len(Pm), bool)
    ids = gc.clone.to_numpy()
    valid = pd.Series(ok).groupby(ids).all()           # clone valid if all its cells valid
    return Pm, Qm, ok, ids, valid


def tvd_of(Pc, Qc):
    return 0.5 * np.abs(Pc.mean(0) - Qc.mean(0)).sum()


def paired(g, var, min_clones, scope):
    rows = []
    nclones = g[g.coupling == "kNN"].groupby("early_cluster").clone.nunique()
    for clu in nclones[nclones >= min_clones].index:
        gcl = g[g.early_cluster == clu]
        labs = sorted({x for col in (f"pred_{var}", f"clonal_{var}") for dct in gcl[col] for x in dct})
        mats = {}
        for cp, gc in gcl.groupby("coupling"):
            Pm, Qm, ok, ids, valid = clone_matrices(gc, var, labs)
            u = valid.index.to_numpy()
            Pc = np.stack([Pm[ok & (ids == c)].mean(0) if (ok & (ids == c)).any() else np.full(len(labs), np.nan) for c in u])
            Qc = np.stack([Qm[ok & (ids == c)].mean(0) if (ok & (ids == c)).any() else np.full(len(labs), np.nan) for c in u])
            mats[cp] = pd.DataFrame(np.hstack([Pc, Qc]), index=u)
        for a, b in PAIRS:
            A, B = mats[a].dropna(), mats[b].dropna()
            common = A.index.intersection(B.index)
            A, B = A.loc[common].to_numpy(), B.loc[common].to_numpy()
            L = len(labs)
            ta, tb = tvd_of(A[:, :L], A[:, L:]), tvd_of(B[:, :L], B[:, L:])
            rng = np.random.default_rng(0)
            idx = rng.integers(0, len(common), size=(N_BOOT, len(common)))
            bd = np.array([tvd_of(A[i, :L], A[i, L:]) - tvd_of(B[i, :L], B[i, L:]) for i in idx])
            rows.append(dict(scope=scope, early_cluster=clu, labels=var, comparison=f"{a} - {b}",
                             n_clones=len(common), tvd_a=ta, tvd_b=tb, diff=ta - tb,
                             diff_lo=np.percentile(bd, 2.5), diff_hi=np.percentile(bd, 97.5),
                             frac_boot_a_better=float((bd < 0).mean())))
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    obs = rq1.load_obs()
    loco, cells, counts = [], [], []
    for tr in rq1.TRANSITIONS:
        s0, s1, _ = rq1.TRANSITIONS[tr]
        for rep in (1, 2, 3):
            tag = f"{CELLSET}_{tr.replace('->', '-')}_r{rep}"
            d, _, _ = rq1.job_cells(obs, CELLSET, tr, rep, cap=True)
            z = np.load(rq1.CACHE / f"{tag}_emb.npz", allow_pickle=True)
            assert (z["ids"] == d.cell_id.to_numpy()).all()
            loco.append(loco_labels(z["X"], d, s1).assign(transition=tr, rep=rep))
            e = d[(d.stage == s0) & d.fg]
            l = d[(d.stage == s1) & d.fg]
            ce = e.groupby("cluster").agg(n_clones=("clone", "nunique"), n_early_fg_cells=("cell_id", "size"))
            lc = l.groupby("clone").size()
            ce["n_later_fg_cells"] = [int(lc.reindex(e[e.cluster == c].clone.unique()).sum()) for c in ce.index]
            counts.append(ce.reset_index().rename(columns={"cluster": "early_cluster"}).assign(transition=tr, rep=rep))
            _, rows = fl.job(obs, tr, rep)
            cells += rows
            print(tr, rep, "done", flush=True)
    loco = pd.concat(loco, ignore_index=True)
    loco.to_csv(OUT / "label_transfer_cells.csv.gz", index=False)
    a1, r1 = acc_rows(loco, ["transition", "rep"])
    a2, r2 = acc_rows(loco, ["transition"])
    pd.concat([a1, a2.assign(rep="pooled")]).to_csv(OUT / "label_transfer_accuracy.csv", index=False)
    pd.concat([r1, r2.assign(rep="pooled")]).to_csv(OUT / "label_transfer_recall.csv", index=False)
    counts = pd.concat(counts, ignore_index=True)
    counts.to_csv(OUT / "counts.csv", index=False)

    cells = pd.DataFrame(cells)
    # reproduction check against A3.5 (unpaired pooled TVDs)
    ref = pd.read_csv(ROOT / "results/15_fate_labels/summary.csv")
    chk = fl.tvd_rows(cells).merge(ref, on=["transition", "early_cluster", "coupling", "labels"], suffixes=("", "_ref"))
    assert len(chk) == len(ref) and np.allclose(chk.tvd, chk.tvd_ref, atol=1e-9), "does not reproduce A3.5"

    rows = []
    for tr, g in cells.groupby("transition"):
        for var in ("k15", "measured"):
            rows += [dict(transition=tr, rep="pooled", **r) for r in paired(g, var, rq1.FATE_MIN_CLONES, "pooled")]
            for rep, gr in g.groupby("rep"):
                rows += [dict(transition=tr, rep=rep, **r) for r in paired(gr, var, 10, "replicate")]
    pd.DataFrame(rows).to_csv(OUT / "paired_tvd.csv", index=False)

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    (OUT / "provenance.txt").write_text("\n".join([
        f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
        "plan: ANALYSIS_PLAN.md amendment A4.3 (2026-11-18, second review round, exploratory)",
        f"python {sys.version.split()[0]} numpy {np.__version__} pandas {pd.__version__}",
        f"input {rq1.OBS.relative_to(ROOT)} sha256 {rq1.sha256(rq1.OBS)}",
        f"input {rq1.QC.relative_to(ROOT)} sha256 {rq1.sha256(rq1.QC)}",
        "input RQ1 embeddings, labels and potentials: data/interim/rq1/primary_*",
        f"check: per-cell fate vectors reproduce results/15_fate_labels/summary.csv ({len(chk)} rows)"]) + "\n")


if __name__ == "__main__":
    main()
