"""A4.7 — CellRank 2 as an established fate-mapping baseline (ANALYSIS_PLAN.md amendment A4.7).

Runs in an isolated environment (venv with system site packages + cellrank 2.4.0; scanpy,
anndata and numpy unchanged): /home/claude/venv_cr/bin/python scripts/25_cellrank.py ...

Per transition x replicate (primary cell set, capped RQ1 pool, RQ1 30-PC embedding):
  * AnnData of the pool cells with obsm X_pca = RQ1 embedding and a categorical time key
    (early stage, later stage); within-stage connectivities from a 30-NN graph in that space
    (scanpy.pp.neighbors, n_neighbors 30, use_rep X_pca, exact kNN);
  * RealTimeKernel with couplings {(early, later): RQ1 balanced OT plan, eps 0.05} (rebuilt from
    the cached potentials), CellRank defaults otherwise, self_transitions = "connectivities",
    conn_weight = 0.2; transition matrix with threshold = "auto" (default);
  * terminal states = the later-stage fate classes of §3 (RQ3): every later cell is labelled with
    its cluster (authors' label, or the RQ1 k = 15 transferred label) mapped to the RQ3 classes of
    that transition (clusters outside them -> "other"); GPCCA estimator, set_terminal_states,
    compute_fate_probabilities for the early cells.
Scored (a) as fate agreement: TVD per early cluster (>= 20 clones) between the mean predicted
class distribution and the clonal class distribution (sisters' classes), clone means then
mean over clones, 1,000-resample clone bootstrap; the OT row (balanced eps 0.05) and the
independent coupling are scored in the same class space for comparison. (b) As RQ3 fate
prediction: a clone's predicted distribution = mean over its early cells; macro-AUROC on the
test clones of each of the 20 RQ3 splits (replayed from scripts/08_rq3_fate.py; no fitting);
compared with the RQ3 representations (results/08_rq3/per_seed.csv, same splits) by paired
bootstrap over seeds. If CellRank cannot run within the workspace limits the item is dropped.

Usage: venv python scripts/25_cellrank.py run --transition T --rep R | summarise
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


rq1 = _load("rq1", "06_rq1_coupling.py")
rq3 = _load("rq3", "08_rq3_fate.py")
CELLSET = "primary"
CACHE = ROOT / "data/interim/cellrank"
OUT = ROOT / "results/25_cellrank"
N_BOOT = 1000


def run_job(obs, transition, rep):
    import anndata as ad
    import cellrank as cr
    import scanpy as sc
    from scipy.spatial.distance import cdist
    tag = f"{CELLSET}_{transition.replace('->', '-')}_r{rep}"
    ck = CACHE / f"{tag}_cells.csv.gz"
    if ck.exists():
        return
    t0 = time.time()
    d, s0, s1 = rq1.job_cells(obs, CELLSET, transition, rep, cap=True)
    z = np.load(rq1.CACHE / f"{tag}_emb.npz", allow_pickle=True)
    assert (z["ids"] == d.cell_id.to_numpy()).all()
    X = z["X"]
    st = d.stage.to_numpy()
    i0, i1 = np.where(st == s0)[0], np.where(st == s1)[0]
    # RQ3 classes for this transition
    obs3, _ = rq3.load(CELLSET)
    early3, fate3, rep3, _ = rq3.clone_table(obs3, transition)
    classes = sorted(set(fate3))
    kept = set(classes) - {"other"}
    late_labels = np.load(rq1.CACHE / f"{tag}_labels.npy", allow_pickle=True)
    late_cls = np.array([c if c in kept else "other" for c in late_labels], dtype=object)
    # OT plan (balanced eps 0.05)
    M = cdist(X[i0], X[i1], "sqeuclidean").astype(np.float32)
    M /= np.float32(np.median(M))
    p = np.load(rq1.CACHE / f"{tag}_bal_0.05.npz")
    P = rq1.plan_from_potentials(M, 0.05, p["log_u"], p["log_v"])
    del M
    order = np.concatenate([i0, i1])          # early cells first, then later cells
    a = ad.AnnData(X=np.zeros((len(order), 1), np.float32))
    a.obs_names = d.cell_id.to_numpy()[order]
    a.obsm["X_pca"] = X[order]
    a.obs["time"] = pd.Categorical([0] * len(i0) + [1] * len(i1), categories=[0, 1], ordered=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # exact kNN (transformer "sklearn"): the default approximate search fails with a numba
        # typing error in this environment for pools > ~8k cells; same 30-NN graph definition
        sc.pp.neighbors(a, n_neighbors=30, use_rep="X_pca", random_state=0, transformer="sklearn")
        k = cr.kernels.RealTimeKernel(a, time_key="time", couplings={(0, 1): P})
        # within-stage connectivities: CellRank recomputes neighbours per time point via
        # scanpy; pass the amendment's 30-NN graph in the RQ1 space (exact kNN, see above)
        k.compute_transition_matrix(self_transitions="connectivities", conn_weight=0.2,
                                    conn_kwargs=dict(n_neighbors=30, use_rep="X_pca",
                                                     transformer="sklearn", random_state=0))
        g = cr.estimators.GPCCA(k)
        term = pd.Series([np.nan] * len(i0) + list(late_cls), index=a.obs_names, dtype=object)
        g.set_terminal_states(states=term.astype("category"))
        g.compute_fate_probabilities()
    fp = g.fate_probabilities
    F = pd.DataFrame(np.asarray(fp), columns=list(fp.names), index=a.obs_names)
    F = F.reindex(columns=classes, fill_value=0.0)
    e = d.iloc[i0].reset_index(drop=True)
    l = d.iloc[i1].reset_index(drop=True)
    efg = np.where(e.fg.to_numpy())[0]
    ce = e.clone.fillna("").astype(str).to_numpy()
    cl = np.where(l.fg.to_numpy(), l.clone.fillna("").astype(str).to_numpy(), "__bg__")
    onehot = pd.get_dummies(pd.Categorical(late_cls, categories=classes)).to_numpy(float)
    W = P[efg].astype(np.float64)
    W /= W.sum(1, keepdims=True)
    rows = []
    for r_, i in enumerate(efg):
        sis = onehot[cl == ce[i]].mean(0)
        rows.append(dict(cell_id=e.cell_id.to_numpy()[i], clone=f"{rep}|{ce[i]}", rep=rep,
                         early_cluster=e.cluster.to_numpy()[i],
                         **{f"cr_{c}": F.loc[e.cell_id.to_numpy()[i], c] for c in classes},
                         **{f"ot_{c}": v for c, v in zip(classes, W[r_] @ onehot)},
                         **{f"ind_{c}": v for c, v in zip(classes, onehot.mean(0))},
                         **{f"clonal_{c}": v for c, v in zip(classes, sis)}))
    CACHE.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).assign(transition=transition).to_csv(ck, index=False)
    print(tag, "done", f"{time.time() - t0:.0f}s", flush=True)


def summarise():
    OUT.mkdir(parents=True, exist_ok=True)
    cells = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_cells.csv.gz"))])
    cells.to_csv(OUT / "cells.csv.gz", index=False)
    frows, arows = [], []
    obs3, _ = rq3.load(CELLSET)
    for tr, g in cells.groupby("transition"):
        classes = [c[3:] for c in g.columns if c.startswith("cr_") and g[c].notna().all()]   # this transition only
        nclo = g.groupby("early_cluster").clone.nunique()
        for clu in nclo[nclo >= rq1.FATE_MIN_CLONES].index:
            gc = g[g.early_cluster == clu]
            Q = gc.groupby("clone")[[f"clonal_{c}" for c in classes]].mean().to_numpy()
            for meth in ("cr", "ot", "ind"):
                Pm = gc.groupby("clone")[[f"{meth}_{c}" for c in classes]].mean().to_numpy()
                rng = np.random.default_rng(0)
                idx = rng.integers(0, len(Q), size=(N_BOOT, len(Q)))
                bt = [0.5 * np.abs(Pm[i].mean(0) - Q[i].mean(0)).sum() for i in idx]
                frows.append(dict(transition=tr, early_cluster=clu, method={"cr": "CellRank2", "ot": "OT_bal_eps0.05",
                                  "ind": "independent"}[meth], n_clones=len(Q),
                                  tvd=0.5 * np.abs(Pm.mean(0) - Q.mean(0)).sum(),
                                  tvd_lo=np.percentile(bt, 2.5), tvd_hi=np.percentile(bt, 97.5)))
        # (b) RQ3 splits
        early3, fate3, rep3, _ = rq3.clone_table(obs3, tr)
        clones = fate3.index.to_numpy()
        Pc = g.groupby("clone")[[f"cr_{c}" for c in classes]].mean()
        Pc.columns = classes
        for seed in rq3.SEEDS:
            rng = np.random.default_rng(seed)
            part = pd.Series(index=clones, dtype=object)
            for r in sorted(rep3.unique()):
                idx = rng.permutation(rep3[rep3 == r].index.to_numpy())
                n_tr, n_va = int(round(0.7 * len(idx))), int(round(0.15 * len(idx)))
                part[idx[:n_tr]], part[idx[n_tr:n_tr + n_va]], part[idx[n_tr + n_va:]] = "train", "val", "test"
            test = part[part == "test"].index
            have = test.intersection(Pc.index)
            yc = fate3[have]
            arows.append(dict(transition=tr, seed=seed, representation="CellRank2",
                              macro_auroc=rq3.macro_auroc(yc.to_numpy(), Pc.loc[have, classes].to_numpy(), classes),
                              n_test_clones=len(test), n_scored=len(have)))
    pd.DataFrame(frows).to_csv(OUT / "fate_agreement.csv", index=False)
    a = pd.DataFrame(arows)
    a.to_csv(OUT / "rq3_per_seed.csv", index=False)
    ref = pd.read_csv(ROOT / "results/08_rq3/per_seed.csv")
    ref = ref[ref.cellset == CELLSET]
    comp = []
    for tr, g in a.groupby("transition"):
        x = g.set_index("seed").macro_auroc
        for rep_name, gr in ref[ref.transition == tr].groupby("representation"):
            y = gr.set_index("seed").macro_auroc.loc[x.index]
            dlt = (x - y).to_numpy()
            rng = np.random.default_rng(0)
            b = rng.choice(dlt, size=(10000, len(dlt))).mean(1)
            comp.append(dict(transition=tr, comparison=f"CellRank2 - {rep_name}", auroc_cellrank=x.mean(),
                             auroc_other=y.mean(), delta=dlt.mean(), ci_lo=np.percentile(b, 2.5),
                             ci_hi=np.percentile(b, 97.5), n_seeds=len(dlt)))
    pd.DataFrame(comp).to_csv(OUT / "rq3_comparison.csv", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    import cellrank
    import scanpy
    (OUT / "provenance.txt").write_text("\n".join([
        f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
        "plan: ANALYSIS_PLAN.md amendment A4.7 (2026-11-18, second review round, exploratory)",
        f"python {sys.version.split()[0]} cellrank {cellrank.__version__} scanpy {scanpy.__version__} numpy {np.__version__}",
        "environment: venv with system site packages + cellrank 2.4.0 (pygpcca 1.1.0, docrep 0.3.2, rich 15.0.0)",
        f"input {rq1.OBS.relative_to(ROOT)} sha256 {rq1.sha256(rq1.OBS)}",
        f"input {rq1.QC.relative_to(ROOT)} sha256 {rq1.sha256(rq1.QC)}",
        "input RQ1 embeddings, labels and balanced eps 0.05 potentials: data/interim/rq1/primary_*",
        "input RQ3 per-seed results: results/08_rq3/per_seed.csv"]) + "\n")


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
