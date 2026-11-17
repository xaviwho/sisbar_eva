"""A3.3 — leave-one-replicate-out for RQ2 ((Pg, D) input) and RQ3 (ANALYSIS_PLAN.md A3.3).

Primary cell set. For each transition, three folds: train on all cross-stage clones of two
replicates (15% of those clones held out for validation, rng seed 0, stratified by
replicate), test on all cross-stage clones of the third. Descriptive only (n = 3 folds), no
p-values.

RQ2 (as scripts/07_rq2_forecast.py except for the split): one early and one later cell
drawn per clone; training/validation targets by coupling (true sister; OT row of the
transductive RQ1 coupling, balanced eps = 0.05, restricted to later cells of training clones
plus background of the same replicate and renormalised; independent: uniform over that
restricted pool); test on the held-out replicate's true pairs. Min-max bounds from the
training clones' true pairs. Models and settings as in 07 (Transformer early stopping on
the validation clones). Reported: MAE per model and delta vs constant median.
Random draws use one np.random.default_rng(0) per (transition, held-out replicate), in the
order split -> sampled cells -> targets (couplings in 07 order, training then validation).

RQ3 (as scripts/08_rq3_fate.py except for the split): fate labels as in 08 (all clones of
the transition); training = early cells of the training clones; 30-PC and LDA fitted on them;
logistic regression per representation; test = clone-mean probabilities for the held-out
replicate's clones. Reported: macro-AUROC and balanced accuracy per representation.

Usage: python scripts/14_loro.py rq2 [--transition ...] | rq3 | summarise
"""
import argparse
import importlib.util
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


rq2 = _load("rq2", "07_rq2_forecast.py")
rq3 = _load("rq3", "08_rq3_fate.py")
pc_mod = _load("a31", "12_rq2_pc_input.py")
CELLSET = "primary"
CACHE = ROOT / "data/interim/loro"
OUT = ROOT / "results/14_loro"


def split_two(rng, ids_by_rep):
    """15% validation within each training replicate (stratified)."""
    part = {}
    for r in sorted(ids_by_rep):
        idx = rng.permutation(ids_by_rep[r])
        n_va = int(round(0.15 * len(idx)))
        for k, c in enumerate(idx):
            part[c] = "val" if k < n_va else "train"
    return part


def run_rq2(transition):
    import torch
    torch.set_num_threads(1)
    obs = pd.read_csv(rq2.OBS, usecols=["cell_id", "stage", "clone"], low_memory=False).set_index("cell_id")
    sc_ = pd.read_csv(rq2.SCORES).set_index("cell_id")
    pools = {r: rq2.load_pool(CELLSET, transition, r, obs) for r in (1, 2, 3)}
    cl = pc_mod.clone_table(pools, obs)
    state = lambda ids: sc_.loc[list(ids), ["Pg", "D"]].to_numpy()
    tt = transition.replace("->", "-")
    CACHE.mkdir(parents=True, exist_ok=True)
    for h in (1, 2, 3):
        ck = CACHE / f"rq2_{tt}_heldout{h}.csv"
        if ck.exists():
            continue
        t0 = time.time()
        rng = np.random.default_rng(0)
        c = cl.copy()
        part = split_two(rng, {r: g.index.to_numpy() for r, g in c[c.rep != h].groupby("rep")})
        c["part"] = [part.get(i, "test") for i in c.index]
        c["ei"] = [rng.choice(x) for x in c.early]
        c["li"] = [rng.choice(x) for x in c.late]
        allowed = {}
        for r, p in pools.items():
            lc = obs.loc[p["l_ids"], "clone"].to_numpy()
            fg_all = set(c[c.rep == r].clone)
            fg_train = set(c[(c.rep == r) & (c.part == "train")].clone)
            allowed[r] = np.where(np.array([(x not in fg_all) or (x in fg_train) for x in lc]))[0]

        def targets(sub, coupling):
            out = []
            for r_, ei, li in zip(sub.rep, sub.ei, sub.li):
                p = pools[r_]
                if coupling == "true":
                    j = li
                elif coupling == "OT":
                    w = rq2.ot_row(p, ei)[allowed[r_]].astype(np.float64)
                    j = allowed[r_][rng.choice(len(w), p=w / w.sum())]
                else:
                    j = rng.choice(allowed[r_])
                out.append(p["l_ids"][j])
            return np.array(out)

        eids = lambda sub: np.array([pools[r_]["e_ids"][i] for r_, i in zip(sub.rep, sub.ei)])
        tr, va, te = c[c.part == "train"], c[c.part == "val"], c[c.part == "test"]
        Btr = np.vstack([state(eids(tr)), state(targets(tr, "true"))])
        lo, hi = Btr.min(0), Btr.max(0)
        scl = lambda A: np.clip((A - lo) / np.where(hi > lo, hi - lo, 1), 0, 1)
        Xte, yte = scl(state(eids(te))), scl(state(targets(te, "true")))
        res = []
        for coupling in rq2.COUPLINGS:
            Xtr, ytr = scl(state(eids(tr))), scl(state(targets(tr, coupling)))
            Xva, yva = scl(state(eids(va))), scl(state(targets(va, coupling)))
            preds, rfp = rq2.fit_models(Xtr, ytr, Xva, yva, 0)
            for m, f in preds.items():
                err = np.abs(f(Xte) - yte)
                res.append(dict(transition=transition, heldout_rep=h, coupling=coupling, model=m,
                                mae=err.mean(), mae_Pg=err[:, 0].mean(), mae_D=err[:, 1].mean(),
                                n_train=len(tr), n_val=len(va), n_test=len(te),
                                rf_params=str(rfp) if m == "RandomForest" else ""))
        pd.DataFrame(res).to_csv(ck, index=False)
        print("rq2", transition, "held-out", h, f"{time.time() - t0:.0f}s", flush=True)


def run_rq3():
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import balanced_accuracy_score
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    CACHE.mkdir(parents=True, exist_ok=True)
    obs, scores = rq3.load(CELLSET)
    for transition in rq3.TRANSITIONS:
        early, fate, rep, late = rq3.clone_table(obs, transition)
        X, _ = rq3.early_counts(early, CELLSET, transition)
        libs = early.library.to_numpy()
        S = scores.loc[early.cell_id, ["P_tcbb", "Pg", "D"]].to_numpy()
        tt = transition.replace("->", "-")
        for h in (1, 2, 3):
            ck = CACHE / f"rq3_{tt}_heldout{h}.csv"
            if ck.exists():
                continue
            rng = np.random.default_rng(0)
            part = split_two(rng, {r: rep[rep == r].index.to_numpy() for r in sorted(rep.unique()) if r != h})
            cp = early.clone.map(lambda x: part.get(x, "test")).to_numpy()
            tr, te = np.where(cp == "train")[0], np.where(cp == "test")[0]
            y = early.fate.to_numpy()
            PC = rq3.pca_train_apply(X, libs, tr)
            lda = LinearDiscriminantAnalysis().fit(PC[tr], y[tr])
            feats = {"P_tcbb_D": S[:, [0, 2]], "Pg_D": S[:, [1, 2]], "PC30": PC, "fate_axis": lda.transform(PC)}
            test_clones = early.clone.to_numpy()[te]
            yc = fate[pd.unique(test_clones)]
            res = []
            for name, F in feats.items():
                clf = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000)).fit(F[tr], y[tr])
                P = pd.DataFrame(clf.predict_proba(F[te]), columns=clf.classes_)
                P["clone"] = test_clones
                Pc = P.groupby("clone").mean().loc[yc.index]
                classes = list(clf.classes_)
                ytrue = yc.to_numpy()
                res.append(dict(transition=transition, heldout_rep=h, representation=name,
                                macro_auroc=rq3.macro_auroc(ytrue, Pc[classes].to_numpy(), classes),
                                balanced_acc=balanced_accuracy_score(ytrue, np.array(classes)[Pc[classes].to_numpy().argmax(1)]),
                                chance_bacc=1 / len(set(ytrue)), n_test_clones=len(yc),
                                n_train_cells=len(tr), hvg_batch_key=rq3.pca_train_apply.batch_hvg))
            pd.DataFrame(res).to_csv(ck, index=False)
            print("rq3", transition, "held-out", h, flush=True)


def summarise():
    OUT.mkdir(parents=True, exist_ok=True)
    a = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("rq2_*.csv"))])
    a.to_csv(OUT / "rq2_per_fold.csv", index=False)
    cm = a[a.model == "ConstMedian"][["transition", "heldout_rep", "coupling", "mae"]].rename(columns={"mae": "const_median_mae"})
    a = a.merge(cm, on=["transition", "heldout_rep", "coupling"])
    a["delta_vs_const_median"] = a.mae - a.const_median_mae
    a["delta_pct_of_const"] = 100 * a.delta_vs_const_median / a.const_median_mae
    a.to_csv(OUT / "rq2_per_fold.csv", index=False)
    a.groupby(["transition", "coupling", "model"]).agg(
        mae_mean=("mae", "mean"), delta_mean=("delta_vs_const_median", "mean"),
        delta_min=("delta_vs_const_median", "min"), delta_max=("delta_vs_const_median", "max"),
        delta_pct_mean=("delta_pct_of_const", "mean"), n_folds=("mae", "size")).reset_index().to_csv(
        OUT / "rq2_summary.csv", index=False)
    b = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("rq3_*.csv"))])
    b.to_csv(OUT / "rq3_per_fold.csv", index=False)
    b.groupby(["transition", "representation"]).agg(
        auroc_mean=("macro_auroc", "mean"), auroc_min=("macro_auroc", "min"), auroc_max=("macro_auroc", "max"),
        bacc_mean=("balanced_acc", "mean"), bacc_min=("balanced_acc", "min"), bacc_max=("balanced_acc", "max"),
        n_folds=("macro_auroc", "size")).reset_index().to_csv(OUT / "rq3_summary.csv", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    (OUT / "provenance.txt").write_text("\n".join([
        f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
        "plan: ANALYSIS_PLAN.md amendment A3.3 (2026-11-17, post-review, exploratory; descriptive, no p-values)",
        f"python {sys.version.split()[0]} numpy {np.__version__} pandas {pd.__version__}",
        f"input {rq2.SCORES.relative_to(ROOT)} sha256 {rq2.sha256(rq2.SCORES)}",
        f"input {rq3.OBS.relative_to(ROOT)} sha256 {rq2.sha256(rq3.OBS)}",
        f"input {rq3.QC.relative_to(ROOT)} sha256 {rq2.sha256(rq3.QC)}",
        "input OT potentials: data/interim/rq1/primary_*_bal_0.05.npz and *_emb.npz (RQ1)",
        "input early counts: data/interim/rq3/primary_*_early_counts.npz (08_rq3)"]) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["rq2", "rq3", "summarise"])
    ap.add_argument("--transition", nargs="*", default=list(rq2.TRANSITIONS))
    a = ap.parse_args()
    if a.cmd == "rq2":
        for t in a.transition:
            run_rq2(t)
    elif a.cmd == "rq3":
        run_rq3()
    else:
        summarise()


if __name__ == "__main__":
    main()
