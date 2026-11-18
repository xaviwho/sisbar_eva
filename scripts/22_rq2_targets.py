"""A4.2 — RQ2 target designs (ANALYSIS_PLAN.md amendment A4.2, with A5).

Primary cell set, all transitions, the 20 splits of 07_rq2: splits, sampled early cells and
true sisters from the 07 replay (data/interim/rq2_draws/, scripts/12_rq2_pc_input.py).
Training pairs for the OT (RQ1 balanced eps 0.05 row) and independent couplings, with the
admissible later set A = later cells of training clones plus background (as §2) or, for the
foreground-only designs, later cells of training clones only:
  (i)   sampled, A as §2: 50 draws, d = 0..49, rng = np.random.default_rng([s, 1000 + d])
        (draws 0-9 are therefore the A3.4 draws; A3.4's (Pg, D) fits are reused for d = 0-9
        after a per-split check that ConstMedian/ConstMean/Ridge reproduce A3.4 exactly);
  (ii)  expected, A as §2: target = sum_{j in A} w_ij y_j / sum_{j in A} w_ij, with y_j the
        min-max scaled (clipped) later state, i.e. the average of the targets the sampled
        design can draw; independent: mean of y over A;
  (iii) foreground-only: A = later cells of training clones; sampled (50 draws,
        rng = np.random.default_rng([s, 5000 + d])) and expected.
Training early cells draw in table order; one draw consumes one rng.choice per cell.
Inputs: (Pg, D) min-max as 07; 30-PC as A3.1 (RQ3 pipeline fitted on all early cells of the
training clones; raw PC scores). Models: ConstMedian, ConstMean, Ridge (alpha 1), CV-tuned RF
(07 grid, random_state = s). Per A5: in sampled designs RF is fitted on d = 0-9 only, constants
and Ridge on all 50 draws; expected designs: one fit per split, all models. Test: held-out true
pairs (07).
Reported per (transition, input, coupling, design, model): delta vs the constant median of the
same design (draw-averaged per split) and vs the true-lineage model (same model and input;
true-lineage MAE from 07 for (Pg, D) and from A3.1 for 30-PC), each with a percentile bootstrap
over the 20 splits (10,000; np.random.default_rng(0)); constants/Ridge summarised on 50 draws
and on draws 0-9. Foreground share of sampled OT targets: per training early cell, the mass of
its A-restricted OT row on training-clone later cells, averaged by transition x replicate.

Usage: python scripts/22_rq2_targets.py run --transition T | summarise
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
a31 = _load("a31", "12_rq2_pc_input.py")
rq3 = a31.rq3
CELLSET = "primary"
DRAWS = ROOT / "data/interim/rq2_draws"
A34 = ROOT / "data/interim/rq2_redraw"
CACHE = ROOT / "data/interim/rq2_targets"
OUT = ROOT / "results/22_rq2_targets"
N_DRAWS, N_RF = 50, 10
MODELS = ["ConstMedian", "ConstMean", "Ridge", "RandomForest"]


def fit(Xtr, ytr, Xte, seed, with_rf):
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import GridSearchCV
    out = {"ConstMedian": np.tile(np.median(ytr, 0), (len(Xte), 1)),
           "ConstMean": np.tile(ytr.mean(0), (len(Xte), 1)),
           "Ridge": Ridge(alpha=1.0).fit(Xtr, ytr).predict(Xte)}
    if with_rf:
        rf = GridSearchCV(RandomForestRegressor(random_state=seed, n_jobs=1), rq2.RF_GRID,
                          scoring="neg_mean_absolute_error", cv=3).fit(Xtr, ytr)
        out["RandomForest"] = rf.predict(Xte)
    return out


def run(transition):
    obs = pd.read_csv(rq2.OBS, usecols=["cell_id", "stage", "clone"], low_memory=False).set_index("cell_id")
    sc_ = pd.read_csv(rq2.SCORES).set_index("cell_id")
    pools = {r: rq2.load_pool(CELLSET, transition, r, obs) for r in (1, 2, 3)}
    cl = a31.clone_table(pools, obs)
    Xc, libs, eidx = a31.early_matrix(transition, pools, cl)
    tt = transition.replace("->", "-")
    CACHE.mkdir(parents=True, exist_ok=True)
    state = lambda ids: sc_.loc[list(ids), ["Pg", "D"]].to_numpy()
    for s in rq2.SEEDS:
        ck = CACHE / f"{tt}_s{s}.csv"
        if ck.exists():
            continue
        t0 = time.time()
        d = pd.read_csv(DRAWS / f"{CELLSET}_{tt}_s{s}.csv")
        tr, te = d[d.part == "train"], d[d.part == "test"]
        Btr = np.vstack([state(tr.early_id), state(tr.target_true)])
        lo, hi = Btr.min(0), Btr.max(0)
        scl = lambda A: np.clip((A - lo) / np.where(hi > lo, hi - lo, 1), 0, 1)
        yte = scl(state(te.target_true))
        # 30-PC input (A3.1 pipeline)
        key = cl.set_index(["rep", "clone"])
        tr_cells = np.concatenate([[pools[r]["e_ids"][i] for i in key.loc[(r, c)].early]
                                   for r, c in zip(tr.rep, tr.clone)])
        PC = rq3.pca_train_apply(Xc, libs, eidx.get_indexer(tr_cells))
        inputs = {"Pg_D": (scl(state(tr.early_id)), scl(state(te.early_id))),
                  "PC30": (PC[eidx.get_indexer(tr.early_id)], PC[eidx.get_indexer(te.early_id)])}
        # admissible sets, scaled later states and OT rows per replicate
        L, rows, fgshare = {}, {}, []
        for r, p in pools.items():
            lc = obs.loc[p["l_ids"], "clone"].to_numpy()
            fg_all = set(d[d.rep == r].clone)
            fg_train = set(d[(d.rep == r) & (d.part == "train")].clone)
            is_tr = np.array([c in fg_train for c in lc])
            A_full = np.where(np.array([(c not in fg_all) or (c in fg_train) for c in lc]))[0]
            A_fg = np.where(is_tr)[0]
            L[r] = dict(A_full=A_full, A_fg=A_fg, y=scl(state(p["l_ids"])))
            for ei in tr[tr.rep == r].ei:
                w = rq2.ot_row(p, ei).astype(np.float64)
                rows[(r, ei)] = w
                fgshare.append(dict(rep=r, share=w[A_fg].sum() / w[A_full].sum()))
        res = []
        designs = [("sampled", "A_full", 1000), ("expected", "A_full", None),
                   ("fg_sampled", "A_fg", 5000), ("fg_expected", "A_fg", None)]
        for coupling in ("OT", "independent"):
            for design, aset, base in designs:
                n_draw = N_DRAWS if base is not None else 1
                for dd in range(n_draw):
                    if base is not None:
                        rng = np.random.default_rng([s, base + dd])
                        tgt = []
                        for r_, ei in zip(tr.rep, tr.ei):
                            A = L[r_][aset]
                            if coupling == "OT":
                                w = rows[(r_, ei)][A]
                                j = A[rng.choice(len(A), p=w / w.sum())]
                            else:
                                j = rng.choice(A)
                            tgt.append(L[r_]["y"][j])
                        ytr = np.array(tgt)
                    else:
                        ytr = []
                        for r_, ei in zip(tr.rep, tr.ei):
                            A = L[r_][aset]
                            if coupling == "OT":
                                w = rows[(r_, ei)][A]
                                ytr.append((w / w.sum()) @ L[r_]["y"][A])
                            else:
                                ytr.append(L[r_]["y"][A].mean(0))
                        ytr = np.array(ytr)
                    for inp, (Xtr, Xte) in inputs.items():
                        reuse = (design == "sampled" and inp == "Pg_D" and dd < N_RF)
                        with_rf = (base is None or dd < N_RF) and not reuse
                        preds = fit(Xtr, ytr, Xte, s, with_rf)
                        mae = {m: float(np.abs(yh - yte).mean()) for m, yh in preds.items()}
                        if reuse:   # A3.4 fits: check, then take its RF
                            old = pd.read_csv(A34 / f"{CELLSET}_{tt}_{coupling}_s{s}_d{dd}.csv").set_index("model").mae
                            for m in ("ConstMedian", "ConstMean", "Ridge"):
                                assert abs(old[m] - mae[m]) < 1e-9, (tt, s, coupling, dd, m)
                            mae["RandomForest"] = float(old["RandomForest"])
                        for m, v in mae.items():
                            res.append(dict(transition=transition, seed=s, coupling=coupling, design=design,
                                            input=inp, draw=dd, model=m, mae=v,
                                            rf_from_a34=bool(reuse and m == "RandomForest")))
        pd.DataFrame(res).to_csv(ck, index=False)
        pd.DataFrame(fgshare).assign(transition=transition, seed=s).to_csv(CACHE / f"fgshare_{tt}_s{s}.csv", index=False)
        print(transition, "seed", s, f"{time.time() - t0:.0f}s", flush=True)


def boot(v, n=10000):
    rng = np.random.default_rng(0)
    b = rng.choice(v, size=(n, len(v))).mean(1)
    return np.percentile(b, [2.5, 97.5])


def summarise():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*-*_s*.csv")) if not p.name.startswith("fgshare")])
    df.to_csv(OUT / "per_draw.csv.gz", index=False)
    true07 = pd.read_csv(ROOT / "results/07_rq2/per_seed.csv")
    true07 = true07[(true07.cellset == CELLSET) & (true07.coupling == "true")].assign(input="Pg_D")
    truePC = pd.read_csv(ROOT / "results/12_rq2_pc/per_seed.csv")
    truePC = truePC[(truePC.coupling == "true") & (truePC.input == "PC30")]
    tru = pd.concat([true07, truePC])[["transition", "seed", "input", "model", "mae"]].rename(columns={"mae": "mae_true"})
    rows = []
    for (tr, inp, cp, des), g in df.groupby(["transition", "input", "coupling", "design"]):
        cm = g[g.model == "ConstMedian"].set_index(["seed", "draw"]).mae
        for m in MODELS:
            gm = g[g.model == m].set_index(["seed", "draw"]).mae
            for subset, sel in (("all_draws", slice(None)), ("draws_0_9", None)):
                x = gm if sel is not None else gm[gm.index.get_level_values("draw") < N_RF]
                if len(x) == 0 or (subset == "all_draws" and m == "RandomForest" and des.endswith("sampled")):
                    continue
                if subset == "draws_0_9" and not des.endswith("sampled"):
                    continue
                dc = (x - cm.loc[x.index]).groupby("seed").mean()
                mt = tru[(tru.transition == tr) & (tru.input == inp) & (tru.model == m)].set_index("seed").mae_true
                dt = (x.groupby("seed").mean() - mt.loc[dc.index])
                lo, hi = boot(dc.to_numpy())
                lo2, hi2 = boot(dt.to_numpy())
                rows.append(dict(transition=tr, input=inp, coupling=cp, design=des, model=m, draws=subset,
                                 n_splits=len(dc), n_draws=int(x.groupby("seed").size().iloc[0]),
                                 mae=x.mean(), const_median_mae=cm.mean(),
                                 delta_vs_const=dc.mean(), ci_lo=lo, ci_hi=hi,
                                 equivalent_to_const=bool(-rq2.DELTA < lo and hi < rq2.DELTA),
                                 delta_pct_of_const=100 * dc.mean() / cm.mean(),
                                 mae_true=mt.mean(), delta_vs_true=dt.mean(), true_ci_lo=lo2, true_ci_hi=hi2,
                                 sd_over_draws=float(x.groupby("seed").std().mean()) if x.groupby("seed").size().iloc[0] > 1 else np.nan))
    pd.DataFrame(rows).to_csv(OUT / "summary.csv", index=False)
    fs = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("fgshare_*.csv"))])
    fs.groupby(["transition", "rep"]).share.agg(["mean", "median", "min", "max", "size"]).reset_index().to_csv(
        OUT / "ot_foreground_share.csv", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    import sklearn
    (OUT / "provenance.txt").write_text("\n".join([
        f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
        "plan: ANALYSIS_PLAN.md amendment A4.2 with A5 (RF on draws 0-9) (2026-11-18, exploratory)",
        f"python {sys.version.split()[0]} sklearn {sklearn.__version__} numpy {np.__version__} pandas {pd.__version__}",
        f"input {rq2.SCORES.relative_to(ROOT)} sha256 {rq2.sha256(rq2.SCORES)}",
        "input splits: data/interim/rq2_draws (07 replay); A3.4 fits for sampled (Pg, D) draws 0-9 "
        "(data/interim/rq2_redraw), reproduced per split for constants and Ridge",
        "input true-lineage MAE: results/07_rq2/per_seed.csv, results/12_rq2_pc/per_seed.csv"]) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "summarise"])
    ap.add_argument("--transition", nargs="*", default=list(rq2.TRANSITIONS))
    a = ap.parse_args()
    if a.cmd == "run":
        for t in a.transition:
            run(t)
    else:
        summarise()


if __name__ == "__main__":
    main()
