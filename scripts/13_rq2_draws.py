"""A3.4 — repeated target draws for RQ2 (ANALYSIS_PLAN.md amendment A3.4; post-review).

Primary cell set, all transitions, the 20 split seeds of 07_rq2. For each split seed s the
split and the sampled early/later cells are those of 07 (replayed and saved by
scripts/12_rq2_pc_input.py in data/interim/rq2_draws/). For the OT and independent couplings,
10 new target draws per split are made, draw seed d = 0..9 nested in s:
rng = np.random.default_rng([s, 1000 + d]). Each draw redraws the training early cells'
targets exactly as in 07 (OT: from the early cell's RQ1 OT row, balanced eps = 0.05,
restricted to later cells of training clones plus background, renormalised; independent:
uniform over the same restricted pool). Constant median, constant mean, Ridge (alpha 1) and
the CV-tuned random forest are refitted per draw; the Transformer is omitted (A3.4).
Test: the held-out true pairs of 07; min-max bounds as in 07 (training clones' true pairs).

Reported per (transition, coupling, model): mean over draws of delta = MAE(model) -
MAE(constant median); variance decomposition of delta_{s,d} into between-split
(variance over s of the draw mean) and within-split (mean over s of the variance over d)
components (ddof = 1); percentile bootstrap 95% CI over splits (10,000) of the draw-averaged
delta, with the delta = 0.01 equivalence check, for context.

Usage: python scripts/13_rq2_draws.py run --transition T [--coupling C] | summarise
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
spec = importlib.util.spec_from_file_location("rq2", ROOT / "scripts/07_rq2_forecast.py")
rq2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rq2)

CELLSET = "primary"
DRAWS = ROOT / "data/interim/rq2_draws"
CACHE = ROOT / "data/interim/rq2_redraw"
OUT = ROOT / "results/13_rq2_draws"
N_DRAWS = 10
MODELS = ["ConstMedian", "ConstMean", "Ridge", "RandomForest"]


def fit_predict(Xtr, ytr, Xte, seed):
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import GridSearchCV
    rf = GridSearchCV(RandomForestRegressor(random_state=seed, n_jobs=1), rq2.RF_GRID,
                      scoring="neg_mean_absolute_error", cv=3).fit(Xtr, ytr)
    return {"ConstMedian": np.tile(np.median(ytr, 0), (len(Xte), 1)),
            "ConstMean": np.tile(ytr.mean(0), (len(Xte), 1)),
            "Ridge": Ridge(alpha=1.0).fit(Xtr, ytr).predict(Xte),
            "RandomForest": rf.predict(Xte)}


def run(transition, couplings):
    obs = pd.read_csv(rq2.OBS, usecols=["cell_id", "stage", "clone"], low_memory=False).set_index("cell_id")
    sc_ = pd.read_csv(rq2.SCORES).set_index("cell_id")
    pools = {r: rq2.load_pool(CELLSET, transition, r, obs) for r in (1, 2, 3)}
    tt = transition.replace("->", "-")
    CACHE.mkdir(parents=True, exist_ok=True)
    state = lambda ids: sc_.loc[list(ids), ["Pg", "D"]].to_numpy()
    for seed in rq2.SEEDS:
        d = pd.read_csv(DRAWS / f"{CELLSET}_{tt}_s{seed}.csv")
        allowed = {}
        for r, p in pools.items():
            lc = obs.loc[p["l_ids"], "clone"].to_numpy()
            fg_all = set(d[d.rep == r].clone)
            fg_train = set(d[(d.rep == r) & (d.part == "train")].clone)
            allowed[r] = np.where(np.array([(c not in fg_all) or (c in fg_train) for c in lc]))[0]
        tr, te = d[d.part == "train"], d[d.part == "test"]
        Btr = np.vstack([state(tr.early_id), state(tr.target_true)])
        lo, hi = Btr.min(0), Btr.max(0)
        scl = lambda A: np.clip((A - lo) / np.where(hi > lo, hi - lo, 1), 0, 1)
        Xtr, Xte, yte = scl(state(tr.early_id)), scl(state(te.early_id)), scl(state(te.target_true))
        # OT rows of training early cells, restricted and normalised once per seed
        rows = {}
        if "OT" in couplings:
            for r_, ei in zip(tr.rep, tr.ei):
                w = rq2.ot_row(pools[r_], ei)[allowed[r_]].astype(np.float64)
                rows[(r_, ei)] = w / w.sum()
        for coupling in couplings:
            for dd in range(N_DRAWS):
                ck = CACHE / f"{CELLSET}_{tt}_{coupling}_s{seed}_d{dd}.csv"
                if ck.exists():
                    continue
                t0 = time.time()
                rng = np.random.default_rng([seed, 1000 + dd])
                tgt = []
                for r_, ei in zip(tr.rep, tr.ei):
                    if coupling == "OT":
                        j = allowed[r_][rng.choice(len(allowed[r_]), p=rows[(r_, ei)])]
                    else:
                        j = rng.choice(allowed[r_])
                    tgt.append(pools[r_]["l_ids"][j])
                preds = fit_predict(Xtr, scl(state(tgt)), Xte, seed)
                pd.DataFrame([dict(cellset=CELLSET, transition=transition, seed=seed, draw=dd,
                                   coupling=coupling, model=m, mae=np.abs(yh - yte).mean())
                              for m, yh in preds.items()]).to_csv(ck, index=False)
                print(transition, coupling, "seed", seed, "draw", dd, f"{time.time() - t0:.0f}s", flush=True)


def summarise():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*.csv"))])
    df.to_csv(OUT / "per_draw.csv.gz", index=False)
    rows = []
    for (tr, cp), g in df.groupby(["transition", "coupling"]):
        w = g.pivot_table(index=["seed", "draw"], columns="model", values="mae")
        for m in MODELS[1:]:
            dl = (w[m] - w["ConstMedian"]).unstack("draw")        # seeds x draws
            per_split = dl.mean(1)
            rng = np.random.default_rng(0)
            boots = rng.choice(per_split.to_numpy(), size=(10000, len(per_split))).mean(1)
            lo_, hi_ = np.percentile(boots, [2.5, 97.5])
            between, within = per_split.var(ddof=1), dl.var(axis=1, ddof=1).mean()
            cm = w["ConstMedian"].mean()
            rows.append(dict(transition=tr, coupling=cp, model=m, n_splits=dl.shape[0], n_draws=dl.shape[1],
                             mean_delta=per_split.mean(), ci_lo=lo_, ci_hi=hi_,
                             equivalent=bool(-rq2.DELTA < lo_ and hi_ < rq2.DELTA),
                             mean_delta_pct_of_const=100 * per_split.mean() / cm,
                             var_between_split=between, var_within_split=within,
                             share_within=within / (between + within) if between + within > 0 else np.nan,
                             frac_draws_model_better=float((dl.to_numpy() < 0).mean()),
                             const_median_mae=cm))
    pd.DataFrame(rows).to_csv(OUT / "summary.csv", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    import sklearn
    (OUT / "provenance.txt").write_text("\n".join([
        f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
        "plan: ANALYSIS_PLAN.md amendment A3.4 (2026-11-17, post-review, exploratory)",
        f"python {sys.version.split()[0]} sklearn {sklearn.__version__} numpy {np.__version__} pandas {pd.__version__}",
        f"input {rq2.SCORES.relative_to(ROOT)} sha256 {rq2.sha256(rq2.SCORES)}",
        "input splits/sampled cells: data/interim/rq2_draws/primary_*_s*.csv (07_rq2 replay, 12_rq2_pc_input)",
        "input OT potentials: data/interim/rq1/primary_*_bal_0.05.npz and *_emb.npz (RQ1)"]) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "summarise"])
    ap.add_argument("--transition", nargs="*", default=list(rq2.TRANSITIONS))
    ap.add_argument("--coupling", nargs="*", default=["OT", "independent"])
    a = ap.parse_args()
    if a.cmd == "run":
        for t in a.transition:
            run(t, a.coupling)
    else:
        summarise()


if __name__ == "__main__":
    main()
