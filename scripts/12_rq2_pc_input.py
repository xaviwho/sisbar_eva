"""A3.1 — RQ2 with a 30-PC input (ANALYSIS_PLAN.md amendment A3.1; post-review, exploratory).

Primary cell set, all three transitions, 20 seeds.

Reuse of the 07_rq2 draws. scripts/07_rq2_forecast.py draws, per seed, from one
np.random.default_rng(seed) in a fixed order: split -> sampled early cell -> sampled later
cell -> (true targets need no draws) -> for each coupling: training targets, validation
targets. This script replays that exact sequence with 07's own functions, so the split, the
sampled cells and every coupling's targets are identical to 07. Reproduction is checked per
seed: ConstMedian, ConstMean and Ridge on the (Pg, D) input (deterministic) must equal the
07 cached MAEs to 1e-9, or the script stops. The replayed draws are saved to
data/interim/rq2_draws/primary_<transition>_s<seed>.csv for reuse (A3.4, A3.3).

30-PC input: the RQ3 pipeline (scripts/08_rq3_fate.py pca_train_apply: HVGs seurat_v3 on raw
counts of training early cells, batch = library with the RQ3 fallback; CP10k, log1p; scale
with training mean/SD, clip 10; PCA 30, arpack, random_state 0), fitted on all early cells of
the seed's training clones and applied to the sampled early cells. Raw PC scores are the
input (no further scaling; none is specified). Target: the later cell's (Pg, D), min-max as
in 07. Models as in 07; Transformer input size 30.

Per-seed MAE is saved for both inputs: (Pg, D) taken from the 07 cache (all five models),
30-PC computed here. Comparison: 30-PC minus (Pg, D) per model x coupling; paired bootstrap
95% CI (10,000), Wilcoxon, Holm within (transition, coupling), d_z, equivalence iff the CI
lies within +/- 0.01.

Usage: python scripts/12_rq2_pc_input.py run [--transition ...] | summarise
"""
import argparse
import functools
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

CELLSET = "primary"
CACHE = ROOT / "data/interim/rq2_pc"
DRAWS = ROOT / "data/interim/rq2_draws"
OUT = ROOT / "results/12_rq2_pc"
_ORIG_MT = rq2.make_transformer


def set_input_size(k):
    rq2.make_transformer = lambda: functools.partial(_ORIG_MT(), input_size=k)


def replay(seed, cl, pools, obs, sc_):
    """Replay 07's draws for one seed. Returns the clone table with sampled ids and targets."""
    rng = np.random.default_rng(seed)
    cl = cl.copy()
    part = np.empty(len(cl), dtype=object)
    for r, g in cl.groupby("rep"):
        idx = rng.permutation(g.index.to_numpy())
        n_tr, n_va = int(round(0.7 * len(idx))), int(round(0.15 * len(idx)))
        part[idx[:n_tr]], part[idx[n_tr:n_tr + n_va]], part[idx[n_tr + n_va:]] = "train", "val", "test"
    cl["part"] = part
    cl["ei"] = [rng.choice(x) for x in cl.early]
    cl["li"] = [rng.choice(x) for x in cl.late]
    allowed = {}
    for r, p in pools.items():
        lc = obs.loc[p["l_ids"], "clone"].to_numpy()
        fg_all = set(cl[cl.rep == r].clone)
        fg_train = set(cl[(cl.rep == r) & (cl.part == "train")].clone)
        ok = np.array([(c not in fg_all) or (c in fg_train) for c in lc])
        allowed[r] = np.where(ok)[0]

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

    cl["early_id"] = [pools[r_]["e_ids"][i] for r_, i in zip(cl.rep, cl.ei)]
    cl["later_id"] = [pools[r_]["l_ids"][i] for r_, i in zip(cl.rep, cl.li)]
    tr, va = cl[cl.part == "train"], cl[cl.part == "val"]
    # same call order as 07: Btr true targets (no draws), test true (no draws), then couplings
    for coupling in rq2.COUPLINGS:
        cl[f"target_{coupling}"] = None
        cl.loc[tr.index, f"target_{coupling}"] = targets(tr, coupling)
        cl.loc[va.index, f"target_{coupling}"] = targets(va, coupling)
    te = cl.part == "test"
    for coupling in rq2.COUPLINGS:
        cl.loc[te, f"target_{coupling}"] = cl.loc[te, "later_id"]
    return cl


def clone_table(pools, obs):
    clones = []
    for r, p in pools.items():
        ce = obs.loc[p["e_ids"], "clone"]
        cl = obs.loc[p["l_ids"], "clone"]
        both = set(ce.dropna()) & set(cl.dropna())
        for c in sorted(both):
            clones.append(dict(rep=r, clone=c,
                               early=np.where(ce.to_numpy() == c)[0],
                               late=np.where(cl.to_numpy() == c)[0]))
    return pd.DataFrame(clones)


def early_matrix(transition, pools, cl):
    """Raw counts for every early cell of the RQ2 clones (from the RQ3 early-count cache)."""
    obs3, _ = rq3.load(CELLSET)
    early3, _, _, _ = rq3.clone_table(obs3, transition)
    X3, _ = rq3.early_counts(early3, CELLSET, transition)
    ids = np.concatenate([[pools[r]["e_ids"][i] for i in e] for r, e in zip(cl.rep, cl.early)])
    pos = pd.Index(early3.cell_id).get_indexer(ids)
    assert (pos >= 0).all(), "RQ2 early cell missing from the RQ3 early-count cache"
    lib = early3.library.to_numpy()[pos]
    return X3[pos], lib, pd.Index(ids)


def run(transition):
    import torch
    torch.set_num_threads(1)
    obs = pd.read_csv(rq2.OBS, usecols=["cell_id", "stage", "clone"], low_memory=False).set_index("cell_id")
    sc_ = pd.read_csv(rq2.SCORES).set_index("cell_id")
    pools = {r: rq2.load_pool(CELLSET, transition, r, obs) for r in (1, 2, 3)}
    cl = clone_table(pools, obs)
    X, libs, eidx = early_matrix(transition, pools, cl)
    CACHE.mkdir(parents=True, exist_ok=True)
    DRAWS.mkdir(parents=True, exist_ok=True)
    tt = transition.replace("->", "-")
    ref = None
    for seed in rq2.SEEDS:
        ck = CACHE / f"{CELLSET}_{tt}_s{seed}.csv"
        if ck.exists():
            continue
        t0 = time.time()
        d = replay(seed, cl, pools, obs, sc_)
        d.drop(columns=["early", "late"]).to_csv(DRAWS / f"{CELLSET}_{tt}_s{seed}.csv", index=False)
        state = lambda ids: sc_.loc[list(ids), ["Pg", "D"]].to_numpy()
        tr, va, te = (d[d.part == p] for p in ("train", "val", "test"))
        Btr = np.vstack([state(tr.early_id), state(tr.target_true)])
        lo, hi = Btr.min(0), Btr.max(0)
        scl = lambda A: np.clip((A - lo) / np.where(hi > lo, hi - lo, 1), 0, 1)
        yte = scl(state(te.target_true))

        # reproduction check against the 07 cache (deterministic models, (Pg, D) input)
        old = pd.read_csv(rq2.CACHE / f"{CELLSET}_{tt}_s{seed}.csv")
        from sklearn.linear_model import Ridge
        for coupling in rq2.COUPLINGS:
            Xtr, ytr = scl(state(tr.early_id)), scl(state(tr[f"target_{coupling}"]))
            Xte = scl(state(te.early_id))
            chk = {"ConstMedian": np.tile(np.median(ytr, 0), (len(Xte), 1)),
                   "ConstMean": np.tile(ytr.mean(0), (len(Xte), 1)),
                   "Ridge": Ridge(alpha=1.0).fit(Xtr, ytr).predict(Xte)}
            for m, yh in chk.items():
                o = old[(old.coupling == coupling) & (old.model == m)].mae.item()
                assert abs(np.abs(yh - yte).mean() - o) < 1e-9, (transition, seed, coupling, m)

        # 30-PC input: fit on all early cells of training clones
        tr_cells = np.concatenate([[pools[r]["e_ids"][i] for i in e] for r, e in zip(tr.rep, tr.early)])
        tr_rows = eidx.get_indexer(tr_cells)
        PC = rq3.pca_train_apply(X, libs, tr_rows)
        pc = lambda ids: PC[eidx.get_indexer(list(ids))]
        set_input_size(PC.shape[1])
        res = []
        for coupling in rq2.COUPLINGS:
            ytr, yva = scl(state(tr[f"target_{coupling}"])), scl(state(va[f"target_{coupling}"]))
            preds, rfp = rq2.fit_models(pc(tr.early_id), ytr, pc(va.early_id), yva, seed)
            for m, f in preds.items():
                err = np.abs(f(pc(te.early_id)) - yte)
                res.append(dict(cellset=CELLSET, transition=transition, seed=seed, coupling=coupling,
                                input="PC30", model=m, mae=err.mean(), mae_Pg=err[:, 0].mean(),
                                mae_D=err[:, 1].mean(), n_train=len(tr), n_val=len(va), n_test=len(te),
                                n_pca_cells=len(tr_rows), hvg_batch_key=rq3.pca_train_apply.batch_hvg,
                                rf_params=str(rfp) if m == "RandomForest" else ""))
        new = pd.DataFrame(res)
        old = old.assign(input="Pg_D")[new.columns.drop(["n_pca_cells", "hvg_batch_key"])]
        pd.concat([old, new], ignore_index=True).to_csv(ck, index=False)
        print(transition, "seed", seed, "reproduced; PC30 done", f"{time.time() - t0:.0f}s", flush=True)


def summarise():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_s*.csv"))])
    df.to_csv(OUT / "per_seed.csv", index=False)
    rows = []
    for (tr, cp), g in df.groupby(["transition", "coupling"]):
        w = g.pivot_table(index="seed", columns=["input", "model"], values="mae")
        fam = []
        for m in rq2.MODELS:
            st = rq2.paired_stats(w[("PC30", m)], w[("Pg_D", m)])
            fam.append(dict(model=m, comparison="PC30 - Pg_D", **st,
                            mae_PC30=w[("PC30", m)].mean(), mae_Pg_D=w[("Pg_D", m)].mean(),
                            const_median_mae=w[("Pg_D", "ConstMedian")].mean()))
        for inp in ("PC30", "Pg_D"):
            for m in rq2.MODELS[1:]:
                st = rq2.paired_stats(w[(inp, m)], w[(inp, "ConstMedian")])
                fam.append(dict(model=m, comparison=f"{inp}: model - ConstMedian", **st,
                                const_median_mae=w[(inp, "ConstMedian")].mean()))
        adj = rq2.holm([f["p"] for f in fam])
        for f, pa in zip(fam, adj):
            f.update(transition=tr, coupling=cp, p_holm=pa, n_seeds=len(w),
                     delta_pct_of_const=100 * f["delta"] / f["const_median_mae"])
            rows.append(f)
    pd.DataFrame(rows).to_csv(OUT / "summary.csv", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    import sklearn
    import torch
    lines = [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
             "plan: ANALYSIS_PLAN.md amendment A3.1 (2026-11-17, post-review, exploratory)",
             f"python {sys.version.split()[0]} torch {torch.__version__} sklearn {sklearn.__version__} "
             f"numpy {np.__version__} pandas {pd.__version__}",
             "environment: Claude cloud workspace (Linux, 2 CPU, 5.8 GB)",
             f"input {rq2.SCORES.relative_to(ROOT)} sha256 {rq2.sha256(rq2.SCORES)}",
             "input (Pg, D) per-seed MAE: data/interim/rq2/primary_*_s*.csv (07_rq2), reproduced per seed",
             "input early counts: data/interim/rq3/primary_*_early_counts.npz (08_rq3)",
             "input OT potentials: data/interim/rq1/primary_*_bal_0.05.npz and *_emb.npz (RQ1)"]
    (OUT / "provenance.txt").write_text("\n".join(lines) + "\n")


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
