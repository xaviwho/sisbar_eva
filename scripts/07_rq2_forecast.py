"""RQ2 — forecasting on true lineage (ANALYSIS_PLAN.md §2, approved at 65cd70d).

State: (Pg, D) = mean log1p(CP10k) over panel genes (Pg = SOX2, NES, HES1, HES5, VIM;
D = TH, DDC, SLC6A3, DRD2, LMX1A, FOXA2). Task: early-cell state -> later sister-cell state.

Per cell set x transition x seed (20 seeds), replicates pooled:
  * Split cross-stage clones 70/15/15, stratified by replicate.
  * One pair per clone: one early and one later cell drawn at random within the clone.
  * Training targets by coupling:
      true         sampled later sister cell;
      OT           later cell drawn from the early cell's RQ1 OT row (balanced, eps=0.05,
                   same cell set / transition / replicate pool), restricted to later cells
                   of TRAINING clones plus background later cells, renormalised;
      independent  uniform draw from the same restricted pool.
  * Validation pairs (model selection: Transformer early stopping) use validation clones'
    early cells with targets drawn by the SAME coupling from the same restricted pool
    (true coupling: their real sisters). No validation/test sister cell is ever a target
    of a non-true coupling.
  * Test pairs are always held-out true clonal pairs.
  * Min-max bounds per score from the early and later cells of training clones' true
    pairs (shared by all couplings within a seed so MAEs are comparable); clipped [0,1].
Models (TCBB settings, copied from xaviwho/leakage-safe-benchmark run_multiseed.py):
  constant median, constant mean, Ridge(alpha=1), RandomForest tuned by 3-fold CV over
  n_estimators=300, max_depth {2,3,5,None}, min_samples_leaf {1,5,10} (MAE scoring),
  Transformer (d_model 128, 8 heads, 4 layers, ff 512, dropout 0.1, GELU, Adam 1e-3,
  wd 1e-5, ReduceLROnPlateau(0.5, 10), MSE, batch 32, 100 epochs, best validation).
  With one input token the Transformer reduces to a feed-forward network (stated in paper).
Statistics per cell set x transition x coupling: paired delta vs constant median and vs
random forest over 20 seeds; percentile paired bootstrap 95% CI (10,000); Wilcoxon,
Holm-corrected within each (cell set, transition, coupling) family; Cohen's dz;
equivalence iff the whole CI lies within +/- 0.01; delta also as % of constant-median MAE.

Usage: python scripts/07_rq2_forecast.py scores | run [--cellset ...] [--transition ...] | summarise
"""
import argparse
import copy
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LIBS = ROOT / "data/processed/libs"
OBS = ROOT / "results/03_objects/cell_obs.csv.gz"
QC = ROOT / "results/04_qc/cell_qc.csv.gz"
RQ1 = ROOT / "data/interim/rq1"
CACHE = ROOT / "data/interim/rq2"
OUT = ROOT / "results/07_rq2"
SCORES = ROOT / "data/interim/cell_scores.csv.gz"

PANELS = {"Pg": ["SOX2", "NES", "HES1", "HES5", "VIM"],
          "D": ["TH", "DDC", "SLC6A3", "DRD2", "LMX1A", "FOXA2"],
          "P_tcbb": ["POU5F1", "NANOG", "SOX2", "UTF1", "TDGF1"]}
TRANSITIONS = {"I->II": ("I", "II", 1), "II->III": ("II", "III", 2), "III->IV": ("III", "IV", 3)}
EPS_OT = 0.05
SEEDS = list(range(20))
COUPLINGS = ["true", "OT", "independent"]
MODELS = ["ConstMedian", "ConstMean", "Ridge", "RandomForest", "Transformer"]
DELTA = 0.01
RF_GRID = {"n_estimators": [300], "max_depth": [2, 3, 5, None], "min_samples_leaf": [1, 5, 10]}


# --------------------------------------------------------------------------- scores
def compute_scores():
    """Per-cell raw panel scores for every primary cell (mean log1p CP10k over panel)."""
    import anndata as ad
    obs = pd.read_csv(OBS, usecols=["cell_id", "library"])
    keep = pd.read_csv(QC, usecols=["cell_id", "keep"])
    ids = set(keep.loc[keep.keep, "cell_id"])
    out = []
    for lib in sorted(obs.library.unique()):
        a = ad.read_h5ad(LIBS / f"{lib}.h5ad")
        a = a[[c in ids for c in a.obs_names]].copy()
        tot = np.asarray(a.X.sum(1)).ravel()
        sym = pd.Index(a.var.symbol)
        rec = {"cell_id": a.obs_names}
        for name, genes in PANELS.items():
            cols = sym.get_indexer(genes)
            assert (cols >= 0).all(), (name, genes)
            X = a.X[:, cols].toarray() / tot[:, None] * 1e4
            rec[name] = np.log1p(X).mean(1)
        out.append(pd.DataFrame(rec))
        print("scored", lib, flush=True)
    pd.concat(out).to_csv(SCORES, index=False)


# --------------------------------------------------------------------------- pools
def load_pool(cellset, transition, rep, obs):
    tag = f"{cellset}_{transition.replace('->', '-')}_r{rep}"
    z = np.load(RQ1 / f"{tag}_emb.npz", allow_pickle=True)
    pot = np.load(RQ1 / f"{tag}_bal_{EPS_OT}.npz")
    assert bool(pot["done"]) and bool(pot["converged"]), tag
    s0, s1, _ = TRANSITIONS[transition]
    ids = z["ids"]
    st = obs.loc[ids, "stage"].to_numpy()
    e_ids, l_ids = ids[st == s0], ids[st == s1]
    X0, X1 = z["X"][st == s0], z["X"][st == s1]
    from scipy.spatial.distance import cdist
    med = np.float32(np.median(cdist(X0, X1, "sqeuclidean").astype(np.float32)))
    return dict(e_ids=e_ids, l_ids=l_ids, X0=X0, X1=X1, med=med,
                lu=pot["log_u"], lv=pot["log_v"], rep=rep)


def ot_row(pool, i):
    from scipy.spatial.distance import cdist
    m = cdist(pool["X0"][i:i + 1], pool["X1"], "sqeuclidean").astype(np.float32)[0] / pool["med"]
    logw = pool["lu"][i] - m / EPS_OT + pool["lv"]
    w = np.exp(logw - logw.max())
    return w


# --------------------------------------------------------------------------- models
def make_transformer():
    import torch
    import torch.nn as nn

    class PositionalEncoding(nn.Module):
        def __init__(self, d_model, max_len=1000, dropout=0.1):
            super().__init__()
            self.dropout = nn.Dropout(p=dropout)
            position = torch.arange(max_len).unsqueeze(1)
            div_term = torch.exp(torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model))
            pe = torch.zeros(max_len, 1, d_model)
            pe[:, 0, 0::2] = torch.sin(position * div_term)
            pe[:, 0, 1::2] = torch.cos(position * div_term)
            self.register_buffer("pe", pe)

        def forward(self, x):
            return self.dropout(x + self.pe[:x.size(0)])

    class TransformerPredictor(nn.Module):   # TCBB src/models/predictors/transformer_predictor.py
        def __init__(self, input_size=2, output_size=2, d_model=128, nhead=8, num_layers=4,
                     dim_feedforward=512, dropout=0.1):
            super().__init__()
            self.input_projection = nn.Linear(input_size, d_model)
            self.pos_encoder = PositionalEncoding(d_model, dropout=dropout)
            layer = nn.TransformerEncoderLayer(d_model, nhead, dim_feedforward, dropout,
                                               activation="gelu", batch_first=True)
            self.transformer_encoder = nn.TransformerEncoder(layer, num_layers=num_layers)
            self.output_projection = nn.Sequential(nn.Linear(d_model, dim_feedforward // 2), nn.GELU(),
                                                   nn.Dropout(dropout), nn.Linear(dim_feedforward // 2, output_size))

        def forward(self, x):
            x = self.input_projection(x).transpose(0, 1)
            x = self.pos_encoder(x).transpose(0, 1)
            return self.output_projection(self.transformer_encoder(x))

    return TransformerPredictor


def fit_transformer(Xtr, ytr, Xva, yva, seed, epochs=100):
    import torch
    import torch.nn as nn
    torch.manual_seed(seed)
    model = make_transformer()()
    t = lambda a: torch.tensor(a, dtype=torch.float32)
    Xtr, ytr, Xva, yva = t(Xtr)[:, None, :], t(ytr)[:, None, :], t(Xva)[:, None, :], t(yva)[:, None, :]
    opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, mode="min", factor=0.5, patience=10)
    crit = nn.MSELoss()
    g = torch.Generator().manual_seed(seed)
    best_val, best_state = float("inf"), None
    for _ in range(epochs):
        model.train()
        perm = torch.randperm(len(Xtr), generator=g)
        for i in range(0, len(Xtr), 32):
            b = perm[i:i + 32]
            opt.zero_grad()
            loss = crit(model(Xtr[b])[:, -1:, :], ytr[b])
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
        model.eval()
        with torch.no_grad():
            vloss = crit(model(Xva)[:, -1:, :], yva).item()
        sched.step(vloss)
        if vloss < best_val:
            best_val, best_state = vloss, copy.deepcopy(model.state_dict())
    model.load_state_dict(best_state)
    model.eval()

    def predict(X):
        with torch.no_grad():
            return model(t(X)[:, None, :])[:, -1, :].numpy()
    return predict


def fit_models(Xtr, ytr, Xva, yva, seed):
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import GridSearchCV
    med, mean = np.median(ytr, 0), ytr.mean(0)
    ridge = Ridge(alpha=1.0).fit(Xtr, ytr)
    rf = GridSearchCV(RandomForestRegressor(random_state=seed, n_jobs=1), RF_GRID,
                      scoring="neg_mean_absolute_error", cv=3).fit(Xtr, ytr)
    tr = fit_transformer(Xtr, ytr, Xva, yva, seed)
    return {"ConstMedian": lambda X: np.tile(med, (len(X), 1)),
            "ConstMean": lambda X: np.tile(mean, (len(X), 1)),
            "Ridge": ridge.predict, "RandomForest": rf.predict, "Transformer": tr}, rf.best_params_


# --------------------------------------------------------------------------- run
def run(cellset, transition):
    import torch
    torch.set_num_threads(1)
    obs = pd.read_csv(OBS, low_memory=False).set_index("cell_id")
    sc_ = pd.read_csv(SCORES).set_index("cell_id")
    s0, s1, e = TRANSITIONS[transition]
    pools = {r: load_pool(cellset, transition, r, obs) for r in (1, 2, 3)}
    # clones per replicate (foreground of the RQ1 pool)
    clones = []
    for r, p in pools.items():
        ce = obs.loc[p["e_ids"], "clone"]
        cl = obs.loc[p["l_ids"], "clone"]
        both = set(ce.dropna()) & set(cl.dropna())
        for c in sorted(both):
            clones.append(dict(rep=r, clone=c,
                               early=np.where(ce.to_numpy() == c)[0],
                               late=np.where(cl.to_numpy() == c)[0]))
    cl = pd.DataFrame(clones)
    CACHE.mkdir(parents=True, exist_ok=True)
    rows = []
    for seed in SEEDS:
        ck = CACHE / f"{cellset}_{transition.replace('->', '-')}_s{seed}.csv"
        if ck.exists():
            rows.append(pd.read_csv(ck))
            continue
        t0 = time.time()
        rng = np.random.default_rng(seed)
        part = np.empty(len(cl), dtype=object)
        for r, g in cl.groupby("rep"):
            idx = rng.permutation(g.index.to_numpy())
            n_tr, n_va = int(round(0.7 * len(idx))), int(round(0.15 * len(idx)))
            part[idx[:n_tr]], part[idx[n_tr:n_tr + n_va]], part[idx[n_tr + n_va:]] = "train", "val", "test"
        cl["part"] = part
        cl["ei"] = [rng.choice(x) for x in cl.early]
        cl["li"] = [rng.choice(x) for x in cl.late]
        # restricted target pools: later cells of training clones + background later cells
        allowed = {}
        for r, p in pools.items():
            lc = obs.loc[p["l_ids"], "clone"].to_numpy()
            fg_all = set(cl[cl.rep == r].clone)
            fg_train = set(cl[(cl.rep == r) & (cl.part == "train")].clone)
            ok = np.array([(c not in fg_all) or (c in fg_train) for c in lc])
            allowed[r] = np.where(ok)[0]
        state = lambda ids: sc_.loc[ids, ["Pg", "D"]].to_numpy()

        def targets(sub, coupling):
            out = []
            for r_, ei, li in zip(sub.rep, sub.ei, sub.li):
                p = pools[r_]
                if coupling == "true":
                    j = li
                elif coupling == "OT":
                    w = ot_row(p, ei)[allowed[r_]].astype(np.float64)
                    j = allowed[r_][rng.choice(len(w), p=w / w.sum())]
                else:
                    j = rng.choice(allowed[r_])
                out.append(p["l_ids"][j])
            return np.array(out)

        def early_ids(sub):
            return np.array([pools[r_]["e_ids"][i] for r_, i in zip(sub.rep, sub.ei)])

        tr, va, te = cl[cl.part == "train"], cl[cl.part == "val"], cl[cl.part == "test"]
        # min-max bounds from training clones' true pairs (shared across couplings)
        Btr = np.vstack([state(early_ids(tr)), state(targets(tr, "true"))])
        lo, hi = Btr.min(0), Btr.max(0)
        sc = lambda A: np.clip((A - lo) / np.where(hi > lo, hi - lo, 1), 0, 1)
        Xte, yte = sc(state(early_ids(te))), sc(state(targets(te, "true")))
        res = []
        for coupling in COUPLINGS:
            Xtr, ytr = sc(state(early_ids(tr))), sc(state(targets(tr, coupling)))
            Xva, yva = sc(state(early_ids(va))), sc(state(targets(va, coupling)))
            preds, rfp = fit_models(Xtr, ytr, Xva, yva, seed)
            for m, f in preds.items():
                yh = f(Xte)
                err = np.abs(yh - yte)
                res.append(dict(cellset=cellset, transition=transition, seed=seed, coupling=coupling,
                                model=m, mae=err.mean(), mae_Pg=err[:, 0].mean(), mae_D=err[:, 1].mean(),
                                n_train=len(tr), n_val=len(va), n_test=len(te),
                                rf_params=json.dumps(rfp) if m == "RandomForest" else ""))
        df = pd.DataFrame(res)
        df.to_csv(ck, index=False)
        rows.append(df)
        print(cellset, transition, "seed", seed, f"{time.time() - t0:.0f}s", flush=True)
    return pd.concat(rows)


# --------------------------------------------------------------------------- stats
def paired_stats(a, b, seed=0):
    from scipy.stats import wilcoxon
    d = np.asarray(a) - np.asarray(b)
    rng = np.random.default_rng(seed)
    boots = rng.choice(d, size=(10000, len(d)), replace=True).mean(1)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0
    sd = d.std(ddof=1)
    return dict(delta=d.mean(), ci_lo=lo, ci_hi=hi, p=p, dz=d.mean() / sd if sd > 0 else np.nan,
                equivalent=bool(-DELTA < lo and hi < DELTA))


def holm(ps):
    ps = np.asarray(ps, float)
    order = np.argsort(ps)
    adj = np.empty_like(ps)
    run_max = 0
    for k, i in enumerate(order):
        run_max = max(run_max, (len(ps) - k) * ps[i])
        adj[i] = min(1.0, run_max)
    return adj


def summarise():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_s*.csv"))])
    df.to_csv(OUT / "per_seed.csv", index=False)
    rows = []
    for (cs, tr, cp), g in df.groupby(["cellset", "transition", "coupling"]):
        w = g.pivot(index="seed", columns="model", values="mae")
        fam = []
        for m in MODELS:
            if m != "ConstMedian":
                fam.append(dict(model=m, reference="ConstMedian", **paired_stats(w[m], w["ConstMedian"])))
            if m not in ("RandomForest",):
                fam.append(dict(model=m, reference="RandomForest", **paired_stats(w[m], w["RandomForest"])))
        adj = holm([f["p"] for f in fam])
        for f, pa in zip(fam, adj):
            f.update(cellset=cs, transition=tr, coupling=cp, p_holm=pa, n_seeds=len(w),
                     mae_model=w[f["model"]].mean(), mae_reference=w[f["reference"]].mean(),
                     const_median_mae=w["ConstMedian"].mean(),
                     delta_pct_of_const=100 * f["delta"] / w["ConstMedian"].mean())
            rows.append(f)
    pd.DataFrame(rows).to_csv(OUT / "comparisons.csv", index=False)
    agg = df.groupby(["cellset", "transition", "coupling", "model"]).mae.agg(["mean", "std", "count"]).reset_index()
    agg.to_csv(OUT / "mae_summary.csv", index=False)
    write_provenance()


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_provenance():
    import scipy
    import sklearn
    import torch
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    lines = [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
             "plan: ANALYSIS_PLAN.md approved at 65cd70d",
             f"python {sys.version.split()[0]} torch {torch.__version__} sklearn {sklearn.__version__} "
             f"numpy {np.__version__} scipy {scipy.__version__} pandas {pd.__version__}",
             "environment: Claude cloud workspace (Linux, 2 CPU, 7 GB)",
             f"input {SCORES.relative_to(ROOT)} sha256 {sha256(SCORES)}",
             "input OT potentials: data/interim/rq1/*_bal_0.05.npz and *_emb.npz (RQ1)"]
    (OUT / "provenance.txt").write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["scores", "run", "summarise"])
    ap.add_argument("--cellset", nargs="*", default=["primary", "sensitivity"])
    ap.add_argument("--transition", nargs="*", default=list(TRANSITIONS))
    a = ap.parse_args()
    if a.cmd == "scores":
        compute_scores()
    elif a.cmd == "run":
        for cs in a.cellset:
            for tr in a.transition:
                run(cs, tr)
    else:
        summarise()


if __name__ == "__main__":
    main()
