"""RQ3 — state adequacy (ANALYSIS_PLAN.md §3, approved at 65cd70d).

Does an early cell's state predict its clone's later fate?

  * Clones: all cross-stage clones of the cell set (all replicates; no OT pool or cap).
  * Label: majority cluster of the clone's later sister cells (ties -> alphabetical first).
    Fates with >= 30 clones per transition kept; the rest pooled as "other".
  * Representations of each early cell:
      (1) P_tcbb, D  (TCBB panel; reduces to SOX2 here)   (2) Pg, D
      (3) 30-PC      (HVG seurat_v3 on raw counts of TRAINING early cells, batch = library;
                      if that fails because a library has too few training cells, HVGs
                      without batch_key, recorded per seed in `hvg_batch_key`;
                      CP10k + log1p; scale with training mean/sd, max 10; PCA fitted on
                      training early cells; applied to validation/test cells)
      (4) fate-bias axis = LDA (K-1 dims) on the training 30-PC, cell-level fate labels.
  * Classifier, identical for all: StandardScaler (fitted on training cells) +
    multinomial logistic regression (L2, C = 1, lbfgs, max_iter 2000). Standardisation is
    an implementation detail fixed here so that only the representation differs.
  * Split: 70/15/15 by clone, stratified by replicate, 20 seeds. Validation is not used
    (nothing is tuned). A test clone's prediction = mean probability over its early cells.
  * Metrics: macro one-vs-rest AUROC over classes present in the test set (primary), and
    balanced accuracy. Chance: AUROC 0.5; balanced accuracy 1/K (K = classes in test).
  * Comparisons: each representation vs chance and every pair of representations, paired
    over seeds (bootstrap 95% CI, 10,000; Wilcoxon; Holm within cell set x transition x
    metric; Cohen's dz).
  * DA lineage (descriptive, III->IV only): clones whose early cells are mostly
    III-NPvMB; measured fate distribution; AUROC of clone-mean Pg, D and LDA component 1
    (LDA fitted on all III->IV clones, no split) for mDA-producing (any later IV-mDAs
    sister) vs other NPvMB clones. No inferential claim.

Usage: python scripts/08_rq3_fate.py run [--cellset ...] | summarise
"""
import argparse
import hashlib
import itertools
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
LIBS = ROOT / "data/processed/libs"
OBS = ROOT / "results/03_objects/cell_obs.csv.gz"
QC = ROOT / "results/04_qc/cell_qc.csv.gz"
SCORES = ROOT / "data/interim/cell_scores.csv.gz"
CACHE = ROOT / "data/interim/rq3"
OUT = ROOT / "results/08_rq3"
TRANSITIONS = {"I->II": ("I", "II", 1), "II->III": ("II", "III", 2), "III->IV": ("III", "IV", 3)}
SEEDS = list(range(20))
MIN_CLONES, N_HVG, N_PC = 30, 2000, 30
REPS = ["P_tcbb_D", "Pg_D", "PC30", "fate_axis"]


def load(cellset):
    obs = pd.read_csv(OBS, low_memory=False)
    qc = pd.read_csv(QC)
    obs = obs.merge(qc[["cell_id", "keep", "doublet_flag"]], on="cell_id")
    obs = obs[obs.keep & ~obs.doublet_flag] if cellset == "sensitivity" else obs[obs.keep]
    sc_ = pd.read_csv(SCORES).set_index("cell_id")
    return obs, sc_


def clone_table(obs, transition):
    s0, s1, e = TRANSITIONS[transition]
    d = obs[obs.in_clone_table & (obs.exp == e) & obs.stage.isin([s0, s1])].copy()
    d["clone"] = d.rep.astype(str) + "|" + d.clone.astype(str)
    st = d.groupby("clone").stage.agg(set)
    cross = st[st.apply(lambda x: {s0, s1} <= x)].index
    d = d[d.clone.isin(cross)]
    late = d[d.stage == s1]
    fate = late.groupby("clone").cluster.agg(lambda s: s.value_counts().sort_index().idxmax())
    vc = fate.value_counts()
    keep = set(vc[vc >= MIN_CLONES].index)
    fate = fate.where(fate.isin(keep), "other")
    early = d[d.stage == s0].copy()
    early["fate"] = early.clone.map(fate)
    rep = d.groupby("clone").rep.first()
    return early.reset_index(drop=True), fate, rep, late


def early_counts(early, cellset, transition):
    import anndata as ad
    ck = CACHE / f"{cellset}_{transition.replace('->', '-')}_early_counts.npz"
    if ck.exists():
        z = np.load(ck, allow_pickle=True)
        X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
        assert (z["ids"] == early.cell_id.to_numpy()).all()
        return X, z["genes"]
    parts, ids = [], []
    for lib, g in early.groupby("library", sort=True):
        a = ad.read_h5ad(LIBS / f"{lib}.h5ad")
        parts.append(a[g.cell_id.to_numpy()].X.tocsr())
        ids.append(g.cell_id.to_numpy())
        genes = a.var_names.to_numpy()
    X = sp.vstack(parts).tocsr()
    order = pd.Index(np.concatenate(ids)).get_indexer(early.cell_id.to_numpy())
    X = X[order]
    np.savez(ck, data=X.data, indices=X.indices, indptr=X.indptr, shape=X.shape,
             ids=early.cell_id.to_numpy(), genes=genes)
    return X, genes


def pca_train_apply(X, libs, tr_idx):
    import anndata as ad
    import scanpy as sc
    a = ad.AnnData(X=X[tr_idx].astype(np.float32))
    a.obs["library"] = libs[tr_idx]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=N_HVG, batch_key="library")
            pca_train_apply.batch_hvg = True
        except Exception:   # pre-specified fallback: too few training cells in a library
            sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=N_HVG)
            pca_train_apply.batch_hvg = False
    hv = a.var.highly_variable.to_numpy()

    def lognorm(M):
        M = M[:, :].astype(np.float32)
        tot = np.asarray(M.sum(1)).ravel()
        M = sp.diags(1e4 / np.maximum(tot, 1)) @ M
        return np.log1p(M[:, hv].toarray())
    Ltr = lognorm(X[tr_idx])
    mu, sd = Ltr.mean(0), Ltr.std(0)
    sd[sd == 0] = 1
    from sklearn.decomposition import PCA
    pca = PCA(n_components=N_PC, svd_solver="arpack", random_state=0)
    pca.fit(np.clip((Ltr - mu) / sd, None, 10))
    Lall = lognorm(X)
    return pca.transform(np.clip((Lall - mu) / sd, None, 10))


def macro_auroc(y, P, classes):
    from sklearn.metrics import roc_auc_score
    vals = []
    for k, c in enumerate(classes):
        pos = y == c
        if 0 < pos.sum() < len(y):
            vals.append(roc_auc_score(pos, P[:, k]))
    return float(np.mean(vals)) if vals else np.nan


def run(cellset):
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import balanced_accuracy_score
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    CACHE.mkdir(parents=True, exist_ok=True)
    obs, scores = load(cellset)
    for transition in TRANSITIONS:
        early, fate, rep, late = clone_table(obs, transition)
        X, genes = early_counts(early, cellset, transition)
        libs = early.library.to_numpy()
        S = scores.loc[early.cell_id, ["P_tcbb", "Pg", "D"]].to_numpy()
        clones = fate.index.to_numpy()
        for seed in SEEDS:
            ck = CACHE / f"{cellset}_{transition.replace('->', '-')}_s{seed}.csv"
            if ck.exists():
                continue
            t0 = time.time()
            rng = np.random.default_rng(seed)
            part = pd.Series(index=clones, dtype=object)
            for r in sorted(rep.unique()):
                idx = rng.permutation(rep[rep == r].index.to_numpy())
                n_tr, n_va = int(round(0.7 * len(idx))), int(round(0.15 * len(idx)))
                part[idx[:n_tr]], part[idx[n_tr:n_tr + n_va]], part[idx[n_tr + n_va:]] = "train", "val", "test"
            cp = early.clone.map(part).to_numpy()
            tr, te = np.where(cp == "train")[0], np.where(cp == "test")[0]
            y = early.fate.to_numpy()
            PC = pca_train_apply(X, libs, tr)
            lda = LinearDiscriminantAnalysis().fit(PC[tr], y[tr])
            feats = {"P_tcbb_D": S[:, [0, 2]], "Pg_D": S[:, [1, 2]], "PC30": PC, "fate_axis": lda.transform(PC)}
            res = []
            test_clones = early.clone.to_numpy()[te]
            yc = fate[pd.unique(test_clones)]
            for name, F in feats.items():
                clf = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000))
                clf.fit(F[tr], y[tr])
                P = pd.DataFrame(clf.predict_proba(F[te]), columns=clf.classes_)
                P["clone"] = test_clones
                Pc = P.groupby("clone").mean().loc[yc.index]
                classes = list(clf.classes_)
                ytrue = yc.to_numpy()
                res.append(dict(cellset=cellset, transition=transition, seed=seed, representation=name,
                                macro_auroc=macro_auroc(ytrue, Pc[classes].to_numpy(), classes),
                                balanced_acc=balanced_accuracy_score(ytrue, np.array(classes)[Pc[classes].to_numpy().argmax(1)]),
                                chance_bacc=1 / len(set(ytrue)), n_test_clones=len(yc),
                                n_train_cells=len(tr), n_classes=len(classes),
                                hvg_batch_key=pca_train_apply.batch_hvg))
            pd.DataFrame(res).to_csv(ck, index=False)
            print(cellset, transition, "seed", seed, f"{time.time() - t0:.0f}s", flush=True)
        if transition == "III->IV":
            da_descriptive(cellset, early, fate, late, X, libs, S)
        pd.DataFrame({"clone": fate.index, "fate": fate.values, "rep": rep[fate.index].values}).to_csv(
            CACHE / f"{cellset}_{transition.replace('->', '-')}_labels.csv", index=False)


def da_descriptive(cellset, early, fate, late, X, libs, S):
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    from sklearn.metrics import roc_auc_score
    PC = pca_train_apply(X, libs, np.arange(len(early)))
    l1 = LinearDiscriminantAnalysis().fit(PC, early.fate.to_numpy()).transform(PC)[:, 0]
    e = early.assign(Pg=S[:, 1], D=S[:, 2], lda1=l1)
    maj = e.groupby("clone").cluster.agg(lambda s: s.value_counts().sort_index().idxmax())
    nb = maj[maj == "III-NPvMB"].index
    late = late.assign(clone=late.rep.astype(str) + "|" + late.clone.astype(str)) if "|" not in str(late.clone.iloc[0]) else late
    mda = late.groupby("clone").cluster.agg(lambda s: (s == "IV-mDAs").any())
    cm = e[e.clone.isin(nb)].groupby("clone")[["Pg", "D", "lda1"]].mean()
    lab = mda.reindex(cm.index).fillna(False).astype(bool)
    rows = [dict(cellset=cellset, quantity="n_NPvMB_clones", value=len(cm)),
            dict(cellset=cellset, quantity="n_mDA_producing", value=int(lab.sum()))]
    for f in ["Pg", "D", "lda1"]:
        v = roc_auc_score(lab, cm[f]) if 0 < lab.sum() < len(lab) else np.nan
        rows.append(dict(cellset=cellset, quantity=f"auroc_{f}_mDA_vs_other", value=v))
    dist = late[late.clone.isin(nb)].cluster.value_counts(normalize=True)
    for k, v in dist.items():
        rows.append(dict(cellset=cellset, quantity=f"fate_share_{k}", value=v))
    pd.DataFrame(rows).to_csv(CACHE / f"{cellset}_III-IV_da_descriptive.csv", index=False)


def paired(a, b):
    from scipy.stats import wilcoxon
    d = np.asarray(a) - np.asarray(b)
    rng = np.random.default_rng(0)
    bt = rng.choice(d, size=(10000, len(d)), replace=True).mean(1)
    sd = d.std(ddof=1)
    return dict(delta=d.mean(), ci_lo=np.percentile(bt, 2.5), ci_hi=np.percentile(bt, 97.5),
                p=wilcoxon(d).pvalue if np.any(d != 0) else 1.0, dz=d.mean() / sd if sd > 0 else np.nan)


def holm(ps):
    ps = np.asarray(ps, float)
    adj, run_max = np.empty_like(ps), 0
    for k, i in enumerate(np.argsort(ps)):
        run_max = max(run_max, (len(ps) - k) * ps[i])
        adj[i] = min(1.0, run_max)
    return adj


def summarise():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.concat([pd.read_csv(p) for p in sorted(CACHE.glob("*_s*.csv"))])
    df.to_csv(OUT / "per_seed.csv", index=False)
    rows = []
    for (cs, tr), g in df.groupby(["cellset", "transition"]):
        for metric, chance_col in [("macro_auroc", None), ("balanced_acc", "chance_bacc")]:
            w = g.pivot(index="seed", columns="representation", values=metric)
            ch = 0.5 if chance_col is None else g.groupby("seed")[chance_col].first().loc[w.index]
            fam = [dict(a=r, b="chance", **paired(w[r], np.broadcast_to(ch, len(w)) if np.isscalar(ch) else ch.to_numpy()))
                   for r in REPS]
            fam += [dict(a=x, b=y, **paired(w[x], w[y])) for x, y in itertools.combinations(REPS, 2)]
            for f, pa in zip(fam, holm([f["p"] for f in fam])):
                f.update(cellset=cs, transition=tr, metric=metric, p_holm=pa,
                         mean_a=w[f["a"]].mean(), mean_b=(w[f["b"]].mean() if f["b"] != "chance" else np.mean(ch)))
                rows.append(f)
    pd.DataFrame(rows).to_csv(OUT / "comparisons.csv", index=False)
    df.groupby(["cellset", "transition", "representation"])[["macro_auroc", "balanced_acc", "chance_bacc"]] \
      .agg(["mean", "std"]).to_csv(OUT / "summary.csv")
    for p in CACHE.glob("*_labels.csv"):
        t = pd.read_csv(p)
        t.groupby("fate").size().rename("n_clones").to_csv(OUT / p.name.replace("_labels", "_class_counts"))
    da = list(CACHE.glob("*_da_descriptive.csv"))
    if da:
        pd.concat([pd.read_csv(p) for p in da]).to_csv(OUT / "da_descriptive.csv", index=False)
    write_provenance()


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_provenance():
    import scanpy as sc
    import sklearn
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    lines = [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
             "plan: ANALYSIS_PLAN.md approved at 65cd70d",
             f"python {sys.version.split()[0]} scanpy {sc.__version__} sklearn {sklearn.__version__} "
             f"numpy {np.__version__} pandas {pd.__version__}",
             "environment: Claude cloud workspace (Linux, 2 CPU, 7 GB)",
             f"input {SCORES.relative_to(ROOT)} sha256 {sha256(SCORES)}",
             f"input {OBS.relative_to(ROOT)} sha256 {sha256(OBS)}",
             f"input {QC.relative_to(ROOT)} sha256 {sha256(QC)}"]
    (OUT / "provenance.txt").write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "summarise"])
    ap.add_argument("--cellset", nargs="*", default=["primary", "sensitivity"])
    a = ap.parse_args()
    if a.cmd == "run":
        for cs in a.cellset:
            run(cs)
    else:
        summarise()


if __name__ == "__main__":
    main()
