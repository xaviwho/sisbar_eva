"""Descriptive checks to back the review of ANALYSIS_PLAN.md (UPDATES.md Q5).

None of these touch an RQ outcome: no coupling is fitted, no model is trained, and no
early-cell feature is related to a later fate. They describe marker expression, score
spread, cell/clone counts and analytic chance levels, to inform settings that must be
fixed before the plan is approved.

A  Marker panels by stage (primary cells): detection rate and mean log1p(CP10k) for the
   TCBB pluripotency panel P, the proposed progenitor panel Pg, and the D panel.
B  Score spread per transition (for the equivalence margin): panel scores as in TCBB
   (mean log1p CP10k over genes), min-max scaled over primary cells of both stages of
   the transition (descriptive stand-in for training-only scaling); on later-stage
   cells of cross-stage clones: SD, IQR, and MAE of the constant median.
C  RQ2 leakage exposure: share of later-stage cells (after the 8,000 cap) that belong
   to cross-stage clones, i.e. potential OT draws from validation/test clones.
D  Effect of the 8,000 cap on chance own-clone share m0 and chance top-k recall.
E  Analytic chance for the independent coupling (AUROC with ties as 0.5; top-k).
F  RQ3 fate classes (>= 30 clones) and the number of III->IV clones ending as mDA.

Panels: P (TCBB dopaminergic) = POU5F1 NANOG SOX2 UTF1 TDGF1; Pg (plan) = SOX2 NES HES1
HES5 VIM; D = TH DDC SLC6A3 DRD2 LMX1A FOXA2. Library size excludes EGFP.
"""
import gzip
import hashlib
import io
import subprocess
import sys
import tarfile
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import scipy.io

ROOT = Path(__file__).resolve().parents[1]
TARS = [ROOT / "data/raw/GSE207921/GSE207921_RAW.tar",
        ROOT / "data/raw/GSE221592/GSE221592_RAW.tar"]
OBS = ROOT / "results/03_objects/cell_obs.csv.gz"
QC = ROOT / "results/04_qc/cell_qc.csv.gz"
OUT = ROOT / "results/05a_plan_checks"
OUT.mkdir(parents=True, exist_ok=True)

P = ["POU5F1", "NANOG", "SOX2", "UTF1", "TDGF1"]
PG = ["SOX2", "NES", "HES1", "HES5", "VIM"]
D = ["TH", "DDC", "SLC6A3", "DRD2", "LMX1A", "FOXA2"]
GENES = sorted(set(P + PG + D))
PAIRS = {1: ("I", "II"), 2: ("II", "III"), 3: ("III", "IV")}
CAP = 8000
KS = (1, 5, 10)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_markers():
    """Per-cell log1p CP10k for GENES, from every library in both tars."""
    frames = []
    for t in TARS:
        with tarfile.open(t) as tf:
            names = tf.getnames()
            libs = sorted({n.rsplit("_", 1)[0] for n in names if n.endswith("_matrix.mtx.gz")})
            for lib in libs:
                rd = lambda suf: gzip.GzipFile(fileobj=io.BytesIO(tf.extractfile(f"{lib}_{suf}").read()))
                feat_name = next(n for n in names if n.startswith(lib + "_") and
                                 (n.endswith("features.tsv.gz") or n.endswith("genes.tsv.gz")))
                feats = pd.read_csv(rd(feat_name[len(lib) + 1:]), sep="\t", header=None)
                bcs = pd.read_csv(rd("barcodes.tsv.gz"), sep="\t", header=None)[0]
                X = scipy.io.mmread(rd("matrix.mtx.gz")).tocsr()  # genes x cells
                sym = feats[1].astype(str).values
                keep_rows = sym != "EGFP"
                libsize = np.asarray(X[keep_rows].sum(axis=0)).ravel()
                idx = {g: np.where(sym == g)[0] for g in GENES}
                assert all(len(v) == 1 for v in idx.values()), {g: len(v) for g, v in idx.items()}
                M = X[[idx[g][0] for g in GENES]].toarray().T
                with np.errstate(divide="ignore", invalid="ignore"):
                    L = np.log1p(M / libsize[:, None] * 1e4)
                df = pd.DataFrame(L, columns=GENES)
                df["cell_id"] = lib + ":" + bcs.str.split("-").str[0]
                frames.append(df)
                print("read", lib, X.shape, flush=True)
                del X
    return pd.concat(frames, ignore_index=True)


def cross_stage_sets(c):
    out = {}
    for exp, (e, l) in PAIRS.items():
        x = c[(c.exp == exp) & c.in_clone_table]
        st = x.groupby(["rep", "clone"])["stage"].agg(set)
        out[exp] = set(st[st.apply(lambda s: {e, l} <= s)].index)
    return out


def main():
    obs = pd.read_csv(OBS, low_memory=False)
    qc = pd.read_csv(QC)
    c = obs.merge(qc[["cell_id", "keep", "doublet_flag"]], on="cell_id", validate="one_to_one")
    c = c[c.keep].copy()
    c["in_clone_table"] = c.in_clone_table.astype(bool)
    mk = read_markers()
    c = c.merge(mk, on="cell_id", how="left", validate="one_to_one")
    assert c[GENES].notna().all().all(), "primary cell missing from matrices"
    c["P"], c["Pg"], c["D"] = c[P].mean(1), c[PG].mean(1), c[D].mean(1)
    c["P_noSOX2"] = c[[g for g in P if g != "SOX2"]].mean(1)

    # ---- A: panels by stage
    rows = []
    for st, g in c.groupby("stage"):
        for gene in GENES:
            rows.append(dict(stage=st, gene=gene,
                             panel="/".join(n for n, pnl in (("P", P), ("Pg", PG), ("D", D)) if gene in pnl),
                             detect_pct=100 * (g[gene] > 0).mean(), mean_log1p=g[gene].mean()))
    A = pd.DataFrame(rows)
    A.to_csv(OUT / "A_marker_detection_by_stage.csv", index=False)
    Ascore = c.groupby("stage")[["P", "P_noSOX2", "Pg", "D"]].agg(["mean", "std"])
    Ascore.to_csv(OUT / "A_panel_scores_by_stage.csv")

    # ---- B: spread of scaled scores on later cells of cross-stage clones
    cs = cross_stage_sets(c)
    key = list(zip(c.rep, c.clone))
    rows = []
    for exp, (e, l) in PAIRS.items():
        t = c[(c.exp == exp) & c.stage.isin([e, l])]
        later = t[(t.stage == l) & pd.Series([k in cs[exp] for k in zip(t.rep, t.clone)], index=t.index)]
        for s in ("P", "Pg", "D"):
            lo, hi = t[s].min(), t[s].max()
            y = ((later[s] - lo) / (hi - lo)).clip(0, 1)
            rows.append(dict(transition=f"{e}->{l}", score=s, n_later=len(y),
                             sd=y.std(), iqr=y.quantile(.75) - y.quantile(.25),
                             const_median_mae=(y - y.median()).abs().mean(),
                             pct_zero=100 * (later[s] == 0).mean()))
    B = pd.DataFrame(rows)
    pair_mae = (B[B.score.isin(["Pg", "D"])].groupby("transition")["const_median_mae"].mean()
                .rename("const_median_mae_PgD_avg"))
    B = B.merge(pair_mae, on="transition")
    B["delta_pct_of_const_mae_PgD"] = 100 * 0.01 / B.const_median_mae_PgD_avg
    B.to_csv(OUT / "B_score_spread.csv", index=False)

    # ---- C/D/E: per transition x replicate, later-stage pool with and without the cap
    rng = np.random.default_rng(0)
    rowsC, rowsD = [], []
    for exp, (e, l) in PAIRS.items():
        for rep in sorted(c.rep.unique()):
            t = c[(c.exp == exp) & (c.rep == rep)]
            inx = pd.Series([k in cs[exp] for k in zip(t.rep, t.clone)], index=t.index)
            later_all = t[t.stage == l]
            in_l = inx[later_all.index]
            other = later_all[~in_l]
            keep_other = other if len(other) <= CAP else other.iloc[rng.choice(len(other), CAP, replace=False)]
            later_cap = pd.concat([later_all[in_l], keep_other])
            rowsC.append(dict(transition=f"{e}->{l}", rep=rep, later_cells_capped=len(later_cap),
                              later_cells_uncapped=len(later_all),
                              crossstage_later_cells=int(in_l.sum()),
                              pct_crossstage_of_capped=100 * in_l.sum() / len(later_cap),
                              barcoded_not_crossstage_later=int((later_all.in_clone_table & ~in_l).sum())))
            # chance per early barcoded cell of a cross-stage clone
            early = t[(t.stage == e) & inx]
            sizes = later_all[in_l].groupby("clone").size()
            for name, n in (("capped", len(later_cap)), ("uncapped", len(later_all))):
                s = early.clone.map(sizes).values
                r = dict(transition=f"{e}->{l}", rep=rep, pool=name, n_later=n,
                         n_early=len(s), m0_mean=float(np.mean(s / n)))
                for k in KS:
                    r[f"chance_top{k}"] = float(np.mean([1 - comb(n - si, k) / comb(n, k) for si in s]))
                rowsD.append(r)
    C = pd.DataFrame(rowsC)
    C.to_csv(OUT / "C_leakage_exposure.csv", index=False)
    Dd = pd.DataFrame(rowsD)
    Dd.to_csv(OUT / "D_chance_vs_cap.csv", index=False)

    # ---- F: fate classes per transition (later-sister majority cluster)
    rowsF = []
    for exp, (e, l) in PAIRS.items():
        t = c[(c.exp == exp) & (c.stage == l) & c.in_clone_table]
        t = t[[k in cs[exp] for k in zip(t.rep, t.clone)]]
        fate = t.groupby(["rep", "clone"])["cluster"].agg(lambda x: x.value_counts().index[0])
        vc = fate.str.replace(r"^[IV]+-", "", regex=True).value_counts()
        for f, nn in vc.items():
            rowsF.append(dict(transition=f"{e}->{l}", fate=f, clones=int(nn), class_ge30=nn >= 30))
    F = pd.DataFrame(rowsF)
    F.to_csv(OUT / "F_fate_classes.csv", index=False)

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    prov = [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
            f"python {sys.version.split()[0]} pandas {pd.__version__} numpy {np.__version__} scipy {scipy.__version__}",
            "environment: Windows (scipy.io.mmread only; no scanpy)"]
    prov += [f"input {p.relative_to(ROOT).as_posix()} sha256 {sha256(p)}" for p in TARS + [OBS, QC]]
    (OUT / "provenance.txt").write_text("\n".join(prov) + "\n")

    pd.set_option("display.width", 220)
    for name, df in (("A detection", A.pivot_table(index=["panel", "gene"], columns="stage", values="detect_pct").round(1)),
                     ("A scores", Ascore.round(3)), ("B", B.round(4)), ("C", C.round(2)),
                     ("D", Dd.round(4)), ("F", F)):
        print(f"\n=== {name}\n{df.to_string()}")


if __name__ == "__main__":
    main()
