"""Results numbers, figures and tables for paper/main.tex (RQ1-RQ4).

Reads only committed result files in results/06_rq1 ... results/09_rq4 (final at
35d89c2) and writes:
  paper/generated/results_numbers.tex   LaTeX macros used in the Results text
  paper/tables/tab_rq2_mae.tex          Table 2 (RQ2 MAE, all couplings and models)
  paper/tables/tab_s*.tex               supplementary tables (sensitivity, cap, replicates)
  paper/figures/fig{2..6}_*.pdf         Figs 2-6
  results/10_results_paper/provenance.txt
No number in the Results is typed by hand; everything below is read from the CSVs.
Primary rows go to the main text; sensitivity rows (cellset == "sensitivity") and the
uncapped I->II replicate-2 rows (uncapped == True) go to supplementary tables.
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"
PAPER = ROOT / "paper"
GEN, TAB, FIG = PAPER / "generated", PAPER / "tables", PAPER / "figures"
OUT = R / "10_results_paper"
for d in (GEN, TAB, FIG, OUT):
    d.mkdir(parents=True, exist_ok=True)

INPUTS = [R / p for p in (
    "06_rq1/summary.csv", "06_rq1/fate_agreement.csv", "06_rq1/solves.csv",
    "06_rq1/blocked_solver_validation.csv",
    "07_rq2/mae_summary.csv", "07_rq2/comparisons.csv", "07_rq2/per_seed.csv",
    "08_rq3/per_seed.csv", "08_rq3/comparisons.csv", "08_rq3/da_descriptive.csv",
    "08_rq3/primary_I-II_class_counts.csv", "08_rq3/primary_II-III_class_counts.csv",
    "08_rq3/primary_III-IV_class_counts.csv",
    "09_rq4/grid_summary.csv", "09_rq4/measured_summary.csv")]

TR = ["I->II", "II->III", "III->IV"]
TN = dict(zip(TR, ["One", "Two", "Three"]))          # macro suffixes
TLAB = {t: t.replace("->", r"$\to$") for t in TR}    # LaTeX labels
TPLOT = {t: t.replace("->", "→") for t in TR}   # figure labels
OT = "OT_bal_eps0.05"
EPS = [0.005, 0.01, 0.05, 0.1]
DELTA = 0.01  # equivalence margin, ANALYSIS_PLAN §2

# Figure style (dataviz reference palette, light; validated: scripts/validate_palette.js)
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, GRID, GREY = "#0b0b0b", "#52514e", "#e4e3df", "#8a8984"
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 7.5, "axes.titlesize": 8,
    "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5, "axes.axisbelow": True,
    "legend.frameon": False, "pdf.fonttype": 42, "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02})
COL1, COL2 = 3.45, 7.0  # IEEE column / page widths (in)

M = {}  # macro name -> LaTeX text


def minus(s):
    """Typeset a leading minus as a true minus sign in text, math and tables."""
    return r"\ensuremath{-}" + s[1:] if s.startswith("-") else s


def f(x, d=2):
    s = f"{x:.{d}f}"
    if float(s) == 0:
        s = s.lstrip("-")
    return minus(s)


def pct(x, d=1):
    return f(100 * x, d)


def put(name, val):
    assert name.isalpha(), name
    assert name not in M, f"duplicate macro {name}"
    M[name] = str(val)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def one(df, **kw):
    q = df
    for k, v in kw.items():
        q = q[q[k] == v]
    assert len(q) == 1, (kw, len(q))
    return q.iloc[0]


# ============================================================ RQ1
def rq1():
    s = pd.read_csv(R / "06_rq1/summary.csv")
    s["rep"] = s.rep.astype(str)
    sv = pd.read_csv(R / "06_rq1/solves.csv")
    val = pd.read_csv(R / "06_rq1/blocked_solver_validation.csv").iloc[0]
    base = dict(cellset="primary", uncapped=False, rep="pooled")
    fw = s[(s.direction == "forward")]
    bw = s[(s.direction == "backward")]

    put("nOTSolves", len(sv))
    put("nOTSolvesConverged", int(sv.converged.sum()))
    put("nRQOneJobs", sv.groupby(["cellset", "transition", "rep", "uncapped"]).ngroups)
    put("valMaxRelPlanDiff", f"{val.max_rel_plan_diff:.0e}".replace("e-0", r"\times10^{-") + "}")
    put("valMeanAbsDiffAUROC", f(val.mean_abs_diff_auroc, 3))

    aur_bal, aur_unb = [], []
    for t in TR:
        n = TN[t]
        o = one(fw, coupling=OT, transition=t, **base)
        ind = one(fw, coupling="independent", transition=t, **base)
        knn = one(fw, coupling="kNN", transition=t, **base)
        put(f"rqOneClones{n}", int(o.n_clones))
        put(f"rqOneCells{n}", int(o.n_cells))
        put(f"rqOneEnr{n}", f(o.enrichment, 1))
        put(f"rqOneEnrLo{n}", f(o.enrichment_lo, 1))
        put(f"rqOneEnrHi{n}", f(o.enrichment_hi, 1))
        put(f"rqOneAUROC{n}", f(o.auroc))
        put(f"rqOneAUROCLo{n}", f(o.auroc_lo))
        put(f"rqOneAUROCHi{n}", f(o.auroc_hi))
        put(f"rqOneMassPct{n}", pct(o.m, 2))
        put(f"rqOneChancePct{n}", pct(o.m0, 3))
        put(f"rqOneTopOnePct{n}", pct(o.top1, 2))
        put(f"rqOneTopTenPct{n}", pct(o.top10, 1))
        put(f"rqOneTopTenChancePct{n}", pct(ind.top10, 2))
        put(f"rqOneKnnAUROC{n}", f(knn.auroc, 3))
        put(f"rqOneKnnEnr{n}", f(knn.enrichment, 1))
        put(f"rqOneKnnEnrLo{n}", f(knn.enrichment_lo, 1))
        put(f"rqOneKnnEnrHi{n}", f(knn.enrichment_hi, 1))
        put(f"rqOneKnnTopTenPct{n}", pct(knn.top10, 1))
        ob = one(bw, coupling=OT, transition=t, **base)
        put(f"rqOneBwdAUROC{n}", f(ob.auroc))
        put(f"rqOneBwdEnr{n}", f(ob.enrichment, 1))
        for e in EPS:
            aur_bal.append(one(fw, coupling=f"OT_bal_eps{e}", transition=t, **base).auroc)
            aur_unb.append(one(fw, coupling=f"OT_unb_eps{e}", transition=t, **base).auroc)
        unb3 = fw[(fw.transition == t) & fw.coupling.str.startswith("OT_unb") &
                  (fw.cellset == "primary") & ~fw.uncapped & (fw.rep == "pooled")]
        put(f"rqOneUnbAUROCMin{n}", f(unb3.auroc.min()))
        put(f"rqOneUnbAUROCMax{n}", f(unb3.auroc.max()))
        # sensitivity
        os_ = one(fw, coupling=OT, transition=t, cellset="sensitivity", uncapped=False, rep="pooled")
        put(f"rqOneSensEnr{n}", f(os_.enrichment, 1))
        put(f"rqOneSensAUROC{n}", f(os_.auroc))
    put("rqOneBalAUROCMin", f(min(aur_bal)))
    put("rqOneBalAUROCMax", f(max(aur_bal)))
    unb3 = fw[(fw.transition == "III->IV") & fw.coupling.str.startswith("OT_unb") &
              (fw.cellset == "primary") & ~fw.uncapped & (fw.rep == "pooled")]
    put("rqOneUnbThreeLoCIBelowHalf", int((unb3.auroc_lo < 0.5).sum()))
    put("rqOneUnbThreeN", len(unb3))
    allknn = fw[(fw.coupling == "kNN")]
    put("rqOneKnnAUROCMinAll", f(allknn.auroc.min(), 3))
    put("rqOneKnnAUROCMaxAll", f(allknn.auroc.max(), 3))
    # per replicate range for the fixed OT coupling
    pr = fw[(fw.coupling == OT) & (fw.cellset == "primary") & ~fw.uncapped & (fw.rep != "pooled")]
    put("rqOneRepAUROCMin", f(pr.auroc.min()))
    put("rqOneRepAUROCMax", f(pr.auroc.max()))
    put("rqOneRepEnrMin", f(pr.enrichment.min(), 1))
    put("rqOneRepEnrMax", f(pr.enrichment.max(), 1))

    # cap robustness (I->II replicate 2)
    cap = {}
    for unc in (False, True):
        for cp in (OT, "kNN"):
            cap[(unc, cp)] = one(fw, coupling=cp, transition="I->II", cellset="primary",
                                 uncapped=unc, rep="2")
    for unc, tag in ((False, "Cap"), (True, "Uncap")):
        o = cap[(unc, OT)]
        put(f"cap{tag}AUROC", f(o.auroc, 3))
        put(f"cap{tag}Enr", f(o.enrichment, 2))
        put(f"cap{tag}TopOnePct", pct(o.top1, 2))
        put(f"cap{tag}ChanceSci", sci(o.m0))
        put(f"cap{tag}KnnAUROC", f(cap[(unc, "kNN")].auroc, 3))
    nlat = sv[(sv.transition == "I->II") & (sv.rep == 2) & (sv.cellset == "primary")]
    put("capNLaterUncapped", f"{int(nlat[nlat.uncapped].n_later.iloc[0]):,}".replace(",", "{,}"))
    put("capNLaterCapped", f"{int(nlat[~nlat.uncapped].n_later.iloc[0]):,}".replace(",", "{,}"))

    tables_rq1(fw, bw, cap, nlat)
    fig2(fw)
    fate()


def sci(x):
    m, e = f"{x:.1e}".split("e")
    return rf"{m}\times10^{{{int(e)}}}"


def tables_rq1(fw, bw, cap, nlat):
    # S1: primary vs sensitivity, pooled forward, key couplings
    rows = []
    lab = {OT: r"OT bal., $\varepsilon{=}0.05$", "OT_unb_eps0.05": r"OT unbal., $\varepsilon{=}0.05$",
           "kNN": r"kNN ($k{=}10$)"}
    for t in TR:
        for i, cp in enumerate(lab):
            cells = []
            for cs in ("primary", "sensitivity"):
                o = one(fw, coupling=cp, transition=t, cellset=cs, uncapped=False, rep="pooled")
                cells += [f"{f(o.enrichment, 1)} [{f(o.enrichment_lo, 1)}, {f(o.enrichment_hi, 1)}]",
                          f"{f(o.auroc)} [{f(o.auroc_lo)}, {f(o.auroc_hi)}]", pct(o.top10, 1)]
            rows.append((TLAB[t] if i == 0 else "") + " & " + lab[cp] + " & " + " & ".join(cells) + r" \\")
        rows.append(r"\addlinespace")
    body = "\n".join(rows[:-1])
    write_table("tab_s1_rq1_sensitivity", rf"""\begin{{table*}}[!t]
\centering
\caption{{RQ1 coupling recovery, primary vs.\ doublet-sensitivity cell set (pooled over
replicates, forward direction). Enrichment = own-clone mass / chance; top-10 = \% of early
cells with an own-clone cell among the 10 heaviest later cells. 95\% bootstrap CIs over clones.}}
\label{{tab:s1}}
\footnotesize
\begin{{tabular}}{{llcccccc}}
\toprule
 & & \multicolumn{{3}}{{c}}{{Primary}} & \multicolumn{{3}}{{c}}{{Sensitivity}} \\
\cmidrule(lr){{3-5}}\cmidrule(lr){{6-8}}
Transition & Coupling & Enrichment & AUROC & Top-10 (\%) & Enrichment & AUROC & Top-10 (\%) \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")

    # S4: cap robustness
    rows = []
    for cp, nm in ((OT, r"OT bal., $\varepsilon{=}0.05$"), ("kNN", r"kNN ($k{=}10$)")):
        for unc, pool in ((False, "capped"), (True, "uncapped")):
            o = cap[(unc, cp)]
            n = int(nlat[nlat.uncapped == unc].n_later.iloc[0])
            rows.append(f"{nm if not unc else ''} & {pool} ({n:,}) & ${sci(o.m0)}$ & ${sci(o.m)}$ & "
                        f"{f(o.enrichment, 2)} & {f(o.auroc, 3)} & {pct(o.top1, 2)} & {pct(o.top10, 2)}" + r" \\")
    write_table("tab_s4_cap", rf"""\begin{{table*}}[!t]
\centering
\caption{{Cap robustness (pre-specified): I$\to$II replicate~2 with the background cap
(8{{,}}000 per stage) and without it. The uncapped OT solve used the block-wise solver of
amendment A1, validated against POT before use. Later-cell pool size in parentheses.}}
\label{{tab:s4}}
\footnotesize
\setlength{{\tabcolsep}}{{3pt}}
\begin{{tabular}}{{llcccccc}}
\toprule
Coupling & Pool & $m_0$ & $m$ & $m/m_0$ & AUROC & Top-1 (\%) & Top-10 (\%) \\
\midrule
{chr(10).join(rows).replace(",", "{,}").replace("{{,}}", "{,}")}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")

    # S5: per replicate, OT bal 0.05, forward and backward, primary
    rows = []
    for t in TR:
        for i, rep in enumerate(["1", "2", "3", "pooled"]):
            o = one(fw, coupling=OT, transition=t, cellset="primary", uncapped=False, rep=rep)
            b = one(bw, coupling=OT, transition=t, cellset="primary", uncapped=False, rep=rep)
            rows.append(f"{TLAB[t] if i == 0 else ''} & {rep} & {int(o.n_clones)} & "
                        f"{f(o.enrichment, 1)} [{f(o.enrichment_lo, 1)}, {f(o.enrichment_hi, 1)}] & "
                        f"{f(o.auroc)} [{f(o.auroc_lo)}, {f(o.auroc_hi)}] & {pct(o.top10, 1)} & "
                        f"{f(b.enrichment, 1)} & {f(b.auroc)}" + r" \\")
        rows.append(r"\addlinespace")
    write_table("tab_s5_rq1_replicates", rf"""\begin{{table*}}[!t]
\centering
\caption{{RQ1 by replicate (primary set, OT balanced, $\varepsilon=0.05$). Replicates are separate
differentiations; couplings are computed within replicate. Pooled = clone-weighted.
Backward = later$\to$early direction (secondary).}}
\label{{tab:s5}}
\footnotesize
\begin{{tabular}}{{llcccccc}}
\toprule
 & & & \multicolumn{{3}}{{c}}{{Forward}} & \multicolumn{{2}}{{c}}{{Backward}} \\
\cmidrule(lr){{4-6}}\cmidrule(lr){{7-8}}
Transition & Replicate & Clones & Enrichment & AUROC & Top-10 (\%) & Enrichment & AUROC \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")


def fig2(fw):
    p = fw[(fw.cellset == "primary") & ~fw.uncapped & (fw.rep == "pooled")]
    fig, axes = plt.subplots(2, 3, figsize=(COL2, 3.4), sharex=True, sharey="row")
    for j, t in enumerate(TR):
        q = p[p.transition == t]
        for i, (met, lo, hi, ylab) in enumerate((("auroc", "auroc_lo", "auroc_hi", "Clone-retrieval AUROC"),
                                                 ("enrichment", "enrichment_lo", "enrichment_hi",
                                                  "Own-clone mass / chance"))):
            ax = axes[i, j]
            for kind, col, mk, nm, dx in (("bal", BLUE, "o", "OT, balanced", 0.94),
                                          ("unb", ORANGE, "s", "OT, unbalanced", 1.06)):
                r = [one(q, coupling=f"OT_{kind}_eps{e}") for e in EPS]
                x = np.array(EPS) * dx
                y = np.array([v[met] for v in r])
                err = np.array([[v[met] - v[lo] for v in r], [v[hi] - v[met] for v in r]])
                ax.errorbar(x, y, yerr=err, color=col, marker=mk, ms=4, lw=1.4, elinewidth=0.8,
                            capsize=0, label=nm, zorder=3)
            knn = one(q, coupling="kNN")
            ax.axhline(knn[met], color=AQUA, lw=1.2, ls="--", label="kNN ($k$=10)", zorder=2)
            ax.axhline(1.0 if met == "enrichment" else 0.5, color=GREY, lw=1.0, ls=":",
                       label="Independent (chance)", zorder=2)
            ax.set_xscale("log")
            ax.set_xticks(EPS)
            ax.set_xticklabels([str(e) for e in EPS])
            ax.minorticks_off()
            if met == "enrichment":
                ax.set_yscale("log")
                ax.set_yticks([0.1, 0.3, 1, 3, 10, 30])
                ax.set_yticklabels(["0.1", "0.3", "1", "3", "10", "30"])
                ax.set_xlabel(r"Entropy $\varepsilon$")
            if i == 0:
                ax.set_title(TPLOT[t])
            if j == 0:
                ax.set_ylabel(ylab)
    h, l = axes[0, 0].get_legend_handles_labels()
    order = [l.index(x) for x in ("OT, balanced", "OT, unbalanced", "kNN ($k$=10)", "Independent (chance)")]
    fig.legend([h[i] for i in order], [l[i] for i in order], loc="upper center", ncol=4,
               bbox_to_anchor=(0.5, 1.04))
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(FIG / "fig2_coupling_recovery.pdf")
    plt.close(fig)


def fate():
    fa = pd.read_csv(R / "06_rq1/fate_agreement.csv")
    p = fa[(fa.cellset == "primary") & ~fa.uncapped]
    ps = fa[(fa.cellset == "sensitivity") & ~fa.uncapped]
    better, total, sagree, sn = 0, 0, 0, 0
    for t in TR:
        q = p[p.transition == t]
        cl = sorted(q.early_cluster.unique())
        n_b = 0
        for c in cl:
            o = one(q, early_cluster=c, coupling=OT)
            i = one(q, early_cluster=c, coupling="independent")
            n_b += int(o.tvd < i.tvd)
            key = c.split("-", 1)[1].replace("+", "")
            put(f"tvdOT{TN[t]}{key}", f(o.tvd))
            put(f"tvdInd{TN[t]}{key}", f(i.tvd))
            put(f"tvdOTLo{TN[t]}{key}", f(o.tvd_lo))
            put(f"tvdOTHi{TN[t]}{key}", f(o.tvd_hi))
            put(f"tvdClones{TN[t]}{key}", int(o.n_clones))
            put(f"tvdKnn{TN[t]}{key}", f(one(q, early_cluster=c, coupling="kNN").tvd))
        put(f"fateNClusters{TN[t]}", len(cl))
        put(f"fateNBetter{TN[t]}", n_b)
        better += n_b
        total += len(cl)
        # sensitivity: same direction (OT < independent) per cluster?
        qs = ps[ps.transition == t]
        agree = 0
        for c in cl:
            if c not in set(qs.early_cluster):
                continue
            so = one(qs, early_cluster=c, coupling=OT).tvd < one(qs, early_cluster=c, coupling="independent").tvd
            agree += int(so == (one(q, early_cluster=c, coupling=OT).tvd < one(q, early_cluster=c, coupling="independent").tvd))
        put(f"fateSensAgree{TN[t]}", agree)
        sagree += agree
        sn += len(set(cl) & set(qs.early_cluster))
        ot_t = q[q.coupling == OT].tvd
        put(f"fateOTMin{TN[t]}", f(ot_t.min()))
        put(f"fateOTMax{TN[t]}", f(ot_t.max()))
        put(f"fateSensN{TN[t]}", len(set(cl) & set(qs.early_cluster)))
    put("fateNBetterAll", better)
    put("fateNClustersAll", total)
    put("fateSensAgreeAll", sagree)
    put("fateSensNAll", sn)

    # Fig 3
    fig, axes = plt.subplots(1, 3, figsize=(COL2, 2.5), sharex=True)
    for ax, t in zip(axes, TR):
        q = p[p.transition == t]
        cl = sorted(q.early_cluster.unique(), key=lambda c: -one(q, early_cluster=c, coupling=OT).n_clones)
        y = np.arange(len(cl))[::-1]
        for cp, col, mk, nm, dy in ((OT, BLUE, "o", r"OT, balanced, $\varepsilon$=0.05", 0.18),
                                    ("kNN", AQUA, "^", "kNN ($k$=10)", 0.0),
                                    ("independent", GREY, "s", "Independent", -0.18)):
            r = [one(q, early_cluster=c, coupling=cp) for c in cl]
            v = np.array([x.tvd for x in r])
            err = np.array([[x.tvd - x.tvd_lo for x in r], [x.tvd_hi - x.tvd for x in r]])
            ax.errorbar(v, y + dy, xerr=err, fmt=mk, color=col, ms=4, elinewidth=0.8, capsize=0,
                        label=nm, zorder=3)
        labs = [f"{c}  (n={int(one(q, early_cluster=c, coupling=OT).n_clones)})" for c in cl]
        ax.set_yticks(y)
        ax.set_yticklabels(labs)
        for tl, c in zip(ax.get_yticklabels(), cl):
            if c.endswith("NPvMB"):
                tl.set_fontweight("bold")
        ax.set_title(TPLOT[t])
        ax.set_xlim(0, 1)
        ax.set_xlabel("TVD to clonal fate distribution")
        ax.grid(axis="y", visible=False)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.06))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIG / "fig3_fate_agreement.pdf")
    plt.close(fig)


# ============================================================ RQ2
MODELS = ["Ridge", "RandomForest", "Transformer"]
MLAB = {"ConstMedian": "Const.\\ median", "ConstMean": "Const.\\ mean", "Ridge": "Ridge",
        "RandomForest": "Random forest", "Transformer": "Transformer"}
CPL = ["true", "OT", "independent"]
CLAB = {"true": "True lineage", "OT": r"OT ($\varepsilon{=}0.05$)", "independent": "Independent"}


def rq2():
    ms = pd.read_csv(R / "07_rq2/mae_summary.csv")
    cmp_ = pd.read_csv(R / "07_rq2/comparisons.csv")
    ps = pd.read_csv(R / "07_rq2/per_seed.csv")
    pc = cmp_[cmp_.cellset == "primary"]
    put("rqTwoSeeds", int(ps.seed.nunique()))
    put("rqTwoDeltaMAE", f(DELTA, 2))
    nt = ps[(ps.cellset == "primary") & (ps.seed == 0) & (ps.coupling == "true") & (ps.model == "ConstMedian")]
    for t in TR:
        r = one(nt, transition=t)
        put(f"rqTwoNTest{TN[t]}", int(r.n_test))
        put(f"rqTwoNTrain{TN[t]}", int(r.n_train))

    # true lineage: models vs constant median
    tv = pc[(pc.coupling == "true") & pc.model.isin(MODELS) & (pc.reference == "ConstMedian")]
    put("rqTwoTrueGainMinPct", f(-tv.delta_pct_of_const.max(), 1))
    put("rqTwoTrueGainMaxPct", f(-tv.delta_pct_of_const.min(), 1))
    put("rqTwoTrueNEquiv", int(tv.equivalent.sum()))
    put("rqTwoTrueN", len(tv))
    put("rqTwoTrueNSig", int((tv.p_holm < 0.05).sum()))
    put("rqTwoTrueNBetter", int((tv.ci_hi < 0).sum()))
    put("rqTwoTrueDzMin", f(-tv.dz.max(), 2))
    put("rqTwoTrueDzMax", f(-tv.dz.min(), 2))
    for t in TR:
        rf = one(tv, transition=t, model="RandomForest")
        put(f"rqTwoRFDelta{TN[t]}", f(1000 * rf.delta, 1))
        put(f"rqTwoRFDeltaLo{TN[t]}", f(1000 * rf.ci_lo, 1))
        put(f"rqTwoRFDeltaHi{TN[t]}", f(1000 * rf.ci_hi, 1))
        put(f"rqTwoRFPct{TN[t]}", f(-rf.delta_pct_of_const, 1))
        put(f"rqTwoConstMAE{TN[t]}", f(one(ms, cellset="primary", transition=t, coupling="true",
                                           model="ConstMedian")["mean"], 3))
        put(f"rqTwoDeltaPct{TN[t]}", f(100 * DELTA / rf.const_median_mae, 1))
    # Transformer vs RF
    trf = pc[(pc.model == "Transformer") & (pc.reference == "RandomForest")]
    put("rqTwoTrfRFNEquiv", int(trf.equivalent.sum()))
    put("rqTwoTrfRFN", len(trf))
    put("rqTwoTrfRFNWorse", int((trf.ci_lo > 0).sum()))
    put("rqTwoTrfRFNBetter", int((trf.ci_hi < 0).sum()))
    # anything beat RF? (any model with CI entirely below 0 vs RF)
    vsrf = pc[(pc.reference == "RandomForest") & pc.model.isin(["Ridge", "Transformer"])]
    put("rqTwoBeatRFN", int((vsrf.ci_hi < 0).sum()))
    put("rqTwoBeatRFTotal", len(vsrf))
    put("rqTwoBeatRFMaxPct", f(-vsrf.delta_pct_of_const.min(), 1))
    put("rqTwoTrueBeatRFN", int((vsrf[vsrf.coupling == "true"].ci_hi < 0).sum()))
    put("rqTwoTrueVsRFN", int((vsrf.coupling == "true").sum()))
    put("rqTwoTrueRidgeWorseRFN", int((vsrf[(vsrf.coupling == "true") & (vsrf.model == "Ridge")].ci_lo > 0).sum()))
    # independent coupling at I->II and II->III: models vs own constant median
    iv = pc[(pc.coupling == "independent") & pc.model.isin(MODELS) & (pc.reference == "ConstMedian") &
            pc.transition.isin(TR[:2])]
    put("rqTwoIndNWorse", int((iv.delta > 0).sum()))
    put("rqTwoIndN", len(iv))
    put("rqTwoIndPctMin", f(iv.delta_pct_of_const.min(), 1))
    put("rqTwoIndPctMax", f(iv.delta_pct_of_const.max(), 1))
    # OT coupling at I->II and II->III
    ov = pc[(pc.coupling == "OT") & pc.model.isin(MODELS) & (pc.reference == "ConstMedian") &
            pc.transition.isin(TR[:2])]
    put("rqTwoOTGainMinPct", f(-ov.delta_pct_of_const.max(), 1))
    put("rqTwoOTGainMaxPct", f(-ov.delta_pct_of_const.min(), 1))
    # III->IV artefact
    t3 = ms[(ms.cellset == "primary") & (ms.transition == "III->IV")]
    for cp in CPL:
        put(f"rqTwoThreeConst{cp.capitalize()}", f(one(t3, coupling=cp, model="ConstMedian")["mean"], 3))
    oi = t3[t3.coupling.isin(["OT", "independent"]) & t3.model.isin(MODELS)]
    put("rqTwoThreeOIMin", f(oi["mean"].min(), 3))
    put("rqTwoThreeOIMax", f(oi["mean"].max(), 3))
    tt = t3[(t3.coupling == "true") & t3.model.isin(MODELS)]
    put("rqTwoThreeTrueMin", f(tt["mean"].min(), 3))
    put("rqTwoThreeTrueMax", f(tt["mean"].max(), 3))
    a3 = pc[(pc.transition == "III->IV") & pc.coupling.isin(["OT", "independent"]) &
            pc.model.isin(MODELS) & (pc.reference == "ConstMedian")]
    put("rqTwoThreeOIGainMinPct", f(-a3.delta_pct_of_const.max(), 1))
    put("rqTwoThreeOIGainMaxPct", f(-a3.delta_pct_of_const.min(), 1))
    # every non-equivalent pre-specified comparison
    ne = pc[~pc.equivalent]
    put("rqTwoNNonEquiv", len(ne))
    put("rqTwoNComparisons", len(pc))
    ne_other = ne[~((ne.transition == "III->IV") & ne.coupling.isin(["OT", "independent"]))]
    assert len(ne_other) == 1, ne_other
    put("rqTwoNNonEquivThree", len(ne) - len(ne_other))
    r = ne_other.iloc[0]
    put("rqTwoNEOtherDesc", f"{MLAB[r.model]} vs.\\ {MLAB[r.reference].lower()} at {TLAB[r.transition]} "
                            f"({CLAB[r.coupling].lower() if r.coupling == 'true' else CLAB[r.coupling]})")
    put("rqTwoNEOtherPct", f(r.delta_pct_of_const, 1))
    # sensitivity summary
    sc = cmp_[cmp_.cellset == "sensitivity"]
    stv = sc[(sc.coupling == "true") & sc.model.isin(MODELS) & (sc.reference == "ConstMedian")]
    put("rqTwoSensTrueGainMinPct", f(-stv.delta_pct_of_const.max(), 1))
    put("rqTwoSensTrueGainMaxPct", f(-stv.delta_pct_of_const.min(), 1))
    put("rqTwoSensTrueNEquiv", int(stv.equivalent.sum()))
    sne = sc[~sc.equivalent]
    same = set(map(tuple, sne[["transition", "coupling", "model", "reference"]].values)) == \
        set(map(tuple, ne[["transition", "coupling", "model", "reference"]].values))
    put("rqTwoSensSameNonEquiv", "the same" if same else "different")
    put("rqTwoSensNNonEquiv", len(sne))

    table_rq2(ms)
    table_s2(cmp_)
    fig4(pc)


def table_rq2(ms):
    p = ms[ms.cellset == "primary"]
    cols = ["ConstMedian", "ConstMean"] + MODELS
    rows = []
    for t in TR:
        for i, cp in enumerate(CPL):
            q = p[(p.transition == t) & (p.coupling == cp)]
            vals = [one(q, model=m) for m in cols]
            best = min(v["mean"] for v in vals[2:])
            cells = []
            for m, v in zip(cols, vals):
                s = f"{f(v['mean'], 3)} ({f(v['std'], 3)})"
                cells.append(rf"\textbf{{{s}}}" if m in MODELS and v["mean"] == best else s)
            rows.append(f"{TLAB[t] if i == 0 else ''} & {CLAB[cp]} & " + " & ".join(cells) + r" \\")
        rows.append(r"\addlinespace")
    write_table("tab_rq2_mae", rf"""\begin{{table*}}[!t]
\centering
\caption{{RQ2 forecasting error on held-out true sister pairs: MAE over (Pg, D), mean (SD) over
{M['rqTwoSeeds']} clone-level splits (primary set). Rows give the coupling used to build
\emph{{training}} pairs; every row is tested on the same held-out true clonal pairs. Each
constant is fitted to its own row's training targets. Bold: lowest model MAE in the row.}}
\label{{tab:rq2}}
\footnotesize
\begin{{tabular}}{{llccccc}}
\toprule
Transition & Training pairs & {' & '.join(MLAB[c] for c in cols)} \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")


def table_s2(cmp_):
    rows = []
    for t in TR:
        for i, cp in enumerate(CPL):
            cells = []
            for cs in ("primary", "sensitivity"):
                q = cmp_[(cmp_.cellset == cs) & (cmp_.transition == t) & (cmp_.coupling == cp)]
                for m in MODELS:
                    r = one(q, model=m, reference="ConstMedian")
                    s = f"{f(r.delta_pct_of_const, 1)}"
                    cells.append(s + ("" if r.equivalent else r"$^\dagger$"))
                r = one(q, model="Transformer", reference="RandomForest")
                cells.append(f"{f(r.delta_pct_of_const, 1)}" + ("" if r.equivalent else r"$^\dagger$"))
            rows.append(f"{TLAB[t] if i == 0 else ''} & {CLAB[cp]} & " + " & ".join(cells) + r" \\")
        rows.append(r"\addlinespace")
    write_table("tab_s2_rq2_sensitivity", rf"""\begin{{table*}}[!t]
\centering
\caption{{RQ2, primary vs.\ doublet-sensitivity set: paired $\Delta$MAE as \% of the
constant-median MAE (negative = model better). Columns 3--5 and 7--9: model vs.\ constant
median of the same training pairs; columns 6 and 10: Transformer vs.\ random forest.
$^\dagger$Not equivalent: the 95\% CI is not inside $\pm\delta$ ($\delta={M['rqTwoDeltaMAE']}$ MAE).}}
\label{{tab:s2}}
\footnotesize
\begin{{tabular}}{{llcccccccc}}
\toprule
 & & \multicolumn{{4}}{{c}}{{Primary}} & \multicolumn{{4}}{{c}}{{Sensitivity}} \\
\cmidrule(lr){{3-6}}\cmidrule(lr){{7-10}}
Transition & Training pairs & Ridge & RF & Transf. & Transf.\,vs\,RF & Ridge & RF & Transf. & Transf.\,vs\,RF \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")


def fig4(pc):
    fig, axes = plt.subplots(1, 3, figsize=(COL2, 2.4), sharey=True)
    ccol = {"true": BLUE, "OT": ORANGE, "independent": GREY}
    cmk = {"true": "o", "OT": "s", "independent": "D"}
    for ax, t in zip(axes, TR):
        q = pc[(pc.transition == t) & (pc.reference == "ConstMedian")]
        ax.axvspan(-DELTA * 1000, DELTA * 1000, color="#f0efec", zorder=0, lw=0)
        ax.axvline(0, color=INK2, lw=0.6, zorder=1)
        yt, yl = [], []
        for gi, m in enumerate(MODELS):
            for ci, cp in enumerate(CPL):
                y = -(gi * 4 + ci)
                r = one(q, model=m, coupling=cp)
                ax.errorbar(1000 * r.delta, y, xerr=[[1000 * (r.delta - r.ci_lo)], [1000 * (r.ci_hi - r.delta)]],
                            fmt=cmk[cp], color=ccol[cp], ms=4, elinewidth=0.9, capsize=0, zorder=3,
                            label=CLAB[cp].replace(r"$\varepsilon{=}0.05$", "ε=0.05") if gi == 0 else None)
            yt.append(-(gi * 4 + 1))
            yl.append(MLAB[m].replace("\\ ", " "))
        ax.set_yticks(yt)
        ax.set_yticklabels(yl)
        ax.set_title(TPLOT[t])
        ax.set_xlabel(r"$\Delta$MAE vs. constant median ($\times10^{-3}$)")
        ax.grid(axis="y", visible=False)
        ax.set_xlim(-25, 12)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, ["Training pairs:  " + l[0]] + l[1:], loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.07))
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIG / "fig4_forecasting.pdf")
    plt.close(fig)


# ============================================================ RQ3
REPS = ["P_tcbb_D", "Pg_D", "PC30", "fate_axis"]
RLAB = {"P_tcbb_D": "(P, D)", "Pg_D": "(Pg, D)", "PC30": "30-PC", "fate_axis": "Fate-bias axis"}
RMAC = {"P_tcbb_D": "PD", "Pg_D": "PgD", "PC30": "PC", "fate_axis": "Fate"}


def rq3():
    ps = pd.read_csv(R / "08_rq3/per_seed.csv")
    cm = pd.read_csv(R / "08_rq3/comparisons.csv")
    da = pd.read_csv(R / "08_rq3/da_descriptive.csv")
    p = ps[ps.cellset == "primary"]
    put("rqThreeSeeds", int(p.seed.nunique()))
    for t in TR:
        n = TN[t]
        q = p[p.transition == t]
        put(f"rqThreeClasses{n}", int(q.n_classes.iloc[0]))
        put(f"rqThreeChanceBacc{n}", f(q.chance_bacc.mean(), 2))
        for r in REPS:
            g = q[q.representation == r]
            put(f"rqThreeAUC{RMAC[r]}{n}", f(g.macro_auroc.mean()))
            put(f"rqThreeAUCSD{RMAC[r]}{n}", f(g.macro_auroc.std(), 2))
            put(f"rqThreeBacc{RMAC[r]}{n}", f(g.balanced_acc.mean(), 2))
        c = cm[(cm.cellset == "primary") & (cm.transition == t) & (cm.metric == "macro_auroc")]
        for a, b, tag in (("Pg_D", "PC30", "PgPC"), ("Pg_D", "fate_axis", "PgFate"),
                          ("PC30", "fate_axis", "PCFate"), ("P_tcbb_D", "Pg_D", "PDPg"),
                          ("Pg_D", "chance", "PgChance")):
            r = one(c, a=a, b=b)
            put(f"rqThreeD{tag}{n}", f(r.delta))
            put(f"rqThreeD{tag}Lo{n}", f(r.ci_lo))
            put(f"rqThreeD{tag}Hi{n}", f(r.ci_hi))
            put(f"rqThreeDz{tag}{n}", f(r.dz, 1))
        cc = pd.read_csv(R / f"08_rq3/primary_{t.replace('->', '-')}_class_counts.csv")
        put(f"rqThreeMinClass{n}", int(cc[cc.fate != "other"].n_clones.min()))
    pm = cm[(cm.cellset == "primary") & (cm.metric == "macro_auroc")]
    mk = pm[pm.b == "chance"]
    put("rqThreeAllAboveChance", "all" if (mk.ci_lo > 0).all() else "not all")
    bt = pm[pm.a.isin(["P_tcbb_D", "Pg_D"]) & pm.b.isin(["PC30", "fate_axis"])]
    put("rqThreeMarkerGapMin", f(-bt.delta.max()))
    put("rqThreeMarkerGapMax", f(-bt.delta.min()))
    put("rqThreeMarkerNSig", int((bt.p_holm < 0.05).sum()))
    put("rqThreeMarkerN", len(bt))
    pf = pm[(pm.a == "PC30") & (pm.b == "fate_axis")].delta.abs()
    put("rqThreePCFateMin", f(pf.min()))
    put("rqThreePCFateMax", f(pf.max()))
    put("rqThreePDPgMaxAbs", f(pm[(pm.a == "P_tcbb_D") & (pm.b == "Pg_D")].delta.abs().max()))
    # sensitivity
    s = ps[ps.cellset == "sensitivity"]
    for t in TR:
        for r in ("Pg_D", "PC30"):
            put(f"rqThreeSensAUC{RMAC[r]}{TN[t]}", f(s[(s.transition == t) & (s.representation == r)].macro_auroc.mean()))
    sm = cm[(cm.cellset == "sensitivity") & (cm.metric == "macro_auroc")]
    sbt = sm[sm.a.isin(["P_tcbb_D", "Pg_D"]) & sm.b.isin(["PC30", "fate_axis"])]
    put("rqThreeSensMarkerNSig", int(((sbt.p_holm < 0.05) & (sbt.delta < 0)).sum()))
    # DA descriptive (primary)
    d = da[da.cellset == "primary"].set_index("quantity").value
    put("daNPvMBClones", int(d["n_NPvMB_clones"]))
    put("daNmDA", int(d["n_mDA_producing"]))
    put("daAUCPg", f(d["auroc_Pg_mDA_vs_other"]))
    put("daAUCD", f(d["auroc_D_mDA_vs_other"]))
    put("daAUCLDA", f(d["auroc_lda1_mDA_vs_other"]))
    fs = d[d.index.str.startswith("fate_share_")].sort_values(ascending=False)
    put("daTopFate", fs.index[0].replace("fate_share_", ""))
    put("daTopFatePct", pct(fs.iloc[0], 0))
    put("daSecondFate", fs.index[1].replace("fate_share_", ""))
    put("daSecondFatePct", pct(fs.iloc[1], 0))
    put("daMDAPct", pct(d["fate_share_IV-mDAs"], 0))
    table_s3(ps)
    fig5(p)


def table_s3(ps):
    rows = []
    for t in TR:
        cells = []
        for cs in ("primary", "sensitivity"):
            q = ps[(ps.cellset == cs) & (ps.transition == t)]
            for r in REPS:
                g = q[q.representation == r].macro_auroc
                cells.append(f"{f(g.mean())} ({f(g.std())})")
        rows.append(f"{TLAB[t]} & " + " & ".join(cells) + r" \\")
    write_table("tab_s3_rq3_sensitivity", rf"""\begin{{table*}}[!t]
\centering
\caption{{RQ3 fate prediction from the early sister cell: macro-AUROC, mean (SD) over
{M['rqThreeSeeds']} clone-level splits, primary vs.\ doublet-sensitivity set. Chance = 0.5.}}
\label{{tab:s3}}
\footnotesize
\begin{{tabular}}{{lcccccccc}}
\toprule
 & \multicolumn{{4}}{{c}}{{Primary}} & \multicolumn{{4}}{{c}}{{Sensitivity}} \\
\cmidrule(lr){{2-5}}\cmidrule(lr){{6-9}}
Transition & {' & '.join(RLAB[r] for r in REPS)} & {' & '.join(RLAB[r] for r in REPS)} \\
\midrule
{chr(10).join(rows)}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")


def fig5(p):
    fig, ax = plt.subplots(figsize=(COL1, 2.3))
    cols = {"P_tcbb_D": YELLOW, "Pg_D": ORANGE, "PC30": BLUE, "fate_axis": AQUA}
    mks = {"P_tcbb_D": "v", "Pg_D": "s", "PC30": "o", "fate_axis": "D"}
    rng = np.random.default_rng(0)  # jitter only
    for j, t in enumerate(TR):
        for k, r in enumerate(REPS):
            g = p[(p.transition == t) & (p.representation == r)].macro_auroc.to_numpy()
            x = j + (k - 1.5) * 0.2
            ax.scatter(x + rng.uniform(-0.05, 0.05, len(g)), g, s=5, color=cols[r], alpha=0.35, lw=0, zorder=2)
            ax.errorbar(x, g.mean(), yerr=g.std(), fmt=mks[r], color=cols[r], ms=4.5, mec=INK, mew=0.4,
                        elinewidth=0.9, capsize=0, zorder=3, label=RLAB[r] if j == 0 else None)
    ax.axhline(0.5, color=GREY, ls=":", lw=1.0, label="Chance", zorder=1)
    ax.set_xticks(range(3))
    ax.set_xticklabels([TPLOT[t] for t in TR])
    ax.set_ylabel("Fate macro-AUROC")
    ax.set_ylim(0.45, 0.92)
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.17), handletextpad=0.2, columnspacing=0.8)
    fig.savefig(FIG / "fig5_state_adequacy.pdf")
    plt.close(fig)


# ============================================================ RQ4
def rq4():
    g = pd.read_csv(R / "09_rq4/grid_summary.csv")
    me = pd.read_csv(R / "09_rq4/measured_summary.csv")
    est = me[me.bound == "estimate"]
    put("rqFourMeasMin", f(est.m.min(), 3))
    put("rqFourMeasMax", f(est.m.max(), 3))
    put("rqFourMeasGap", f(1 - est.m.min(), 3))
    for t in TR:
        put(f"rqFourMass{TN[t]}", f(one(est, transition=t, rate_kind="1-own_clone_mass", rep=1).m, 4))
        put(f"rqFourTop{TN[t]}", f(one(est, transition=t, rate_kind="1-top1_recall", rep=1).m, 4))
    for rep in (1, 2, 3):
        put(f"rqFourRMSEOne{['', 'A', 'B', 'C'][rep]}", f(one(g, m=1.0, rep=rep).rmse, 2))
        put(f"rqFourRMSEHalf{['', 'A', 'B', 'C'][rep]}", f(one(g, m=0.5, rep=rep).rmse, 2))
        put(f"rqFourRMSETenth{['', 'A', 'B', 'C'][rep]}", f(one(g, m=0.1, rep=rep).rmse, 3))
    ratio = [one(est, transition=t, rate_kind=k, rep=rep).rmse / one(g, m=1.0, rep=rep).rmse
             for t in TR for k in ("1-own_clone_mass", "1-top1_recall") for rep in (1, 2, 3)]
    put("rqFourRatioMinPct", f(100 * min(ratio), 1))
    put("rqFourRatioMaxPct", f(100 * max(ratio), 1))
    put("rqFourROneMin", f(g[g.m == 1.0].r.min(), 2))
    put("rqFourROneMax", f(g[g.m == 1.0].r.max(), 2))
    put("rqFourMaxAbsBias", sci(max(g.bias.abs().max(), me.bias.abs().max())))
    put("rqFourMaxAbsBiasCI", sci(max(g.bias_lo.abs().max(), g.bias_hi.abs().max())))
    put("rqFourEffectMean", "0.5")
    fig6(g, est)


def fig6(g, est):
    fig, axes = plt.subplots(1, 2, figsize=(COL1, 1.75), sharex=True)
    cols = {1: BLUE, 2: ORANGE, 3: AQUA}
    mks = {1: "o", 2: "s", 3: "^"}
    for ax, (met, ylab) in zip(axes, (("rmse", "Per-cell effect RMSE"), ("r", "Pearson $r$ with true effect"))):
        for rep in (1, 2, 3):
            q = g[g.rep == rep].sort_values("m")
            ax.fill_between(q.m, q[f"{met}_lo"], q[f"{met}_hi"], color=cols[rep], alpha=0.2, lw=0)
            ax.plot(q.m, q[met], color=cols[rep], marker=mks[rep], ms=3, lw=1.2, label=f"Replicate {rep}")
            mm = est[(est.rep == rep) & (est.rate_kind == "1-own_clone_mass")]
            ax.plot(mm.m, mm[met], ls="none", marker="*", ms=7, mfc="none", mec=INK, mew=0.6, zorder=4,
                    label="RQ1-measured OT rate" if rep == 1 else None)
        ax.set_xlabel("Mismatch rate $m$")
        ax.set_ylabel(ylab)
        ax.set_xticks([0, 0.5, 1])
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.2))
    fig.tight_layout(rect=(0, 0, 1, 0.93), w_pad=1.0)
    fig.savefig(FIG / "fig6_perturbation.pdf")
    plt.close(fig)


# ============================================================ output
def write_table(name, text):
    (TAB / f"{name}.tex").write_text(
        f"% Generated by scripts/10_results_paper.py from results/06-09. Do not edit.\n{text}",
        encoding="utf-8")


def check_wording():
    """Results text uses categorical words ("all", "every", "none") for these; fail if untrue."""
    eq = [("rqTwoTrueNBetter", "rqTwoTrueN"), ("rqTwoTrueNEquiv", "rqTwoTrueN"),
          ("rqTwoTrfRFNEquiv", "rqTwoTrfRFN"), ("rqTwoIndNWorse", "rqTwoIndN"),
          ("rqThreeMarkerNSig", "rqThreeMarkerN"), ("nOTSolvesConverged", "nOTSolves")]
    for a, b in eq:
        assert M[a] == M[b], (a, M[a], b, M[b])
    assert M["rqTwoTrueBeatRFN"] == "0"
    assert M["rqThreeAllAboveChance"] == "all"


def main():
    rq1()
    rq2()
    rq3()
    rq4()
    check_wording()
    lines = ["% Generated by scripts/10_results_paper.py from results/06_rq1 ... 09_rq4.",
             "% Do not edit by hand."]
    lines += [rf"\newcommand{{\{k}}}{{{v}}}" for k, v in M.items()]
    (GEN / "results_numbers.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
    pd.Series(M, name="value").rename_axis("macro").to_csv(OUT / "results_numbers.csv")

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT,
                           capture_output=True, text=True).stdout.strip()
    prov = [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
            "plan: ANALYSIS_PLAN.md approved at 65cd70d (amendment A1)",
            f"python {sys.version.split()[0]} pandas {pd.__version__} numpy {np.__version__} "
            f"matplotlib {matplotlib.__version__}",
            "environment: Windows (VS Code session)"]
    prov += [f"input {p.relative_to(ROOT).as_posix()} sha256 {sha256(p)}" for p in INPUTS]
    (OUT / "provenance.txt").write_text("\n".join(prov) + "\n")
    print(f"{len(M)} macros; tables and figures written")
    print(json.dumps(M, indent=0)[:6000])


if __name__ == "__main__":
    main()
