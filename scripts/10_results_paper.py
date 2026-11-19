"""Results numbers, figures and tables for paper/main.tex (RQ1-RQ4).

Reads only committed result files in results/06_rq1 ... results/09_rq4 (final at
35d89c2) and the post-review results 11-17 (amendment A3), and writes:
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

# Figure style: IEEE compsoc (sans captions) -> Arial; drawn at final print size so LaTeX
# never rescales (text width 7.17 in, column 3.50 in). Okabe-Ito colours, checked for
# colour-vision deficiency; every series also has its own marker. Fixed meaning across
# figures: OT = blue, unbalanced OT = vermilion, kNN = green, independent/chance = grey,
# true lineage = black.
C_OT, C_UNB, C_KNN, C_IND, C_TRUE = "#0072B2", "#D55E00", "#009E73", "#7F7F7F", "#000000"
C_ORANGE, C_PURPLE = "#E69F00", "#CC79A7"
INK, INK2, BAND = "#000000", "#4D4D4D", "#EBEBEB"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "font.size": 7, "axes.titlesize": 7.5, "axes.titleweight": "normal", "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "text.color": INK, "axes.linewidth": 0.6, "lines.linewidth": 1.0,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "legend.frameon": False, "legend.handlelength": 1.6, "legend.columnspacing": 1.2,
    "legend.handletextpad": 0.4, "errorbar.capsize": 1.5,
    "pdf.fonttype": 42, "figure.dpi": 300})
COL1, COL2 = 3.50, 7.17  # IEEE compsoc column / text width (in), measured from IEEEtran


def panel(ax, letter, x=-0.16, y=1.02):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left")


def save(fig, name):
    # no bbox cropping: the PDF has exactly the figure size; no timestamp, so reruns on
    # unchanged results give byte-identical files
    fig.savefig(FIG / name, metadata={"CreationDate": None})
    plt.close(fig)


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
    # cross-transition ranges for the fixed coupling (Discussion)
    fx = fw[(fw.coupling == OT) & (fw.cellset == "primary") & ~fw.uncapped & (fw.rep == "pooled")]
    kx = fw[(fw.coupling == "kNN") & (fw.cellset == "primary") & ~fw.uncapped & (fw.rep == "pooled")]
    put("rqOneAUROCMinT", f(fx.auroc.min()))
    put("rqOneAUROCMaxT", f(fx.auroc.max()))
    put("rqOneEnrMinT", f(fx.enrichment.min(), 1))
    put("rqOneEnrMaxT", f(fx.enrichment.max(), 1))
    put("rqOneMassPctMaxT", pct(fx.m.max(), 2))
    put("rqOneTopOnePctMaxT", pct(fx.top1.max(), 2))
    put("rqOneKnnEnrMinT", f(kx.enrichment.min(), 1))
    put("rqOneKnnEnrMaxT", f(kx.enrichment.max(), 1))
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
(8{{,}}000 per stage) and without it. The uncapped OT solve used the memory-efficient
block-wise Sinkhorn solver, validated against POT before use. Later-cell pool size in parentheses.}}
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


def fig1():
    """Overview schematic (no result values): lineage ground truth and analysis flow."""
    from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch
    W, H, YB = COL2, 2.75, 0.10  # drawing coordinates in inches; bottom YB cropped
    fig = plt.figure(figsize=(W, H - YB))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W); ax.set_ylim(YB, H); ax.axis("off")
    C1, C2 = C_ORANGE, C_PURPLE  # two clones

    def box(x0, y0, x1, y1, title, body, ec=INK, fc="white", lw=0.7, ts=7, bs=6.3):
        ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0,rounding_size=0.05",
                                    ec=ec, fc=fc, lw=lw))
        ax.text((x0 + x1) / 2, y1 - 0.07, title, ha="center", va="top", fontsize=ts, fontweight="bold")
        if body:
            ax.text((x0 + x1) / 2, y1 - 0.25, body, ha="center", va="top", fontsize=bs, color=INK2,
                    linespacing=1.3)

    def arrow(p, q, color=INK, lw=0.8, ls="-", cs="arc3"):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=7, color=color, lw=lw,
                                     ls=ls, connectionstyle=cs, shrinkA=0, shrinkB=0))

    def cell(x, y, c, r=0.065, later=False):
        ax.add_patch(Circle((x, y), r, fc=c if not later else "white", ec=c, lw=1.2 if later else 0.5))
        if later:
            ax.add_patch(Circle((x, y), r * 0.45, fc=c, ec="none"))

    # ---- (a) lineage ground truth
    ax.text(0.04, H - 0.06, "a", fontsize=9, fontweight="bold", va="top")
    ax.text(1.0, H - 0.08, "Lineage ground truth", ha="center", va="top", fontsize=7, fontweight="bold")
    cell(0.65, 2.30, C1, 0.08); cell(1.35, 2.30, C2, 0.08)
    ax.text(1.0, 2.30, "barcoded\nprogenitors", ha="center", va="center", fontsize=5.8, color=INK2)
    for xc, c in ((0.65, C1), (1.35, C2)):
        arrow((xc, 2.20), (xc, 1.98), color=INK2, lw=0.6)
        for dx in (-0.21, -0.07, 0.07, 0.21):
            cell(xc + dx, 1.90, c, 0.055)
    ax.text(1.0, 2.09, "3–4 divisions", ha="center", va="center", fontsize=5.8, color=INK2)
    ax.plot([1.0, 1.0], [1.70, 0.62], color=INK2, lw=0.6, ls=(0, (2, 2)))
    ax.text(1.0, 1.72, "split", ha="center", va="bottom", fontsize=5.8, color=INK2,
            bbox=dict(fc="white", ec="none", pad=0.5))
    for x0, lab in ((0.12, "profiled at $t$"), (1.08, "profiled at $t{+}1$")):
        ax.add_patch(FancyBboxPatch((x0, 0.68), 0.80, 0.86, boxstyle="round,pad=0,rounding_size=0.05",
                                    ec=INK2, fc="#F7F7F7", lw=0.5))
        ax.text(x0 + 0.40, 1.57, lab, ha="center", va="bottom", fontsize=6)
    for (x, y, c) in ((0.32, 1.12, C1), (0.55, 0.88, C1), (0.70, 1.20, C2), (0.38, 0.84, C2)):
        cell(x, y, c)
    for (x, y, c) in ((1.30, 1.15, C1), (1.62, 0.90, C1), (1.70, 1.22, C2), (1.36, 0.86, C2)):
        cell(x, y, c, later=True)
    arrow((0.40, 1.16), (1.23, 1.17), color=C1, lw=0.7, cs="arc3,rad=-0.25")
    ax.text(0.82, 1.38, "sisters", ha="center", va="bottom", fontsize=5.8, color=C1)
    ax.text(1.0, 0.55, "I→II, II→III, III→IV\n3 biological replicates each", ha="center", va="top",
            fontsize=6, color=INK2, linespacing=1.3)

    # ---- (b) analysis flow
    X0, X1 = 2.15, 7.10
    ax.text(2.05, H - 0.06, "b", fontsize=9, fontweight="bold", va="top")
    yb0, yb1 = 1.98, 2.62
    box(X0 + 0.10, yb0, 3.35, yb1, "scRNA-seq", "both stages,\nper replicate")
    box(3.60, yb0, 5.05, yb1, "Expression space", "CP10k, log(1+x), 2,000 HVGs,\n30 PCs; panel scores Pg, D")
    box(5.30, yb0, X1, yb1, "Couplings $W$", "entropic OT (balanced, unbalanced;\n$\\varepsilon$ grid), kNN, independent",
        ec=C_OT, lw=1.0)
    arrow((3.35, 2.30), (3.60, 2.30)); arrow((5.05, 2.30), (5.30, 2.30))
    gap = 0.10
    wq = (X1 - X0 - 3 * gap) / 4
    xs = [X0 + i * (wq + gap) for i in range(4)]
    yq0, yq1 = 0.86, 1.56
    rq = [("RQ1 Coupling recovery", "own-clone mass vs chance,\nAUROC, top-$k$,\nfate agreement (TVD)"),
          ("RQ2 Forecasting", "train on true / OT /\nindependent pairs;\ntest on held-out true pairs"),
          ("RQ3 Fate prediction", "(P, D), (Pg, D), 30-PC,\nfate-bias axis $\\to$\nclone's later fate"),
          ("RQ4 Mismatch cost", "simulated effect;\nper-cell RMSE, bias\nat mismatch rate $m$")]
    for x, (t, b) in zip(xs, rq):
        box(x, yq0, x + wq, yq1, t, b, ts=6.6, bs=6.0)
    # bus from expression space and couplings to all four questions
    ybus = 1.78
    ax.plot([4.325, 4.325], [yb0, ybus], color=INK, lw=0.8)
    ax.plot([6.20, 6.20], [yb0, ybus], color=C_OT, lw=0.8)
    ax.plot([xs[0] + wq / 2, xs[3] + wq / 2], [ybus, ybus], color=INK, lw=0.8)
    for x in xs:
        arrow((x + wq / 2, ybus), (x + wq / 2, yq1))
    # ground truth into the questions
    arrow((1.88, 1.11), (X0, 1.11), color=C_TRUE, lw=1.0)
    ax.text((1.88 + X0) / 2, 1.16, "sister\nlinks", ha="center", va="bottom", fontsize=5.8)
    # RQ1 measured mismatch feeds RQ4
    yarc = 0.74
    ax.plot([xs[0] + wq / 2, xs[0] + wq / 2, xs[3] + wq / 2], [yq0, yarc, yarc], color=C_OT, lw=0.8)
    arrow((xs[3] + wq / 2, yarc), (xs[3] + wq / 2, yq0), color=C_OT)
    ax.text((xs[1] + xs[2] + wq) / 2, yarc - 0.03, "measured mismatch rate $m$ = 1 $-$ own-clone mass",
            ha="center", va="top", fontsize=5.8, color=C_OT)
    # protocol strip
    ax.add_patch(FancyBboxPatch((X0, 0.16), X1 - X0, 0.24, boxstyle="round,pad=0,rounding_size=0.04",
                                ec="none", fc=BAND))
    ax.text((X0 + X1) / 2, 0.28, "Protocol: clone-level splits (70/15/15)  ·  training-only scaling and PCA  ·  "
            "20 seeds  ·  paired intervals, Holm  ·  $\\delta$ fixed in advance",
            ha="center", va="center", fontsize=6)
    save(fig, "fig1_overview.pdf")


def fig2(fw):
    p = fw[(fw.cellset == "primary") & ~fw.uncapped & (fw.rep == "pooled")]
    fig, axes = plt.subplots(2, 3, figsize=(COL2, 3.15), sharex=True, sharey="row",
                             gridspec_kw=dict(hspace=0.22, wspace=0.12))
    fig.subplots_adjust(left=0.075, right=0.995, top=0.86, bottom=0.11)
    rows = (("auroc", "auroc_lo", "auroc_hi", "Clone-retrieval AUROC"),
            ("enrichment", "enrichment_lo", "enrichment_hi", "Own-clone mass / chance"))
    for j, t in enumerate(TR):
        q = p[p.transition == t]
        for i, (met, lo, hi, ylab) in enumerate(rows):
            ax = axes[i, j]
            ref = 1.0 if met == "enrichment" else 0.5
            ax.axhline(ref, color=C_IND, lw=0.8, ls=(0, (1, 1.5)), zorder=1,
                       label="Independent (= chance)")
            ax.axhline(one(q, coupling="kNN")[met], color=C_KNN, lw=0.9, ls=(0, (4, 2)), zorder=1,
                       label="kNN (k = 10)")
            for kind, col, mk, mfc, nm, dx in (("bal", C_OT, "o", C_OT, "OT, balanced", 0.95),
                                               ("unb", C_UNB, "s", "white", "OT, unbalanced (\u03c4 = 0.9)", 1.05)):
                r = [one(q, coupling=f"OT_{kind}_eps{e}") for e in EPS]
                y = np.array([v[met] for v in r])
                err = np.array([[v[met] - v[lo] for v in r], [v[hi] - v[met] for v in r]])
                ax.errorbar(np.array(EPS) * dx, y, yerr=err, color=col, marker=mk, ms=3.6, mfc=mfc,
                            mec=col, mew=0.9, lw=1.0, elinewidth=0.7, zorder=3, label=nm)
            ax.set_xscale("log")
            ax.set_xticks(EPS)
            ax.set_xticklabels(["0.005", "0.01", "0.05", "0.1"])
            ax.minorticks_off()
            ax.set_xlim(0.0036, 0.14)
            if met == "enrichment":
                ax.set_yscale("log")
                ax.set_yticks([0.1, 0.3, 1, 3, 10, 30])
                ax.set_yticklabels(["0.1", "0.3", "1", "3", "10", "30"])
                ax.minorticks_off()
                ax.set_xlabel("Entropy \u03b5")
            else:
                ax.set_ylim(0.44, 0.68)
                ax.set_title(TPLOT[t], pad=3)
            if j == 0:
                ax.set_ylabel(ylab)
            panel(ax, "abcdef"[i * 3 + j], x=-0.2 if j == 0 else -0.08)
    h, l = axes[0, 0].get_legend_handles_labels()
    order = [2, 3, 1, 0]
    fig.legend([h[k] for k in order], [l[k] for k in order], loc="upper center", ncol=4,
               bbox_to_anchor=(0.53, 1.0))
    save(fig, "fig2_coupling_recovery.pdf")


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

    fig3(p)


def fig3(p):
    fig, axes = plt.subplots(1, 3, figsize=(COL2, 2.15), sharex=True,
                             gridspec_kw=dict(wspace=0.55))
    fig.subplots_adjust(left=0.115, right=0.99, top=0.84, bottom=0.17)
    for k, (ax, t) in enumerate(zip(axes, TR)):
        q = p[p.transition == t]
        cl = sorted(q.early_cluster.unique(), key=lambda c: -one(q, early_cluster=c, coupling=OT).n_clones)
        y = np.arange(len(cl))[::-1].astype(float)
        o = [one(q, early_cluster=c, coupling=OT) for c in cl]
        ind = [one(q, early_cluster=c, coupling="independent") for c in cl]
        kn = [one(q, early_cluster=c, coupling="kNN") for c in cl]
        for yy, a, b in zip(y, ind, o):   # dumbbell: independent -> OT
            ax.plot([a.tvd, b.tvd], [yy, yy], color="#BDBDBD", lw=1.6, zorder=1, solid_capstyle="round")
        ax.errorbar([x.tvd for x in o], y, xerr=[[x.tvd - x.tvd_lo for x in o], [x.tvd_hi - x.tvd for x in o]],
                    fmt="o", ms=3.8, color=C_OT, elinewidth=0.7, zorder=4, label="OT, balanced, \u03b5 = 0.05")
        ax.plot([x.tvd for x in ind], y, "s", ms=3.4, mfc="white", mec=C_IND, mew=0.9, zorder=3,
                label="Independent")
        ax.plot([x.tvd for x in kn], y, "^", ms=3.4, mfc="white", mec=C_KNN, mew=0.9, zorder=3,
                label="kNN (k = 10)")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{c.split('-', 1)[1]} ({int(x.n_clones)})" for c, x in zip(cl, o)])
        for tl, c in zip(ax.get_yticklabels(), cl):
            if c.endswith("NPvMB"):
                tl.set_fontweight("bold")
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.set_ylim(-0.6, len(cl) - 0.4)
        ax.set_xlim(0, 0.85)
        ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
        ax.set_title(f"{TPLOT[t]}  (early cluster, clones)", pad=3)
        ax.set_xlabel("TVD to clonal fate distribution")
        panel(ax, "abc"[k], x=-0.36)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend([h[2], h[1], h[0]], [l[2], l[1], l[0]], loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.0))
    save(fig, "fig3_fate_agreement.pdf")


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
    fig, axes = plt.subplots(1, 3, figsize=(COL2, 2.0), sharey=True, sharex=True,
                             gridspec_kw=dict(wspace=0.08))
    fig.subplots_adjust(left=0.1, right=0.99, top=0.83, bottom=0.19)
    ccol = {"true": C_TRUE, "OT": C_OT, "independent": C_IND}
    cmk = {"true": "o", "OT": "s", "independent": "D"}
    cnm = {"true": "True lineage", "OT": "OT (\u03b5 = 0.05)", "independent": "Independent"}
    for k, (ax, t) in enumerate(zip(axes, TR)):
        q = pc[(pc.transition == t) & (pc.reference == "ConstMedian")]
        ax.axvspan(-DELTA * 1000, DELTA * 1000, color=BAND, zorder=0, lw=0)
        ax.axvline(0, color=INK, lw=0.6, zorder=1)
        yt, yl = [], []
        for gi, m in enumerate(MODELS):
            for ci, cp in enumerate(CPL):
                y = -(gi * 4 + ci)
                r = one(q, model=m, coupling=cp)
                ax.errorbar(1000 * r.delta, y, xerr=[[1000 * (r.delta - r.ci_lo)], [1000 * (r.ci_hi - r.delta)]],
                            fmt=cmk[cp], color=ccol[cp], ms=3.4, mfc=ccol[cp] if cp != "independent" else "white",
                            mew=0.9, elinewidth=0.8, zorder=3, label=cnm[cp] if gi == 0 else None)
            yt.append(-(gi * 4 + 1))
            yl.append(MLAB[m].replace("\\ ", " "))
        ax.set_yticks(yt)
        ax.set_yticklabels(yl)
        ax.tick_params(axis="y", length=0)
        ax.set_ylim(-11, 1)
        ax.set_xlim(-25, 10)
        ax.set_xticks([-20, -10, 0, 10])
        ax.set_title(TPLOT[t], pad=3)
        ax.set_xlabel(r"$\Delta$MAE vs. constant median ($\times10^{-3}$)")
        if k == 0:
            ax.text(-9.4, -10.8, "\u00b1\u03b4", ha="left", va="bottom", fontsize=6.5, color=INK2)
        panel(ax, "abc"[k], x=-0.27 if k == 0 else -0.05)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, ["Training pairs: " + l[0]] + l[1:], loc="upper center", ncol=3,
               bbox_to_anchor=(0.55, 1.0))
    save(fig, "fig4_forecasting.pdf")


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
    mean_auc = p.groupby(["transition", "representation"]).macro_auroc.mean()
    for r in REPS:
        v = mean_auc.xs(r, level="representation")
        put(f"rqThreeAUC{RMAC[r]}MinT", f(v.min()))
        put(f"rqThreeAUC{RMAC[r]}MaxT", f(v.max()))
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
    fig, ax = plt.subplots(figsize=(COL1, 2.05))
    fig.subplots_adjust(left=0.15, right=0.99, top=0.84, bottom=0.12)
    cols = {"P_tcbb_D": C_ORANGE, "Pg_D": C_UNB, "PC30": C_KNN, "fate_axis": C_PURPLE}
    mks = {"P_tcbb_D": "v", "Pg_D": "s", "PC30": "o", "fate_axis": "D"}
    rng = np.random.default_rng(0)  # horizontal jitter only
    for j, t in enumerate(TR):
        for k, r in enumerate(REPS):
            g = p[(p.transition == t) & (p.representation == r)].macro_auroc.to_numpy()
            x = j + (k - 1.5) * 0.19
            ax.scatter(x + rng.uniform(-0.045, 0.045, len(g)), g, s=3, color=cols[r], alpha=0.45,
                       lw=0, zorder=2)
            ax.errorbar(x, g.mean(), yerr=g.std(ddof=1), fmt=mks[r], color=cols[r], ms=3.8,
                        mec=INK, mew=0.5, elinewidth=0.9, capsize=1.5, zorder=3,
                        label=RLAB[r] if j == 0 else None)
    ax.axhline(0.5, color=C_IND, ls=(0, (1, 1.5)), lw=0.8, zorder=1)
    ax.text(2.48, 0.505, "chance", ha="right", va="bottom", fontsize=6, color=INK2)
    ax.set_xticks(range(3))
    ax.set_xticklabels([TPLOT[t] for t in TR])
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylabel("Fate macro-AUROC")
    ax.set_ylim(0.45, 0.92)
    ax.legend(loc="upper center", ncol=4, bbox_to_anchor=(0.45, 1.2), handletextpad=0.2,
              columnspacing=0.8)
    save(fig, "fig5_state_adequacy.pdf")


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
    # (a) scalar mismatch (random wrong partners); (b) structured partners (A4.1(b), linear effect)
    st = pd.read_csv(R / "21_rq4_structured/summary.csv")
    st = st[st.effect == "linear"]
    fig, axes = plt.subplots(1, 2, figsize=(COL1, 1.85), sharey=True,
                             gridspec_kw=dict(wspace=0.12, width_ratios=[1, 1.25]))
    fig.subplots_adjust(left=0.13, right=0.99, top=0.8, bottom=0.30)
    cols = {1: C_OT, 2: C_UNB, 3: C_KNN}
    mks = {1: "o", 2: "s", 3: "^"}
    ax = axes[0]
    for rep in (1, 2, 3):
        q = g[g.rep == rep].sort_values("m")
        ax.fill_between(q.m, q.rmse_lo, q.rmse_hi, color=cols[rep], alpha=0.2, lw=0)
        ax.plot(q.m, q.rmse, color=cols[rep], marker=mks[rep], ms=2.6, mfc="white", mew=0.8,
                lw=0.9, label=f"Replicate {rep}")
        mm = est[(est.rep == rep) & (est.rate_kind == "1-own_clone_mass")]
        ax.plot(mm.m, mm.rmse, ls="none", marker="*", ms=6.5, color=INK, mew=0, zorder=5,
                label="RQ1-measured OT rate" if rep == 1 else None)
    ax.set_xlabel("Mismatch rate m (random partners)")
    ax.set_ylabel("Per-cell effect RMSE")
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(["0", "0.5", "1"])
    panel(ax, "a", x=-0.38)
    ax = axes[1]
    order = ["random", "kNN_sampled", "OT_sampled", "kNN_averaged", "OT_barycentric"]
    labs = ["Random", "kNN\nsamp.", "OT\nsamp.", "kNN\navg.", "OT\nbary."]
    for rep, dx in ((1, -0.18), (2, 0.0), (3, 0.18)):
        for i, pt in enumerate(order):
            r = one(st, rep=rep, partner=pt)
            ax.errorbar(i + dx, r.rmse, yerr=[[r.rmse - r.rmse_lo], [r.rmse_hi - r.rmse]], color=cols[rep],
                        marker=mks[rep], ms=2.6, mfc="white", mew=0.8, lw=0.8, capsize=0)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(labs, fontsize=6)
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_xlabel("Partner from control half", labelpad=1)
    panel(ax, "b", x=-0.08)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend([h[0], h[2], h[3], h[1]], [l[0], l[2], l[3], l[1]], loc="upper center", ncol=4, bbox_to_anchor=(0.53, 1.0), handlelength=1.2,
               columnspacing=0.7)
    save(fig, "fig6_perturbation.pdf")


# ============================================================ post-review (amendment A3)
# Exploratory additions made after peer review; settings fixed in ANALYSIS_PLAN.md A3
# (2026-11-17) before any was run. Main text cites them as post-review additions; full
# results go to supplementary Tables S6-S13.
LEARNED = ["Ridge", "RandomForest", "Transformer"]
MLAB = {"ConstMedian": "Const.\\ median", "ConstMean": "Const.\\ mean", "Ridge": "Ridge",
        "RandomForest": "Random forest", "Transformer": "Transformer"}
CLAB = {"true": "True", "OT": "OT", "independent": "Indep."}
REPLAB = {"P_tcbb_D": "(P, D)", "Pg_D": "(Pg, D)", "PC30": "30-PC", "fate_axis": "Fate-bias axis"}
EFFLAB = {"linear": "Linear (primary)", "nonlinear": "Nonlinear", "multi_axis": "Multi-axis",
          "reversed": "Reversed", "sparse": "Sparse (20 genes)"}
INPUTS += [R / p for p in (
    "11_single_barcode/summary.csv", "12_rq2_pc/summary.csv", "12_rq2_pc/per_seed.csv",
    "13_rq2_draws/summary.csv", "14_loro/rq2_per_fold.csv", "14_loro/rq3_per_fold.csv",
    "14_loro/rq3_summary.csv", "15_fate_labels/summary.csv",
    "15_fate_labels/pool_label_shares.csv", "16_rq4_sens/summary.csv",
    "17_rq2_inductive/comparisons.csv")]


def ci(x, lo, hi, d=4):
    return f"{f(x, d)} [{f(lo, d)}, {f(hi, d)}]"


def post_review():
    # ---- A3.1 RQ2 with 30-PC input
    s = pd.read_csv(R / "12_rq2_pc/summary.csv")
    pc = s[(s.comparison == "PC30 - Pg_D") & s.model.isin(LEARNED)]
    tt = pc[(pc.coupling == "true") & (pc.transition == "II->III")]
    assert (tt.ci_hi < -DELTA).all()  # beyond the margin for every learned model
    put("pcTrueTwoGainMinPct", f(-tt.delta_pct_of_const.max(), 1))
    put("pcTrueTwoGainMaxPct", f(-tt.delta_pct_of_const.min(), 1))
    put("pcTrueTwoDeltaMin", f(-tt.delta.max(), 3))
    put("pcTrueTwoDeltaMax", f(-tt.delta.min(), 3))
    to = pc[(pc.coupling == "true") & (pc.transition != "II->III")]
    assert to.equivalent.all()
    put("pcTrueOtherMaxPct", f(-to.delta_pct_of_const.min(), 1))
    oi = pc[pc.coupling != "true"]
    assert oi.equivalent.all()
    put("pcOIMaxGainPct", f(-oi.delta_pct_of_const.min(), 1))
    put("pcOTThreeWorseMinPct", f(pc[(pc.coupling == "OT") & (pc.transition == "III->IV")].delta_pct_of_const.min(), 1))
    put("pcOTThreeWorseMaxPct", f(pc[(pc.coupling == "OT") & (pc.transition == "III->IV")].delta_pct_of_const.max(), 1))
    vc = s[(s.comparison == "PC30: model - ConstMedian") & s.model.isin(LEARNED)
           & (s.coupling == "true") & (s.transition == "II->III")]
    put("pcTrueTwoVsConstMinPct", f(-vc.delta_pct_of_const.max(), 1))
    put("pcTrueTwoVsConstMaxPct", f(-vc.delta_pct_of_const.min(), 1))
    # descriptive (not a pre-specified contrast): true minus OT training pairs, same input
    ps = pd.read_csv(R / "12_rq2_pc/per_seed.csv")
    w = ps[ps.model.isin(LEARNED)].pivot_table(index=["transition", "input", "model", "seed"],
                                               columns="coupling", values="mae").reset_index()
    w["d"] = w["true"] - w["OT"]
    md = w.groupby(["transition", "input", "model"]).d.mean()
    put("pcTrueMinusOTTwoMin", f(-md.loc[("II->III", "PC30")].max(), 3))
    put("pcTrueMinusOTTwoMax", f(-md.loc[("II->III", "PC30")].min(), 3))
    put("pgTrueMinusOTTwoMaxAbs", f(md.loc[("II->III", "Pg_D")].abs().max(), 3))
    rows = []
    for t in TR:
        for c in ["true", "OT", "independent"]:
            for i, m in enumerate(LEARNED):
                r = one(pc, transition=t, coupling=c, model=m)
                v = one(s, transition=t, coupling=c, model=m, comparison="PC30: model - ConstMedian")
                rows.append(f"{TLAB[t] if (c == 'true' and i == 0) else ''} & {CLAB[c] if i == 0 else ''} & "
                            f"{MLAB[m]} & {f(r.mae_Pg_D, 4)} & {f(r.mae_PC30, 4)} & "
                            f"{ci(r.delta, r.ci_lo, r.ci_hi)} & {f(r.delta_pct_of_const, 1)} & "
                            f"{f(v.delta_pct_of_const, 1)}" + r" \\")
        rows.append(r"\addlinespace")
    write_table("tab_s6_rq2_pc", rf"""\begin{{table*}}[!t]
\centering
\caption{{Post-review addition: RQ2 with a 30-PC input (primary set, 20 splits). Splits,
sampled cells and training targets are those of the main analysis; only the input changes.
$\Delta$ = MAE(30-PC) $-$ MAE((Pg, D)), mean over splits with 95\% split-stability interval;
\%: as a percentage of the constant-median MAE. Last column: 30-PC model minus constant
median, \% of the constant-median MAE. Equivalence margin $\delta=0.01$.}}
\label{{tab:s6}}
\footnotesize
\begin{{tabular}}{{lllccccc}}
\toprule
Transition & Pairs & Model & MAE (Pg, D) & MAE 30-PC & $\Delta$ [95\% interval] & $\Delta$ (\%) & vs.\ const.\ (\%) \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")

    # ---- A3.2 inductive OT arm
    c = pd.read_csv(R / "17_rq2_inductive/comparisons.csv")
    tr = c[c.reference == "transductive_OT_same_model"]
    assert ((tr.ci_lo > -DELTA) & (tr.ci_hi < DELTA)).all()
    put("indMaxAbsDelta", f(tr.delta.abs().max(), 4))
    put("indMaxAbsPct", f(tr.delta_pct_of_const.abs().max(), 1))
    put("indSeeds", int(tr.n_seeds.max()))
    rows = []
    for t in TR:
        for i, m in enumerate(MLAB):
            r = one(tr, transition=t, model=m)
            vc_ = c[(c.transition == t) & (c.model == m) & (c.reference == "ConstMedian")]
            vtxt = f(vc_.iloc[0].delta_pct_of_const, 1) if len(vc_) else "--"
            rows.append(f"{TLAB[t] if i == 0 else ''} & {MLAB[m]} & {f(r.mae_model, 4)} & "
                        f"{f(r.mae_reference, 4)} & {ci(r.delta, r.ci_lo, r.ci_hi)} & "
                        f"{f(r.delta_pct_of_const, 1)} & {vtxt}" + r" \\")
        rows.append(r"\addlinespace")
    write_table("tab_s7_rq2_inductive", rf"""\begin{{table*}}[!t]
\centering
\caption{{Post-review addition: inductive OT training pairs (primary set, splits 0--{int(tr.n_seeds.max()) - 1}).
Gene selection, PCA, cost scale and OT plan are refitted per split without validation and
test clones; validation cells are projected and their rows use the fitted later-cell
potential. $\Delta$ = inductive $-$ transductive MAE for the same model and splits, with 95\%
split-stability interval; \%: of the constant-median MAE. Last column: inductive model minus
inductive constant median (\%).}}
\label{{tab:s7}}
\footnotesize
\begin{{tabular}}{{llccccc}}
\toprule
Transition & Model & MAE inductive & MAE transductive & $\Delta$ [95\% interval] & $\Delta$ (\%) & vs.\ const.\ (\%) \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")

    # ---- A3.3 leave-one-replicate-out
    q = pd.read_csv(R / "14_loro/rq3_per_fold.csv")
    qw = q.pivot_table(index=["transition", "heldout_rep"], columns="representation", values="macro_auroc")
    put("loroNFolds", len(qw))
    put("loroNPCAbove", int((qw.PC30 > qw.Pg_D).sum()))
    put("loroNFateAbove", int((qw.fate_axis > qw.Pg_D).sum()))
    qs = pd.read_csv(R / "14_loro/rq3_summary.csv")
    for t in TR:
        for rep, nm in (("Pg_D", "PgD"), ("PC30", "PC"), ("fate_axis", "Fate")):
            put(f"loro{nm}{TN[t]}", f(one(qs, transition=t, representation=rep).auroc_mean))
    l2 = pd.read_csv(R / "14_loro/rq2_per_fold.csv")
    lt = l2[(l2.coupling == "true") & l2.model.isin(LEARNED)]
    put("loroTrueMinPct", f(lt.delta_pct_of_const.min(), 1))
    put("loroTrueMaxPct", f(lt.delta_pct_of_const.max(), 1))
    put("loroTrueNBeyond", int((lt.delta_vs_const_median < -DELTA).sum()))
    put("loroTrueN", len(lt))
    rows = []
    for t in TR:
        for i, rep in enumerate(REPLAB):
            r = one(qs, transition=t, representation=rep)
            rows.append(f"{TLAB[t] if i == 0 else ''} & {REPLAB[rep]} & {f(r.auroc_mean)} "
                        f"[{f(r.auroc_min)}, {f(r.auroc_max)}] & {f(r.bacc_mean)} "
                        f"[{f(r.bacc_min)}, {f(r.bacc_max)}]" + r" \\")
        rows.append(r"\addlinespace")
    rows2 = []
    for t in TR:
        for cpl in ["true", "OT", "independent"]:
            for i, m in enumerate(LEARNED):
                g = l2[(l2.transition == t) & (l2.coupling == cpl) & (l2.model == m)].sort_values("heldout_rep")
                rows2.append(f"{TLAB[t] if (cpl == 'true' and i == 0) else ''} & {CLAB[cpl] if i == 0 else ''} & {MLAB[m]} & "
                             + " & ".join(f(x, 1) for x in g.delta_pct_of_const) + r" \\")
        rows2.append(r"\addlinespace")
    write_table("tab_s8_loro", rf"""\begin{{table}}[!t]
\centering
\caption{{Post-review addition: leave-one-replicate-out RQ3 (primary set; descriptive).
Train on two replicates, test on the third; mean [min, max] over the three held-out
replicates.}}
\label{{tab:s8}}
\footnotesize
\setlength{{\tabcolsep}}{{3pt}}
\begin{{tabular}}{{llcc}}
\toprule
Transition & Representation & Macro-AUROC & Balanced accuracy \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table}}

\begin{{table}}[!t]
\centering
\caption{{Post-review addition: leave-one-replicate-out RQ2 (primary set, (Pg, D) input;
descriptive). Model minus constant median trained on the same pairs, as \% of the
constant-median MAE, for held-out replicates 1, 2 and 3.}}
\label{{tab:s9}}
\footnotesize
\setlength{{\tabcolsep}}{{4pt}}
\begin{{tabular}}{{lllccc}}
\toprule
Transition & Pairs & Model & Rep.~1 & Rep.~2 & Rep.~3 \\
\midrule
{chr(10).join(rows2[:-1])}
\bottomrule
\end{{tabular}}
\end{{table}}
""")

    # ---- A3.4 repeated target draws
    d = pd.read_csv(R / "13_rq2_draws/summary.csv")
    dl = d[d.model.isin(LEARNED)]
    put("drawShareMinPct", f(100 * dl.share_within.min(), 0))
    put("drawShareMaxPct", f(100 * dl.share_within.max(), 0))
    put("drawNDraws", int(d.n_draws.max()))
    e12 = dl[dl.transition != "III->IV"]
    assert e12.equivalent.all()
    put("drawEarlyMaxPct", f(e12.mean_delta_pct_of_const.abs().max(), 1))
    rows = []
    for t in TR:
        for cpl in ["OT", "independent"]:
            for i, m in enumerate(["ConstMean"] + LEARNED[:2]):
                r = one(d, transition=t, coupling=cpl, model=m)
                rows.append(f"{TLAB[t] if (cpl == 'OT' and i == 0) else ''} & {CLAB[cpl] if i == 0 else ''} & {MLAB[m]} & "
                            f"{ci(r.mean_delta, r.ci_lo, r.ci_hi)} & {f(r.mean_delta_pct_of_const, 1)} & "
                            f"{pct(r.share_within, 0)} & {pct(r.frac_draws_model_better, 0)}" + r" \\")
        rows.append(r"\addlinespace")
    write_table("tab_s10_rq2_draws", rf"""\begin{{table*}}[!t]
\centering
\caption{{Post-review addition: repeated training-target draws (primary set; 20 splits
$\times$ {int(d.n_draws.max())} draws; (Pg, D) input). $\Delta$ = model minus constant median,
mean over splits and draws, with 95\% interval; \%: of the constant-median MAE. Within-split
share: fraction of the variance of $\Delta$ due to the target draw. Last column: share of
draws in which the model beat the constant.}}
\label{{tab:s10}}
\footnotesize
\begin{{tabular}}{{lllcccc}}
\toprule
Transition & Pairs & Model & $\Delta$ [95\% interval] & $\Delta$ (\%) & Within-split share (\%) & Model better (\%) \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")

    # ---- A3.5 fate agreement on measured labels
    fl = pd.read_csv(R / "15_fate_labels/summary.csv")
    sh = pd.read_csv(R / "15_fate_labels/pool_label_shares.csv")
    for t in TR:
        put(f"fateTransferredPct{TN[t]}", pct(one(sh, transition=t, rep="pooled").share_later_transferred, 0))
    otk = fl[(fl.coupling == OT) & (fl.labels == "k15")]
    o3 = otk[otk.transition == "III->IV"].mass_on_transferred
    put("fateOTMassTransThreeMinPct", pct(o3.min(), 0))
    put("fateOTMassTransThreeMaxPct", pct(o3.max(), 0))
    kk = fl[(fl.coupling == OT) & fl.labels.isin(["k5", "k15", "k30"])]
    put("fateKRangeMaxOT", f(kk.groupby(["transition", "early_cluster"]).tvd.agg(lambda x: x.max() - x.min()).max()))
    pw = fl.pivot_table(index=["transition", "early_cluster", "labels"], columns="coupling", values="tvd").reset_index()
    pw["better"] = pw[OT] < pw["independent"]
    for t in TR:
        mm = pw[(pw.transition == t) & (pw.labels == "measured")]
        kq = pw[(pw.transition == t) & (pw.labels == "k15")]
        put(f"fateMeasNBetter{TN[t]}", int(mm.better.sum()))
        put(f"fateMeasNClusters{TN[t]}", len(mm))
        assert int(kq.better.sum()) == int(M[f"fateNBetter{TN[t]}"])  # k15 reproduces RQ1
    put("fateMeasNBetterAll", sum(int(M[f"fateMeasNBetter{TN[t]}"]) for t in TR))
    # I->II clusters that switch from OT-better (all labels) to OT-worse (measured labels)
    k1 = pw[(pw.transition == "I->II") & (pw.labels == "k15")].set_index("early_cluster")
    m1 = pw[(pw.transition == "I->II") & (pw.labels == "measured")].set_index("early_cluster")
    sw = [cl for cl in k1.index if k1.better[cl] and not m1.better[cl]]
    assert sw == ["I-NPpHB", "I-VLMCs"], sw
    put("fateMeasSwitchMaxDiffOne", f(max(m1[OT][cl] - m1["independent"][cl] for cl in sw), 3))
    for cl in ["III-GABAs", "III-Gluts", "III-NPvMB", "III-NPaHB", "III-NPpHB", "III-Motors"]:
        nm = cl.split("-")[1]
        for lab, tag in (("measured", "Meas"), ("k15", "Trans")):
            r = pw[(pw.early_cluster == cl) & (pw.labels == lab)].iloc[0]
            put(f"tvd{tag}OT{nm}", f(r[OT]))
            put(f"tvd{tag}Ind{nm}", f(r["independent"]))
    rows = []
    for t in TR:
        cls = sorted(pw[pw.transition == t].early_cluster.unique())
        for i, cl in enumerate(cls):
            a = pw[(pw.early_cluster == cl) & (pw.labels == "k15")].iloc[0]
            b = pw[(pw.early_cluster == cl) & (pw.labels == "measured")].iloc[0]
            o = one(fl, early_cluster=cl, coupling=OT, labels="k15")
            kr = kk[kk.early_cluster == cl].tvd
            rows.append(f"{TLAB[t] if i == 0 else ''} & {cl} ({int(o.n_clones)}) & {f(a[OT])} & {f(a['independent'])} & "
                        f"{f(b[OT])} & {f(b['independent'])} & {pct(o.mass_on_transferred, 0)} & {f(kr.max() - kr.min())}" + r" \\")
        rows.append(r"\addlinespace")
    write_table("tab_s11_fate_labels", rf"""\begin{{table*}}[!t]
\centering
\caption{{Post-review addition: fate agreement (TVD, lower is better) with all later cells
(authors' labels plus $k=15$ transferred labels; as in the main text) and with the authors'
measured labels only (predicted distribution renormalised over measured cells). OT: balanced,
$\varepsilon=0.05$. Transferred mass: share of the cluster's OT mass on cells with transferred
labels. $k$ range: largest minus smallest OT TVD over $k\in\{{5,15,30\}}$. Clones in
parentheses. Later-pool share with transferred labels: {M['fateTransferredPctOne']}\%,
{M['fateTransferredPctTwo']}\% and {M['fateTransferredPctThree']}\% for the three transitions.}}
\label{{tab:s11}}
\footnotesize
\begin{{tabular}}{{llcccccc}}
\toprule
 & & \multicolumn{{2}}{{c}}{{All labels}} & \multicolumn{{2}}{{c}}{{Measured labels}} & & \\
\cmidrule(lr){{3-4}}\cmidrule(lr){{5-6}}
Transition & Early cluster & OT & Indep. & OT & Indep. & Transferred mass (\%) & $k$ range \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")

    # ---- A3.6 RQ4 effect models
    e = pd.read_csv(R / "16_rq4_sens/summary.csv")
    em = e[(e.point == "rq1_measured") & (e.bound == "estimate") & (e.rate_kind == "1-own_clone_mass")]
    dense = em[em.effect_model != "sparse"]
    spread = dense.groupby(["transition", "rep"]).rmse.agg(lambda x: x.max() - x.min()).max()
    assert spread < 1e-9, spread  # identical RMSE for every dense effect shape
    put("effRMSEMin", f(dense.rmse.min()))
    put("effRMSEMax", f(dense.rmse.max()))
    put("effRMin", f(dense.r.min()))
    put("effRMax", f(dense.r.max()))
    sp = em[em.effect_model == "sparse"]
    put("effSparseRMSEMin", f(sp.rmse.min()))
    put("effSparseRMSEMax", f(sp.rmse.max()))
    put("effSparseRMin", f(sp.r.min()))
    put("effSparseRMax", f(sp.r.max()))
    put("effMaxAbsBiasCI", f(e[["bias_lo", "bias_hi"]].abs().max().max(), 3))
    g = e[e.point == "grid"]
    rows = []
    for k, lab in EFFLAB.items():
        cells = [lab]
        for rep in (1, 2, 3):
            r5 = one(g, effect_model=k, m=0.5, rep=rep)
            r1 = one(g, effect_model=k, m=1.0, rep=rep)
            cells.append(f"{f(r5.rmse)} / {f(r1.rmse)}")
        x = em[em.effect_model == k]
        cells += [f"{f(x.rmse.min())}--{f(x.rmse.max())}", f"{f(x.r.min())}--{f(x.r.max())}"]
        rows.append(" & ".join(cells) + r" \\")
    write_table("tab_s12_rq4_effects", rf"""\begin{{table*}}[!t]
\centering
\caption{{Post-review addition: RQ4 with other effect shapes (II$\to$III stage-III cells, 100
seeds). Nonlinear: $e_i=0.5+0.25(\tilde z_i^2-1)/\sqrt2$; multi-axis: mean of the first three
standardised PCs, scaled by $1/\sqrt3$; reversed: $e_i=-0.5+0.25\tilde z_i$; sparse: linear
effect on 20 of the 200 genes, estimated over those 20. Per-cell RMSE at $m=0.5$ / $m=1$ by
replicate, and RMSE and Pearson $r$ at the mismatch rates measured in RQ1 (1 $-$ own-clone
mass; all transitions and replicates). Bias: every 95\% interval within $\pm${M['effMaxAbsBiasCI']}.}}
\label{{tab:s12}}
\footnotesize
\begin{{tabular}}{{lccccc}}
\toprule
 & \multicolumn{{3}}{{c}}{{RMSE, $m=0.5$ / $m=1$}} & \multicolumn{{2}}{{c}}{{At measured rates}} \\
\cmidrule(lr){{2-4}}\cmidrule(lr){{5-6}}
Effect & Rep.~1 & Rep.~2 & Rep.~3 & RMSE & $r$ \\
\midrule
{chr(10).join(rows)}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")

    # ---- A3.7 single-barcode clones
    b = pd.read_csv(R / "11_single_barcode/summary.csv")
    bo = b[b.coupling == OT]
    for t in TR:
        a = one(bo, transition=t, subset="all")
        assert abs(a.auroc - float(M[f"rqOneAUROC{TN[t]}"].replace(r"\ensuremath{-}", "-"))) < 0.006
        sb = one(bo, transition=t, subset="single_barcode")
        put(f"sbAUROC{TN[t]}", f(sb.auroc))
        put(f"sbNClones{TN[t]}", int(sb.n_clones))
        put(f"sbPctMulti{TN[t]}", f(a.pct_cells_multi, 1))
    put("sbMaxAbsAUROCDiff", f(max(abs(one(bo, transition=t, subset="single_barcode").auroc
                                      - one(bo, transition=t, subset="all").auroc) for t in TR), 3))
    rows = []
    for t in TR:
        for i, cpl in enumerate([OT, "kNN", "independent"]):
            for j, sub in enumerate(["all", "single_barcode"]):
                r = one(b, transition=t, coupling=cpl, subset=sub)
                rows.append(f"{TLAB[t] if (i == 0 and j == 0) else ''} & "
                            f"{({OT: 'OT', 'kNN': 'kNN', 'independent': 'Indep.'})[cpl] if j == 0 else ''} & "
                            f"{'All' if sub == 'all' else 'Single'} & {int(r.n_clones)} & "
                            f"{f(r.auroc)} [{f(r.auroc_lo)}, {f(r.auroc_hi)}] & "
                            f"{f(r.enrichment, 1)} [{f(r.enrichment_lo, 1)}, {f(r.enrichment_hi, 1)}] & "
                            f"{pct(r.top10, 1)}" + r" \\")
        rows.append(r"\addlinespace")
    write_table("tab_s13_single_barcode", rf"""\begin{{table*}}[!t]
\centering
\caption{{Post-review addition: RQ1 restricted to clones marked by a single viral barcode
(primary set, forward, pooled over replicates; clone-weighted means with 95\% clone-bootstrap
CIs). Cells of cross-stage clones carrying more than one barcode: {M['sbPctMultiOne']}\%,
{M['sbPctMultiTwo']}\% and {M['sbPctMultiThree']}\% for the three transitions.}}
\label{{tab:s13}}
\footnotesize
\begin{{tabular}}{{lllcccc}}
\toprule
Transition & Coupling & Clones & $n$ & AUROC & Enrichment & Top-10 (\%) \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")


# ============================================================ second review round (A4, A5, A6)
# Exploratory additions after the second review; settings fixed in ANALYSIS_PLAN.md A4
# (2026-11-18), with compute reductions A5 (A4.2) and A6 (A4.1(b)) written before those runs.
# Main text cites them as post-review additions; full results in Tables S14-S22.
INPUTS += [R / p for p in (
    "18_fate_unc/label_transfer_accuracy.csv", "18_fate_unc/paired_tvd.csv",
    "19_rq1_baselines/summary.csv", "20_partner_error/summary.csv",
    "21_rq4_structured/summary.csv", "22_rq2_targets/summary.csv",
    "22_rq2_targets/ot_foreground_share.csv", "23_rq1_generality/summary.csv",
    "23_rq1_generality/solves.csv", "24_larry/summary.csv", "24_larry/solves.csv",
    "24_larry/fate_agreement.csv", "25_cellrank/rq3_comparison.csv",
    "25_cellrank/fate_agreement.csv")]
LT = ["d2->d4", "d4->d6"]
LN = dict(zip(LT, ["A", "B"]))
LLAB = {"d2->d4": r"Day 2$\to$4", "d4->d6": r"Day 4$\to$6"}
BASE = [(OT, "OT balanced"), ("OT_unb_eps0.05", "OT unbalanced"), ("Gauss_eps0.05", "Gaussian kernel"),
        ("kNN5", "kNN, $k=5$"), ("kNN10", "kNN, $k=10$"), ("kNN20", "kNN, $k=20$"),
        ("kNN50", "kNN, $k=50$"), ("MNN10", "MNN, $k=10$"), ("independent", "Independent")]
PARTNERS = [("random", "Random"), ("kNN_sampled", "kNN, sampled"), ("OT_sampled", "OT, sampled"),
            ("kNN_averaged", "kNN, averaged"), ("OT_barycentric", "OT, barycentric")]


def sg(x, d=2):
    """Signed number: explicit plus sign for positive values."""
    s = f(x, d)
    return "$+$" + s if (x > 0 and float(f"{x:.{d}f}") != 0) else s


def rng_(x, d=2, scale=1.0):
    """'a--b' (or 'a' when equal at this precision) for a series."""
    lo, hi = f(scale * min(x), d), f(scale * max(x), d)
    return lo if lo == hi else f"{lo}--{hi}"


def second_review():
    # ---- A4.4 proximity baselines (SISBAR)
    b = pd.read_csv(R / "19_rq1_baselines/summary.csv")
    bp = b[b.rep == "pooled"].copy()
    bp["apr"] = bp.ap / bp.ap_chance
    ga = bp[bp.coupling.str.startswith("Gauss")]
    assert (ga.groupby("transition").auroc.agg(lambda x: x.max() - x.min()) < 1e-6).all()  # rank = distance (float ties)
    g5 = ga[ga.coupling == "Gauss_eps0.05"]
    put("gaussAUROCMin", f(g5.auroc.min()))
    put("gaussAUROCMax", f(g5.auroc.max()))
    for t in TR:
        put(f"gaussAUROC{TN[t]}", f(one(g5, transition=t).auroc))
    o = bp[bp.coupling == OT]
    for t in TR:  # reproduces RQ1
        assert abs(one(o, transition=t).auroc - float(M[f"rqOneAUROC{TN[t]}"].replace(r"\ensuremath{-}", "-"))) < 0.006
    put("apOTRatioMin", f(o.apr.min(), 0))
    put("apOTRatioMax", f(o.apr.max(), 0))
    kn = bp[bp.coupling.isin(["kNN5", "kNN10", "kNN20", "kNN50"])]
    put("knnAllAUROCMin", f(kn.auroc.min()))
    put("knnAllAUROCMax", f(kn.auroc.max()))
    put("knnFiveEnr", rng_(bp[bp.coupling == "kNN5"].enrichment, 1))
    put("knnFiftyEnr", rng_(bp[bp.coupling == "kNN50"].enrichment, 1))
    put("mnnCov", rng_(bp[bp.coupling == "MNN10"].coverage_cells, 0, 100))
    put("otNeffMin", f"{o.n_eff.min():,.0f}".replace(",", "{,}"))
    put("otNeffMax", f"{o.n_eff.max():,.0f}".replace(",", "{,}"))
    put("otWnnTenPct", rng_(o.w_nn10, 0, 100))
    rows = []
    for t in TR:
        for i, (cp, lab) in enumerate(BASE):
            r = one(bp, transition=t, coupling=cp)
            rows.append(f"{TLAB[t] if i == 0 else ''} & {lab} & {f(r.auroc)} [{f(r.auroc_lo)}, {f(r.auroc_hi)}] & "
                        f"{f(r.apr, 1)} & {f(r.enrichment, 1)} [{f(r.enrichment_lo, 1)}, {f(r.enrichment_hi, 1)}] & "
                        f"{f(r.n_eff, 0)} & {pct(r.w_nn10, 0)} & {pct(r.coverage_cells, 0)}" + r" \\")
        rows.append(r"\addlinespace")
    write_table("tab_s14_baselines", rf"""\begin{{table*}}[!t]
\centering
\caption{{Post-review addition: proximity baselines (primary set, forward, pooled over
replicates; clone-bootstrap 95\% CIs). OT: $\varepsilon=0.05$. Gaussian kernel: row-normalised
$\exp(-C_{{ij}}/\varepsilon)$, $\varepsilon=0.05$; its AUROC and AP are identical for
$\varepsilon\in\{{0.01,0.05,0.1\}}$ because it ranks by distance alone. MNN: uniform over mutual
nearest neighbours; cells without one are excluded (coverage). AP ratio: average precision
divided by its chance level. $n_{{\mathrm{{eff}}}}$: effective number of partners
$\exp H(\mathbf{{w}}_i)$. NN10: weight on the ten nearest later cells.}}
\label{{tab:s14}}
\footnotesize
\begin{{tabular}}{{llcccccc}}
\toprule
Transition & Coupling & AUROC & AP ratio & Enrichment & $n_{{\mathrm{{eff}}}}$ & NN10 (\%) & Coverage (\%) \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")

    # ---- A4.1(a) lineage-grounded partner error
    p = pd.read_csv(R / "20_partner_error/summary.csv")
    for cp, nm in ((OT, "OT"), ("OT_unb_eps0.05", "Unb"), ("kNN10", "Knn")):
        for t in TR:
            r = one(p, transition=t, rep="pooled", coupling=cp)
            put(f"pe{nm}{TN[t]}", f(r.R))
            put(f"pe{nm}CI{TN[t]}", f"[{f(r.R_lo)}, {f(r.R_hi)}]")
    rows = []
    for t in TR:
        for i, rep in enumerate(["pooled", "1", "2", "3"]):
            cells = []
            for cp in (OT, "OT_unb_eps0.05", "kNN10"):
                r = one(p, transition=t, rep=rep, coupling=cp)
                cells.append(f"{f(r.R)} [{f(r.R_lo)}, {f(r.R_hi)}]")
            r0 = one(p, transition=t, rep=rep, coupling=OT)
            rows.append(f"{TLAB[t] if i == 0 else ''} & {'Pooled' if rep == 'pooled' else 'Rep.~' + rep} & "
                        f"{int(r0.n_clones_floor)} & " + " & ".join(cells) + r" \\")
        rows.append(r"\addlinespace")
    write_table("tab_s15_partner_error", rf"""\begin{{table*}}[!t]
\centering
\caption{{Post-review addition: lineage-grounded partner error (primary set, RQ1 30-PC space).
$D(i)=\sum_j w_{{ij}}\,\overline{{\|\mathbf{{z}}_j-\mathbf{{z}}_s\|^2}}$, averaged over the
sisters $s$ of $i$; $R=(D-D_{{\mathrm{{sis}}}})/(D_{{\mathrm{{indep}}}}-D_{{\mathrm{{sis}}}})$, where
$D_{{\mathrm{{sis}}}}$ is the distance between two sisters of the same clone (0: partners as close
as a sister; 1: as far as a random later cell). Clone-weighted means with clone-bootstrap 95\%
CIs over clones with at least two later sisters (clones).}}
\label{{tab:s15}}
\footnotesize
\begin{{tabular}}{{llcccc}}
\toprule
Transition & Replicate & Clones & OT balanced & OT unbalanced & kNN ($k=10$) \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")

    # ---- A4.1(b) RQ4 with structured partners
    s = pd.read_csv(R / "21_rq4_structured/summary.csv")
    assert s.all_ot_converged.all() and (s.n_seeds == 20).all()
    put("strSeeds", int(s.n_seeds.max()))
    for eff, tag in (("linear", ""), ("sparse", "Sp")):
        e = s[s.effect == eff]
        for pt, nm in (("OT_sampled", "OTSamp"), ("OT_barycentric", "OTBary"), ("kNN_sampled", "KnnSamp"),
                       ("kNN_averaged", "KnnAvg"), ("random", "Rand")):
            x = e[e.partner == pt]
            put(f"str{tag}{nm}RMSE", rng_(x.rmse, 2))
            put(f"str{tag}{nm}R", rng_(x.r, 2))
        oo = e[e.partner.str.startswith("OT")]
        put(f"str{tag}OTRMSEMax", f(oo.rmse.max()))
        put(f"str{tag}OTR", rng_(oo.r, 2))
    ot_s = s[s.partner.str.startswith("OT")]
    put("strOTMaxAbsBias", f(ot_s.bias.abs().max(), 3))
    kl = s[(s.partner.str.startswith("kNN")) & (s.effect == "linear")]
    assert (kl.bias_hi < 0).all()
    put("strKnnBias", f"{f(kl.bias.min(), 3)} to {f(kl.bias.max(), 3)}")
    put("strNeff", rng_(s.ot_neff_mean, 0))
    put("strWnnTenPct", rng_(s.ot_wnn10_mean, 0, 100))
    rows = []
    for eff, elab in (("linear", "Linear"), ("sparse", "Sparse (20 genes)")):
        for i, (pt, plab) in enumerate(PARTNERS):
            cells = []
            for rep in (1, 2, 3):
                r = one(s, effect=eff, partner=pt, rep=rep)
                cells.append(f"{f(r.rmse, 3)} [{f(r.rmse_lo, 3)}, {f(r.rmse_hi, 3)}] & {f(r.r)}")
            rows.append(f"{elab if i == 0 else ''} & {plab} & " + " & ".join(cells) + r" \\")
        rows.append(r"\addlinespace")
    r1 = s[s.effect == "linear"].groupby("rep").first()
    write_table("tab_s16_rq4_structured", rf"""\begin{{table*}}[!t]
\centering
\caption{{Post-review addition: RQ4 with structured partners (II$\to$III stage-III cells; per
seed, a random control half and a treated half, each capped at 4{{,}}000 cells; {M['strSeeds']}
seeds; 95\% intervals over seeds). Treated cells carry the effect of Section~2.6; partners are
drawn from the control half: a random cell, a cell sampled from the kNN ($k=10$) or OT
($\varepsilon=0.05$, fitted between the halves in their own 30-PC space) row, or the row
average (kNN averaged, OT barycentric). The oracle partner (the cell's own untreated profile)
has error 0. OT rows: effective number of partners {M['strNeff']}; weight on the ten nearest
controls {M['strWnnTenPct']}\%. Bias: OT within $\pm${M['strOTMaxAbsBias']}; kNN, linear effect,
{M['strKnnBias']}.}}
\label{{tab:s16}}
\footnotesize
\begin{{tabular}}{{llcccccc}}
\toprule
 & & \multicolumn{{2}}{{c}}{{Replicate 1}} & \multicolumn{{2}}{{c}}{{Replicate 2}} & \multicolumn{{2}}{{c}}{{Replicate 3}} \\
\cmidrule(lr){{3-4}}\cmidrule(lr){{5-6}}\cmidrule(lr){{7-8}}
Effect & Partner & RMSE & $r$ & RMSE & $r$ & RMSE & $r$ \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")

    # ---- A4.2 RQ2 target designs
    d = pd.read_csv(R / "22_rq2_targets/summary.csv")
    sh = pd.read_csv(R / "22_rq2_targets/ot_foreground_share.csv")
    for t in TR:
        put(f"fgShare{TN[t]}", rng_(sh[sh.transition == t]["mean"], 0, 100))
    # one row per design x model: Ridge on all draws (50), RF on draws 0-9 (A5)
    d = d[((d.model == "Ridge") & (d.draws == "all_draws")) | ((d.model == "RandomForest") & (d.draws != "all_draws"))
          | ((d.model == "RandomForest") & d.design.str.contains("expected"))].copy()
    d = d.drop_duplicates(["transition", "input", "coupling", "design", "model"])
    assert len(d) == 3 * 2 * 2 * 4 * 2
    q = lambda **kw: one(d, **kw)
    put("tdThreeSampOTPct", f(q(transition="III->IV", input="Pg_D", coupling="OT", design="sampled", model="Ridge").delta_pct_of_const, 1))
    put("tdThreeSampIndPct", f(q(transition="III->IV", input="Pg_D", coupling="independent", design="sampled", model="Ridge").delta_pct_of_const, 1))
    put("tdThreeExpOTPct", f(q(transition="III->IV", input="Pg_D", coupling="OT", design="expected", model="Ridge").delta_pct_of_const, 1))
    put("tdThreeFgExpOTPct", f(q(transition="III->IV", input="Pg_D", coupling="OT", design="fg_expected", model="Ridge").delta_pct_of_const, 1))
    put("tdThreeFgExpIndPct", f(q(transition="III->IV", input="Pg_D", coupling="independent", design="fg_expected", model="Ridge").delta_pct_of_const, 1))
    fg3 = d[(d.transition == "III->IV") & d.design.str.startswith("fg")]
    assert fg3.equivalent_to_const.all()  # both inputs, both couplings, both models
    put("tdThreeFgPct", f"{sg(fg3.delta_pct_of_const.min(), 1)} to {sg(fg3.delta_pct_of_const.max(), 1)}")
    put("tdThreeVsTrueSamp", f(q(transition="III->IV", input="Pg_D", coupling="OT", design="sampled", model="Ridge").delta_vs_true, 3))
    put("tdThreeVsTrueFg", f(q(transition="III->IV", input="Pg_D", coupling="OT", design="fg_expected", model="Ridge").delta_vs_true, 3))
    early = d[d.transition != "III->IV"]
    pg_e = early[early.input == "Pg_D"]
    assert pg_e.equivalent_to_const.all()
    po = pg_e[pg_e.coupling == "OT"]
    put("tdEarlyOTPct", f"{f(po.delta_pct_of_const.min(), 1)} to {f(po.delta_pct_of_const.max(), 1)}")
    put("tdEarlyVsTrue", f"{sg(po.delta_vs_true.min(), 3)} to {sg(po.delta_vs_true.max(), 3)}")
    pc2 = d[(d.transition == "II->III") & (d.input == "PC30") & (d.coupling == "OT")]
    assert ((pc2.ci_lo < -DELTA) & (pc2.ci_hi > -DELTA)).all()  # neither equivalent nor beyond the margin
    put("tdTwoPCOTPct", f"{f(pc2.delta_pct_of_const.min(), 1)} to {f(pc2.delta_pct_of_const.max(), 1)}")
    put("tdTwoPCOTVsTrue", f"{f(pc2.delta_vs_true.min(), 3)}--{f(pc2.delta_vs_true.max(), 3)}")
    ind = d[d.coupling == "independent"]
    ne = ind[~ind.equivalent_to_const]
    assert set(zip(ne.transition, ne.design)) == {("III->IV", "sampled")}, ne
    DES = [("sampled", "Sampled"), ("expected", "Expected"), ("fg_sampled", "Training clones, sampled"),
           ("fg_expected", "Training clones, expected")]
    rows = []
    for t in TR:
        for inp, ilab in (("Pg_D", "(Pg, D)"), ("PC30", "30-PC")):
            for cpl in ("OT", "independent"):
                for i, (des, dlab) in enumerate(DES):
                    a = q(transition=t, input=inp, coupling=cpl, design=des, model="Ridge")
                    rf = q(transition=t, input=inp, coupling=cpl, design=des, model="RandomForest")
                    rows.append(f"{TLAB[t] if (inp == 'Pg_D' and cpl == 'OT' and i == 0) else ''} & "
                                f"{ilab if (cpl == 'OT' and i == 0) else ''} & {CLAB[cpl] if i == 0 else ''} & {dlab} & "
                                f"{f(a.delta_pct_of_const, 1)} [{f(100 * a.ci_lo / a.const_median_mae, 1)}, {f(100 * a.ci_hi / a.const_median_mae, 1)}] & "
                                f"{f(rf.delta_pct_of_const, 1)} [{f(100 * rf.ci_lo / rf.const_median_mae, 1)}, {f(100 * rf.ci_hi / rf.const_median_mae, 1)}] & "
                                f"{f(a.delta_vs_true, 3)} & {f(rf.delta_vs_true, 3)}" + r" \\")
            rows.append(r"\addlinespace")
    write_table("tab_s17_rq2_targets", rf"""\begin{{table*}}[p]
\centering
\caption{{Post-review addition: RQ2 training-target designs (primary set, 20 splits). Sampled:
one target drawn from the coupling row restricted to $\mathcal{{A}}$ (training-clone and
background later cells), as in the main analysis; expected: the row-weighted mean target;
training clones: $\mathcal{{A}}$ restricted to later cells of training clones. Sampled designs:
Ridge over 50 target draws, random forest over draws 0--9. Model minus constant median trained
on the same targets, \% of the constant-median MAE, with 95\% split-stability interval;
equivalence margin $\delta=0.01$ MAE. Last two columns: MAE minus that of the same model
trained on true sister pairs. Mean share of the OT target mass on training-clone cells, by
replicate: {M['fgShareOne']}\%, {M['fgShareTwo']}\% and {M['fgShareThree']}\% for the three
transitions.}}
\label{{tab:s17}}
\scriptsize
\setlength{{\tabcolsep}}{{4pt}}
\begin{{tabular}}{{llllcccc}}
\toprule
 & & & & \multicolumn{{2}}{{c}}{{vs.\ constant (\%)}} & \multicolumn{{2}}{{c}}{{vs.\ true pairs (MAE)}} \\
\cmidrule(lr){{5-6}}\cmidrule(lr){{7-8}}
Transition & Input & Pairs & Target design & Ridge & Random forest & Ridge & RF \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")

    # ---- A4.3 fate labels: leave-one-clone-out transfer accuracy and paired differences
    a = pd.read_csv(R / "18_fate_unc/label_transfer_accuracy.csv")
    for t in TR:
        r = one(a, transition=t, rep="pooled")
        put(f"ltAcc{TN[t]}", f(r.accuracy))
        put(f"ltBacc{TN[t]}", f(r.balanced_accuracy))
    ar = a[a.rep != "pooled"]
    put("ltBaccMin", f(ar.balanced_accuracy.min()))
    pt = pd.read_csv(R / "18_fate_unc/paired_tvd.csv")
    pi = pt[(pt.rep == "pooled") & (pt.comparison == f"{OT} - independent")]
    w = pi.pivot_table(index=["transition", "early_cluster"], columns="labels", values=["diff_lo", "diff_hi"])
    better = (w["diff_hi"]["k15"] < 0) & (w["diff_hi"]["measured"] < 0)
    worse = (w["diff_lo"]["k15"] > 0) & (w["diff_lo"]["measured"] > 0)
    put("fpNClusters", len(w))
    put("fpBothBetter", int(better.sum()))
    put("fpBothWorse", int(worse.sum()))
    assert sorted(w.index[worse].get_level_values(1)) == ["I-NbpHB", "II-Motors", "III-NPpHB"]
    for cl in ("III-GABAs", "III-Gluts"):  # the A3.5 'reversals' are not distinguishable from zero
        x = w.xs(cl, level=1).iloc[0]
        assert x["diff_lo"]["k15"] < 0 < x["diff_hi"]["k15"] and x["diff_lo"]["measured"] < 0 < x["diff_hi"]["measured"]
    for cl, nm in (("III-NPaHB", "NPaHB"), ("III-NPvMB", "NPvMB")):
        for lab, tg in (("k15", "Trans"), ("measured", "Meas")):
            r = one(pi, early_cluster=cl, labels=lab)
            put(f"fp{tg}{nm}", f"{sg(r['diff'])} [{sg(r.diff_lo)}, {sg(r.diff_hi)}]")
    pk = pt[(pt.rep == "pooled") & (pt.comparison == f"{OT} - kNN") & (pt.labels == "k15")]
    put("fpKnnBetter", f"{int((pk.diff_hi < 0).sum())}")
    put("fpKnnN", f"{len(pk)}")
    rows = []
    for t in TR:
        cls = sorted(pi[pi.transition == t].early_cluster.unique())
        for i, cl in enumerate(cls):
            cells = []
            for comp in ("independent", "kNN"):
                for lab in ("k15", "measured"):
                    r = one(pt, rep="pooled", early_cluster=cl, labels=lab, comparison=f"{OT} - {comp}")
                    cells.append(f"{f(r['diff'])} [{f(r.diff_lo)}, {f(r.diff_hi)}]")
            n = int(one(pi, early_cluster=cl, labels="k15").n_clones)
            rows.append(f"{TLAB[t] if i == 0 else ''} & {cl} ({n}) & " + " & ".join(cells) + r" \\")
        rows.append(r"\addlinespace")
    acc = "; ".join(f"{TLAB[t]} {M['ltAcc' + TN[t]]} ({M['ltBacc' + TN[t]]})" for t in TR)
    write_table("tab_s18_fate_paired", rf"""\begin{{table*}}[!t]
\centering
\caption{{Post-review addition: paired fate-agreement differences (primary set, pooled over
replicates). $\Delta$TVD = TVD(OT) $-$ TVD(reference) for the same clones, with 95\% paired
clone-bootstrap intervals (negative: OT closer to the clonal fates); all labels (authors' plus
$k=15$ transferred) and measured labels only. Clones in parentheses. Leave-one-clone-out
accuracy of the label transfer (balanced accuracy in parentheses), pooled: {acc}; lowest
per-replicate balanced accuracy {M['ltBaccMin']}.}}
\label{{tab:s18}}
\footnotesize
\setlength{{\tabcolsep}}{{4pt}}
\begin{{tabular}}{{llcccc}}
\toprule
 & & \multicolumn{{2}}{{c}}{{OT $-$ independent}} & \multicolumn{{2}}{{c}}{{OT $-$ kNN}} \\
\cmidrule(lr){{3-4}}\cmidrule(lr){{5-6}}
Transition & Early cluster & All labels & Measured & All labels & Measured \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")

    # ---- A4.5 within-dataset generality
    g = pd.read_csv(R / "23_rq1_generality/summary.csv")
    assert pd.read_csv(R / "23_rq1_generality/solves.csv").converged.all()
    gp = g[g.rep == "pooled"]
    sp = gp[gp.space != "rq1"]
    put("genOTAUROC", rng_(sp[sp.coupling == OT].auroc, 2))
    put("genKnnAUROC", rng_(sp[sp.coupling == "kNN"].auroc, 2))
    put("genOTRThree", rng_(sp[(sp.coupling == OT) & (sp.transition == "III->IV")].R, 2))
    gr = gp[gp.space == "rq1"]
    for lam, nm in (("0.5", "Half"), ("1.0", "One")):
        for t in TR:
            put(f"grAUROC{nm}{TN[t]}", f(one(gr, coupling=f"OT_growth_bal_lam{lam}", transition=t).auroc))
            put(f"grR{nm}{TN[t]}", f(one(gr, coupling=f"OT_growth_bal_lam{lam}", transition=t).R))
    put("grUnbAUROC", rng_(gr[gr.coupling.str.startswith("OT_growth_unb")].auroc, 2))
    SPL = [("hvg2000_pc15", "2{,}000 HVGs, 15 PCs"), ("hvg2000_pc50", "2{,}000 HVGs, 50 PCs"),
           ("hvg4000_pc30", "4{,}000 HVGs, 30 PCs")]
    rows = []
    for t in TR:
        k = 0
        for spc, slab in SPL:
            for cp, clab in ((OT, "OT balanced"), ("kNN", "kNN")):
                r = one(gp, transition=t, space=spc, coupling=cp)
                rows.append(f"{TLAB[t] if k == 0 else ''} & {slab if cp == OT else ''} & {clab} & "
                            f"{f(r.auroc)} [{f(r.auroc_lo)}, {f(r.auroc_hi)}] & {f(r.enrichment, 1)} & {f(r.R)} [{f(r.R_lo)}, {f(r.R_hi)}]" + r" \\")
                k += 1
        for cp, clab in (("OT_growth_bal_lam0.5", r"Balanced, $\lambda=0.5$"), ("OT_growth_bal_lam1.0", r"Balanced, $\lambda=1$"),
                         ("OT_growth_unb_lam0.5", r"Unbalanced, $\lambda=0.5$"), ("OT_growth_unb_lam1.0", r"Unbalanced, $\lambda=1$")):
            r = one(gp, transition=t, space="rq1", coupling=cp)
            rows.append(f" & {'Growth-informed OT' if 'bal_lam0.5' in cp and 'unb' not in cp else ''} & {clab} & "
                        f"{f(r.auroc)} [{f(r.auroc_lo)}, {f(r.auroc_hi)}] & {f(r.enrichment, 1)} & {f(r.R)} [{f(r.R_lo)}, {f(r.R_hi)}]" + r" \\")
        rows.append(r"\addlinespace")
    write_table("tab_s19_generality", rf"""\begin{{table*}}[!t]
\centering
\caption{{Post-review addition: coupling recovery in other expression spaces and with
growth-informed OT (primary set, forward, pooled; clone-bootstrap 95\% CIs). OT: $\varepsilon=0.05$;
kNN: $k=10$. Growth-informed OT uses the RQ1 space with row marginal $a_i\propto\exp(\lambda
\tilde g_i)$, $\tilde g_i$ the standardised mean expression of the S-phase and G2/M genes
of~\cite{{Tirosh_2016}}; unbalanced: $\tau=0.9$. Partner error $R$ as in Table~S15. RQ1 space
(2{{,}}000 HVGs, 30 PCs), for reference: OT AUROC {M['rqOneAUROCOne']}, {M['rqOneAUROCTwo']},
{M['rqOneAUROCThree']}; $R$ {M['peOTOne']}, {M['peOTTwo']}, {M['peOTThree']}.}}
\label{{tab:s19}}
\footnotesize
\begin{{tabular}}{{llcccc}}
\toprule
Transition & Space & Coupling & AUROC & Enrichment & $R$ \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")

    # ---- A4.6 LARRY
    l = pd.read_csv(R / "24_larry/summary.csv")
    ls = pd.read_csv(R / "24_larry/solves.csv")
    assert ls.converged.all()
    l["apr"] = l.ap / l.ap_chance
    l["pool"] = l["pool"].astype(str)
    lpo = l[l.pool != "all_pools"]
    lal = l[l.pool == "all_pools"]
    put("larryNClonesA", f"{int(ls[ls.transition == 'd2->d4'].n_clones.max()):,}".replace(",", "{,}"))
    put("larryNClonesB", rng_(ls[ls.transition == "d4->d6"].n_clones, 0))
    for t in LT:
        r = one(lal, transition=t, coupling=OT)
        put(f"larryAUROC{LN[t]}", f(r.auroc))
        put(f"larryAPR{LN[t]}", f(r.apr, 0))
        put(f"larryR{LN[t]}", f(r.R))
        put(f"larryEnr{LN[t]}", f(r.enrichment, 0))
        put(f"larryKnnR{LN[t]}", f(one(lal, transition=t, coupling="kNN10").R))
        put(f"larryGaussAUROC{LN[t]}", f(one(lal, transition=t, coupling="Gauss_eps0.05").auroc))
        put(f"larryAUROCPools{LN[t]}", rng_(lpo[(lpo.transition == t) & (lpo.coupling == OT)].auroc, 2))
    bal = lal[lal.coupling.str.startswith("OT_bal")]
    put("larryOTEpsAUROC", rng_(bal.auroc, 2))
    put("larryKnnAUROC", rng_(lal[lal.coupling.str.startswith("kNN")].auroc, 2))
    put("larryPoolMaxDev", f(max((lpo[(lpo.transition == t) & (lpo.coupling == OT)].auroc
                                  - one(lal, transition=t, coupling=OT).auroc).abs().max() for t in LT), 3))
    lf = pd.read_csv(R / "24_larry/fate_agreement.csv")
    lw = lf.pivot_table(index=["transition", "early_cluster", "n_clones"], columns="coupling", values="tvd").reset_index()
    for t in LT:
        x = lw[lw.transition == t]
        put(f"larryFateNBetter{LN[t]}", int((x[OT] < x["independent"]).sum()))
        put(f"larryFateN{LN[t]}", len(x))
    nr = one(lw, transition="d2->d4", early_cluster="Neutrophil")
    put("larryTvdOTNeu", f(nr[OT]))
    put("larryTvdIndNeu", f(nr["independent"]))
    ur = one(lw, transition="d4->d6", early_cluster="Undifferentiated")
    put("larryTvdOTUnd", f(ur[OT]))
    put("larryTvdIndUnd", f(ur["independent"]))
    rows = []
    for t in LT:
        for i, (cp, lab) in enumerate(BASE):
            r = one(lal, transition=t, coupling=cp)
            pr = lpo[(lpo.transition == t) & (lpo.coupling == cp)]
            rows.append(f"{LLAB[t] if i == 0 else ''} & {lab} & {f(r.auroc)} ({rng_(pr.auroc, 2)}) & {f(r.apr, 1)} & "
                        f"{f(r.enrichment, 1)} & {f(r.R)} ({rng_(pr.R, 2)}) & {f(r.n_eff, 0)} & {pct(r.coverage_cells, 0)}" + r" \\")
        rows.append(r"\addlinespace")
    rows2 = []
    for t in LT:
        x = lw[lw.transition == t].sort_values("early_cluster")
        for i, (_, r) in enumerate(x.iterrows()):
            rows2.append(f"{LLAB[t] if i == 0 else ''} & {r.early_cluster} ({int(r.n_clones)}) & {f(r[OT])} & "
                         f"{f(r['kNN10'])} & {f(r['independent'])}" + r" \\")
        rows2.append(r"\addlinespace")
    nfg = ls.groupby("transition")[["n_early_fg", "n_later_fg"]].agg(["min", "max"])
    write_table("tab_s20_larry", rf"""\begin{{table*}}[!t]
\centering
\caption{{Post-review addition: second lineage dataset, LARRY in vitro haematopoiesis
\cite{{Weinreb_2020}} (GEO GSE140802). Per transition, three pools, each a subsample of whole
cross-stage clones (at most 6{{,}}000 foreground cells per time point) plus background cells
(8{{,}}000 cells per time point); pools are clone subsamples, not biological replicates. At day
2$\to$4 all {M['larryNClonesA']} cross-stage clones fit in every pool, so the pools differ only
in background cells; at day 4$\to$6 they hold {M['larryNClonesB']} clones. Metrics as in
Tables~S14 and S15, computed over all pools; per-pool ranges in parentheses.}}
\label{{tab:s20}}
\footnotesize
\begin{{tabular}}{{llcccccc}}
\toprule
Transition & Coupling & AUROC (pools) & AP ratio & Enrichment & $R$ (pools) & $n_{{\mathrm{{eff}}}}$ & Coverage (\%) \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}

\vspace{{6pt}}
\begin{{tabular}}{{llccc}}
\toprule
Transition & Early cell type (clones) & OT & kNN & Independent \\
\midrule
{chr(10).join(rows2[:-1])}
\bottomrule
\end{{tabular}}
\end{{table*}}
""")

    # ---- A4.7 CellRank 2
    c = pd.read_csv(R / "25_cellrank/rq3_comparison.csv")
    for t in TR:
        put(f"crAUROC{TN[t]}", f(one(c, transition=t, comparison="CellRank2 - Pg_D").auroc_cellrank))
        r = one(c, transition=t, comparison="CellRank2 - Pg_D")
        put(f"crVsPgD{TN[t]}", f"{f(r.delta)} [{f(r.ci_lo)}, {f(r.ci_hi)}]")
    lo = c[c.comparison.isin(["CellRank2 - PC30", "CellRank2 - fate_axis"])]
    assert (lo.ci_hi < 0).all()
    put("crVsPCFate", f"{f(lo.delta.max())} to {f(lo.delta.min())}")
    cf = pd.read_csv(R / "25_cellrank/fate_agreement.csv")
    cw = cf.pivot_table(index=["transition", "early_cluster", "n_clones"], columns="method", values="tvd").reset_index()
    put("crFateNBetterInd", int((cw.CellRank2 < cw.independent).sum()))
    put("crFateNBetterOT", int((cw.CellRank2 < cw[OT]).sum()))
    put("crFateN", len(cw))
    v = one(cw, early_cluster="III-NPvMB")
    put("crTvdNPvMB", f(v.CellRank2))
    rows = []
    for t in TR:
        x = cw[cw.transition == t].sort_values("early_cluster")
        for i, (_, r) in enumerate(x.iterrows()):
            rows.append(f"{TLAB[t] if i == 0 else ''} & {r.early_cluster} ({int(r.n_clones)}) & {f(r.CellRank2)} & "
                        f"{f(r[OT])} & {f(r.independent)}" + r" \\")
        rows.append(r"\addlinespace")
    rows2 = []
    for t in TR:
        for i, (cmp_, lab) in enumerate((("P_tcbb_D", "(P, D)"), ("Pg_D", "(Pg, D)"), ("PC30", "30-PC"), ("fate_axis", "Fate-bias axis"))):
            r = one(c, transition=t, comparison=f"CellRank2 - {cmp_}")
            rows2.append(f"{TLAB[t] if i == 0 else ''} & {f(r.auroc_cellrank) if i == 0 else ''} & {lab} & {f(r.auroc_other)} & "
                         f"{f(r.delta)} [{f(r.ci_lo)}, {f(r.ci_hi)}]" + r" \\")
        rows2.append(r"\addlinespace")
    write_table("tab_s21_cellrank", rf"""\begin{{table}}[!t]
\centering
\caption{{Post-review addition: CellRank~2~\cite{{Weiler_2024}} fate probabilities (real-time
kernel from the balanced OT plan, $\varepsilon=0.05$, with within-stage connectivity kernels,
weight 0.2; terminal states = later fate classes). Top: fate agreement (TVD, lower is better;
clones in parentheses). Bottom: as an RQ3 fate predictor on the test clones of the 20 splits
(no fitting); difference in macro-AUROC from each RQ3 representation, with 95\%
split-stability interval.}}
\label{{tab:s21}}
\footnotesize
\setlength{{\tabcolsep}}{{3pt}}
\begin{{tabular}}{{llccc}}
\toprule
Transition & Early cluster & CellRank 2 & OT & Indep. \\
\midrule
{chr(10).join(rows[:-1])}
\bottomrule
\end{{tabular}}

\vspace{{6pt}}
\begin{{tabular}}{{lclcc}}
\toprule
Transition & CellRank 2 & Representation & AUROC & Difference \\
\midrule
{chr(10).join(rows2[:-1])}
\bottomrule
\end{{tabular}}
\end{{table}}
""")


# ============================================================ output
def write_table(name, text):
    (TAB / f"{name}.tex").write_text(
        f"% Generated by scripts/10_results_paper.py from results/06-25. Do not edit.\n{text}",
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
    assert M["loroNPCAbove"] == M["loroNFateAbove"] == M["loroNFolds"]  # "each ... in all"


def main():
    fig1()
    rq1()
    rq2()
    rq3()
    rq4()
    post_review()
    second_review()
    check_wording()
    lines = ["% Generated by scripts/10_results_paper.py from results/06_rq1 ... 25_cellrank.",
             "% Do not edit by hand."]
    lines += [rf"\newcommand{{\{k}}}{{{v}}}" for k, v in M.items()]
    (GEN / "results_numbers.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
    pd.Series(M, name="value").rename_axis("macro").to_csv(OUT / "results_numbers.csv")

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT,
                           capture_output=True, text=True).stdout.strip()
    prov = [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
            "plan: ANALYSIS_PLAN.md approved at 65cd70d (amendments A1-A6; A3, A4 = post-review additions; A5, A6 = compute reductions for A4)",
            f"python {sys.version.split()[0]} pandas {pd.__version__} numpy {np.__version__} "
            f"matplotlib {matplotlib.__version__}",
            "environment: Windows (VS Code session)"]
    prov += [f"input {p.relative_to(ROOT).as_posix()} sha256 {sha256(p)}" for p in INPUTS]
    (OUT / "provenance.txt").write_text("\n".join(prov) + "\n")
    print(f"{len(M)} macros; tables and figures written")
    print(json.dumps(M, indent=0)[:6000])


if __name__ == "__main__":
    main()
