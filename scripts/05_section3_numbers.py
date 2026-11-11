"""Generate every number used in paper Section 3 (Data and setup) from committed results.

Outputs
  paper/generated/data_numbers.tex    one \\newcommand per number used in the prose
  paper/tables/tab_data_counts.tex    cells and cross-stage clones per transition
  results/05_section3/provenance.txt  inputs (SHA-256) and git commit

Inputs (all committed): results/01_clone_counts, results/02_library_map,
results/03_objects, results/04_qc. No number in paper/main.tex is typed by hand; the
only constants here are the published counts' source label and the stage days
(from GEO sample metadata), which are design facts, not results.

Consistency checks: cross-stage counts recomputed from per-cell tables must equal
results/04_qc/clone_counts.csv; otherwise the script stops.
"""
import hashlib
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"
IN = {
    "clone_counts_01": R / "01_clone_counts/clone_counts.csv",
    "suffix_map": R / "02_library_map/suffix_library_map.csv",
    "exp3_gap": R / "02_library_map/exp3_gap_check.csv",
    "obs": R / "03_objects/cell_obs.csv.gz",
    "obj_summary": R / "03_objects/summary.csv",
    "qc_after_paper": R / "03_objects/clone_counts_after_qc.csv",
    "cell_qc": R / "04_qc/cell_qc.csv.gz",
    "qc_summary": R / "04_qc/summary.csv",
    "qc_lib": R / "04_qc/qc_by_library.csv",
    "qc_clones": R / "04_qc/clone_counts.csv",
}
SEC = ROOT / "paper/generated"
TAB = ROOT / "paper/tables"
OUT = R / "05_section3"
for d in (SEC, TAB, OUT):
    d.mkdir(parents=True, exist_ok=True)

PAIRS = {1: ("I", "II"), 2: ("II", "III"), 3: ("III", "IV")}
TR = {1: "I->II", 2: "II->III", 3: "III->IV"}
DAYS = {"I": 13, "II": 21, "III": 30, "IV": 45}  # GEO sample characteristics


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def n(x):
    return f"{int(x):,}"


def pct(a, b, d=1):
    return f"{100 * a / b:.{d}f}"


def summary(p):
    s = pd.read_csv(p, index_col=0, dtype=str)["value"]
    num = pd.to_numeric(s, errors="coerce")
    return pd.Series([v if pd.isna(x) else x for v, x in zip(s, num)], index=s.index,
                     dtype=object)


def cross_stage(cells):
    """cells: DataFrame with exp, stage, clone. Returns {exp: (clones, cells, shared_set)}."""
    out = {}
    for exp, (e, l) in PAIRS.items():
        c = cells[cells.exp == exp]
        st = c.groupby("clone")["stage"].agg(set)
        shared = set(st[st.apply(lambda s: {e, l} <= s)].index)
        out[exp] = (len(shared), int(c.clone.isin(shared).sum()), shared)
    return out


def main():
    obs = pd.read_csv(IN["obs"], low_memory=False)
    qc = pd.read_csv(IN["cell_qc"])
    obs = obs.merge(qc[["cell_id", "keep", "doublet_flag"]], on="cell_id", how="left",
                    validate="one_to_one")
    ct = obs[obs.in_clone_table.astype(bool)].copy()
    objs, qcs = summary(IN["obj_summary"]), summary(IN["qc_summary"])
    qcc = pd.read_csv(IN["qc_clones"]).set_index(["transition", "cell_set"])
    pqc = pd.read_csv(IN["qc_after_paper"]).set_index(["transition", "cells_used"])
    lib = pd.read_csv(IN["qc_lib"])
    smap = pd.read_csv(IN["suffix_map"])
    gap = summary(IN["exp3_gap"])

    # ---- consistency checks against committed clone counts
    sets = {"all_clone_table_cells": ct,
            "primary": ct[ct.keep],
            "sensitivity_doublets_excluded": ct[ct.keep & ~ct.doublet_flag]}
    cs = {k: cross_stage(v) for k, v in sets.items()}
    for k, res in cs.items():
        for exp, (cl, ce, _) in res.items():
            ref = qcc.loc[(TR[exp], k)]
            assert (cl, ce) == (ref.clones, ref.cells), (k, exp, cl, ce, ref.to_dict())

    pub = {exp: (qcc.loc[(TR[exp], "primary")].published_clones,
                 qcc.loc[(TR[exp], "primary")].published_cells) for exp in PAIRS}

    # ---- derived numbers
    multi = ct.clone.str.contains(";")
    by_sheet = ct.assign(m=multi).groupby(["stage", "exp"])["m"].mean()
    umi_cap = 30000
    thr = objs["thresholds"]
    assert f"UMI<={umi_cap}" in thr, thr
    ct_iv = ct[ct.stage == "IV"]
    swap = smap[~smap.rep_consistent.astype(bool)]
    # viral replicate -> set of GEO sample names it actually maps to
    swap_txt = {r: sorted(set(g.gsm)) for r, g in swap.groupby("rep_from_virus")}
    da_iv = ct[(ct.exp == 3) & (ct.stage == "IV") & ct.keep
               & ct.clone.isin(cs["primary"][3][2]) & (ct.cluster == "IV-mDAs")]
    ratio = lib.median_umi_flagged / lib.median_umi_keep
    e3_cl_all, e3_ce_all, _ = cs["all_clone_table_cells"][3]

    M = {
        "nLibraries": n(objs["libraries"]),
        "nGenes": n(objs["genes"]),
        "nBarcodes": n(objs["barcodes"]),
        "nPublishedCells": n(objs["published_hq_cells"]),
        "nCloneTableCells": n(len(ct)),
        "pctMultiBarcode": pct(multi.sum(), len(ct), 0),
        "pctMultiBarcodeMin": f"{100 * by_sheet.min():.0f}",
        "pctMultiBarcodeMax": f"{100 * by_sheet.max():.0f}",
        "nStageIVCells": n(gap["stageIV_cells"]),
        "nStageIVCollisions": n(gap["stageIV_bc_in_gt1_IV_library"]),
        "nExpThreeClones": n(e3_cl_all),
        "nExpThreeCells": n(e3_ce_all),
        "nExpThreePubClones": n(pub[3][0]),
        "nExpThreePubCells": n(pub[3][1]),
        "pctExpThreeCloneGap": pct(pub[3][0] - e3_cl_all, pub[3][0]),
        "pctExpThreeCellGap": pct(pub[3][1] - e3_ce_all, pub[3][1]),
        "swapRepTwoGSM": ", ".join(swap_txt[2]),
        "swapRepThreeGSM": ", ".join(swap_txt[3]),
        # paper QC
        "nPaperQCPass": n(objs["qc_pass"]),
        "nCloneTableFailPaperQC": n(objs["clone_table_cells_failing_paper_qc"]),
        "pctCloneTableFailPaperQC": pct(objs["clone_table_cells_failing_paper_qc"], len(ct)),
        "nCloneTableOverUMICap": n((ct.total_counts > umi_cap).sum()),
        "pctStageIVOverUMICap": pct((ct_iv.total_counts > umi_cap).sum(), len(ct_iv), 0),
        "nExpThreeClonesPaperQC": n(pqc.loc[("III->IV", "qc_pass_only")].clones),
        "umiCap": n(umi_cap),
        # primary rule
        "nPrimaryCells": n(qcs["primary_cells"]),
        "nCloneTablePrimary": n(qcs["clone_table_primary"]),
        "nCloneTableDropped": n(len(ct) - qcs["clone_table_primary"]),
        "pctCloneTableKept": pct(qcs["clone_table_primary"], len(ct)),
        # doublets
        "nScrubletAutoTotal": n(qcs["scrublet_auto_predicted_total"]),
        "nScrubletAutoMax": n(lib.scrublet_n_predicted.max()),
        "nLibsScrubletZero": n((lib.scrublet_n_predicted == 0).sum()),
        "rhoDoubletUMI": f"{qcs['doublet_score_vs_umi_spearman']:.2f}",
        "umiRatioFlagMin": f"{ratio.min():.1f}",
        "umiRatioFlagMax": f"{ratio.max():.1f}",
        "nSensitivityCells": n(qcs["sensitivity_cells"]),
        "nCloneTableSensitivity": n(qcs["clone_table_sensitivity"]),
        "nDACellsCrossStage": n(len(da_iv)),
    }
    # rule constants, parsed from the committed rule strings (not typed)
    m = re.search(r"genes>=(\d+); mt%<=([\d.]+)", qcs["primary_rule"])
    M["minGenes"], M["maxMt"] = n(m.group(1)), f"{float(m.group(2)):g}"
    m = re.search(r"top ([\d.]+)%-per-([\d,]+)-barcodes", qcs["sensitivity_rule"])
    M["dblRatePct"], M["dblRatePer"] = m.group(1), m.group(2)
    assert re.search(rf"genes>={m and M['minGenes'].replace(',', '')}", thr)
    for s, d in DAYS.items():
        M[f"day{s.replace('IV', 'Four').replace('III', 'Three').replace('II', 'Two').replace('I', 'One')}"] = str(d)

    lines = ["% AUTO-GENERATED by scripts/05_section3_numbers.py -- do not edit by hand."]
    lines += [f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in M.items()]
    (SEC / "data_numbers.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # ---- table: rows = cell set, columns = transition (clones / cells)
    rows = [("Published~\\cite{You_2023}", {e: pub[e] for e in PAIRS}),
            ("All clone-table cells", {e: cs["all_clone_table_cells"][e][:2] for e in PAIRS}),
            ("Paper thresholds (not used)",
             {e: tuple(pqc.loc[(TR[e], "qc_pass_only")][["clones", "cells"]]) for e in PAIRS}),
            ("Primary (this study)", {e: cs["primary"][e][:2] for e in PAIRS}),
            ("Doublet sensitivity", {e: cs["sensitivity_doublets_excluded"][e][:2] for e in PAIRS})]
    hdr = (" & ".join(f"{a}$\\to${b}" for a, b in PAIRS.values()) + " \\\\\n"
           + " & " + " & ".join(f"d{DAYS[a]}$\\to$d{DAYS[b]}" for a, b in PAIRS.values()))
    t = ["% AUTO-GENERATED by scripts/05_section3_numbers.py -- do not edit by hand.",
         "\\begin{table}[!t]",
         "\\centering",
         "\\caption{Cross-stage clones / cells in those clones, per transition and cell set.}",
         "\\label{tab:data_counts}",
         "\\footnotesize",
         "\\setlength{\\tabcolsep}{3pt}",
         "\\begin{tabular}{lrrr}",
         "\\hline",
         f"Cell set & {hdr} \\\\",
         "\\hline"]
    for name, vals in rows:
        t.append(name + " & " + " & ".join(f"{n(vals[e][0])} / {n(vals[e][1])}" for e in PAIRS)
                 + " \\\\")
    t += ["\\hline",
          "\\end{tabular}",
          "\\end{table}"]
    (TAB / "tab_data_counts.tex").write_text("\n".join(t) + "\n", encoding="utf-8")

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                            capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT,
                           capture_output=True, text=True).stdout.strip()
    prov = [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
            f"python {sys.version.split()[0]} pandas {pd.__version__}"]
    prov += [f"input {p.relative_to(ROOT).as_posix()} sha256 {sha256(p)}" for p in IN.values()]
    (OUT / "provenance.txt").write_text("\n".join(prov) + "\n")
    pd.Series(M).to_csv(OUT / "section3_numbers.csv", header=["value"])
    print(pd.Series(M).to_string())
    print("\n".join(t))


if __name__ == "__main__":
    main()
