# CLAUDE.md — Lineage-grounded evaluation (SISBAR) project

Context for Claude Code sessions on this repository. Author: Victor Ikenna Kanu,
Network Systems Lab, Kumoh National Institute of Technology (PI: Prof. Dong-Seong Kim).

## Session workflow — read first

1. `git pull`, then read **UPDATES.md**. Work its *Open instructions* in order.
2. At the end of each work session, add a dated entry to UPDATES.md's *Log*, move
   finished instructions to *Done* (with commit hash), and list blocking decisions
   under *Open questions*. Commit and push.
3. Copy durable facts and decisions into this file too.

## What this project is

An evaluation paper (target: IEEE TCBB, fallback IEEE JBHI) that uses measured
clonal lineage to test two assumptions that differentiation models rely on:

1. that inferred couplings (optimal transport, OT) link cells to their true descendants;
2. that marker-based states such as (P, D) capture cell fate.

It is NOT a new-model paper. Do not add new architectures unless asked.

It builds on the author's earlier work:
- **TCBB benchmark paper** ("Benchmarking Single-Cell Marker-State Forecasting Under
  Independent and Coupled Cell Pairings"): leakage-safe protocol, Proposition 1
  (independent pairing ⇒ constant optimal predictor), OT positive controls.
  Code: https://github.com/xaviwho/leakage-safe-benchmark — reuse it.
- **DT1** (proteomics-integrated 3D state space; under review at CMPB).
- **DT2** (control-paired gate framework; under review at CMPB).

## Data

You et al., "Mapping of clonal lineages across developmental stages in human neural
differentiation", Cell Stem Cell 2023, DOI 10.1016/j.stem.2023.02.007 (PMID 36933556).
Method: SISBAR (single-cell split barcoding). One hESC line (H9), ventral
midbrain–hindbrain differentiation.

| Transition | Days     | Cross-stage clones | Cells |
|------------|----------|--------------------|-------|
| I → II     | 13 → 21  | 768                | 2,762 |
| II → III   | 21 → 30  | 1,213              | 5,856 |
| III → IV   | 30 → 45  | 559                | 2,816 |

Counts are pooled over all three biological replicates. GEO: **GSE207921 = replicate 1,
all four stages; GSE221592 = replicates 2–3, all four stages** (plus the clone table).
Viral barcodes end in `-R-E` (replicate, experiment); the cell-barcode `_N` suffix is
the 10x library, except in the stage-IV sheet where it is missing — use the replicate
from the viral barcode there.

- 173,564 cells total; 26% barcoded, in 10,702 clones; 2,540 clones span two stages.
- **No clone spans three stages.** All forecasting is one-step (t → t+1).
- Ground truth is **sister cells**: clones split after ~3–4 divisions, one half
  profiled early, the other later. Not the same cell.
- Stage IV contains midbrain DA neurons (NR4A2+/PITX3+/EN1+). DA-specific clones are
  few (72 from the stage-III ventral midbrain progenitor) → DA analyses are secondary.
- Supporting data: Jerber et al. 2021 (D11/D30/D52) for transfer of the fate axis.

### Data check results (2026-11-11)
- Clone table: `GSE221592_Virusbar_processed_data_combined.xlsx` (GSE221592 suppl;
  SHA-256 ec64e23ff7ccf2c0e3d245bdb14e16444a948b70d52918892554e56a3961eb2b).
  Six sheets (Clone_stage_I_Exp#1 … Clone_stage_IV_Exp#3); columns Cell.Barcode,
  Virus.Barcode (`;`-separated when a cell carries several), Cluster.Id. Covers all
  three replicates; the `_N` suffix on Cell.Barcode is the library index.
- **Clone definition that reproduces the paper:** a clone = the exact set of viral
  barcodes on a cell (the full Virus.Barcode string), NOT each barcode separately.
  ~36–53% of cells carry more than one viral barcode.
- Reproduced exactly: Exp1 I→II 768 clones / 2,762 cells; Exp2 II→III 1,213 / 5,856.
- Exp3 III→IV gives 555 / 2,692 vs published 559 / 2,816. Stage-IV cell barcodes have
  no library suffix; ~18 barcodes appear in more than one of the three stage-IV
  libraries (GSM6322995, GSM6889988, GSM6889994). Resolve before RQ1.
- **Only 8 stage-IV mDA cells sit in cross-stage clones** (and 21 mGluts). DA-endpoint
  analyses are not feasible; treat DA as progenitor fate bias (III-NPvMB clones), not
  as a forecast target.
- Expression: per-library 10x matrices (barcodes/features/matrix) on GEO; replicate 1
  in GSE207921, replicates 2–3 in GSE221592. 22 libraries, 215,960 barcodes in total
  (unfiltered vs the paper's 173,564 post-QC cells).
- **Library map done** (`scripts/02_library_map.py` → `results/02_library_map/`): every
  clone-table cell is found in exactly the library its suffix/replicate predicts (100%
  overlap; runner-up ≤0.6%). Per-cell map: `results/02_library_map/cell_library_map.csv.gz`.
- **Exp3 replicate labels 2↔3 are swapped** between the viral barcode and the GEO sample
  name, consistently in stage III and IV (viral rep 2 → GSM6889993 "Rep3", GSM6889994
  "Rep3"; viral rep 3 → GSM6889987 "Rep2", GSM6889988 "Rep2"). Clone pairing is
  unaffected. Use the viral-barcode replicate as the replicate ID everywhere; never
  the GEO sample name.
- **Exp3 gap: not closed — accepted.** All 5,329 stage-IV cells are in their own matrix;
  matching changes nothing (still 555 / 2,692). 20 stage-IV cell barcodes occur in
  more than one stage-IV library (chance collisions; resolved by replicate). Report
  555 / 2,692 vs 559 / 2,816 (0.7% clones, 4.4% cells) in the Data section.

### Decisions (2026-11-11)
- **Exp3 gap: Option 2, time-boxed.** Match stage-IV cells to the replicate 1–3
  expression matrices using the replicate from the viral barcode. This can resolve
  ambiguous cells but cannot recover cells absent from the deposited table. If the
  gap is not closed after this one pass, stop: accept 555 / 2,692 as the deposited-data
  count, report the 0.7% clone / 4.4% cell discrepancy in the Data section, and (optional)
  ask the corresponding author. Do not tune filters to hit 559.
- Provenance: commit the script first, then run it, so the recorded commit is the one
  that produced the output.
- Scanpy is blocked on the Windows machine (Application Control blocks a numcodecs DLL);
  do not try to bypass the security policy. **Expression/OT work runs in the Claude
  cloud workspace** (Linux), results committed to this repo. Only public SISBAR data
  go there — never FOUNDIN-PD/PPMI data. Large outputs (.h5ad) stay in git-ignored
  `data/`; scripts and small result tables are committed.
- QC: the paper's stated thresholds (≥1,000 genes, ≤30,000 UMIs, ≤15% mt) drop 9.9% of
  the authors' own clone-table cells, nearly all via the UMI cap (40% of stage-IV clone
  cells). **Decided (author, 2026-11-11, Q2 option b + Q3 option i):**
  - Primary cell set: ≥1,000 genes and ≤15% mt; no UMI cap; no doublet removal.
    `results/04_qc/cell_qc.csv.gz` column `keep`.
  - Scrublet is unstable here (0–18% flagged depending on settings), so it is NOT a
    filter. Scores are recorded (`mean_center=False`, `random_state=0`, per library);
    `doublet_flag` = top expected-rate fraction (0.8% per 1,000 barcodes) per library.
  - **Every RQ1–RQ3 analysis is re-run excluding `doublet_flag` cells as a
    pre-specified sensitivity analysis**, reported alongside the primary result.
  - Never tune QC to hit a target count.

### Expression objects (2026-11-11, `scripts/03_build_objects.py`)
- `data/processed/libs/<library>.h5ad` (cloud workspace only, git-ignored): raw counts,
  33,694 genes; EGFP removed from genes and stored as `obs.egfp_counts`.
- `results/03_objects/cell_obs.csv.gz`: all 215,960 barcodes with library, stage, exp,
  rep (viral), QC metrics, `qc_pass` (paper thresholds), clone, cluster. Cell ID =
  `<library>:<core barcode>`.

### Next data tasks
1. ~~Map each `_N` suffix (and stage-IV viral-barcode replicate) to its GSM library.~~ Done.
2. ~~Resolve or document the Exp3 discrepancy (rule above).~~ Accepted as 555 / 2,692.
3. ~~Build the per-clone expression objects (cloud workspace).~~ Done (`d45f505`).
4. Apply the QC rule chosen in Q2; add doublet scores if option (b).

## Research questions

- **RQ1 Coupling recovery:** own-clone OT mass vs chance and vs nearest neighbour;
  clone-retrieval AUROC / top-k recall; fate-level agreement for stage-III progenitors.
- **RQ2 Forecasting on true lineage:** true vs OT vs independent coupling; do models
  beat the constant baseline; does anything beat the random forest.
- **RQ3 State adequacy:** does the earlier sister cell's (P, D) state predict its
  clone's later fate? Compare (P, D), 30-PC PCA, and a lineage-derived fate-bias axis.
- **RQ4 Perturbation-model consequence:** simulate control→treated response with the
  pair-mismatch rate measured in RQ1; report added effect error.

## Non-negotiable protocol rules

- **Split by clone, never by cell.** Sister cells in train and test = leakage.
- Training-partition-only normalisation (min–max bounds, scalers, PCA fits).
- One shared split per seed, used identically by every model.
- 20 seeds; paired bootstrap 95% CIs; Wilcoxon with Holm correction; Cohen's dz.
- Equivalence margin fixed BEFORE analysis; a CI containing zero is not equivalence.
- Always report constant baselines: training median (matched to MAE) and mean.
- Models: Ridge, complexity-matched random forest, Transformer. No new models.
- OT: entropic, in 30-PC PCA of 2,000 HVGs (as in TCBB); fit on all cells, score on
  the barcoded subset only. Sweep entropy and unbalanced growth settings.
- Controls: TCBB simulation (positive); independent pairing and clone-label shuffle
  within stage (negative). Do NOT use the stage-order shuffle (shown insensitive).
- Every number in figures/tables is generated programmatically from saved per-seed
  outputs. Never hand-enter a result.
- Record SHA-256 of every input file and the git commit of each analysis run.

## Paper outline (≈10–12 IEEE double-column pages, 6 figures, 3 tables)

1. Introduction · 2. Related work · 3. Data & setup (Fig 1) · 4. Methods (Table 1)
5. Results: coupling recovery (Fig 2), fate agreement (Fig 3), forecasting (Fig 4,
   Table 2), state adequacy + Jerber transfer (Fig 5), perturbation consequence (Fig 6)
6. Discussion (Table 3 checklist) · 7. Limitations & conclusion

## Open decisions (ask the author, don't assume)

- Journal: TCBB vs JBHI.
- P/D marker panels for a midbrain–hindbrain system.
- Equivalence margin value.
- Whether to pre-register thresholds before RQ2–RQ3.

## Related, separate project — do not mix

FOUNDIN-PD (PPMI Tier 3, Project 150) cell-to-clinic study is a separate,
pre-registered project awaiting data access. Do not pull PPMI clinical outcomes for
FOUNDIN-PD participants before the OSF registration is filed, and never put PPMI
individual-level data in this repository.

## Style

- Writing: IEEE style, plain and specific, effect sizes with intervals before p-values.
- Report null results as fully as positive ones.
