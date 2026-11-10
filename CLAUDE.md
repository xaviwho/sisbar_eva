# CLAUDE.md — Lineage-grounded evaluation (SISBAR) project

Context for Claude Code sessions on this repository. Author: Victor Ikenna Kanu,
Network Systems Lab, Kumoh National Institute of Technology (PI: Prof. Dong-Seong Kim).

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

| Transition | Days     | Cross-stage clones | Cells | GEO                      |
|------------|----------|--------------------|-------|--------------------------|
| I → II     | 13 → 21  | 768                | 2,762 | GSE207921 (replicate 1)  |
| II → III   | 21 → 30  | 1,213              | 5,856 | GSE221592 (replicates 2–3) |
| III → IV   | 30 → 45  | 559                | 2,816 | GSE221592                |

- 173,564 cells total; 26% barcoded, in 10,702 clones; 2,540 clones span two stages.
- **No clone spans three stages.** All forecasting is one-step (t → t+1).
- Ground truth is **sister cells**: clones split after ~3–4 divisions, one half
  profiled early, the other later. Not the same cell.
- Stage IV contains midbrain DA neurons (NR4A2+/PITX3+/EN1+). DA-specific clones are
  few (72 from the stage-III ventral midbrain progenitor) → DA analyses are secondary.
- Supporting data: Jerber et al. 2021 (D11/D30/D52) for transfer of the fate axis.

### First data tasks (do these before any modelling)
1. Download both GEO series; confirm whether barcode-to-cell tables are deposited
   and which samples map to which experiment.
2. Reproduce the published clone counts (768 / 1,213 / 559) as a sanity check.
   If they don't match, stop and report — don't tune until they do.

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
