# ANALYSIS_PLAN.md — pre-specified experiments (RQ1–RQ4)

Status: **APPROVED 2026-11-11 by the author** (Victor Ikenna Kanu), including the review
edits in UPDATES.md Q5 (items 1–9), which were merged into the text before approval. The
commit that sets this line is the pre-registration record: every setting below is fixed
before any RQ result exists. Changes after approval are logged as dated amendments at the
bottom, with the reason, and never silently.

Applies CLAUDE.md protocol rules throughout. All runs in the cloud workspace. Scripts are
committed before they run; each writes `results/NN_*/provenance.txt`.

---

## 0. Common settings

- **Cell sets.** Primary = `results/04_qc/cell_qc.csv.gz` `keep`. Sensitivity = `keep` and
  not `doublet_flag`. Every RQ1–RQ3 analysis runs on both; primary is reported first.
- **Transitions.** I→II (Exp 1), II→III (Exp 2), III→IV (Exp 3). Cells of a transition =
  primary cells of both stages from that experiment.
- **Replicates are separate differentiations.** Clones are replicate-specific (viral
  barcode `-R-E`), so every coupling is computed **within replicate**; results are
  reported per replicate and pooled (clone-weighted).
- **Expression space (per transition × replicate).** Library-size normalise to 10,000,
  log1p; 2,000 HVGs (`scanpy.pp.highly_variable_genes`, flavor `seurat_v3` on raw counts,
  `batch_key` = library); scale (max 10); PCA, 30 components; fitted on both stages pooled,
  barcoded and unbarcoded cells alike. Fixed `random_state = 0`.
  *For RQ2–RQ3 the PCA is refit on the training partition only (leakage rule).*
- **Cell-count cap (memory).** Background cells are subsampled to at most 8,000 per stage
  per replicate (seed 0); every cell of a cross-stage clone is always kept. **Background**
  = unbarcoded cells plus barcoded cells whose clone is not cross-stage; background cells
  are never positives. (The cap binds only in I→II and in II→III replicates 2–3;
  `results/05a_plan_checks/C_leakage_exposure.csv`.)
- **Fate labels.** The authors' `Cluster.Id` exists only for clone-table cells. Unbarcoded
  cells get a label by kNN vote (k = 15) from clone-table cells of the same stage and
  replicate in the 30-PC space; used only for RQ1 fate-level agreement.
- **Uncertainty.** 95% percentile bootstrap over clones (1,000 resamples) for RQ1 and RQ4;
  20 clone-level split seeds for RQ2–RQ3 (CLAUDE.md statistics).

## 1. RQ1 — coupling recovery

**Question.** How much of an early cell's coupling mass lands on its own clone's later
sister cells?

**Couplings** (per transition × replicate, cost = squared Euclidean distance in 30-PC
space divided by its median):
1. Entropic OT, balanced, ε ∈ {0.005, 0.01, 0.05, 0.1} (POT, log-domain Sinkhorn).
2. Entropic OT, unbalanced, same ε, marginal relaxation τ = 0.9 (`reg_m = ε·τ/(1−τ)`).
3. kNN baseline: each early cell spreads its mass uniformly over its k = 10 nearest later
   cells.
4. Independent (chance): mass proportional to cell counts.

The full ε grid is reported. **No setting is chosen by clone recovery.** For RQ2 and RQ4
the OT coupling used is fixed now: **balanced, ε = 0.05.**

**Metrics** (barcoded early cells whose clone has later sister cells; direction forward,
row-normalised early → later; backward reported as secondary):
- **Own-clone mass** m = mass on own-clone later cells; reported with chance m₀ (own-clone
  share of later cells) as enrichment m/m₀ and difference m − m₀.
- **Clone-retrieval AUROC:** later cells ranked by coupling weight; positives = own clone.
- **Top-k recall,** k ∈ {1, 5, 10}: an own-clone cell among the k heaviest later cells.
- **Fate-level agreement** (Fig 3): for each early-stage cluster with ≥ 20 cross-stage
  clones, total-variation distance between the OT-predicted later-cluster distribution and
  the clonal one (measured sister-cell clusters). Priority cluster: III-NPvMB (the
  midbrain DA progenitor).

**Chance and ties.** Chance is computed analytically per early cell from the pool actually
used (capped), per replicate: m₀ = s/n; top-k chance = 1 − C(n−s, k)/C(n, k) (s = own-clone
later cells, n = later cells). AUROC is the Mann–Whitney statistic with ties counted as ½,
so the independent coupling has AUROC = 0.5 exactly. For kNN, top-k (k ≤ 10) is the
hypergeometric expectation within the 10 equally weighted neighbours. No random
tie-breaking.

**Cap robustness (pre-specified).** I→II replicate 2, where the cap changes chance most
(24,391 later cells uncapped), is rerun without the cap for the balanced ε = 0.05 coupling
and kNN; reported as a supplementary comparison.

**Expected** (stated before running): OT above chance but far from perfect; lowest at
III→IV, where fates branch; OT not clearly better than kNN.

**Outputs:** `results/06_rq1/` — per-cell metrics, per replicate × transition × coupling
summary with CIs, fate-agreement table, provenance.

## 2. RQ2 — forecasting on true lineage

**Question.** With real lineage, do models beat the constant baseline, and does anything
beat the random forest?

- **State** (2-D, TCBB-comparable): **D** = mean log1p of {TH, DDC, SLC6A3, DRD2, LMX1A,
  FOXA2} (unchanged from DT1/TCBB). In this system D is carried by FOXA2, DDC and LMX1A
  (SLC6A3 ≤0.5%, DRD2 ≤4.5% detected) and is zero in 30–37% of later cells, so the paper
  describes it as a floor-plate/ventral-midbrain identity score, not a DA-neuron score.
  **P** is replaced by a **neural-progenitor score Pg** = mean log1p of {SOX2, NES, HES1,
  HES5, VIM}, because the pluripotency genes other than SOX2 are silent from day 13.
  Scores = mean log1p(CP10k) over panel genes, as in TCBB. Min–max scaled on training
  cells only, clipped to [0, 1]. (Evidence: `results/05a_plan_checks/A_*`, `B_*`.)
- **Task.** Early-cell state → later sister-cell state (one step).
- **Pairs per seed.** One pair per cross-stage clone (one early and one later cell drawn
  at random within the clone).
- **Training couplings.** (a) true clonal pairs; (b) OT pairs: for each training early cell,
  a later cell drawn from its OT row (balanced, ε = 0.05); (c) independent: a random later
  cell from the same replicate. **Draws in (b) and (c) are restricted to later cells of
  training clones plus background cells** (OT row renormalised over them), so no
  validation/test sister cell is ever a training target. **All are tested on held-out true
  clonal pairs.** The coupling's PCA is fitted on all cells (unsupervised, as in RQ1); this
  is disclosed.
- **Split.** 70/15/15 by clone, stratified by replicate; 20 seeds; one split per seed shared
  by all models and couplings.
- **Models.** Constant median, constant mean, Ridge (α = 1), complexity-matched random
  forest (3-fold CV over depth and leaf size), Transformer as in TCBB. With one input
  token the Transformer reduces to a feed-forward network; it is kept for continuity and
  this is stated in the paper.
- **Metric.** MAE averaged over the two state components.
- **Statistics.** Paired Δ vs constant median and vs random forest; percentile paired
  bootstrap 95% CI; Wilcoxon, Holm-corrected; Cohen's dz. **Equivalence margin
  δ = 0.01 MAE** (same as TCBB). Equivalence only if the whole CI lies within ±δ. δ is
  8.6–12.4% of the constant-median MAE of (Pg, D) on later cells of cross-stage clones
  (`B_score_spread.csv`); every Δ is also reported as a % of the constant-median MAE.
- **Expected:** true and OT pairs beat the constant; independent pairs do not
  (Proposition 1); no Transformer advantage.

**Outputs:** `results/07_rq2/`.

## 3. RQ3 — state adequacy

**Question.** Does an early cell's state predict its clone's later fate?

- **Label.** A clone's fate = majority cluster of its later sister cells. Fates with
  ≥ 30 clones per transition only; others pooled as "other". Classes per transition
  (primary set, current counts): I→II NPaHB, NPvMB, Motors, NPMHB, NPpHB, Gluts + other;
  II→III NPaHB, NPvMB, Motors, NPpHB, Gluts + other; III→IV hGluts, Seros, LUM+VLMCs,
  Motors, OTX2+VLMCs, GABAs + other.
- **Representations** of the early cell: (1) original (P, D) with the TCBB pluripotency
  panel P = {POU5F1, NANOG, SOX2, UTF1, TDGF1}, kept for continuity; here it reduces to a
  SOX2 score (the other four genes ≤0.8% detected at every stage), and SOX2 is shared with
  Pg; (2) (Pg, D); (3) 30-PC; (4) **fate-bias axis** =
  LDA on 30-PC fitted on training clones only (K−1 dims).
- **Classifier.** Multinomial logistic regression (L2, C = 1), identical for every
  representation, so only the representation differs.
- **Split.** 70/15/15 by clone, 20 seeds. Early cells of one clone share the clone's label;
  a test prediction is the mean probability over the clone's early cells.
- **Metrics.** Macro-AUROC (primary) and balanced accuracy, vs a stratified chance
  baseline. Paired comparisons between representations, Holm-corrected.
- **DA lineage (secondary, descriptive).** Only 6 III→IV clones end as mDA. Instead, for
  III-NPvMB clones, report the measured fate distribution and whether (Pg, D) or the
  fate-bias axis separates mDA-producing from other NPvMB clones; no inferential claim.
- **Jerber transfer.** Separate later instruction (needs the Jerber data downloaded).
- **Expected:** (P, D) near chance; 30-PC and the fate-bias axis well above it.

**Outputs:** `results/08_rq3/`.

## 4. RQ4 — consequence for perturbation models

**Question.** If a perturbation model pairs control and treated cells with OT-like
accuracy, how much error does mismatching add to the estimated response?

- **Simulation** on II→III later-stage cells (largest set), per replicate: "control" =
  real cells; "treated" = the same cells with a known heterogeneous effect added to 200
  HVGs, effect size linear in PC1 (mean shift 0.5 log units, SD across cells 0.25).
- **Pairing with mismatch.** Each treated cell is paired with its true control, except a
  fraction m reassigned to a random control. m ∈ {0, 0.1, …, 1.0}, plus the mismatch rates
  measured in RQ1 (OT balanced ε = 0.05, per transition): **primary 1 − own-clone mass**,
  the probability that a partner drawn from the OT row (as in RQ2) is not a sister;
  secondary 1 − top-1 recall (argmax pairing). Sisters are not the same cell, so these are
  upper bounds on the pairing error relevant to perturbation models.
- **Outputs measured.** Per-cell effect RMSE and correlation with the true effect;
  population mean-effect bias. 100 simulation seeds.
- **Expected:** mean effect robust to mismatch; per-cell (heterogeneous) effects degrade
  steeply; RQ1's measured rate sits in the degraded range.

**Outputs:** `results/09_rq4/`.

## 5. Sensitivity and reporting

- RQ1–RQ3 rerun on the sensitivity cell set; reported as a supplementary table, and any
  conclusion that changes is stated in the main text.
- Every figure and table is generated from saved per-seed outputs by scripts in
  `scripts/`; numbers in the paper come from generated `.tex` macros.

## Amendments

_None._
