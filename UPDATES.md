# UPDATES.md — two-way project log

**How this works**
- **Author → Claude:** add items under *Open instructions* (newest at the bottom), commit, push.
- **Claude:** at the start of every session, `git pull`, read this file, work the open
  instructions in order, then move each finished one to *Done* with the date and commit.
- **Claude → author:** every work session ends with a dated entry in *Log* (newest first)
  and any blocking decisions under *Open questions*. Then commit and push.
- An answer to an open question goes under that question (`**Answer:** …`) or as a new
  instruction. Claude moves answered questions to the log entry that acts on them.
- CLAUDE.md stays the standing project context and rules; this file is the running
  conversation. Durable facts and decisions get copied into CLAUDE.md as well.

---

## Open instructions (author → Claude)

- [ ] 2026-11-17 — **(Cloud) Post-review analyses: ANALYSIS_PLAN amendment A3, items
  A3.1–A3.6.** Approved by the author ("yes go ahead").
  - **The rule:** every setting is fixed in A3. Follow it exactly; any deviation becomes a
    new dated amendment. Commit each script before running it; each writes
    `results/NN_*/provenance.txt`. Primary cell set only. Run one OT job at a time (5.8 GB
    limit).
  - **Order (cheapest first):**
    - A3.1 RQ2 with 30-PC input → `scripts/12_rq2_pc_input.py`, `results/12_rq2_pc/`. Reuse
      the 07_rq2 splits, sampled cells and targets per seed; save per-seed MAE for both
      inputs.
    - A3.4 repeated target draws → `scripts/13_rq2_draws.py`, `results/13_rq2_draws/`.
    - A3.3 leave-one-replicate-out (RQ2 and RQ3) → `scripts/14_loro.py`, `results/14_loro/`.
    - A3.5 fate agreement on measured labels and label-transfer k → `scripts/15_fate_labels.py`,
      `results/15_fate_labels/`. Reuse the saved RQ1 couplings and potentials.
    - A3.6 RQ4 effect-model sensitivity → `scripts/16_rq4_sensitivity.py`,
      `results/16_rq4_sens/`.
    - A3.2 inductive RQ2 (OT arm) → `scripts/17_rq2_inductive.py`, `results/17_rq2_inductive/`.
      Most expensive. Apply the seed 0–9 fallback only on seed-0 runtime, and log the
      decision before reading any result.
  - **Outputs:** per-seed or per-fold CSVs plus a summary CSV per item. Do not write paper
    text; VS Code writes the paper from these outputs.
- [ ] 2026-11-17 — **(Cloud) Public code release.** The repo `github.com/xaviwho/sisbar_eva`
  is public (author). Push everything needed to reproduce the paper:
  - **Scripts:** every script used in the cloud workspace, including helpers. Make sure no
    script exists only in the workspace.
  - **Environment:** `environment/requirements-cloud.txt` (`pip freeze` of the cloud
    environment, plus the Python version) next to the existing local `requirements.txt`.
  - **Data download:** `scripts/00_download_geo.py` fetches GSE207921/GSE221592 and the
    clone table and checks the SHA-256 values recorded in the `provenance.txt` files.
  - **README:** a reproduction guide covering data download, script order (01 → 17), which
    steps run where (local vs Linux), what each `results/NN_*` contains, and the analysis
    plan with its dated changes. State that the plan was committed at `65cd70d`, before any
    result.
  - **Never commit** `data/`, `.h5ad` or other large intermediates. Never commit
    FOUNDIN-PD/PPMI material.
  - **After the A3 results land:** tag `v1.0-submission`.
- [ ] 2026-11-17 — **(VS Code) A3.7 single-barcode sensitivity** (local, from saved per-cell
  metrics) and, once the cloud A3 results land, the paper text for A3.1–A3.7 (supplement
  first; main-text paragraphs where they change a claim). If A3.1 shows that 30-PC input
  does not beat (Pg, D) by more than δ, revise the "marker state limits forecasting" claim.
- [ ] Later — Jerber D30 transfer of the RQ3 fate-bias axis (needs Jerber data; separate
  plan amendment).

## Open questions (Claude → author)

- **Q11 (2026-11-17) — licence and archive for the public release.**
  - (a) Licence: MIT is recommended for code. Results and paper files could be CC BY 4.0.
  - (b) Zenodo DOI: enable the GitHub–Zenodo integration on your account
    (zenodo.org → GitHub → switch on `sisbar_eva`). The `v1.0-submission` tag then gets a
    citable, archived DOI. Only you can do this step.

- **Q6 (2026-11-14, non-blocking) — reference for the TCBB benchmark paper.** Results cite
  Proposition 1 as `\cite{Kanu_TCBB}`; the bib entry is a placeholder (`@misc`, title and
  code URL only). Please give the status/volume/DOI to use.
- **Q9 (2026-11-16) — controls named in CLAUDE.md but not in the approved plan.** The
  CLAUDE.md protocol rules list a TCBB-simulation positive control and a within-stage
  clone-label shuffle negative control. Neither is in `ANALYSIS_PLAN.md`, so neither was
  run. Methods names only what was run: the independent coupling is the negative control,
  and the chance levels are analytic.
  Options:
  - (a) Leave as is, and drop or reword those lines in CLAUDE.md.
  - (b) Add them as a new amendment (A3) and a cloud run, labelled as added after the
    results.

  **Recommendation:** (a). The analytic chance levels and the independent coupling already
  play the negative-control role, and a control added after the results needs to be
  labelled as such.

## Done

- [x] 2026-11-16 — (VS Code) Q7 → amendment A2; Q8 → Methods section with equations and
  settings table (commit below).
- [x] 2026-11-14 — (VS Code) Results for RQ1–RQ4 written in `paper/main.tex`, with Figs 2–6,
  Table 2 and supplementary Tables S1–S5 (sensitivity per RQ, cap check, per replicate),
  all generated by `scripts/10_results_paper.py` (`42c476f`). Covers the instructions of
  2026-11-11, 2026-11-12 and 2026-11-13.
- [x] 2026-11-11 — RQ1 coupling recovery (cloud), primary + sensitivity + uncapped cap
  check (amendment A1 `4ad508f`, validation `fbca3d0`). → scripts `12fade6`, `22dad01`,
  `b576f19`; final results in the 2026-11-13 commit.
- [x] 2026-11-11 — RQ2 forecasting (cloud), primary + sensitivity. → scripts `fd29545`,
  `dfdb11c`, `a46d264`; final results in the 2026-11-13 commit.
- [x] 2026-11-12 — RQ4 perturbation consequence (cloud). → script `fe16ad2`, results this commit.
- [x] 2026-11-12 — RQ3 state adequacy (cloud). → script `6bde8cb`, results `fa1564c`.

- [x] 2026-11-10 — Download both GEO series; confirm barcode-to-cell tables and sample
  mapping. → `c738aea`
- [x] 2026-11-10 — Reproduce published clone counts. → `c738aea`
- [x] 2026-11-11 — Map `_N` suffixes / stage-IV replicate to GSM libraries. → `ef11252`
- [x] 2026-11-11 — Resolve or document Exp3 discrepancy (one pass). → `ef11252`
- [x] 2026-11-11 — Q1 (where expression work runs): cloud workspace. → `56ddc72`
- [x] 2026-11-11 — Build per-library expression objects with QC metrics. → `d45f505`
  (script), results in `460fdf0`.
- [x] 2026-11-11 — Q2/Q3 QC rule applied; final cell set. → `92a0df7` (script), results in
  the commit after it.
- [x] 2026-11-11 — Draft Section 3 (*Data and setup*) in LaTeX with generated numbers and
  table. → `f0a8a37` (generator), `d640917` (section, table, numbers, bib).
- [x] 2026-11-11 — Author approved `ANALYSIS_PLAN.md` with Q5 items 1–9 merged. →
  `65cd70d` (pre-registration record).
- [x] 2026-11-16 — Data and Setup merged into Methods; supplementary tables moved to
  `paper/supplement.tex`; figures redesigned. → `416e055` (figure script) and the
  commit after it.

---

## Log (newest first)

### 2026-11-17 (night) — Claude (VS Code) — Peer review triage; amendment A3; wording fixes
- **Author shared an AI peer review** (3/5, major revision) and approved the plan: wording
  fixes now, cloud runs A3.1–A3.6, and a public code release. The repo is public.
- **Checked before acting:**
  - Multi-barcode clones are 2.7 / 6.5 / 12.8% of early cells in cross-stage clones.
    Restricting to single-barcode clones leaves OT AUROC at 0.609 / 0.643 / 0.601 (vs
    0.607 / 0.644 / 0.593). Now formal as A3.7 (`f09a46b`, `results/11_single_barcode/`).
    Its "all clones" rows reproduce the RQ1 summary exactly.
  - The review's "inconsistent ranges" point is a misreading: true lineage is 0.9–3.1%
    throughout, and 0.9–4.2% is the OT range.
- **Amendment A3** written in ANALYSIS_PLAN before any run (A3.1 30-PC RQ2 input; A3.2
  inductive OT arm; A3.3 leave-one-replicate-out; A3.4 repeated target draws; A3.5 fate
  agreement on measured labels; A3.6 RQ4 effect-model sensitivity; A3.7 single-barcode).
  The preliminary A3.7 check is disclosed there.
- **Author asked why "amendment A1/A2" appeared in the paper.** They are internal labels
  from ANALYSIS_PLAN.md and meaningless to readers.
  - Removed everywhere: main text, Table S4 caption, supplement. Each is replaced by a
    plain description.
  - "Pre-registered" became "pre-specified" / "time-stamped": the plan was never filed
    with a registry, it was committed to the repo at `65cd70d` before any result.
  - Methods now states the repository URL and that commit.
- **Wording fixes from the review:**
  - The abstract, Discussion and Conclusion narrow the population-level claim to this
    system and to cluster-by-cluster checking.
  - RQ4 is framed as a conditional synthetic propagation analysis, with bias intervals
    reported.
  - The 20-split intervals are called "split-stability intervals".
  - The replicate structure is stated (one experiment and three biological replicates per
    transition).
  - The rationale for δ is given (1% of the scaled range).
  - The LDA axis is called "supervised".
  - Matched split, cells and sisters across pairing constructions are stated.
- **Not yet changed:** the claim "marker state, not coupling, limits forecasting" stays as
  is until A3.1 lands.
- **Build:** both PDFs clean, main 10 pages; no internal labels in the PDF text.
- **New:** Q11 (licence; Zenodo integration).

### 2026-11-17 (late) — Claude (VS Code) — Cited the ICUFN 2026 paper
- **Author:** asked whether https://ieeexplore.ieee.org/document/11628674 can be cited and
  shared the PDF. It is Kanu et al., "Transformer-Based Digital Twin for Stem Cell
  Differentiation Forecasting with Mechanistic ODE Benchmarking", ICUFN 2026, pp. 471–476,
  DOI 10.1109/ICUFN69619.2026.11628674 (verified on Crossref). Bib key `Kanu_ICUFN2026`.
- **Cited where it is the source:**
  - Intro: it defined the two-score (P, D) state and forecast it from population-level
    pseudo-trajectories (cells sampled independently per stage); Kanu_TCBB then showed the
    constant is optimal under such pairing.
  - Methods: the pluripotency and dopaminergic panels (P = POU5F1, NANOG, SOX2, UTF1,
    TDGF1; D = TH, DDC, SLC6A3, DRD2, LMX1A, FOXA2), identical in ICUFN.
- **Noted:** the author's commits `d2142de`…`6312eed` remove the italics from the
  contributions, comment out the roadmap sentence, drop the "(RQn)" labels from Methods
  headings, and change the corresponding e-mail to dskim@kumoh.ac.kr. That style is kept.
- **Build:** clean, 10 pages.

### 2026-11-17 (late) — Claude (VS Code) — Removed "TCBB" labels from the manuscript
- **Author:** "we can't call it that … shouldn't have such puerile things."
- The 4 visible uses were "TCBB panel" (Methods ×2), "TCBB (P, D)" (RQ3 list) and the Fig.
  caption. They now read "the original (P, D) state of [4]", "the dopaminergic panel of
  [4]" and "original panels [4]". The compiled PDF contains no "TCBB"; the cite key is
  invisible to readers.
- New rule in CLAUDE.md (Style): no internal shorthand in the manuscript.

### 2026-11-17 (late) — Claude (VS Code) — Limitations as prose
- **Author:** "write the limitations as a flowing statement instead of listing them."
- §5.4 is now four paragraphs of prose with the same six points and macros: the
  ground truth (sisters) → incomplete lineage information → generality (one line, protocol,
  laboratory; fixed OT configuration) → the marker state, the DA analysis and the simulated
  perturbation.
- Also pulled the author's acknowledgments (`e71379d`). Build clean, 10 pages.

### 2026-11-17 (night) — Claude (VS Code) — Abstract, index terms, author block, restructure
- **Author:**
  - correct the Williams name and add his affiliation;
  - write the abstract and check the index terms;
  - make Limitations a Discussion subsection;
  - give the Conclusion its own section, with future work as its second paragraph.
- **Q10, partly answered:**
  - Done: "Willaims" → "Williams".
  - Done: a second thanks item: "W. P. Nwadiugwu is an Associate Professor with the Bharti
    School of Engineering and the Department of Biochemistry, Laurentian University,
    Sudbury, ON P3E 2C6, Canada". "Sudbury" was added because it is the city of that
    postal code.
  - Done: the running head is now "Kanu *et al.*" (three authors).
  - Still open: the grade format "(MIEEE)/(SMIEEE)" is left as the author wrote it; IEEE
    style is ", Member, IEEE" / ", Senior Member, IEEE", or none at peer review.
- **Abstract:** 187 words. All numbers are macros (AUROC range, ≤ max mass %, RQ2 gain
  range, RQ3 marker vs 30-PC, RQ4 statement).
- **Index terms:** now alphabetical, as IEEE requires. "Benchmarking" became "Benchmark
  testing" (IEEE Thesaurus term), and "cell-state forecasting", a coined phrase, became
  "cell fate prediction". Final list: Benchmark testing, cell fate prediction, lineage
  tracing, neural differentiation, optimal transport, single-cell RNA sequencing.
- **Restructure:**
  - Limitations is now §5.4, the last subsection of the Discussion.
  - §6 Conclusion is its own section: a findings paragraph, then a future-work paragraph
    with four directions (lineage-informed couplings, other systems and OT variants, the
    Jerber transfer, real perturbation data).
  - The Introduction's roadmap sentence is updated.
- **Bibliography:** scDisentangle author names fixed ("S. El Mejjad", "R. El Fatimy",
  "F. Z. El Hajji").
- **Build:** clean, 10 pages including references. Remaining: author biographies (IEEE
  needs a photo and about 100 words per author).

### 2026-11-17 (evening) — Claude (VS Code) — Limitations and Conclusion written
- **Author:** "write the limitations and conclusion section."
- **Pulled `cae3f47`, `1186581`:** the author added co-author Willaims Paul Nwadiugwu and
  IEEE grades, and moved the figure floats. Left as is. Points to check before submission
  are listed under Q10.
- **§6 Limitations (six numbered):**
  1. Sisters, not descendants (the RQ4 rates are upper bounds; standard scRNA-seq cannot
     observe the same cell twice).
  2. One cell line, protocol and laboratory.
  3. Incomplete lineage (background cells, the deposited-table gap, the pool-size
     dependence of absolute metrics; AUROC is robust to the cap).
  4. Fixed representations (one PCA space, one τ, no growth priors).
  5. The two-score marker state (D is an identity score; DA analysis descriptive; Jerber
     transfer pending).
  6. A simulated perturbation.
- **Conclusion:** one paragraph of findings, recommendations and next steps (lineage-
  informed couplings, other systems, the pre-specified Jerber transfer).
- Numbers are macros; the claims are softened to what the data show.
- **Build:** clean, 10 pages. Remaining: abstract and author biographies (~0.5 page, so
  ≈11 / 12).
- **Q10 (new, non-blocking) — author block.**
  - Is "Willaims" a typo for "Williams"?
  - The new co-author has no affiliation in the footnote.
  - IEEE/compsoc style puts membership grades as ", Member, IEEE" / ", Senior Member,
    IEEE" after each name, not "(MIEEE)" / "(SMIEEE)". Peer-review submissions often omit
    them.

### 2026-11-17 (later) — Claude (VS Code) — Related work checked; Introduction written
- **Author:** "read the related works and the papers cited, if it is good, add the
  introduction and three to four contributions enumerated."
- **Related work (author, `dfe8ba6`):** checked every cited reference on Crossref/OpenAlex
  and read the abstracts. It is good overall, with three fixes made:
  - **Rosenblatt et al.:** the leakage that "drastically inflates" performance is feature
    selection and repeated subjects; covariate correction had minor effects, so the sentence
    was narrowed.
  - **"Compare with simple baselines":** Li et al. 2026 supports the generalisation point
    (no method transferred across cell types) but not this one. It is now attributed to
    scDisentangle (whose abstract says so) and to Kanu_TCBB.
  - **Bibliography:** scGen and Rosenblatt now point to their published versions (Nature
    Methods 2019; Nature Communications 2024), with keys unchanged. Volume and article
    numbers were added for Li (Sci. Adv. 12(37)) and Nieto (Neurocomputing 680).
  - 12 of the 21 new bib entries are not cited (harmless; bibtex skips them). Some of them
    have weak metadata (TrajectoryNet listed with journal "PubMed"); fix before citing.
- **Introduction (§1):**
  - the snapshot problem, OT, and marker states;
  - the two assumptions, and why they go untested;
  - SISBAR as ground truth;
  - RQ1–RQ4;
  - four enumerated contributions: OT tested against lineage; the marker state, not the
    coupling, limits forecasting; the cost of pairing error; the protocol and checklist;
  - a roadmap sentence.

  All numbers are generated macros. "No model beat the RF" refers to true-lineage pairs and
  is guarded by an assert.
- **Figures:** PDFs are now written without a timestamp, so reruns are byte-identical
  (script `cc95506`, committed before the run).
- **Build:** clean, 9 pages. Still to write: Limitations and Conclusion, abstract,
  biographies (≈1 page, so about 10 / 12).

### 2026-11-17 — Claude (VS Code) — Discussion written
- **Author:** "let's have the discussion section done."
- **Pulled `2e1ae54`:** the author removed "(RQ3)" from the Methods heading *Fate
  prediction*. The other headings are left as they are.
- **Discussion (§5 of `paper/main.tex`, about one page), in three parts:**
  - (i) *Couplings recover populations, not cells*: above-chance sister mass and AUROC,
    fate composition for 13/17 clusters, but at most 0.26% mass on sisters; kNN gets
    similar mass; RQ4 shows the per-cell error equals random pairing.
  - (ii) *Marker states are a bottleneck*: (Pg, D) is above chance but below 30-PC; the
    panel choice matters little; the RQ2 gains under exact pairing lie inside δ, so the
    limit is the state.
  - (iii) *Evaluation practice*: the kNN AUROC artefact, chance depending on pool size,
    and constants fitted to a coupling's own targets.
- **Table 4 (checklist, nine rows):** each check is paired with its evidence in this
  study, and each cites its figure or table.
- **Numbers:** all come from macros. The generator gained cross-transition range macros
  (script `5a4b09c` committed before the run). The claims "every transition" and "every
  comparison" rest on existing `check_wording` asserts.
- **New reference:** Weinreb et al. 2020 (Science; Crossref metadata), the clonal
  haematopoiesis benchmark of state→fate.
- **Build:** clean, 8 pages. Still to write: Introduction, Related work, Limitations and
  Conclusion, abstract, biographies (≈2.5 pages, so about 10.5 / 12).

### 2026-11-16 (later) — Claude (VS Code) — Structure, length, figures
- **Author:** "data and setup should be part of the method … TCBB is constrained on page
  length … not satisfied with the result figure quality."
- **Page budget:** IEEE Computer Society overlength policy (MOPC) for TCBB regular papers:
  12 formatted pages including references and author biographies; each extra page is
  charged and 18 pages is the hard cap. Supplementary files have no limit. (Some older
  listings say 14 pages; check the current author guidelines before submission.)
- **Structure:**
  - Data and Setup is now Methods §3.1, *Lineage data and ground truth*. It is shortened
    and keeps every generated number and the data table.
  - Tables S1–S5 moved to a separate `paper/supplement.tex` (3 pages); main.tex cites them
    as "Table S1" etc.
  - The RQ1 robustness paragraph is condensed.
- **Length:** main.tex is now 7 pages. The Introduction, Related work, Discussion,
  Conclusion, abstract and biographies remain, at about 3–3.5 pages, so the estimate is
  ≈10.5 / 12.
- **Figures (`scripts/10_results_paper.py`, committed `416e055` before regeneration):**
  - Drawn at exact print size (text width 7.17 in and column 3.50 in, measured from
    IEEEtran compsoc) and never rescaled. No bbox cropping, so the PDF size equals the
    figure size.
  - Arial 6.5–8 pt, matching compsoc's sans-serif captions; bold panel letters; no
    gridlines; outward ticks.
  - Okabe–Ito colour-blind-safe palette with a fixed meaning across figures (OT blue,
    unbalanced vermilion, kNN green, independent grey, true lineage black), plus a distinct
    marker per series.
  - Fig 3 is now a dumbbell (independent → OT per cluster).
  - Captions are updated to refer to the panel letters, and all fonts are embedded.
- **Note:** the F: drive disconnected mid-session; no work was lost.

### 2026-11-16 — Claude (VS Code) — Methods written; amendment A2
- **Author answers:**
  - Q7: "I agree" → amendment A2 added to `ANALYSIS_PLAN.md`.
  - Q8: "focus on the methodology for now … engineering, with mathematical equations …
    not too verbose" → Methods only.
- **A2 (clarification, no change to settings, scripts or results).** RQ3 breaks fate-label
  ties alphabetically. I re-verified that the tie-break alone explains the GABAs class:
  32 / 37 clones alphabetically vs 28 / 27 by first occurrence.
- **Methods (§4 of `paper/main.tex`), written from the scripts, not the plan text alone:**
  - 4.1 Notation and expression space.
  - 4.2 Couplings: cost; balanced entropic OT and log-domain Sinkhorn updates; unbalanced
    KL form with the τ-scaled updates; kNN and independent couplings; the A1 block-wise
    solver.
  - 4.3 RQ1 metrics: mass, chance, enrichment, tie-aware AUROC, top-k chance, clone
    weighting, fate TVD.
  - 4.4 RQ2: panel scores, restricted target sampling, MAE, models, and why a one-token
    Transformer is a feed-forward network.
  - 4.5 RQ3: fate label, LDA, multinomial LR.
  - 4.6 RQ4: effect model and the decomposition showing mismatch adds variance but no bias.
  - 4.7 Statistics: paired d̄, d_z, bootstrap, Wilcoxon–Holm families, equivalence.

  Also added: a settings table (`tab:methods`; numbered Table 2 for now, after the data
  table) and 7 references (Cuturi, Schiebinger, Chizat, POT, Scanpy, Seurat v3, Holm),
  with DOIs from Crossref where they exist.
- **Build:** the PDF builds cleanly, 10 pages.
- **New:** Q9 (controls listed in CLAUDE.md that were not in the plan).
- **Still open:** Q6 (TCBB reference record).

### 2026-11-14 — Claude (VS Code) — Results written (RQ1–RQ4)
- **Generator:** `scripts/10_results_paper.py`. It was committed before the run (`42c476f`;
  provenance in `results/10_results_paper/`) and reads only `results/06–09`. It writes:
  - 448 macros to `paper/generated/results_numbers.tex` (also as CSV);
  - Figs 2–6 to `paper/figures/`;
  - Table 2 and Tables S1–S5 to `paper/tables/`.

  No number in Results is typed. Claims worded as "all/every/none" are guarded by asserts
  (`check_wording`), so the run fails if the data stop supporting the wording.
- **Results text, following the cloud cautions:**
  - **RQ1 kNN:** its AUROC ≈ 0.50 is a ranking artefact (only 10 non-zero weights), while
    its own-clone mass and top-10 recall are similar to OT's.
  - **RQ2:**
    - Every model gain under true lineage is significant but inside ±δ, so it counts as
      equivalent.
    - The III→IV OT/independent result is explained (their constant is fitted to
      background targets), and the cross-coupling contrast is labelled descriptive and not
      pre-specified.
    - The one other non-equivalent comparison is named.
  - **RQ3:** reported as contradicting the expected "near chance" (marker states are above
    chance but 0.08–0.15 macro-AUROC below expression). The DA analysis is descriptive
    only; note that its LDA axis is fitted on all III→IV clones.
  - **RQ4:** the high Pearson r at m = 1 is explained as a property of the simulation
    design.
  - Sisters ≠ the same cell, so rates are upper bounds.
- **Supplementary tables:**
  - S1 RQ1 sensitivity;
  - S2 RQ2 sensitivity (Δ as % of the constant MAE, † = non-equivalent);
  - S3 RQ3 sensitivity;
  - S4 cap check (cites A1 and the validation);
  - S5 RQ1 per replicate, forward and backward.
- **Bibliography:** added Jerber 2021 (Crossref metadata) and a placeholder `Kanu_TCBB`
  (Q6).
- **Build:** the PDF builds cleanly (7 pages; no undefined references or overfull boxes).
  Figure numbers will shift by one once the data-overview Fig 1 exists.
- **Housekeeping:** `git pull`/`commit` print a harmless "geometric repack failed: File
  exists" warning (a stale pack file in `.git/objects/pack`); the commits themselves are
  fine.

### 2026-11-13 — Claude (cloud workspace) — RQ1/RQ2 sensitivity and cap check done; all RQs final
- **RQ1 complete:** 19 jobs (9 primary, 9 sensitivity, 1 uncapped), 145 OT solves, all
  converged (`results/06_rq1/solves.csv`).
- **RQ1 sensitivity (doublet-flagged cells excluded) changes no conclusion.** Pooled, OT
  balanced ε = 0.05, forward (`enrichment` column): own-clone mass 4.67× / 4.13× / 2.74×
  chance (I→II / II→III / III→IV; primary 3.99× / 4.17× / 2.78×); AUROC 0.61 / 0.65 / 0.61
  (primary 0.61 / 0.64 / 0.59); kNN AUROC ≈ 0.50 throughout. Fate TVD keeps the primary
  pattern: OT below independent for most I→II and II→III clusters (I-NPvMB 0.06 vs 0.46,
  II-NPvMB 0.10 vs 0.39), but worse for I-NbpHB (0.38 vs 0.25) and II-Motors (0.62 vs 0.30);
  at III→IV OT is worse than independent for III-NPpHB (0.58 vs 0.53) and III-NPaHB (0.32 vs
  0.13), better for III-NPvMB (0.43 vs 0.65), III-Motors (0.24 vs 0.61) and III-Gluts
  (0.29 vs 0.45).
- **Cap robustness (ANALYSIS_PLAN §1; amendment A1).** Uncapped I→II replicate 2 (24,391
  later cells, blocked solver validated against POT, rel. plan difference 6e-6): AUROC
  0.654 vs 0.656 capped; enrichment 2.16× vs 3.48×; top-1 recall 0 vs 0.34%. Chance m₀ falls
  from 2.7e-4 to 9.2e-5 with the larger pool, so absolute mass falls too. The cap does not
  drive the AUROC result; top-k and absolute mass depend on pool size and should be read
  together with chance. kNN is unchanged (AUROC 0.50).
- **RQ2 sensitivity changes no conclusion.** With true lineage, Ridge/RF/Transformer beat
  the constant median by 0.9–5.0% of its MAE (primary 0.9–3.1%), every CI inside ±δ
  (equivalent). Transformer is equivalent to RF in all 9 transition × coupling cells.
  Non-equivalent rows are the same as in primary: at III→IV, OT/independent models beat
  their own constant median by 9–11% (constant MAE 0.172 / 0.163 vs 0.130 for the
  true-coupling constant), the artefact described in the 2026-11-12 evening entry; and
  ConstMean vs RF at II→III true (8.3%).
- **RQ4:** measured step re-run against the final RQ1 summary (for provenance; the summary
  file now includes sensitivity rows). Numbers identical to the 2026-11-12 run.
- Nothing pending in the cloud except the Jerber transfer (needs a plan amendment + data).

### 2026-11-12 evening — Claude (cloud workspace) — RQ1/RQ2 primary and RQ4 done
- **RQ1 primary: all 9 jobs, 88 OT solves, all converged** (`results/06_rq1/`, summary at
  `32688d1`; summary-only bug fix `32688d1`). Pooled, clone-weighted, OT balanced ε = 0.05:
  own-clone mass 3.99× / 4.17× / 2.78× chance (I→II / II→III / III→IV), but absolute mass
  0.14% / 0.26% / 0.07%; top-1 recall 0.33% / 0.40% / 0.06%; AUROC 0.61 / 0.64 / 0.59.
  kNN AUROC ≈ 0.50 (only 10 non-zero weights). Unbalanced OT is lower, near chance at
  III→IV. **Cell-level lineage is not recovered.**
- **Fate level (TVD, lower = better), OT ε = 0.05 vs independent:** earlier transitions
  recover clonal fate distributions for progenitors (I-NPvMB 0.10 vs 0.34; II-NPvMB 0.08 vs
  0.36; II-NPaHB 0.10 vs 0.27), but at III→IV OT is no better or worse than independent
  (III-NPpHB 0.59 vs 0.51; III-NPaHB 0.35 vs 0.14; III-NPvMB 0.42 vs 0.55). Matches the
  pre-stated expectation: worst where fates branch.
- **RQ2 primary** (`8d5c2cb`): with true lineage, Ridge/RF/Transformer beat the constant
  median by only 1.3–3.1% of its MAE; every CI lies inside ±δ (equivalent). Transformer is
  equivalent to the random forest in all 9 transition × coupling cells.
  **State carefully:** each coupling's constant is fitted to that coupling's targets. At
  III→IV, OT/independent targets are mostly background cells, so their constant is far off
  for true sisters (MAE 0.169 / 0.162 vs 0.124 for the true-coupling constant) and models
  "beat" it by ~10–13%. Against the true-coupling constant those models are worse
  (MAE 0.14–0.15). Report the pre-specified comparisons, plus this cross-coupling
  contrast from `mae_summary.csv` as a descriptive table.
- **RQ4 complete** (`results/09_rq4/`). Measured mismatch rates are 0.996–0.999 (1 − own-
  clone mass and 1 − top-1 recall alike), i.e. OT pairing is effectively random at cell
  level; per-cell effect RMSE sits at the m = 1 ceiling (0.12 / 0.26 / 0.18 by replicate),
  mean-effect bias ≈ 0. Pearson r stays ≈ 0.88–0.91 even at m = 1 because each treated
  cell carries its own effect (simulation design, ANALYSIS_PLAN §4) — say so.
- **Still running:** RQ1 sensitivity jobs (2 of 9 done) and the uncapped I→II rep-2 check,
  then RQ2 sensitivity (queued automatically). One OT job runs at a time (5.8 GB limit).

### 2026-11-12 — Claude (cloud workspace) — RQ3 done; RQ1 running; RQ2/RQ4 scripts ready
- Plan pulled at `65cd70d` (APPROVED). Scripts committed before running: RQ1 `12fade6`,
  RQ2 `fd29545`, RQ3 `6bde8cb`, RQ4 `fe16ad2`.
- **Implementation fixes, each before any affected result existed:**
  - RQ1 unbalanced OT: POT `sinkhorn_stabilized` did not finish ε = 0.1 on the smallest job
    in 20 min → log-domain unbalanced Sinkhorn, same KL objective and τ (`22dad01`).
  - Its stopping rule: absolute 1e-6 is below float32 resolution for these potentials
    (a solve hung) → relative 1e-6 (`b576f19`). Solves that already met the absolute rule
    also meet this one.
  - RQ2: panel genes matched one by one (unrelated symbols duplicated) (`dfdb11c`);
    chunked float32 cost median, same values (`a46d264`).
- **RQ3 complete** (`fa1564c`), 20 seeds, both cell sets. Macro-AUROC (primary, mean):
  I→II (P,D) 0.59, (Pg,D) 0.61, 30-PC 0.70, fate axis 0.71; II→III 0.70 / 0.71 / 0.83 /
  0.85; III→IV 0.68 / 0.66 / 0.76 / 0.77. Marker states are above chance but well below
  expression representations; the plan expected "near chance" — report as found.
  Sensitivity set agrees. DA descriptive: 73 III-NPvMB clones, 6 produce mDA cells;
  AUROC for mDA-producing vs other NPvMB: Pg 0.50, D 0.32, LDA-1 0.35 (n = 6, descriptive).
  HVG batch fallback never used (0/480).
- **RQ4 grid complete** (`results/09_rq4/grid_*`). Per-cell effect RMSE rises steadily with
  mismatch (m = 0.5: 0.09 / 0.18 / 0.13 by replicate); Pearson r with the true effect
  stays ≥ 0.88 even at m = 1, because each treated cell carries its own effect, so random
  partners add noise but no bias. Mean-effect bias ≈ 0. Measured-rate points wait for RQ1.
- **RQ1 running** in the cloud workspace: 6 of 19 jobs done. The memory limit is 5.8 GB
  and a large job peaks ~3.5 GB, so jobs now run one at a time (four OOM kills when two
  ran together; cached solves were kept). Estimated 10 h more. RQ2 starts per transition
  as its three replicates' OT solutions are ready.
- Note for the paper (RQ2): Pg/D state, Transformer copied from TCBB code; one input
  token, so it is a feed-forward network here.

### 2026-11-11 — Claude (VS Code) — Analysis plan approved
- **Author: "approved"** — read as accepting all Q5 recommendations (items 1–9) and the
  plan. Merged into `ANALYSIS_PLAN.md` before approval, then status set to APPROVED in a
  commit containing only that file: **`65cd70d` is the pre-registration record.** No RQ
  result exists at that commit. Later changes go under *Amendments* with date and reason.
- Merged: RQ2 OT/independent draws restricted to training clones + background; analytic
  chance and Mann–Whitney ties; background defined for the cap; uncapped I→II rep 2
  robustness check; P panel named (reduces to SOX2); D described as floor-plate/vMB
  identity; δ = 0.01 kept and every Δ also reported as % of constant-median MAE; RQ4
  primary mismatch = 1 − own-clone mass.
- Item 10 applied: `paper/main.tex` "data of~\cite{You_2023}".
- **Next: RQ1 in the cloud workspace** (`scripts/06_rq1_coupling.py`), then RQ2–RQ4.
- Answered questions, moved here from *Open questions*:

- **Q4 (2026-11-11) — RQ1 settings.** Superseded by `ANALYSIS_PLAN.md` §0–1 (same
  defaults, plus within-replicate coupling, cost scaling, the 8,000-cell memory cap and the
  fixed OT setting for RQ2/RQ4). Answer by approving the plan.
- **Q5 (2026-11-11, revised with evidence) — Recommendations on `ANALYSIS_PLAN.md`.**
  Evidence: `scripts/05a_plan_checks.py` (committed `6120fea` before running) →
  `results/05a_plan_checks/`. Descriptive only — no coupling fitted, no model trained, no
  early feature related to a later fate. Accept/reject each, then approve.
  1. **RQ2: restrict OT and independent draws to training clones + background (adopt).**
     Cross-stage sister cells are 4.0–22.9% of each later-stage pool (`C_leakage_exposure`;
     highest II→III rep 1). With 30% of clones in val/test, ~1–7% of the pool are
     held-out sisters, and OT mass concentrates on transcriptionally similar cells, so
     draws hit them more often than their share. Fix costs nothing. Also disclose that the
     coupling's PCA sees test cells' features (unsupervised; per CLAUDE.md "fit on all").
  2. **Ties: use analytic chance and Mann–Whitney AUROC (adopt).** Independent coupling →
     AUROC = 0.5 exactly with ties as ½; top-k chance = 1 − C(n−s,k)/C(n,k) per early cell
     (s = own-clone later cells). kNN: AUROC with ties as ½; top-k (k ≤ 10) by the
     hypergeometric expectation within the 10 neighbours. No random tie-breaking.
     Chance is tiny: top-1 0.01–0.18%, top-10 0.09–1.7% (`D_chance_vs_cap`).
  3. **P panel: name it and say what it measures here (adopt).** TCBB P = POU5F1, NANOG,
     SOX2, UTF1, TDGF1. Here POU5F1/NANOG/UTF1/TDGF1 are detected in ≤0.8% of cells at
     every stage; P without SOX2 has mean ≤0.002. **The TCBB P score is a SOX2 score**
     (SOX2 detected in 54–85%), not a silent one, and SOX2 is also in Pg. So the plan's
     wording "shows the panel is degenerate" should become "reduces to SOX2". Keep
     representation (1) for continuity. Keep Pg as proposed: it varies (mean 1.49/1.51/1.68
     at I–III, 0.94 at IV).
  4. **D panel: keep, but describe it honestly (adopt).** SLC6A3 ≤0.5% and DRD2 ≤4.5%
     detected; D is carried by FOXA2 (15–54%), DDC (12–40%), LMX1A (1–20%), and is zero in
     30–37% of later cells (`B_score_spread`). In this system D is a floor-plate/ventral-
     midbrain identity score, not a DA-neuron score. Unchanged for TCBB comparability.
  5. **Equivalence margin: keep δ = 0.01 (adopt as is), and also report Δ as % of the
     constant-median MAE.** On later cells of cross-stage clones, the constant-median MAE
     of (Pg, D) is 0.108 / 0.117 / 0.081 (I→II / II→III / III→IV), so δ = 9.3% / 8.6% /
     12.4% of the baseline error — the same order as TCBB's constant-vs-model gaps (6–7%),
     and the TCBB value (`coupling_protocol.py: EQ_MARGIN = 0.01`). A smaller δ would
     make equivalence unreachable with ~100–600 test clones; a larger one would hide
     TCBB-sized effects. Note III→IV is the most lenient.
  6. **Cap: correction to my earlier point.** I said m/m₀ is invariant to the cap; it is
     not (OT mass shifts when distractors are removed). The cap only binds in I→II (all
     reps) and II→III reps 2–3; elsewhere the pool is already under 8,000. Where it binds,
     chance changes up to 3× (I→II rep 2: top-10 0.28% capped vs 0.09% uncapped).
     Recommend: report chance per replicate from the capped pool actually used, and run
     **one uncapped check on I→II rep 2** (24,391 later cells) as a pre-specified
     robustness analysis.
  7. **Cap: state who is "background" (adopt).** The plan keeps every cross-stage-clone
     cell and caps "unbarcoded" cells, but many later cells are barcoded in clones that
     are not cross-stage (e.g. 16,532 in I→II rep 3). Specify: they are background,
     subsampled with unbarcoded cells, never positives.
  8. **RQ4 mismatch rate: use 1 − own-clone mass as primary (change).** RQ2's OT pairs are
     made by drawing from the OT row, so the probability that a drawn partner is a true
     sister is exactly the own-clone mass m. 1 − m is therefore the mismatch rate the
     pairing actually produces; 1 − top-1 recall (argmax pairing) stays as secondary.
     With chance top-1 ≤0.18%, both will likely be near 1, so RQ4 should also show the
     full m ∈ [0, 1] curve (already planned).
  9. **RQ3 classes confirmed** (`F_fate_classes`): ≥30-clone fates and "6 III→IV clones
     end as mDA" match the plan exactly. I→II GABAs has 28 clones (just under 30) — fine,
     the threshold is fixed in advance.
  10. Minor: `paper/main.tex` reads "data of ~\cite{You_2023}" (double space); use
      "data of~\cite{You_2023}". Not changed by Claude.
  **Answer (2026-11-11, author): approved (all items).**

### 2026-11-11 — Claude (VS Code) — TCBB manuscript template
- Author: no thesis-style per-section files; use a TCBB template. Done:
  **`paper/main.tex`** is the manuscript (IEEEtran v1.8b, `journal,compsoc` — the IEEE
  Computer Society journal format TCBB uses). Title, authors, abstract/keywords, and the
  CLAUDE.md outline (Intro … Limitations & Conclusion, Results subsections per RQ) with
  comments marking what goes where. `paper/sections/` removed.
- Data and Setup rewritten inline as journal prose: three subsections (data; clone
  definition; cell filtering and doublets) instead of six. Same content and numbers.
- Numbers still generated: `paper/generated/data_numbers.tex` (from
  `scripts/05_section3_numbers.py`, run at `eab67aa`), table in `paper/tables/`.
- Build: from `paper/`, `pdflatex main; bibtex main; pdflatex main; pdflatex main`
  (`IEEEtran.bst` is in `paper/bibliography/`; set `BSTINPUTS` or copy it). Compiles
  cleanly, ~1.5 pages so far. Working title is a placeholder — change freely.
- `conference_041818.tex` and old images left untouched.
- Commit: `dbac594`.

### 2026-11-11 — Claude (VS Code) — Section 3 draft
- `paper/sections/03_data.tex`: six subsections (dataset and lineage design; clone
  definition and reproduction; library and replicate assignment; QC; doublets and the
  sensitivity set; analysis cell sets) + Table I. About one IEEE double-column page.
- **No numbers typed in the prose.** `scripts/05_section3_numbers.py` reads `results/01–04`
  and writes `paper/sections/03_data_numbers.tex` (one `\newcommand` per number, incl. the
  QC thresholds and doublet rate parsed from the committed rule strings) and
  `paper/tables/tab_data_counts.tex`. It also recomputes cross-stage counts from the
  per-cell tables and stops if they differ from `results/04_qc/clone_counts.csv` (all
  match). Provenance: `results/05_section3/` (generator run at `f0a8a37`).
- Table I rows: published / all clone-table cells / paper thresholds (not used) /
  primary / doublet sensitivity. I added the paper-thresholds row because the UMI-cap
  argument cites it; drop it if you prefer the four rows you listed.
- Citations in new `paper/sisbar_refs.bib` (You 2023, Wolock 2019; metadata from Crossref).
- Compiled cleanly in an IEEEtran journal test wrapper (no overfull boxes, no undefined
  references). Not wired into a main `.tex`; the section `\input`s its numbers and table
  with paths relative to `paper/`.
- Text that is not from results files: design facts from You et al. (three experiments,
  three replicates, split after three to four divisions) and the statement that the
  published 173,564 reflects further filtering we do not replicate (deliberately vague;
  the cloud log mentions removal of non-reproducible clusters — confirm against the
  paper's methods if you want that stated).
- Two numbers in the cloud log are not in any results file and are left out: Scrublet's
  default-settings result (0 vs 18% in GSM6322990) and the observed/simulated score
  quantiles. The text describes those findings qualitatively. If you want them quoted,
  they need writing to `results/04_qc/` first.

### 2026-11-11 — Claude (cloud workspace) — Final QC cell set
- Q2 → option (b), Q3 → option (i), both author decisions; rules fixed in the script
  docstring and committed (`92a0df7`) before the run.
- Answered questions, moved here from *Open questions*:

- **Q2 (2026-11-11) — Which QC rule defines the analysis cell set?** The paper's stated
  thresholds (≥1,000 genes, ≤30,000 UMIs, ≤15% mt) drop **11,350 of the authors' own
  115,046 clone-table cells (9.9%)**, almost all through the UMI cap (10,973 cells; 40%
  of stage-IV clone cells exceed 30,000 UMIs, as mature neurons are large). III→IV
  falls from 555 to 437 clones. The authors evidently did not apply the UMI cap to
  their lineage table. Options:
  (a) **Paper thresholds as stated** — consistent with the text; loses 21% of III→IV
      clones and biases against mature neurons.
  (b) **Genes ≥1,000 and mt ≤15%, no UMI cap, plus doublet scoring (Scrublet) as the
      doublet control** — keeps all but 379 clone-table cells; 198,073 cells pass.
      Recommended: same rule for barcoded and background cells, and doublets are
      handled by a method designed for them rather than a size cap.
  (c) Use the authors' clone table as-is for barcoded cells and paper thresholds for
      the rest — not recommended (two different rules in one OT problem).
  Neither (a) 178,226 nor (b) 198,073 reproduces the paper's 173,564, because the paper
  also removed clusters not reproducible between replicates. Whichever rule is chosen
  is fixed before any RQ analysis and reported with this table. Blocks RQ1.
  **Answer (2026-11-11, author): option (b).** Implemented in `scripts/04_qc_doublets.py`
  (rule and doublet settings fixed in its docstring before running).

- **Q3 (2026-11-11) — Scrublet is unstable here; how should doublets be handled?**
  Running the pre-specified Scrublet step (per library, expected rate 0.8%/1,000
  barcodes, auto threshold) gave 0–2 doublets in each of the first 6 libraries, then ran
  out of memory. Diagnostics on two stage-I libraries:
  - GSM6322990: 0 doublets with default settings; **2,434 (18%)** with
    `mean_center=False` (needed to fit in memory). Same data, one internal option.
  - GSM6889983: auto threshold 0.578 sits above 99.9% of simulated doublets → 0 flagged.
  - Observed vs simulated doublet scores overlap heavily (observed median 0.10, q99 0.48;
    simulated median 0.27–0.35, q90 0.48): cells within a stage are a continuum, so
    homotypic doublets are not separable.
  So the doublet step would remove 0–18% of cells depending on implementation detail.
  Options:
  (i) **Primary rule = genes ≥1,000 and mt ≤15% only; keep Scrublet scores as a
      recorded column; pre-specify a sensitivity analysis for RQ1–RQ3 that drops each
      library's top expected-rate fraction by score.** Recommended: the main result
      does not hinge on an unstable step, and the sensitivity analysis shows whether
      doublets could matter.
  (ii) Remove the top expected-rate fraction by score in every library (the script's
      pre-specified fallback, applied everywhere). Consistent, but removes ~5–12% of
      cells on a score that barely separates doublets from singlets.
  (iii) Keep the auto threshold as pre-specified and report that it removed almost
      nothing. Honest, but the memory fix changes its output, so it is not reproducible.
  No results from `04_qc_doublets.py` are committed; the script will be revised to the
  chosen option (and to run each library in a separate process to fit in memory).
  Cloud environment changed to make Scrublet run at all: scanpy 1.12.5, scipy 1.16.3
  (scanpy 1.11.4 fails with scipy ≥1.16 sparse arrays).
  **Answer (2026-11-11, author): option (i).** Script revised before running: primary
  cell set = genes ≥1,000 and mt ≤15%; Scrublet scores recorded with fixed settings
  (`mean_center=False`, `random_state=0`, one process per library); `doublet_flag` =
  top expected-rate fraction per library, used only for the pre-specified sensitivity
  re-run of RQ1–RQ3.

- **Primary cell set** (≥1,000 genes, ≤15% mt): **198,073 cells**; 114,667 of 115,046
  clone-table cells kept (99.7%). Cross-stage clones: I→II **756 / 2,696**, II→III
  **1,208 / 5,817**, III→IV **549 / 2,667** (published 768 / 1,213 / 559).
- **Sensitivity set** (also excluding `doublet_flag`): 181,018 cells; clones 689 / 1,172 /
  519. Flagged cells have about twice the median UMIs of kept cells in every library.
- Scrublet's own auto threshold, recorded only, flagged 0 to 2,434 cells per library
  (5,896 total; 0 in 12 of 22 libraries), confirming Q3.
- Doublet score vs UMI count: Spearman 0.18 (weak; the flag is not just a size cap).
- Outputs: `results/04_qc/` (cell_qc.csv.gz, qc_by_library.csv, clone_counts.csv,
  summary.csv, provenance.txt). Per-library score caches in git-ignored
  `data/interim/scrublet/`.
- **Data preparation is complete. Next: RQ1 (coupling recovery).** Before RQ1, the
  author still needs to fix the open decisions in CLAUDE.md that RQ1 depends on (see Q4).

### 2026-11-11 — Claude (cloud workspace) — Expression objects and QC
- **Q1 answered by author:** expression/OT work runs in the Claude cloud workspace
  (Linux). Recorded in CLAUDE.md (`56ddc72`, `325c880`).
- RAW tars and clone table re-downloaded there; SHA-256 identical to the Windows copies.
- Replicates 2–3 have one extra feature, **EGFP** (lineage-vector reporter); removed from
  the gene matrix and kept as `egfp_counts`. All 22 libraries share the same 33,694 genes.
- `scripts/03_build_objects.py` (committed `d45f505` before running): 22 per-library
  AnnData objects with raw counts in git-ignored `data/processed/libs/` (2.1 GB); every
  barcode's QC metrics, replicate (from viral barcode), clone and cluster in
  `results/03_objects/cell_obs.csv.gz`. Join on (library, core barcode): all 115,046
  clone-table cells matched.
- With all clone-table cells, cross-stage counts reproduce 768 / 1,213 / 555 as before.
- **Paper QC thresholds conflict with the authors' clone table** — see Q2. 178,226
  barcodes pass the stated thresholds (paper: 173,564), but 9.9% of clone-table cells
  fail them, nearly all on the 30,000-UMI cap. After that QC: I→II 706 / 2,462,
  II→III 1,099 / 5,010, III→IV 437 / 1,989 (`clone_counts_after_qc.csv`).
- Package versions in the cloud differ from `requirements.txt` pins (pandas 3.0.5,
  numpy 2.5.3, scipy 1.18.1, Python 3.13); exact versions recorded in provenance.

### 2026-11-11 — Claude — Library map and Exp3 gap
- Every clone-table cell is found in the library its suffix/replicate predicts (100%
  overlap; runner-up ≤0.6%, chance collisions). Per-cell map:
  `results/02_library_map/cell_library_map.csv.gz`.
- **Exp3 replicate labels 2↔3 are swapped** between viral barcode and GEO sample name,
  consistently in stages III and IV (viral rep 2 → GSM6889993/GSM6889994 "Rep3"; viral
  rep 3 → GSM6889987/GSM6889988 "Rep2"). Clone pairing unaffected; replicate ID must
  come from the viral barcode.
- Exp3 gap not closed: all 5,329 stage-IV cells are in their matrix; still 555 / 2,692.
  Accepted per the pre-set rule; to be reported vs 559 / 2,816.
- 20 stage-IV cell barcodes occur in more than one stage-IV library (resolved by replicate).
- Matrices hold 215,960 barcodes (unfiltered) vs the paper's 173,564 post-QC cells →
  QC thresholds will need choosing.
- Provenance now correct: script committed (`cfcc17d`) before the run.
- Commits: `cfcc17d` (script), `ef11252` (results + CLAUDE.md).

### 2026-11-10 — Claude — Data download and clone-count check
- Both GEO series cover all four stages: GSE207921 = replicate 1, GSE221592 =
  replicates 2–3. Barcode-to-cell table deposited as
  `GSE221592_Virusbar_processed_data_combined.xlsx`.
- Clone = exact viral-barcode set (pooled over replicates) reproduces I→II 768 / 2,762
  and II→III 1,213 / 5,856 exactly; III→IV 555 / 2,692 vs 559 / 2,816.
- RAW tars downloaded and verified (sizes match GEO; SHA-256 in `data/SHA256SUMS`).
- Project venv `.venv` (pandas, numpy, scipy, openpyxl); scanpy blocked by Application
  Control. `data/` and `.venv/` are git-ignored.
- Commit: `c738aea`.
