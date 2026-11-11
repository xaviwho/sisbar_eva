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

- [ ] 2026-11-11 — **Author: approve `ANALYSIS_PLAN.md`.** It fixes every setting for
  RQ1–RQ4 before any result exists (answers Q4 and the CLAUDE.md open decisions on the
  P/D panel, the equivalence margin and pre-registration). Reply "approved" (or edit the
  file). Claude then changes its status line to APPROVED and commits: that commit is the
  pre-registration record. **Instructions below wait for this.**
- [ ] 2026-11-11 — **(Cloud) RQ1 coupling recovery** → `scripts/06_rq1_coupling.py`,
  `results/06_rq1/`. ANALYSIS_PLAN §0–1. Primary set, then sensitivity set.
- [ ] 2026-11-11 — **(Cloud) RQ2 forecasting on true lineage** → `scripts/07_rq2_forecast.py`,
  `results/07_rq2/`. ANALYSIS_PLAN §2. Uses the RQ1 OT coupling (balanced, ε = 0.05).
- [ ] 2026-11-11 — **(Cloud) RQ3 state adequacy** → `scripts/08_rq3_fate.py`,
  `results/08_rq3/`. ANALYSIS_PLAN §3. Can run in parallel with RQ2.
- [ ] 2026-11-11 — **(Cloud) RQ4 perturbation consequence** → `scripts/09_rq4_mismatch.py`,
  `results/09_rq4/`. ANALYSIS_PLAN §4. Needs RQ1's measured mismatch rates.
- [ ] 2026-11-11 — **(VS Code) Results sections, as each RQ lands:** write the matching
  Results subsection in `paper/main.tex`, with figures and numbers generated from
  `results/06–09` by scripts (no typed numbers). Report nulls as fully as positives.
- [ ] Later — Jerber D30 transfer of the RQ3 fate-bias axis (needs Jerber data; separate
  plan amendment).

## Open questions (Claude → author)

- **Q4 (2026-11-11) — RQ1 settings.** Superseded by `ANALYSIS_PLAN.md` §0–1 (same
  defaults, plus within-replicate coupling, cost scaling, the 8,000-cell memory cap and the
  fixed OT setting for RQ2/RQ4). Answer by approving the plan.

## Done

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

---

## Log (newest first)

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
