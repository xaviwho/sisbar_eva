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

_None yet. Add items as `- [ ] YYYY-MM-DD — instruction`._

## Open questions (Claude → author)

- **Q1 (2026-11-11) — Where should scanpy / expression work run?** Scanpy is blocked on
  the Windows machine by Application Control. Options: (a) WSL2 on this machine, if
  policy allows; (b) a Linux lab server (needs access details). Blocks next data task 3
  (per-clone expression objects).
  **Answer (2026-11-11, author):** run it in the Claude cloud workspace (Linux; public
  SISBAR data only). Being acted on — see next log entry.

## Done

- [x] 2026-11-10 — Download both GEO series; confirm barcode-to-cell tables and sample
  mapping. → `c738aea`
- [x] 2026-11-10 — Reproduce published clone counts. → `c738aea`
- [x] 2026-11-11 — Map `_N` suffixes / stage-IV replicate to GSM libraries. → `ef11252`
- [x] 2026-11-11 — Resolve or document Exp3 discrepancy (one pass). → `ef11252`

---

## Log (newest first)

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
