# Prompt for Codex, round 2: audit the OLD tracks claim by claim (nothing sampled)

Paste everything below the line into Codex. The work list is fixed in `paper/codex_audit2/chunks.json` (74 chunks). Progress is kept in `paper/codex_audit2/progress.md`, so you can stop and resume at any time (the user switches tools when a limit hits).

---

You are an independent auditor of the OLD work in the repo `C:\Users\vihan\20206 science fair`: Track A (scripts 01-41), Track B (42-48), Track C and C-Extended (49-65, `src/neural_pde`, `src/rl`, `src/sensing`, `src/hil`), the forecasting and PredictGBM work (scripts 66-139), the old ledger entries, and all Markdown files. A first audit (`paper/codex_audit_report.md`) only sampled this work. This round must cover EVERY item in the work list. Nothing is exempt because it is old, archived, or already reported.

## Rules
1. Use the cheapest, smallest model for sub-agents. Run up to 4 sub-agents in parallel, one chunk each. Use a larger model only for `SUMMARY.md`.
2. READ-ONLY on the repo. You may write only inside `paper/codex_audit2/`. Do not edit, commit or push anything else.
3. Never run heavy jobs. Never run a script in place if it writes to `output/`. If you run a script, copy the repo to a scratch folder first. Light checks only: `python -m py_compile`, reading files, counting rows, opening JSON.
4. Every verdict needs evidence: file path, line or JSON key, the value you saw. Write "could not verify" instead of guessing. Do not praise.
5. No chunk may be skipped, merged or marked done without its file. If a chunk is too big, split it into `chunk_NNNa`, `chunk_NNNb` and say so in `progress.md`.

## Procedure
1. Run `python paper/codex_audit2/make_chunks.py` once (it keeps existing progress). Read `chunks.json` and `progress.md`.
2. For every chunk with status TODO: do the chunk task below, write `paper/codex_audit2/chunk_NNN.md`, then change its row in `progress.md` to DONE (or PARTIAL with the reason). Update `progress.md` after EVERY chunk.
3. When all chunks are DONE or PARTIAL, write `paper/codex_audit2/SUMMARY.md` (see end).

## Chunk task by kind
### kind = ledger (entries of `vault/PAPER_FINDINGS_LEDGER.md`)
Split each listed entry into its individual claims (every number, comparison, p-value, dataset size, "supported/not supported", date, script reference). For EACH claim make one table row:
`entry | claim (short quote) | cited file/key | value found | status | evidence | one-line correction`.
Status is exactly one of: CORRECT, WRONG_NUMBER, NO_SOURCE (no file cited or file missing), STALE (later entry, amendment or output contradicts or supersedes it; name it), OVERCLAIM (simulation or assumed parameter or post hoc or spent set presented as validated; "first", "novel"), LIMIT_MISSING (the weakest point is not stated), UNVERIFIABLE (needs data or compute you do not have).
Also check: is the entry marked superseded when a later entry replaces it? Does the script named exist? Does the script write the file named?

### kind = markdown (Markdown files)
For each file: list every factual or numeric claim and check it as above against `output/`, the ledger and the code. Also record: is the file current, historical or wrong for the current goal (`START_HERE.md` says which); does it carry the HISTORICAL banner if it should; does it call anything "first", "novel", "revolutionary", "validated", "clinically viable" without a source; does it describe simulated or assumed results as real. Same table and same status labels.

### kind = script (Python files in `src/`)
For each script: (a) `python -m py_compile` result; (b) imports that do not exist in `Requirements.txt` or the repo; (c) output paths: relative to the working directory (violates `CLAUDE.md`) or built from `__file__`; (d) every output file the script writes: does it exist in `output/` (or `archive_old_gbm/output/` if moved); (e) does the header or docstring claim something the code does not do; (f) hard-coded numbers that look like results typed in by hand; (g) random seeds present or missing; (h) any use of test or locked data before a freeze (`analysis_plan_v1.md` amendments); (i) any assumed constant (kill rate, growth rate, threshold) that the ledger presents as measured. One table row per script and per problem found: `script | check | finding | evidence | severity (high, medium, low)`.

### kind = tests
For each test file: what it really asserts; any test that cannot fail (asserts True, only checks shapes, mocks the thing under test); any hard-coded expected value copied from a result; which scripts it covers. State which claims in the repo cite "N tests passed" and whether that number is right (run `pytest tests/ -q` ONLY in a scratch copy).

## Final files
- `paper/codex_audit2/SUMMARY.md`: (1) counts by status per kind and per track; (2) the 30 worst findings, each with file and evidence; (3) a list of every ledger entry and Markdown file that needs a correction, with the one-line correction; (4) scripts that should be excluded from the runnable inventory; (5) "Could not verify" list; (6) questions for the author. Keep it plain.
- `paper/codex_audit2/progress.md`: all chunks DONE or PARTIAL with reasons.

Do not fix anything. Report only.
