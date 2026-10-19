# Prompt for Codex: audit the ENTIRE repo, old tracks included (read-only)

Paste everything below the line into Codex.

---

You are an independent auditor. Audit EVERYTHING: the new probabilistic-RANO work AND the older Track A (MSOS, scripts 01-41), Track B (PDE cohort, 42-48), Track C and C-Extended (49-65 and later), the H-1 to H-4 forecasting work, the PredictGBM scripts, the viewer, the vault notes, the whole ledger, every README, every MD file, and every JSON in `output/` that a document cites. Nothing is exempt because it is old. Work in the repo at `C:\Users\vihan\20206 science fair`. Do NOT trust any summary. Check files on disk. Use the cheapest, smallest model available to you for sub-agents (the smallest "mini" or "nano" model in your current family). Use sub-agents to split the work by area: (1) new RANO work, (2) Track A, (3) Track B, (4) Track C and C-Extended, (5) forecasting and PredictGBM scripts 80-139, (6) ledger entries before 2br, (7) all MD files, (8) `output/` JSON and CSV files, (9) tests and environment files (`Requirements.txt`, `pyproject.toml`, `Makefile`, `run_all.*`, `Dockerfile`). Each sub-agent writes its own part to `paper/codex_audit/<area>.md`. You may use a larger model only to write the final report.

## Read first
`START_HERE.md`, then `AGENTS.md`, `CLAUDE.md`, `analysis_plan_v1.md` (Amendments 14 to 21 only), `vault/PAPER_FINDINGS_LEDGER.md` (ALL entries), `paper/chatgpt_review_response.md` (Rounds 3 and 4).

## Project in one line
"Probabilistic RANO": how much does segmentation noise change RANO 2.0 volumetric progression calls and simulated trial PFS and power in glioblastoma. Old GBM work is in `archive_old_gbm/` and is background only.

## Hard rules for you
1. READ-ONLY on everything except new files in `paper/codex_audit/` and the final `paper/codex_audit_report.md`. Do not edit code, docs, ledger, outputs or amendments. Do not commit. Do not push.
2. Do NOT re-run scripts 140 to 145 in place (they overwrite `output/`). If you want to re-run one, first copy the repo to a scratch folder and run there. Report whether the re-run matches the saved output to 3 decimals.
3. Do not fix anything. Report only.
4. Every finding must name the file and line, quote the number or text, and say how you checked.
5. Do not give praise. Do not guess. Write "could not verify" when you could not.

## Tasks
A. Number check. For every number in ledger entries 2br, 2bs, 2bt, 2bu, 2bv, 2bw, 2bx, open the cited file in `output/` and compare. List every mismatch, rounding error, or number with no file.
B. Order check. Use `git log` to confirm that for each of Amendments 15, 16, 17, 18, 19, 20, 21 the declaration commit is earlier than the code commit, and the code commit is earlier than the result commit. List any violation.
C. Declaration vs code. For scripts `140` to `145`, `rano_calc.py` and `tests/test_rano_calc.py`: does the code do exactly what its amendment says (formulas, grid, seeds, thresholds, estimands, pass rules)? List every deviation, even small ones.
D. Prediction check. For each declared prediction (H-S1 to H-S4, PA1, PA2, PB1, PD1, PX1 to PX4, PU1 to PU4), recompute supported or not from the saved output files and compare with what the ledger says.
E. Tests. Run `python -m pytest tests/test_rano_calc.py -q` and report the result. Run the whole `pytest tests/` in a scratch copy and report failures (some old tests may fail because data is missing; say which).
F. Claims in docs. Read `README.md`, `START_HERE.md`, `paper/overview_for_review_rano.md`, `paper/phd_level_plan_probabilistic_rano.md`, `paper/problem_first_rano_noise.md`. List every sentence that is false, stale, overclaims ("first", "novel", "validated", "proves"), or hides a limit. Check that the 0.5 cm3 floor is never called a RANO 2.0 rule.
G. Science check. Look for errors in the method: Is the Kaplan-Meier median, log-rank, and Peto log HR implemented correctly in `src/141_trial_sim.py`? Is the nadir rule right? Is the AR(1) noise right in `src/142_rano_sensitivity.py`? Is the patient-cluster bootstrap right in `src/140_rano_tool_dependence.py`? Is the conversion table in `src/145_conversion_audit.py` correct arithmetic? Do `src/rano_calc.py` formulas match the normal-error model?
H. Leakage and reuse. Confirm no script reads outcome data before its declaration. Confirm the H-2 and H-3 sets are not used as locked tests.
I. Repo hygiene. List files in the repo root and `paper/` that are stale, duplicated, or in conflict with `START_HERE.md`. List large files that should not be tracked.
K. OLD WORK, wrong and outdated. For EVERY ledger entry before 2br and every claim in `README_STATUS.md`, `vault/*.md`, `paper/*.md` and `archive_old_gbm/**`: (a) does the cited `output/` file exist and hold that number; (b) is the claim later contradicted by a newer ledger entry, amendment or output (stale); (c) does it overclaim (real data vs simulation, assumed kill rates, "validated", "first", "novel", locked test vs spent set); (d) is it labelled with its weakest point. Make a table: claim | where | status (correct / wrong number / stale / overclaim / no source) | evidence | suggested correction (one line).
L. OLD CODE. For each script in `src/` numbered 01-139 and each script in `src/neural_pde`, `src/rl`, `src/sensing`, `src/hil`: does it import and parse (`python -m py_compile`), does it use relative output paths instead of `__file__` paths (CLAUDE.md rule), does its declared output file exist in `output/`, does it write something different from what the ledger says it wrote. Report broken or unreproducible scripts. Do not run heavy ones; say "not run".
M. Cross-document consistency. List every number or fact that appears in two or more documents with different values (script counts, patient counts, Dice or gain values, p-values, dataset sizes, dates).
N. Dead files. List files that nothing references and that look abandoned (old logs, duplicate outputs, `.err` files, copies), and tracked files over 10 MB.
J. Anything else you find that would make a judge or reviewer lose trust.

## Output
Write `paper/codex_audit_report.md` with these sections: Summary (max 15 lines, worst problems first, with counts of wrong / stale / overclaim items per area), then one section per task A to N, each a table: finding | file:line | evidence | severity (high, medium, low). End with a list "Cannot verify" and "Questions for the author". Keep it plain and short.
