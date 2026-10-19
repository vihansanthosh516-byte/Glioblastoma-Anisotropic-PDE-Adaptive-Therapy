# Prompt for Codex: audit the whole repo (read-only)

Paste everything below the line into Codex.

---

You are an independent auditor. Work in the repo at `C:\Users\vihan\20206 science fair`. Do NOT trust any summary. Check files on disk. Use the cheapest, smallest model available to you for sub-agents (the smallest "mini" or "nano" model in your current family). Use sub-agents to split the work by area. You may use a larger model only to write the final report.

## Read first
`START_HERE.md`, then `AGENTS.md`, `CLAUDE.md`, `analysis_plan_v1.md` (Amendments 14 to 21 only), `vault/PAPER_FINDINGS_LEDGER.md` (entries 2br to 2bx only), `paper/chatgpt_review_response.md` (Rounds 3 and 4).

## Project in one line
"Probabilistic RANO": how much does segmentation noise change RANO 2.0 volumetric progression calls and simulated trial PFS and power in glioblastoma. Old GBM work is in `archive_old_gbm/` and is background only.

## Hard rules for you
1. READ-ONLY on everything except one new file: `paper/codex_audit_report.md`. Do not edit code, docs, ledger, outputs or amendments. Do not commit. Do not push.
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
J. Anything else you find that would make a judge or reviewer lose trust.

## Output
Write `paper/codex_audit_report.md` with these sections: Summary (max 10 lines, worst problems first), then one section per task A to J, each a table: finding | file:line | evidence | severity (high, medium, low). End with a list "Cannot verify" and "Questions for the author". Keep it plain and short.
