# Recurrence reliability: execution status (dated handoff)

Last update: 2026-11-16, Claude. Executing `paper/CLAUDE_RESEARCH_PAPER_HANDOFF_2026-11-16.md`.

| Phase | State | Record |
|---|---|---|
| A novelty | Done: passes narrowly (methods/benchmark, not biology) | `paper/recurrence_reliability_prior_art.md` |
| B eligibility | Done: 23 of 27 pairs | `paper/recurrence_reliability_eligibility.md`, `output/recurrence_reliability/metadata_manifest.json` |
| C QC | Done: 21 patients (pairs 19, 20 fail malignant support after fixed filters); feature files are symbol-only (recorded deviation) | `paper/recurrence_reliability_qc_plan.md`, `output/recurrence_reliability/qc_report.json` |
| D protocol | LOCKED in commit 642f9ac with `configs/recurrence_reliability.json` | `paper/recurrence_reliability_protocol.md` |
| E code | Written; 19 synthetic tests pass in a scratch copy | `src/recurrence_reliability/`, `tests/test_recurrence_reliability.py` |
| F run | Started once 2026-11-16 (`run_20261116_locked`); stdout `output/recurrence_reliability/run_20261116_stdout.txt`. If interrupted, do not rerun silently: record the interruption first, delete the partial run dir with a dated note, then rerun. | `output/recurrence_reliability/run_20261116_locked/` |

## Exact next actions

1. Wait for the run to finish; read `result.json` and `patient_changes.tsv` back (counts, finite values, denominators, hashes).
2. Write `paper/recurrence_reliability_result.md`, update START_HERE, draft manuscript.
