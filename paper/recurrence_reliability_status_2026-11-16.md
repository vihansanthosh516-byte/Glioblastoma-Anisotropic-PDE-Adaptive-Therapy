# Recurrence reliability: execution status (dated handoff)

Last update: 2026-11-16, Claude. Executing `paper/CLAUDE_RESEARCH_PAPER_HANDOFF_2026-11-16.md`.

| Phase | State | Record |
|---|---|---|
| A novelty | Done: passes narrowly (methods/benchmark, not biology) | `paper/recurrence_reliability_prior_art.md` |
| B eligibility | Done: 23 of 27 pairs | `paper/recurrence_reliability_eligibility.md`, `output/recurrence_reliability/metadata_manifest.json` |
| C QC | Plan frozen (commit cce5406); matrices downloading/downloaded to `data/recurrence_reliability/raw/gsm/` (git-ignored) | `paper/recurrence_reliability_qc_plan.md` |
| D protocol | Drafted, not locked | `paper/recurrence_reliability_protocol.md` |
| E code | Written; 19 synthetic tests pass in a scratch copy | `src/recurrence_reliability/`, `tests/test_recurrence_reliability.py` |
| F run | Not run. No program score has been computed on real data. | — |

## Exact next actions

1. Confirm all 98 matrix/feature files downloaded (`data/recurrence_reliability/matrix_urls.txt`; rerun `tr -d '\r' < data/recurrence_reliability/matrix_urls.txt | xargs -n1 -P8 curl -sS -m 900 --retry 3 -O` inside `raw/gsm/` if any is short).
2. `py -3 src/recurrence_reliability/qc.py` → `output/recurrence_reliability/qc_report.json`. Read it back.
3. Write `configs/recurrence_reliability.json` with `locked: true`, the QC report SHA-256 and the final cohort; update the protocol n; commit and push (the freeze).
4. Run once: `py -3 src/recurrence_reliability/run.py --config configs/recurrence_reliability.json`.
