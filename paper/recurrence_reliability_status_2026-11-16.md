# Recurrence reliability: execution status (dated handoff)

Last update: 2026-11-16, Claude. Executing `paper/CLAUDE_RESEARCH_PAPER_HANDOFF_2026-11-16.md`.

| Phase | State | Record |
|---|---|---|
| A novelty | Done: passes narrowly (methods/benchmark, not biology) | `paper/recurrence_reliability_prior_art.md` |
| B eligibility | Done: 23 of 27 pairs | `paper/recurrence_reliability_eligibility.md`, `output/recurrence_reliability/metadata_manifest.json` |
| C QC | Done: 21 patients (pairs 19, 20 fail malignant support after fixed filters); feature files are symbol-only (recorded deviation) | `paper/recurrence_reliability_qc_plan.md`, `output/recurrence_reliability/qc_report.json` |
| D protocol | LOCKED in commit 642f9ac with `configs/recurrence_reliability.json` | `paper/recurrence_reliability_protocol.md` |
| E code | Written; 19 synthetic tests pass in a scratch copy (rerun after each code change) | `src/recurrence_reliability/`, `tests/test_recurrence_reliability.py` |
| F run | Done. First attempt stopped before scoring (spec conflict) → amendment 1 → `run_20261116_locked_a1` completed in 122 s. Discordance 3/21 (95% CI 3.0–36.3%), inconclusive. | `paper/recurrence_reliability_result.md`, `output/recurrence_reliability/run_20261116_locked_a1/` |
| Manuscript | Draft written | `paper/recurrence_reliability_manuscript_draft.md` |
| G external | Metadata audit started; CARE remains conditional, no expression/scoring opened | `paper/recurrence_reliability_care_metadata_audit_2026-11-16.md` |

## Exact next actions

1. Do not rerun `run_20261116_locked_a1`; it is the reported result.
2. Recommended next: Phase G CARE audit (GSE274546/GSE274548): paired IDs, CNV annotation linkage (GBM-CARE-WT repository), overlap with Wang, terms, prior exposure. Metadata only; reserve before any scoring; apply the frozen pipeline once. This is the only route to a narrower interval.
3. Alternative: prespecified multi-scorer extension (AUCell, JASMINE, ANS, z-score) on the same 21 patients, labelled exploratory.
4. Figures for the manuscript from `patient_changes.tsv` (patients with SE bars), not cells.

## Post-run review (Codex, 2026-11-16)

Saved counts, exact interval arithmetic, config hash, and all recorded source hashes were verified without rerunning the study. See `paper/recurrence_reliability_postrun_review_2026-11-16.md`. CARE is the preferred independent-cohort lead, not the only possible source of additional patients. A 59-patient cohort cannot confirm a rate below 5% using the locked two-sided interval even with zero events (upper limit about 6.06%). The 59/72 sample-size calculations assume zero events and do not promise reliability at the observed rate.

Before any CARE outcome analysis, a separate prospective protocol must resolve patient overlap, annotation/barcode linkage, cohort/QC, shared-reference inference assumptions, and whether the Wang fitted reference is transported or the scoring algorithm is refitted on CARE. The original Wang run, protocol, and result remain unchanged.
