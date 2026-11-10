# Track A — expression and regional biology

## Current role

Track A is a completed exploratory single-cell analysis. Its frozen A5 candidate did not meet its five-of-six donor success rule. The integrated phase uses Track A as conceptual motivation for a molecular-like channel. It uses a declared synthetic signal-quality grid, not observed A5 metrics, as q_signal settings. Neither the expression labels nor the synthetic channel are drug-resistance measurements.

## Primary records

- [Frozen Yu et al. A5 protocol](../../paper/trackA_yu2020/PROTOCOL.md)
- [A5 execution audit](../../paper/trackA_execution_audit_20261031.md)
- [A5 null-result record](../../paper/trackA_null_result.md)
- [Novelty decision](../../paper/trackA_novelty_decision_20261101.md)
- [Alternate-cohort screen](../../paper/trackA_alternate_dataset_plan_20261101.md)
- [Integrated plan (active, revised 2026-11-08)](../../paper/trackABC_integrated_research_plan_revised_2026-11-08.md)

## Main implementation and checks

- `src/189_trackA_yu2020_lodo.py` — nested donor-held-out A5 analysis.
- `src/151_trackA_donor_table.py`, `src/157_trackA_composition_audit.py`, `src/trackA_identity.py` — donor and composition checks.
- `tests/test_trackA_identity.py`, `tests/test_trackA_donor_table.py`, `tests/test_trackA_composition_audit.py`, `tests/test_trackA_runner_helpers_synthetic.py` — focused checks.
- `paper/trackA_yu2020/` — frozen protocol and donor list.
- `output/trackA_yu2020/` and `output/trackA_composition_audit/` — saved records; read the protocol and audit before interpreting them.

## Interpretation boundary

The A5 labels are author-supplied MRI-marked sampling-point labels and the analyzed expression is mixed-cell. A5 failed its frozen donor-level success rule. It does not support malignant-cell-specific, resistance, causal, or independent external-confirmation claims. Deferred analyses and alternate-cohort leads remain deferred unless a dated amendment changes scope. The separate 119-of-121 gene score is an exploratory published-response proxy (paired differences split 3 and 3), not a drug-resistance test. A further metadata-only cohort search on 2026-11-09 found no dataset that passes every gate; see the [evidence table](../../paper/gap_status_2026-11-09.md).
