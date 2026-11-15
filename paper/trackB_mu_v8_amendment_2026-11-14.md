# Track B v8 (two-part collapse model) — amendment, locked before v8 scoring (2026-11-14)

## Unchanged
Same question, data, pairs, folds, target, primary score, 6 baselines, bootstrap seed (`paper/trackB_mu_tcks_protocol_2026-11-14.md`). Multiple-try rule: try 8 passes only if the upper bound of the 99.375% (1 - 0.05/8) CI vs best baseline is < 0.

## v8 method (only change)
- Collapse label (training only): core_vox_out < 0.05 * core_vox_in.
- Part 1: logistic regression (sklearn LogisticRegression, default settings, max_iter 1000) on the v5 feature set, standardized on training folds, gives p = P(collapse).
- Part 2: v5 median regression refit on NON-collapse training pairs only.
- Point forecast (median of the two-part mixture, simple rule): if p > 0.5, predict the training median y_out of collapse pairs; else the part-2 forecast.
- "Multiple surgeries" is NOT used (it is known only after the fact).
No tuning.

## Honesty note
Designed after seeing that collapse pairs hold a large share of error (post hoc, try 8, reused cohort). Exploratory only. Two-part (hurdle) models are standard; the novelty is only the use here.
