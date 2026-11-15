# Track B hybrid-method search — locked protocol (2026-11-14)

Written before any hybrid candidate is built or scored.

## Same question, same score (unchanged)
Predict the next tumor-core volume. Target y = log(core voxels + 1), core = labels {1, 3}. Primary score: patient-averaged MAE of y. Data: MU primary pairs (`data/manifests/forecast_pairs_mu.csv`), 5 patient folds (`data/manifests/split_mu.csv`). Baselines: the 6 from `paper/trackB_mu_tcks_protocol_2026-11-14.md`, plus v8 as the best earlier method.

## What counts as a hybrid
Each candidate must join at least two different model families into ONE new method (for example: mechanistic growth model + statistical learner; state-space filter + classifier; gated mixture of experts). A plain known model alone is not allowed.

## Rules for every candidate
- Uses only `src/track_b_hybrid_harness.py` for data, folds and scoring. The harness is not edited by candidates.
- Implements `fit_predict(tr, te, comp_tr, comp_te) -> np.ndarray`. Fitting uses training folds only. Any setting choice uses inner CV inside the training folds only.
- Never reads UCSD-PTGBM or any other cohort.
- One fixed method per candidate; result written once to `output/track_b_hybrid_<name>/result.json`.

## Selection and correction
- 6 candidates. Dev pass rule for each: upper bound of the (1 - 0.05/14) patient-bootstrap CI vs the best of the 6 baselines < 0 (14 = 8 earlier tries + 6 new).
- Finalist = lowest dev MAE among candidates. Ties: simpler method.
- Finalist confirmation: one frozen run on UCSD-PTGBM with `src/track_b_v8_ucsd_external.py` rules (fit on all MU; comparator pop_rate; 95% CI).

## Honesty (stated in advance)
- UCSD was already used once to test v8. A second use for the finalist is a second exposure; it will be labelled so.
- A positive result cannot be promised. If no candidate passes, that is reported.
- Novelty of the finalist will be checked by a prior-art search before any "new" claim.
