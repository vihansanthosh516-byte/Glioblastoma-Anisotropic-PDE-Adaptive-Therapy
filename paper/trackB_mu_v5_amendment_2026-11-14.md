# Track B v5 (input-scan tumor-composition latent model) — amendment, locked before v5 scoring (2026-11-14)

## Unchanged
Same question, data, pairs, folds, target y = log(core_vox_out + 1), primary score, 6 baselines, bootstrap seed (`paper/trackB_mu_tcks_protocol_2026-11-14.md`). Multiple-try rule: try 5 passes only if the upper bound of the 99% (1 - 0.05/5) CI vs best baseline is < 0.

## v5 method (only change)
v1-v4 used only sizes and dates. v5 adds what the INPUT scan's tumor mask shows (labels from the dataset mask at tp_in only; the target scan is never read for features):
- log edema volume (label 2), necrosis fraction of core (label 1 / core), enhancing fraction of core (label 3 / core), log whole-tumor volume, core compactness (core voxels / bounding-box voxels).
- Plus the v3 features (v2 latent x_hat, y_in, dt, dt-by-phase, has_prior).
- Model: median (L1) linear regression, sklearn QuantileRegressor(quantile=0.5, alpha=0, solver="highs"), features standardized on training folds. No tuning.

## Honesty note
Post hoc, try 5, reused cohort. Exploratory only. Composition-based growth prediction has prior art in radiomics; the result file will state this.
