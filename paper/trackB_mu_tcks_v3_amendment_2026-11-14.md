# Track B v3 (median-trained latent-size regression) — amendment, locked before v3 scoring (2026-11-14)

## Unchanged (same question and focus)
Same question, data, pairs, folds, target y = log(core_vox_out + 1), primary score (patient MAE of y), 6 baselines, bootstrap seed 20261114, secondary outputs. See `paper/trackB_mu_tcks_protocol_2026-11-14.md`.

## Multiple-try correction (new, applies from v3 on)
Each new method on this data is one more try. Try k passes only if the upper bound of the two-sided (1 - 0.05/k) bootstrap CI vs the best baseline is < 0. v3 is try 3: 98.33% CI. v1 and v2 results are not changed.

## v3 method (only change)
- The score is an absolute error, so the best point forecast is a median. v1/v2 used mean-type fits.
- v3 fits a median (L1) linear regression on training folds: y_out ~ x_hat + y_in + dt + dt*phase dummies + has_prior.
- x_hat = v2 latent size (same formula, training-fold parameters).
- sklearn QuantileRegressor(quantile=0.5, alpha=0, solver="highs"). No tuning.

## Honesty note
Post hoc (designed after v1, v2 results) on a reused cohort. Exploratory only.
