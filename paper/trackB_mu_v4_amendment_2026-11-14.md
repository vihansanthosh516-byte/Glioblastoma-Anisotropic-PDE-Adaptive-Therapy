# Track B v4 (median-loss boosted latent model) — amendment, locked before v4 scoring (2026-11-14)

## Unchanged
Same question, data, pairs, folds, target, primary score, 6 baselines, bootstrap seed, as in `paper/trackB_mu_tcks_protocol_2026-11-14.md`. Multiple-try rule from v3: try 4 passes only if the upper bound of the 98.75% (1 - 0.05/4) CI vs best baseline is < 0.

## v4 method (only change)
- Same v3 features (v2 latent size x_hat, y_in, dt, dt-by-phase, has_prior) plus phase dummies.
- Nonlinear median model: sklearn HistGradientBoostingRegressor(loss="absolute_error", random_state=0), default settings. No tuning.

## Honesty note
After v3, I looked at the 32 pairs with log error > 3 (they hold 36% of the total error; many are tumor collapse to 0 or sudden jumps, likely surgery or segmentation drop-out). That look used outcomes, so v4 is post hoc. Exploratory only.
