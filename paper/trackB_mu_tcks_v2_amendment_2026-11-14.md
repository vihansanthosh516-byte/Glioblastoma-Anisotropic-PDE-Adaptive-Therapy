# Track B TCKS v2 (noise-aware latent size) — amendment, locked before v2 scoring (2026-11-14)

## What does NOT change (same question and focus)
- Question: does the new method predict the next tumor-core volume better than all standard volume-forecast methods?
- Same data, pairs, 5 patient folds, target y = log(core_vox_out + 1).
- Same primary score (patient-averaged MAE of y), same bootstrap (5000, seed 20261114), same success rule (CI upper bound < 0 vs best baseline).
- Same 6 scored baselines, unchanged code. Same secondary outputs.

## What changes (method only)
v1 TCKS failed (`output/track_b_mu_tcks_20261114/result.json`). v2 treats each observed size as noisy:
- Prior estimate of the true size at day_in from the previous scan: x_prior = y_prev + r_phase dt_prior, with variance P = s2 + tau2 dt_prior^2.
- Update with the observed y_in (noise variance s2): K = P / (P + s2); x_hat = x_prior + K (y_in - x_prior).
- No prior scan: x_hat = y_in.
- Forecast: y_hat = x_hat + r_phase dt.
- r_phase, tau2, s2: estimated on training folds exactly as in v1.

## Honesty note
v2 was designed after seeing v1 results (v1 showed past rate is noise-dominated). So v2 is post hoc on a reused cohort: exploratory only. The v1 failure stays reported.
