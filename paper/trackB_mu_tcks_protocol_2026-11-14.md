# Track B: Treatment-Clock Kinetic Shrinkage (TCKS) on MU-Glioma-Post — locked protocol (2026-11-14)

Status: written and committed BEFORE any TCKS or baseline score was computed.
Label: EXPLORATORY. MU-Glioma-Post was used before in this project (e.g. `output/mu_fitted_response/results.json`, PDE vs persistence Dice). It is not an untouched test set.

## Question
Does a forecast that shrinks each patient's own growth rate toward a treatment-phase population rate predict the next tumor-core volume better than all standard volume-forecast methods?

## Data (fixed)
- Pairs: `data/manifests/forecast_pairs_mu.csv`, rows with `in_primary` and `input_core_nonempty` true.
- Folds: `data/manifests/split_mu.csv` (5 patient folds, seed 20261004, made 2026-10). Not changed.
- Treatment dates: `data/tcia/MU-Glioma-Post_ClinicalData-July2025.xlsx`, sheet "MU Glioma Post".
- Volume = core voxel count from the manifest. No image is read.

## Target and score
- y = log(core_vox_out + 1).
- Primary score: mean absolute error (MAE) of y, averaged per patient first, then over patients.
- Out-of-fold: every model is fit on 4 folds and scored on the 5th.
- Primary test: TCKS patient MAE minus best baseline patient MAE (best = lowest pooled MAE). 95% patient bootstrap CI, 5000 draws, seed 20261114.
- Success rule: CI upper bound < 0. Otherwise "not better".
- Secondary (reported, not used for success): MAE in pairs with a prior scan; 25%-growth call AUC.

## Baselines (all run on the same pairs and score)
1. Persistence: y_hat = log v_in.
2. Two-point exponential extrapolation (prior scan rate; persistence if no prior scan).
3. Population constant rate (training median log-rate × dt).
4. Empirical-Bayes random-slope (patient rate shrunk to one population rate, no treatment).
5. Gompertz with population parameters fit on training folds.
6. Gradient-boosted regression on features (log v_in, dt, prior rate, has prior, treatment phase, age, MGMT, GBM yes/no). Fixed settings: sklearn HistGradientBoostingRegressor defaults, random_state 0.
7. Fisher-KPP PDE: not rerun. Its earlier MU mask result did not beat persistence (`output/mu_fitted_response/results.json`). Named, not scored here.

## TCKS (fixed before scoring)
- Phase of the forecast interval, from its midpoint day m:
  - post_progression: first-progression day <= day_in (known at forecast time).
  - chemoradiation: RT start <= m <= RT end.
  - adjuvant: additional-therapy start <= m <= its end.
  - other: none of these.
- Assumption: the treatment schedule is known at forecast time.
- r_phase = training median of observed log-rate (y_out - y_in)/dt in that phase.
- r_hist = prior-scan rate for the patient (if any).
- Shrinkage weight w = tau2 / (tau2 + s2 / dt_prior^2); tau2 and s2 estimated on training folds by method of moments (variance of (r_obs - r_phase) across pairs, and residual noise of log-volume from pairs with dt < 60 days). w = 0 if no prior scan.
- y_hat = y_in + (w r_hist + (1 - w) r_phase) dt.
- No tuning after scoring. Any change is a new, labelled version.

## Weak points (stated in advance)
- Reused cohort; exploratory only.
- Volume is from the dataset's automatic segmentation, not checked here.
- Treatment schedule known in advance is an assumption.
- Novelty: phase-conditioned rate shrinkage overlaps with tumor-growth-inhibition mixed-effects models (Ribba 2012, Claret 2009). Prior-art check is in the result file.
