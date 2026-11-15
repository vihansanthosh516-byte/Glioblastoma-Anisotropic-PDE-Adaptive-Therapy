# Track B: visit-informed growth hybrid (VIG) — protocol, locked before scoring (2026-11-14)

Question, target and score are unchanged: from a post-treatment scan plus scan dates, forecast y = log(core voxels + 1) at the next scan; core = necrosis + enhancing; score = mean absolute error, averaged per patient first.

## Status of the data (read this first)

- Development: MU-Glioma-Post (334 pairs, 135 patients), UCSD-PTGBM (46 pairs, 33 patients), UCSF-ALPTDG (186 pairs, 186 patients). All three are spent. They were used freely to design this method. Nothing from them is a test result.
- Outside set: LUMIERE (79 patients, 400 eligible pairs; `py -3.14 src/track_b_vig_lumiere.py --inventory`, input-side counts only).
- **LUMIERE is not untouched.** Its volumes were used before: script 91 (enhancing-volume forecast), script 104 (noise floor, ledger 2ab), script 130 (H-3v frozen test of a different model on enhancing volume, ledger 2ba), H-4 training (61 patients, ledger 2bh), and RANO scripts 140/142. To my knowledge the target of this protocol (necrosis + enhancing from DeepBraTumIA, next-scan, all consecutive pairs) was not scored before; this was not verified exhaustively. A pass here is therefore an "exposed outside set" result, not the untouched confirmation the project rule asks for.
- Before this lock, only input-side LUMIERE facts were looked at for this protocol: table layout, voxel spacing, number of scans, pairs and patients, scan-interval quantiles, share of pairs with an earlier scan, input-size quantiles. No outcome statistic was computed.

## Method (one method, fixed)

Hybrid of three model families: an exponential growth law, a visit-process (scan scheduling) model, and a history term for the last measured change.

For a pair with input log size `y_in`, interval `dt` days, and last observed change `c = y_in - y_prev` (0 if there is no earlier post-operative scan with core > 0):

1. Visit weight: `w(dt) = exp(-dt / tau)`. It stands for the chance that growth was suspected, given how soon the next scan was done.
2. Size change: `g(dt) = w * rho * dt + (1 - w) * drift`.
3. History: `h(c) = lam_rise * max(c, 0) + lam_fall * min(c, 0)`.
4. Forecast: `y_hat = y_in + g(dt) - h(c)`.

Fitting, on all 566 development pairs, absolute loss: (`rho`, `drift`, `tau`) by Nelder-Mead from three fixed starts, lowest loss kept; then (`lam_rise`, `lam_fall`) on the 213 pairs with an earlier scan, coordinate grid 0 to 0.6 in steps of 0.02. The fit is deterministic. Development fit seen before the lock: rho 0.0481/day, drift -0.045, tau 21.3 days, lam_rise 0.06, lam_fall 0.34. Curve `g`: +0.35 at 20 days, +0.22 at 45, +0.13 at 60, +0.02 at 90, -0.04 at 180.

No LUMIERE outcome is used to fit or choose anything.

## Development evidence (exploratory, post hoc; scripts and outputs in `output/track_b_loco_explore_20261114/`)

- Scan interval is the only input found whose effect has the same sign in all three cohorts (`explore1_stdout.txt`, `explore2_stdout.txt`). Input size, composition and time since surgery help on MU and UCSD and hurt on UCSF.
- Interval part alone, leave-one-cohort-out, patient MAE minus persistence: MU -0.018 (95% CI -0.038 to +0.002), UCSD -0.012 (-0.037 to +0.012), UCSF -0.026 (-0.047 to -0.004); pooled -0.022 (-0.035 to -0.008); pooled vs population rate -0.027 (-0.045 to -0.010) (`explore4_stdout.txt`, row `exp`).
- History part with separate rise and fall terms, 5-fold patient CV repeated 5 times, patients with an earlier scan, minus the interval part alone: MU -0.031, UCSD -0.102 (12 patients) (`explore5_stdout.txt`, column `h2-base`). In one repeat the pooled 95% CI was -0.095 to +0.011, so this part is not established.
- About 20 model variants were looked at in `explore2` to `explore5`. The final combination as a whole was not scored leave-one-cohort-out. All of this is selection on spent data.

## Outside test (one run)

- Source: `data/external/lumiere/LUMIERE-pyradiomics-deepbratumia-features.csv` (rows with Sequence = CT1; voxel count per label; 1 mm isotropic, asserted). Core voxels = Necrosis + Contrast-enhancing. A missing voxel count is taken as 0 (assumption: empty mask).
- Scans rated "Pre-Op" in `LUMIERE-ExpertRating-v202211.csv` are removed. Time = week number times 7 days.
- Pairs: consecutive remaining scans of one patient. Eligible: interval 14 to 365 days and input core > 0 (the rule of the UCSD and UCSF protocols).
- Baselines get: age; MGMT missing; GBM = 1; treatment phase "other" (LUMIERE has no treatment dates in these tables).
- Comparators, all fitted on the same 566 development pairs and frozen: persistence, two-point exponential, population rate, empirical-Bayes random slope, Gompertz, gradient boosting (the six locked baselines, same code as before), plus Gompertz and gradient boosting refitted with absolute loss.
- Primary success rule: for EACH of the eight comparators, the upper limit of the 95% CI of (VIG minus comparator), per-patient paired bootstrap, 50,000 resamples, seed 20261114, is below 0. At least 25 patients. All eight must pass.
- Reported with it: paired sign-flip test; VIG without the history term; error by interval bin.
- One run. The script refuses to overwrite `output/track_b_vig_lumiere_20261114/result.json`. The result is reported unchanged, pass or fail. No rerun, no change of rule after the run.
- Multiple testing: one method, one run, one outside set; "all eight must pass" needs no split of the 5% level.

## Weakest points

- The expected gain is small: about 0.02 log units (1.5% to 3% of the error) in development.
- The interval signal is clinical behaviour (a sooner scan when growth is suspected), not tumor biology. In real use the date of the next scan is often not known at the time of the forecast. Inside the locked question the scan dates are given inputs, and every baseline that uses a rate also uses the interval.
- `rho` = 0.048/day is an effective fitted value, not a measured growth rate.
- LUMIERE median interval is 91 days, where the interval part is near 0. The run therefore mostly tests the history part, which is the weaker part.
- LUMIERE masks are fully automatic (DeepBraTumIA); MU masks were radiologist-refined; label conventions differ in number but are matched by name.
- LUMIERE is exposed (see above). 79 patients.
- Novelty is not verified. No prior-art search was done for "informative scan interval" in tumor-volume forecasting.
