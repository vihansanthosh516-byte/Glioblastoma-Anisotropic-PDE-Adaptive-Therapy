# Track B: development of a second hybrid ("VIG-2") — written before any run (2026-11-19)

Question, target, score unchanged (next-scan log core volume, patient-level MAE). This is **development**, not a test. No fresh outside cohort exists yet (`paper/trackB_next_data_plan_2026-11-19.md`).

## Data (development pool)
All six cohorts we hold: MU, UCSD, UCSF, LUMIERE (`output/track_b_code_audit_20261114/pooled4_pairs.csv`, 976 pairs) and CFB, WAW (`output/track_b_vig_cfb_waw_20261115/test_pairs.csv`, 295 pairs; spent as a test, now development data). 1,271 pairs.

## Honest limit, stated first
The idea below was chosen after seeing that Gompertz with absolute loss was the best method on CFB+WAW (`result.json`). So the CFB and WAW folds are **not clean**: the candidate family was picked because it fit them. The only clean test is a future untouched cohort.

## Known facts that motivate the candidates
- Earlier four-cohort leave-one-cohort-out run (`fix_and_loco4_stdout.txt`): Gompertz-L1 (pull toward a population level K) helped on LUMIERE and MU, hurt on UCSF (+0.126) and UCSD. VIG helped on MU and UCSF.
- CFB+WAW: Gompertz-L1 best (1.119), VIG 1.151, persistence 1.154.
- So level-reversion pays where the input size is noisy and costs where it is not. A hybrid that has both an interval/history part (VIG) and a level-reversion part, fitted together, is the untried object. (Blend 0.5·(VIG + Gompertz-L1) was tried earlier and did not help pooled; joint fitting is different.)

## Candidates (all absolute loss, fitted on the five training cohorts of each fold)
All use `y_hat = y_in + change`. `g` = interval curve of VIG (`vig_curve`, params rho, drift, tau). `r(y) = (1 - exp(-a·dt)) · (K - y_in)` = Gompertz level-reversion (params K, a ≥ 0). `h` = history pull-back of VIG (lam_rise, lam_fall; grid search as in `fit_vig`).
- **H1 VIG+G:** `y_in + g(dt) + r − h`. Five smooth parameters fitted jointly by Nelder-Mead on absolute loss (several starts), then lam by grid.
- **H2 G+hist:** `y_in + r − h` (no interval curve).
- **H3 VIG+Gpeak:** like H1 but K is the patient's earlier peak when one exists (`max(y_peak_prev, y_in)`), else the population K.

References, same folds: persistence, population rate, Gompertz-L1 (`fit_gomp_l1`), VIG (`fit_vig`, unchanged), pk_gomp.

## Fold design and scoring
Leave one cohort out, six folds. Per fold: patient-level MAE of each method (error averaged per patient first). Difference to persistence per cohort with a patient-level bootstrap (20,000 resamples, seed 20261119).
Summaries: (a) cohort-equal mean of the six differences; (b) number of cohorts better than persistence; (c) same on the four earlier cohorts only; (d) vs VIG and vs Gompertz-L1.

## Selection rule (fixed now)
Pick the candidate with the lowest cohort-equal mean difference to persistence, if it also (1) beats persistence in at least 4 of 6 cohorts and (2) has a lower mean than both VIG and Gompertz-L1. If no candidate does this: **nothing is frozen; the result is "no development gain"**. No candidate is added, dropped or retuned after the run. Grid values and starts are those in the code written for this run.

## What a pass would mean
A pass means only: a candidate to freeze and test once on a future untouched cohort (power: `paper/trackB_vig_cfb_waw_protocol_2026-11-15.md`, section 8, redone for the new effect size). It is not a claim of improvement.
