# Track B: noise-aware level-reversion, development result (2026-11-19)

Pre-registration (committed first, `392396f`): `paper/trackB_vig2_noise_prereg_2026-11-19.md`. Script `src/track_b_vig2_noise_dev.py`; numbers from `output/track_b_vig3_noise_dev_20261119/noise_dev_result.json`. Development only. Six cohorts, leave one cohort out.

## Result: **no candidate passes. Nothing is frozen.**

Difference to persistence (patient-level MAE, negative = better):

| Method | MU | UCSD | UCSF | LUMIERE | CFB | WAW | Mean (6) | Mean (4 earlier) | Better in |
|---|---|---|---|---|---|---|---|---|---|
| VIG | -0.034 | -0.021 | -0.024 | -0.022 | +0.010 | -0.043 | -0.022 | -0.025 | 5 of 6 |
| H1 | -0.066 | +0.018 | +0.083 | -0.072 | -0.009 | -0.149 | -0.032 | -0.009 | 4 of 6 |
| N0 (H1, K bounded) | -0.067 | +0.020 | +0.092 | -0.067 | -0.018 | -0.191 | -0.039 | -0.006 | 4 of 6 |
| N1 (size-dependent pull) | -0.067 | +0.021 | +0.103 | -0.076 | -0.014 | -0.209 | -0.040 | -0.005 | 4 of 6 |
| N2 (N1 + history) | -0.067 | +0.021 | +0.166 | -0.075 | -0.013 | -0.211 | -0.030 | +0.011 | 4 of 6 |

N1 and N2 had a lower six-cohort mean than VIG and H1 (N1), but failed three of the four rules: the four earlier cohorts are worse than VIG, UCSF is worse than persistence by +0.10 / +0.17 (limit +0.03), and only 4 of 6 cohorts beat persistence (need 5).

## What this tells us
- Making the pull depend on input size did **not** fix UCSF. UCSF is hurt by any level-reversion, whether the size is small or not. So "small cores are noisier" is not the reason UCSF differs.
- Bounding K (N0) did not change the picture. The 15.3 value in the earlier run was not the cause.
- The gain in the six-cohort mean comes from WAW (23 patients, -0.2) and LUMIERE. It is not stable.
- **The simplest hybrid, VIG, is the most consistent:** better than persistence in 5 of 6 cohorts, worst cohort +0.010. But its gain is small (-0.022 cohort-equal mean), and on the locked outside test it was -0.003 with an interval crossing 0.

## Weakest points
- The whole family was chosen after seeing the CFB/WAW outcome; folds are partly informed.
- Six cohorts is too few to learn any cohort-level feature; I dropped the cohort-noise idea for that reason before running.
- WAW has 23 patients and an error near 3.0; it can swing the mean.

## Power note (new, for the plan)
If the true gain is about 0.022 (VIG's development mean) with a patient-level SD about 0.15 (`result.json` scale), a 5% two-sided test with 80% power needs about 365 patients against persistence alone; beating all eight comparators needs more. This is my calculation from those two numbers; it assumes the 0.022 holds, which the outside test did not support.

## Decision for the project
Stop adding level-reversion variants. Two honest paths remain on the same question:
1. **Fresh large cohort (350+ patients) with VIG unchanged** as the frozen method. Power is the issue, not the method.
2. **A different kind of signal**, not another curve on size: use the images (shape, location, enhancement fraction) as inputs. This is a new pre-registration and needs the same fresh cohort for a clean test.
