# Track B: VIG-2 development result (2026-11-19)

Pre-registration (committed before the run, `0ba4dd4`): `paper/trackB_vig2_prereg_2026-11-19.md`. Script `src/track_b_vig2_dev.py`; every number below is from `output/track_b_vig2_dev_20261119/dev_result.json`. Development only, not a test. Six cohorts, leave one cohort out.

## Result by the pre-registered rule
Selected: **H1 (VIG interval curve + Gompertz level-reversion + history pull-back, all fitted jointly with absolute loss).** H2 also passed; H3 did not.

Difference to persistence (patient-level MAE, negative = better), held-out cohort, mean and 95% patient bootstrap:

| Cohort (patients) | VIG | H1 | H2 |
|---|---|---|---|
| MU (135) | -0.034 [-0.056, -0.013] | -0.066 [-0.130, -0.001] | -0.048 [-0.114, 0.017] |
| UCSD (42) | -0.021 [-0.063, 0.018] | +0.018 [-0.064, 0.104] | +0.016 [-0.070, 0.115] |
| UCSF (186) | -0.024 [-0.041, -0.008] | **+0.083 [0.018, 0.152]** | +0.108 [0.042, 0.178] |
| LUMIERE (79) | -0.022 [-0.050, 0.006] | -0.072 [-0.108, -0.036] | -0.062 [-0.092, -0.033] |
| CFB (158) | +0.010 [-0.012, 0.034] | -0.009 [-0.049, 0.033] | -0.008 [-0.046, 0.030] |
| WAW (23) | -0.043 [-0.251, 0.084] | -0.149 [-0.455, 0.085] | -0.197 [-0.516, 0.057] |
| **Cohort-equal mean, six** | **-0.0224** | **-0.0324** | -0.0320 |
| Cohorts better than persistence | 5 of 6 | 4 of 6 | 4 of 6 |
| Cohort-equal mean, first four (clean) | **-0.0253** | -0.0091 | +0.0035 |

## Four questions
1. **Positive?** Barely. H1 passes the rule, with a six-cohort mean 0.010 lower than VIG. That gap has no interval and is small.
2. **How much?** About 0.03 log units better than persistence on average (about 3% of a 1.0 error). Same size as the old exploratory VIG effect.
3. **More unique?** H1 is a new joint model (interval curve + level reversion + history), but it is not clearly better than VIG.
4. **Fix?** See "what this means".

## Weakest points (read these first)
- **The gain over VIG is not clean.** On the four cohorts that were not used to choose the idea, VIG is better than H1 (-0.025 vs -0.009). H1's edge over VIG comes from CFB, WAW and LUMIERE, and WAW has only 23 patients.
- H1 is clearly worse than persistence on UCSF (+0.083, interval above 0). The level-reversion hurts where input sizes are accurate.
- CFB and WAW folds are contaminated: the candidate family was chosen after seeing Gompertz-L1 win there.
- The fitted level K is not stable: in the CFB fold K = 15.3 in log units (a = 0.0017 per day), outside the data range. The reversion term then acts as a linear drift, not a pull toward a size.
- Selection was by a six-cohort mean in which cohort differences are large; no formal test of H1 vs VIG was in the pre-registration.

## What this means
There is no development evidence that a joint hybrid beats VIG on cohorts where VIG works, and the pattern across cohorts says: **level-reversion helps when input size is noisy, hurts when it is accurate.** CFB (tool-made targets) and WAW (resection-cavity effects) are the noisy cohorts; UCSF is the accurate one.

That points to the real next idea: let the pull strength depend on how noisy the input size is (a noise-aware shrinkage), instead of one global rate. It must be written down as a new pre-registration with its own candidate list before it runs, and it needs a clean cohort to test.

Nothing is frozen. No test is spent. A fresh cohort is still required (`paper/trackB_next_data_plan_2026-11-19.md`).
