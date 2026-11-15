# Track B MU method tries log (2026-11-14) — exploratory, reused cohort

Same question, data (334 pairs, 135 patients), folds, target and score for every try. Best baseline in all tries: population growth rate, patient MAE 1.528 (log volume).

| Try | Method | Patient MAE | Minus best (CI) | CI level | Pass | 25%-growth AUC | Result file |
|---|---|---|---|---|---|---|---|
| 1 | TCKS (phase-rate shrinkage) | 1.640 | +0.112 (0.063, 0.163) | 95% | no | 0.438 | output/track_b_mu_tcks_20261114/result.json |
| 2 | TCKS v2 (noise-aware latent size) | 1.555 | +0.027 (-0.004, 0.059) | 95% | no | 0.531 | output/track_b_mu_tcks_v2_20261114/result.json |
| 3 | v3 median latent regression | 1.529 | +0.001 (-0.071, 0.079) | 98.33% | no | 0.587 | output/track_b_mu_v3_20261114/result.json |
| 4 | v4 median-loss boosted | 1.582 | +0.054 (-0.115, 0.221) | 98.75% | no | 0.567 | output/track_b_mu_v4_20261114/result.json |
| 5 | v5 input-mask composition + latent | 1.497 | -0.031 (-0.148, 0.100) | 99% | no | 0.641 | output/track_b_mu_v5_20261114/result.json |
| 6 | v6 composition + edge shape | 1.505 | -0.023 (-0.155, 0.117) | 99.17% | no | 0.658 | output/track_b_mu_v6_20261114/result.json |
| 7 | v7 MRI brightness + composition | 1.508 | -0.020 (-0.141, 0.102) | 99.29% | no | 0.636 | output/track_b_mu_v7_20261114/result.json |
| 8 | v8 two-part collapse model | 1.449 | -0.079 (-0.287, 0.121) | 99.375% | no | 0.650 | output/track_b_mu_v8_20261114/result.json |

Weak points: every try after 1 is post hoc; the cohort was used before; segmentations are the dataset's own and not checked here; best baseline AUC for comparison: Gompertz 0.580.

## External test of frozen v8 (UCSD-PTGBM, one run)
Protocol: paper/trackB_v8_ucsd_external_protocol_2026-11-14.md. Fit on all MU, applied once. 46 pairs, 33 patients.
Patient MAE: v8 1.219; persistence 1.342; Gompertz 1.438; pop_rate 1.441; EB slope 1.488; two-point 1.582; GBR 1.687.
v8 minus pop_rate (pre-set comparator): -0.222, 95% CI (-0.700, 0.066). Pass: no (CI crosses 0). v8 was lowest of all 7 methods.
Weak points: small n (33 patients); UCSD exposed before for a different target; phase windows assumed; only 1 of 4 collapse pairs flagged.
Source: output/track_b_v8_ucsd_external_20261114/result.json

## Hybrid search (6 subagent candidates, protocol paper/trackB_hybrid_search_protocol_2026-11-14.md)
Dev patient MAE (99.64% CI rule; none passed): moe 1.453, compartment 1.463, pkpd 1.473, stack 1.474, deconv 1.480, switch 1.590 (worse than pop_rate). v8 dev 1.449 stays lowest.
Finalist moe on UCSD (second exposure): 1.256 vs pop_rate 1.441; diff -0.184, 95% CI (-0.643, 0.079). Not passed. v8 on UCSD stays lowest (1.219). Source: output/track_b_hybrid_moe_ucsd_20261114/result.json

## Frozen v8 on UCSF-ALPTDG (protocol paper/trackB_v8_ucsf_alptdg_protocol_2026-11-14.md, one run)
Eligible pairs 186 of 298 (dt 14-365 days, non-empty time1 core). PRIMARY untouched patients: 137.
Primary patient MAE: persistence 0.832, pop_rate 0.851, Gompertz 0.955, v8 1.072, GBR 1.205.
v8 minus pop_rate: +0.221, 95% CI (0.102, 0.340). v8 is significantly WORSE. Primary success: no.
Secondary (49 used-before patients): v8 0.880 vs pop_rate 0.743; diff +0.137, CI (-0.070, 0.351).
Reading: with expert masks (baseline error 0.83 vs 1.53 on MU), v8's MU-learned corrections hurt. The MU gains look tied to MU's automatic-mask noise, not general tumor biology.
Source: output/track_b_v8_ucsf_alptdg_20261114/result.json

## VIG hybrid on LUMIERE (protocol paper/trackB_vig_lumiere_protocol_2026-11-14.md, locked in commit 5a68658, one run)
Method: scan-interval growth weight + rebound-after-fall history term. Fit on all 566 development pairs (MU + UCSD + UCSF). LUMIERE: 400 pairs, 79 patients, DeepBraTumIA masks. LUMIERE was exposed before for other targets; it is not untouched.
Patient MAE: Gompertz 1.220, VIG 1.239, Gompertz (absolute loss) 1.239, persistence 1.261, GBR (absolute loss) 1.264, pop_rate 1.289, GBR 1.317, EB slope 1.497, two-point 2.041.
VIG minus persistence -0.022, 95% CI (-0.056, +0.012). VIG minus Gompertz +0.019 (-0.039, +0.080). VIG minus pop_rate -0.050 (-0.088, -0.011).
Primary success (all eight CI upper limits below 0): NO. Passed 3 of 8 (two-point, pop_rate, EB slope).
Secondary: without the history term VIG is +0.007 vs persistence (-0.009, +0.024), so the scan-interval signal seen in MU, UCSD and UCSF did not show on LUMIERE. The history term gave about -0.029.
Weak points: LUMIERE exposed; 79 patients; about 20 variants were looked at on development data before the lock; LUMIERE has no treatment dates in the tables used.
Source: output/track_b_vig_lumiere_20261114/result.json
Development exploration (post hoc): output/track_b_loco_explore_20261114/. Problem review: paper/trackB_problem_review_2026-11-14.md.

## Code and data audit + four-cohort leave-one-cohort-out (exploratory, post hoc)
Record: paper/trackB_code_audit_2026-11-14.md. No large code bug found. Two small data fixes (UCSD 46 -> 55 pairs; one LUMIERE pre-operative scan).
Pooled over 442 patients, each cohort held out in turn: VIG 1.121, persistence 1.148, pop_rate 1.154, Gompertz (abs. loss) 1.189, Gompertz 1.217. VIG minus persistence -0.027, 95% CI (-0.039, -0.014); every comparator CI upper limit is below 0.
NOT a confirmed win: the pooled rule was chosen after the locked LUMIERE run failed; all four cohorts are spent; the VIG design was chosen with three of the four cohorts in view; gain is about 2.4%.
Source: output/track_b_code_audit_20261114/loco4_all_baselines_stdout.txt

## Prior-art check for VIG (2026-11-15, literature only, no new run)
Record: paper/trackB_vig_prior_art_2026-11-14.md.
The scan-interval effect (shorter gap, more enhancing growth) is already published on a UCSF-ALPTDG subset: Lindner & Fiehler 2026, PMC13421742, odds ratio 0.91 per day (0.84 to 0.99). Informative visit times are a known idea.
Not found: a next-scan volume forecast that uses the interval as an informative-visit weight in a growth law with a two-slope reversal term, tested against all six baselines.
Weak points: one paper index, about 16 queries, most papers abstract only. "Not found" is not "proven absent".
