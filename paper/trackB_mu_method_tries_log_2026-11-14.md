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
