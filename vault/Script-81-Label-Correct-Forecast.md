# Script 81 — Label-Correct Forecast (Core Target)

**Status:** Committed (bc223ff)

## What Changed
Previous target was `mask > 0`, which includes resection cavity (label 4) and edema (label 2) — both change from surgery/treatment, not growth. Script 81 scores only the enhancing + necrotic core (labels 1, 3).

## Result (133 patients)
| Arm | Out-of-fold Dice |
|---|---|
| anisotropic | 0.259 |
| iso_same | 0.260 |
| iso_homog | 0.261 |
| **no_change** | **0.233** |

- Anisotropic vs no_change: +0.026 (95% CI 0.008–0.044, Holm p = 0.0015)
- Anisotropic vs iso_same: −0.001 (p = 0.97)
- Model beats no-change in 55% of patients

## Limitations
- Absolute Dice is low (0.26)
- Modest gain (+11% relative)
- Model still loses on shrinking tumors (−0.039 Dice)

## Source
- output/forecast_labels_core/results.json

## Related
- [[Script-78-Horizon-Crossover]]
- [[NEGATIVES_REVISITED]] (B1c row)