> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Script 81 — Label-Correct Forecast (Core Target)

**Status:** Committed (bc223ff)

## What Changed
Previous target was `mask > 0`, which includes resection cavity (label 4) and edema (label 2) — both change from surgery/treatment, not growth. Script 81 scores only the enhancing + necrotic core (labels 1, 3).

## Retrospective internal out-of-fold result (133 scored;19 skipped)

Not external/prospective validation. There are126 cases with nonempty target in both scans; denominators differ. Target selection followed the earlier volume-label probe (post hoc); later paired manifest comparisons use different folds/windows and are not a direct replication.
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
- Model loses on shrinking/same-size tumors: -0.0390814 Dice,n=50 (`output/forecast_labels_core/results.json:descriptive_shrank_or_same.tests.anisotropic_vs_no_change.mean_diff`)

## Source
- output/forecast_labels_core/results.json

## Related
- [[Script-78-Horizon-Crossover]]
- [[NEGATIVES_REVISITED]] (B1c row)