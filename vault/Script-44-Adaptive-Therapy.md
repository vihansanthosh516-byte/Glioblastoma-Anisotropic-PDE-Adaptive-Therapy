# Script 44 — Adaptive Therapy Engine

**Status:** Committed 17bfade. Threshold test running.

## Purpose
Simulate adaptive vs MTD dosing on real MU-Glioma patients with 
per-patient rho and E_max.

## Cohort
- 61 positive-rho patients (rho > 0 from fitted longitudinal data)
- 42 negative-rho responders excluded (41% of cohort, matches TMZ response 
  rate in MGMT-methylated GBM)

## Key parameters
- E_MAX_RATIO = 1000 (per-patient E_max = rho * 1000)
- THRESHOLD_OFF = 0.80 (drug holiday when tumor < 80% of baseline)
- THRESHOLD_ON = 1.0 (resume drug when tumor > 100%)
- Horizon: 500 days (5000 steps at dt=0.1 day)

## Results (threshold=0.80)
| Metric | MTD | Adaptive |
|---|---|---|
| Mean drug use | 100% | 32% |
| Median resistant fraction | 99.2% | 30.6% |
| Median VR ratio | — | 1.99 |
| Patients with adaptive earlier progression | — | 10/61 |

## Resistance selection
- 54/61 patients: MTD selected more resistance
- 0/61 patients: adaptive selected more resistance
- Pearson r (rho vs drug sparing) = −0.84

## Files
- src/44_adaptive_therapy.py
- output/adaptive_cohort_summary.json
- output/adaptive_geometry_metrics.json (61 records, 30 fields each)
- output/adaptive_therapy_comparison.png

## Threshold test
Rerunning with THRESHOLD_OFF=0.50 to test high-rho hypothesis. Result pending.