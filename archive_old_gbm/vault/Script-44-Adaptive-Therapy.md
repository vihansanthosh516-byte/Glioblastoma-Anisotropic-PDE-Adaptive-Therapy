> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Script 44 — Adaptive Therapy Engine

**Status:** Committed 17bfade. Historical threshold status superseded by ledger 2bn; no rerun authorized.

## Purpose
Simulate adaptive vs MTD dosing with fitted MU-Glioma growth inputs and assumed per-patient E_max=rho*1000. Outcomes are model simulations, not measured treatment responses.

## Cohort
- 61 positive-rho patients (rho > 0 from fitted longitudinal data)
- 42 negative-rho fits excluded from 103 records (40.8%). Calling them treatment responders is an assumption; no MGMT-response-prevalence comparison is established (`output/adaptive_cohort_summary.json:cohort_filter`).

## Key parameters
Source: `src/44_adaptive_therapy.py`, ADAPTIVE_DT and THRESHOLD_OFF/ON defaults. Selected simulation constants, not measurements. Root summary keys e_max_ratio/horizon_days/horizon_steps match the stated ratio/horizon.
- E_MAX_RATIO = 1000 (per-patient E_max = rho * 1000)
- THRESHOLD_OFF = 0.80 (drug holiday when tumor < 80% of baseline)
- THRESHOLD_ON = 1.0 (resume drug when tumor > 100%)
- Horizon: 500 days (5000 steps at dt=0.1 day)

## Results (threshold=0.80)
| Metric | MTD | Adaptive |
|---|---|---|
| Mean drug use | 100% | 32.36% (`mean_drug_percent_of_mtd`) |
| Median modeled resistant fraction, ratio1000 sweep | about 99% | about 31%; `output/emax_sweep_summary.json`, not the old vault-only exact 99.2/30.6 |
| Mean final-mass ratio adaptive/MTD | — | 2.3913 (`adaptive_cohort_summary.json:mean_vr_ratio`); old median1.99 unverified |
| Patients with adaptive earlier progression | — | 10/61 |

## Resistance selection
- Source `output/emax_sweep_summary.json`, ratio1000: 54/61 modeled cases have more resistance under MTD; seven ties
- 0/61 patients: adaptive selected more resistance
- Historical rho/sparing r=-0.84 is unverified without an exact output key; no inferential claim retained.

## Files
- src/44_adaptive_therapy.py
- output/adaptive_cohort_summary.json
- output/adaptive_geometry_metrics.json (61 records, 30 fields each)
- output/adaptive_therapy_comparison.png

## Threshold test
The old .50 test attribution to baseline adaptive_cohort_summary.json is insufficient: its .80 copy repeats baseline values. Exact provenance for 5.8 holidays at .50 is unverified. Later declared eight-cell sweep is `vault/audit_2026-10-18/P4_threshold_sweep_results.csv` and report, ledger2bn: assumed-model drug use 27.4–43.3% across cells, progression pattern sensitive to ON/OFF. This is not clinical rescue or authority to retune.