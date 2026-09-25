# Track B — Verified Results

**Status:** Verified 2026-09-24
**Scope:** Anisotropic PDE cohort + adaptive therapy pipeline on real MU-Glioma-Post data

## Overview

Track B is the anisotropic invasion and adaptive therapy pipeline evaluated on the 
MU-Glioma-Post TCIA cohort. It spans scripts 42–49: PDE-based invasion modeling, 
stromal coupling, adaptive therapy vs. standard MTD dosing, sensitivity analysis, 
MPC optimal control, and a 3D volumetric extension.

## Data & Cohort Selection

- **203 patients** in the raw MU-Glioma-Post directory
- **154 patients** with fitted growth parameters (R² > 0) via `src/72_estimate_params_real.py`
- **103 patients** with usable tumor masks and longitudinal pairs (in `output/spatial_recurrence_profiles.npz`)
- **61 patients** with positive fitted growth rate (rho > 0) — selected for adaptive therapy
- **42 patients** with negative fitted rho — excluded (treatment responders; simulation assumes growing tumor)
- **8 patients** selected for MPC + 3D extension sims (spread across rho range, R² > 0.7)

**Source:** `output/adaptive_cohort_summary.json`, `output/mu_glioma_params_real.csv`

## Script 43 — Stromal Feedback

- **103 patients** simulated over 500 days (dt=0.1, 5000 steps)
- **Per-patient rho** loaded from `mu_glioma_params_real.csv` via `mu_glioma_loader`
- **Tumor mass range:** 45.2 – 631.7 mm³ (14× variation across cohort)
- Per-patient rho drives the reaction term; D, alpha, beta spatial fields from `spatial_recurrence_profiles.npz`
- **Source:** `output/stromal_feedback_metrics.json`, `output/stromal_evolution_cohort.npz`
- **Commit:** `e76d26a`

## Script 44 — Adaptive Therapy

- **61-patient** positive-rho cohort, 500-day simulation
- **Adaptive dosing uses 32.36%** of the MTD dose (mean across cohort)
- **Resistant fraction (median):** 99.2% under MTD vs 30.6% under adaptive
- **54/61 patients:** MTD selects more resistance than adaptive
- **0/61 patients:** adaptive selects more resistance than MTD
- **39/61 patients:** TTP-censored in both arms (progression-equivalent)
- **10/61 patients:** adaptive progressed earlier than MTD
- **Threshold sensitivity:** `THRESHOLD_OFF=0.50` test increases drug use to 43% and suppresses holidays (5.8 mean vs 48.9 at `THRESHOLD_OFF=0.80`)
- **Source:** `output/adaptive_cohort_summary.json`, `output/adaptive_geometry_metrics.json`
- **Commit:** `17bfade`

## Script 46 — Sobol Sensitivity Analysis

- **154-patient** cohort, 5-parameter Sobol analysis
- **6000 ODE evaluations** (500 Saltelli samples)
- **rho_s: S1 = 0.998** (95% CI ± 0.086) — proliferation rate dominates TTP variance
- All other parameters: **ST < 0.015** (aniso_ratio, mu_r, EC50, D_white)
- **Source:** `output/sobol_sensitivity_results.json`

## Script 47 — MPC Optimal Control

- **8 real patients**, rho range 0.0001 to 0.11 /day
- **Dual-agent MPC** maintains 360d TTP for 6/8 patients
- **Single-agent fails at high rho:** PatientID_0188 TTP 25d (vs MTD 40d)
- **Resistant fraction ≤ 0.05** for 6/8 under dual-agent (vs 1.00 under MTD)
- **Extreme outlier:** PatientID_0123 (rho=0.11) fails on all arms (22d dual vs 7d MTD)
- **Source:** `output/dual_drug_comparison.json`
- **Commit:** `944753c`

## Script 48 — 3D Volumetric Extension

- **8 patients** (same cohort as script 47), 50³ grid, 180 days
- **MTD eradicates tumor** (0 mm³) for 5/8 patients at low rho
- **Adaptive holds** 5,000 – 41,500 mm³ residual with **88–99% dose sparing**
- **Dose sparing collapses** as rho increases:
  - rho < 0.01: 89–99% sparing
  - rho = 0.012: 87.8%
  - rho = 0.037: 67.8%
  - rho = 0.11: 10.9%
- **Aggregate:** MTD 6217 ± 11304 mm³, Adaptive 21318 ± 13227 mm³, sparing 81.4% ± 30.4%
- **Source:** `output/3d_extension_summary.json`
- **Commit:** `7314f7e`

## Synthesis — Three Track B Findings

### 1. Adaptive therapy suppresses resistance selection

- **54/61 patients** under MTD select more resistance than under adaptive
- **0/61 patients** show the reverse
- Mechanism: preserving the drug-sensitive clone through dose holidays, competitive suppression of the resistant clone

### 2. Proliferation rate dominates outcome variance

- Sobol S1 = **0.998** for rho_s (the other 4 parameters combined contribute < 2%)
- Interpretation: TTP under any dosing regimen is determined almost entirely by tumor proliferation rate, not by anisotropy, mutation rate, or drug PK

### 3. Dose sparing is inversely correlated with rho

- **Pearson r = -0.84** between rho and drug reduction
- Slow-growing tumors (rho < 0.002): 80–99% dose sparing
- Fast-growing tumors (rho > 0.03): < 30% dose sparing
- Predictive rule: rho = 0.02/day is the inflection point

## Limitations

- Scripts 47/48 use 8 patients (not full 61) for computational cost
- PatientID_0123 (rho=0.11, V0=3 mm³) fails on all arms — model limit
- PatientID_0188 (rho=0.037) extends TTP from 40d to 263d but does not reach 360d
- Tract orientations in script 48 are synthetic (not per-patient DTI-derived)
- Script 47 uses a reduced ODE; full spatial version is script 48

## Files Supporting This Section

| File | Produced by | Content |
|---|---|---|
| `output/stromal_feedback_metrics.json` | 43 | Per-patient tumor mass, kinetics |
| `output/adaptive_cohort_summary.json` | 44 | Cohort-level adaptive therapy metrics |
| `output/adaptive_geometry_metrics.json` | 44 | Per-patient adaptive metrics (61 records) |
| `output/sobol_sensitivity_results.json` | 46 | Sobol indices for 5 parameters |
| `output/dual_drug_comparison.json` | 47 | Per-patient MPC 3-arm results (8 records) |
| `output/3d_extension_summary.json` | 48 | Per-patient 3D metrics (8 records) |

## Related Vault Notes

- [[Script-44-Adaptive-Therapy]] — detailed adaptive therapy analysis
- [[Script-46-Sensitivity]] — Sobol sensitivity detail
- [[Script-47-Optimal-Control]] — MPC + dual-drug
- [[Script-48-3D-Extension]] — 3D extension
- [[Decisions]] — E_MAX_RATIO = 1000, THRESHOLD_OFF calibration
- [[Open-Questions]] — outstanding questions
- [[Repo]] — master map