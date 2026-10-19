> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Script 48 — 3D Volumetric Extension

**Status:** Committed [hash]

## Cohort
- Eight selected fitted MU-Glioma inputs in an assumed-treatment numerical simulation (same selection as script 47)
- rho: 0.0001 to 0.11 /day
- Grid: 50³ voxels, 180 days, dt=0.05

## Key findings

### Simulated thresholded volume reaches zero
Five of eight modeled cases reach thresholded voxel volume0 by day180; not biological eradication.
Model dose remains100%; no clinical full-response conclusion follows.

### Adaptive holds residual, saves drug  
Adaptive preserves 8,339-41,503 mm³ residual with 10.9-98.9% dose sparing (corrected 2026-10-04 from `output/3d_extension_summary.json`)

### Legacy-run dose sparing by selected rho (not a general rule)
| rho bucket | Dose sparing |
|---|---|
| < 0.01 | 89-99% |
| 0.012 | 87.8% |
| 0.037 | 67.8% |
| 0.11 | 10.9% |

Same stratification as scripts 44 and 47.

## Legacy aggregate statistics (n=8, mean±SD; simulated grid outcomes)
- MTD: 6217 ± 11304 mm³
- Adaptive: 21318 ± 13227 mm³  
- Dose sparing: 81.4% ± 30.4%

## Limitations
- 8-patient cohort (not full 61)
- Tract orientations synthetic
- **Legacy diffusivity:** D_parallel = 0.013, D_perp = 0.0013 mm²/d, about 10x below Swanson-type values (see [[Script-78-Horizon-Crossover]]). All 3D numbers inherit this. A separate alternate-regime rerun exists below; do not pool the two regimes. Comparable literature parameter sourcing is not established here.
- PatientID_0123 V0=3 mm³ seeds at clip floor 6235 mm³

## Files
- src/48_3d_extension.py
- output/3d_master_cohort_volumes.npz (15 MB)
- output/3d_extension_summary.json
## Literature-D re-run (2026-10-04)
Run with `GBM_D_REGIME=swanson` (D 0.13 / 0.013). Outputs carry a `_swanson` suffix: `output/3d_extension_summary_swanson.json`, `3d_master_cohort_volumes_swanson.npz`. Sparing falls in three of the four higher-rho selected cases; one rises (0.887->0.837, 0.877->0.887, 0.678->0.346, 0.109->0.051); the four lowest-rho selected cases stay about0.989. MTD thresholded volume is zero in five of eight, not every case. n=8. Default (no env var) still runs legacy D.

Source keys in root 3d_extension_summary{,_swanson}.json: cohort_size,grid_size,sim_days,dt_days,patients[].adaptive.{final_volume_mm3,dose_sparing_fraction}, aggregate. Tracts and treatment effects are assumed; legacy residual-volume and sparing ranges are eight-case model ranges with substantial variance. PatientID_0123 does not start at fitted V0: its initial_volume_mm3=6235. The single rho trend is sensitive to regime and cannot guide patient therapy.
