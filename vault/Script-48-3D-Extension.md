# Script 48 — 3D Volumetric Extension

**Status:** Committed [hash]

## Cohort
- 8 real MU-Glioma patients (same selection as script 47)
- rho: 0.0001 to 0.11 /day
- Grid: 50³ voxels, 180 days, dt=0.05

## Key findings

### MTD eradicates tumor at low rho
5 of 8 patients: MTD drives tumor to 0 mm³ by day 180
Full response but at 100% drug exposure

### Adaptive holds residual, saves drug  
Adaptive preserves 8,339-41,503 mm³ residual with 10.9-98.9% dose sparing (corrected 2026-10-04 from `output/3d_extension_summary.json`)

### Dose sparing collapses with rho
| rho bucket | Dose sparing |
|---|---|
| < 0.01 | 89-99% |
| 0.012 | 87.8% |
| 0.037 | 67.8% |
| 0.11 | 10.9% |

Same stratification as scripts 44 and 47.

## Aggregate statistics
- MTD: 6217 ± 11304 mm³
- Adaptive: 21318 ± 13227 mm³  
- Dose sparing: 81.4% ± 30.4%

## Limitations
- 8-patient cohort (not full 61)
- Tract orientations synthetic
- **Legacy diffusivity:** D_parallel = 0.013, D_perp = 0.0013 mm²/d, about 10x below Swanson-type values (see [[Script-78-Horizon-Crossover]]). All 3D numbers inherit this. Scale and rho-dependence of sparing were not re-run at literature D.
- PatientID_0123 V0=3 mm³ seeds at clip floor 6235 mm³

## Files
- src/48_3d_extension.py
- output/3d_master_cohort_volumes.npz (15 MB)
- output/3d_extension_summary.json
## Literature-D re-run (2026-10-04)
Run with `GBM_D_REGIME=swanson` (D 0.13 / 0.013). Outputs carry a `_swanson` suffix: `output/3d_extension_summary_swanson.json`, `3d_master_cohort_volumes_swanson.npz`. Sparing for the four high-rho patients falls (0.887->0.837, 0.877->0.887, 0.678->0.346, 0.109->0.051); the four low-rho patients stay at 0.989. MTD still reaches 0 mm3. n=8. Default (no env var) still runs legacy D.
