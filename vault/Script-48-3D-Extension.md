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
Adaptive preserves 5k-40k mm³ residual with 88-99% dose sparing

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
- PatientID_0123 V0=3 mm³ seeds at clip floor 6235 mm³

## Files
- src/48_3d_extension.py
- output/3d_master_cohort_volumes.npz (15 MB)
- output/3d_extension_summary.json