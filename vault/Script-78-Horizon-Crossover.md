# Script 78: Horizon-crossover ablation (B2 revisited)

**Date:** 2026-10-02
**Script / output:** `src/78_ablation_horizon_crossover.py` -> `output/ablation_horizon_crossover.json` (16 cells, 76,490 s compute, all values read back from disk).
**Design:** fixed in the script header before the run; no change after seeing data. The full pre-specified grid completed, so there is no deviation. (An earlier note in this session called the run "killed" and proposed cutting the 1 mm arm. That was wrong: the harness stopped only its shell wrapper, the python process kept running, and all 16 cells finished.)

## Question
Script 77 found the ablation null at 2 mm / 90 d and called it structural. Script 60 uses D_white = 0.013, D_gray = 0.0013 mm²/d. Is the null a property of the model, or of those D values?

## Critical finding: script 60's D values are ~10x below the literature
| | D_white (mm²/d) | D_gray (mm²/d) |
|---|---|---|
| Script 60 | 0.013 | 0.0013 |
| Swanson-type (arXiv 2402.02273 and others) | 0.13 | 0.013 |
| Ratio | 0.10 | 0.10 |

The script-60 numbers match the Swanson cm²/d values read as mm²/d. The unit slip is a likely cause, inferred from the numbers alone. Not confirmed in the git history. The user's note quoted gray as 0.001; the literature gray value is 0.013.

## Method
- Same script-60 tensor and solver, natural growth, seed peak 0.08 K, rho 0.02 /d.
- Arms: `full` (DTI tensor) vs `iso_matched` (trace(D)/3 · I, same mean diffusivity; the primary comparison) and `x_aligned` (script 60's "No DTI").
- D_gray = D_white / 10. Horizons 90, 180, 365, 600, 900 d. Grids 2 mm and 1 mm. Endpoint: mask Dice at u ≥ 0.16 K. Crossover = first horizon with Dice < 0.95.

## Result 1: 2 mm sweep (rho 0.02), Dice(full vs iso_matched), u ≥ 0.16 K
| D_white | 90 d | 180 d | 365 d | 600 d | 900 d | Crossover |
|---|---|---|---|---|---|---|
| 0.01 | 0.997 | 0.988 | 0.989 | 0.992 | 0.992 | none |
| 0.03 | 0.985 | 0.970 | 0.977 | 0.980 | 0.977 | none |
| 0.1 | 0.923 | 0.931 | 0.949 | 0.950 | 0.946 | 90 d |
| 0.3 | 0.838 | 0.874 | 0.917 | 0.921 | 0.918 | 90 d |

- Effect size rises with D and with L_diff/dx. At D = 0.01 and 0.03 (script 60's range), L_diff/dx = 0.67 and 1.16 at 90 d and the effect is under 2% Dice at 90 d.
- At D = 0.1 (L_diff/dx = 2.1) the full model differs from the matched isotropic one by 8% Dice at 90 d. At D = 0.3 the difference is 16%.
- The effect is **largest early and shrinks at long horizons**, because the tumour grows to fill the domain. There is no monotonic rise with horizon.
- Script 60's own ablation arm (`x_aligned`, "No DTI") at 90 d: Dice 0.997 (D 0.01), 0.947 (0.03), 0.828 (0.1), 0.607 (0.3). Script 77 reported 0.994 at its D range.

## Result 2: 1 mm grid at all four D values (rho 0.02)
| D_white | Dice at 90 d, 2 mm | 1 mm | Dice at 900 d, 2 mm | 1 mm |
|---|---|---|---|---|
| 0.01 | 0.997 | 0.998 | 0.992 | 0.994 |
| 0.03 | 0.985 | 0.980 | 0.977 | 0.981 |
| 0.1 | 0.923 | 0.943 | 0.946 | 0.951 |
| 0.3 | 0.838 | 0.848 | 0.918 | 0.920 |

The effect is nearly the same at 1 mm as at 2 mm. It is a property of the physics, not of the grid. The crossover horizon is identical (90 d at D 0.1 and 0.3, none at 0.01 and 0.03).

Grid convergence of the full arm (1 mm downsampled vs 2 mm): Dice 0.83-0.91 at 90 d (masks are only 85-170 voxels at 2 mm), rising to 0.95-0.98 at 365 d and later. Report 90-d values with that resolution error in mind.

## Result 3: rho 0.01 and 0.04, 2 mm (secondary)
- rho 0.04: crossover none at D 0.01, 0.03, 0.1; 90 d at D 0.3. Same D ordering.
- rho 0.01: **censored at 90 d.** The tumour is 2 voxels (D 0.03) or 0 voxels (D 0.1, 0.3) above the 0.16 K threshold at 90 d. The reported crossovers for these cells (90 d, 180 d) come from an undefined or degenerate Dice and must not be quoted. From 180 d on, the Dice is 0.90-0.96 at D 0.1 and 0.2-0.90 at D 0.3 (0.200 at 180 d where the mask is still tiny).
- rho 0.04, D 0.3, 900 d: the mask touches the domain boundary (censored).

## Caveats
- Dice 0.92-0.95 at D = 0.1 is a modest effect. No comparison to imaging or segmentation noise was run, so whether it is detectable on real scans is untested.
- Elongation differences change sign across horizons (e.g. D 0.3: -0.41 at 180 d, +0.24 at 600 d), so elongation is not a clean monotone readout here. Script 74's elongation result (600 d, 103 patients) still stands for the patient-tensor model.
- D = 0.1 is close to, not equal to, the literature 0.13 mm²/d.
- Deterministic model, one seed, no statistical test.
- The comparison to script 77 is approximate: script 77 used 30 LHS scenarios, script 78 uses a fixed grid.

## Reframing
**The ablation null is a parameter artifact, not structural.** At script 60's D (0.01-0.03) diffusion cannot cross a 2 mm voxel in 90 d, so the null follows. At literature diffusivities (D ≈ 0.1) the DTI orientation changes the day-90 tumour mask by ~8% Dice, and this holds on a 1 mm grid. Volume is still insensitive (mass conservation, script 77), so the effect is a shape/extent effect. Treatment-planning conclusions drawn from script 60 (including its Sobol S1(rho) = 0.998) used a diffusivity 10x too low and should be re-checked at D = 0.1-0.13.
