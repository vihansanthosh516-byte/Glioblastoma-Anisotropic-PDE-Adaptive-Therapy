# Script 80 — PDE Sobol (Swanson vs Legacy)

**Status:** Committed

## Method
- Real PDE Sobol analysis at Swanson D (D_white=0.13, D_gray=0.013 mm²/d)
- Comparison against same Sobol at legacy D (0.013/0.0013)
- N=128 base samples (wide CIs)
- Replaces the reduced-ODE script 46 result (S1=0.998 was an artifact)

## Result (volume, Swanson regime)
| Parameter | S1 | ST |
|---|---|---|
| α_sens | **0.43 ± 0.14** | 0.76 |
| ρ_s | 0.20 ± 0.18 | 0.39 |
| D_w | 0.000 | 0.000 |

- D_w enters extent weakly (ST = 0.020 ± 0.012)
- The archived S1(rho_s)=0.998 is retracted

## Source
- output/sobol_pde_swanson_vs_legacy.json

## Related
- [[Script-46-Sensitivity]] (superseded reduced-ODE)
- [[Script-60-66-Swanson-D]]