# Script 46 — Sobol Sensitivity Analysis

**Status:** Committed

## Method
- 5-parameter reduced-ODE Sobol analysis
- Real MU-Glioma parameter distributions (rho, D_white from 154 patients)
- N=500 Saltelli samples → 6000 ODE evaluations
- E_MAX calibrated to scripts 44/47/48 (`E_MAX_RATIO = 1000`)

## Key finding

**Proliferation rate dominates TTP variance:**
- S1 for `rho_s` = **0.998** (95% CI ±0.086)
- ST for `rho_s` = 1.015 (within Monte Carlo noise of 1.0)
- All other parameters ST < 0.015

| Parameter | S1 | ST |
|---|---|---|
| rho_s | 0.998 | 1.015 |
| aniso_ratio | 0.002 | 0.014 |
| mu_r | -0.001 | 0.001 |
| EC50 | 0.003 | 0.005 |
| D_white | 0.002 | 0.005 |

**Interpretation:** Treatment outcome is determined almost entirely by 
tumor proliferation rate. The other parameters (anisotropy, mutation rate, 
drug sensitivity, diffusivity) contribute negligibly to TTP variance 
under the calibrated model.

## Files
- `src/46_sensitivity_analysis.py`
- `output/sobol_sensitivity_results.json`
- `output/sobol_tornado_plot.png`

## Related
- [[Script-44-Adaptive-Therapy]] — same real MU-Glioma cohort
- [[Script-47-Optimal-Control]] — MPC on the same parameter ranges
- [[Decisions]] — E_MAX calibration