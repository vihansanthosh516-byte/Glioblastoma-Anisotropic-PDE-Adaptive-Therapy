# Script 46 — Sobol Sensitivity Analysis

**Status:** Committed

## Method
- 5-parameter reduced-ODE Sobol analysis
- Real MU-Glioma parameter distributions (rho, D_white)
- N=500 Saltelli samples → 6000 evaluations

## Key finding
- rho_s dominates TTP variance: S1 = 0.998 (95% CI ±0.086)
- All other parameters ST < 0.015
- Interpretation: proliferation rate alone determines treatment outcome

## Files
- src/46_sensitivity_analysis.py
- output/sobol_sensitivity_results.json
- output/sobol_tornado_plot.png