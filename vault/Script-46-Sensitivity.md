# Script 46 — Sobol Sensitivity Analysis

> **SUPERSEDED (2026-10-03):** the S1=0.998 result came from a reduced ODE, not the PDE. Real PDE Sobol: α_sens S1=0.43, ρ_s S1=0.20, D_w S1=0.000. See [[Script-60-66-Swanson-D]].

**Status:** Committed

## Method
- 5-parameter reduced-ODE Sobol analysis
- Real MU-Glioma parameter distributions (rho, D_white)
- N=500 Saltelli samples → 6000 evaluations

## Key finding
- Real PDE Sobol: α_sens S1 = 0.43, ρ_s S1 = 0.20, D_w S1 = 0.000
  (95% CI: see [[Script-60-66-Swanson-D]])
- Interpretation: drug sensitivity (α_sens) is the dominant driver of 
  tumor volume at literature diffusivity

## Superseded (reduced-ODE, kept for history)
- Original archived claim: rho_s dominated TTP variance with S1 = 0.998.
  This came from a reduced ODE with arbitrary coupling, not the PDE. 
  The archived value is retracted.
## Files
- src/46_sensitivity_analysis.py
- output/sobol_sensitivity_results.json
- output/sobol_tornado_plot.png