> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Script 46 — Sobol Sensitivity Analysis

> **SUPERSEDED (2026-10-03):** the S1=0.998 result came from a reduced ODE, not the PDE. Real PDE Sobol: α_sens S1=0.43, ρ_s S1=0.20, D_w S1≈0.0000062 (legacy; small nonzero). See [[Script-60-66-Swanson-D]].

**Status:** Committed

## Method
- 5-parameter reduced-ODE Sobol analysis
- Historical claim of MU-fitted sampling distributions needs an exact run/code manifest; these are not observed treatment effects.
- Historical N=500 /6000-evaluation assertion needs exact sampling metadata; do not transfer it to the distinct PDE comparison.

## Key finding
- Real PDE Sobol: α_sens S1 = 0.43, ρ_s S1 = 0.20, D_w S1≈0.0000062 legacy, 0.0000749 alternate regime
  (95% CI: see [[Script-60-66-Swanson-D]])
- Interpretation: alpha_sens has a larger point-estimate S1 in this assumed model, with wide overlapping uncertainty (rough confidence half-widths alpha .14, rho .18). No stable causal/clinical ranking or matched literature calibration is established. Source `output/sobol_pde_swanson_vs_legacy.json:legacy.V_total`, `swanson.V_total`; original ODE .998390 remains superseded.

## Superseded (reduced-ODE, kept for history)
- Original archived claim: rho_s dominated TTP variance with S1 = 0.998.
  This came from a reduced ODE with arbitrary coupling, not the PDE. 
  The archived value is retracted.
## Files
- src/46_sensitivity_analysis.py
- output/sobol_sensitivity_results.json
- output/sobol_tornado_plot.png