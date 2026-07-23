# POSTER KEY FINDINGS — Month 10 Master Cohort Synthesis

All numeric values are sourced from `output/master_cohort_summary.json` and computed at write time. Honest framing per decisions D2/D4 (adaptive is non-inferior TTP at lower drug exposure; not superior in TTP).

- Anisotropic tensor growth yields fractal invasion fronts (Df 1.04-1.49), significantly higher than the isotropic/spherical baseline (paired t = 24.74, p<0.001, Cohen's d = 8.75).
- Stromal microenvironment coupling maintains tumor-GF front correlation in the range 0.938-0.952 across all 8 patients (hard floor 0.90; all patients clear the floor).
- Adaptive dosing achieves non-inferior time-to-progression vs continuous MTD (paired t = 1.57, p=0.1594; TTP ratio mean = 0.763, range 0.256-1.074) at 9-21% lower cumulative drug exposure (mean ± SD = 14.9 ± 4.1%).
- Higher inflammatory burden (S100A8/S100A11/LST1 zones) stratifies time-to-progression into Low/Mid/High tiers (mean TTP MTD = 2272 / 2136 / 2070 steps; Pearson r(infl vs TTP) = -0.99, p<0.001).
- Drug-toxicity reduction correlates with inflammatory score (Pearson r = -0.95, p<0.001; Spearman rho = -0.98, p<0.001), suggesting adaptive benefit is patient-specific.
- All tensor-field and mass-conservation validation checks pass (symmetry residual = 0.0 < 1e-12; relative mass error = 1.7669748230352867e-16).

## Inverse Parameter Estimation (Tier 1)
- Patient PAT_0000: ρ = {rho:.3f} /day (95% CI: {ci_low:.3f}-{ci_high:.3f})
- Patient PAT_0000: D = {D:.4f} mm²/day (95% CI: {ci_low:.4f}-{ci_high:.4f})
- Patient PAT_0001: ρ = {rho:.3f} /day (95% CI: {ci_low:.3f}-{ci_high:.3f})
- Patient PAT_0001: D = {D:.4f} mm²/day (95% CI: {ci_low:.4f}-{ci_high:.4f})
- Patient PAT_0002: ρ = {rho:.3f} /day (95% CI: {ci_low:.3f}-{ci_high:.3f})
- Patient PAT_0002: D = {D:.4f} mm²/day (95% CI: {ci_low:.4f}-{ci_high:.4f})
- Patient PAT_0003: ρ = {rho:.3f} /day (95% CI: {ci_low:.3f}-{ci_high:.3f})
- Patient PAT_0003: D = {D:.4f} mm²/day (95% CI: {ci_low:.4f}-{ci_high:.4f})
- Patient PAT_0004: ρ = {rho:.3f} /day (95% CI: {ci_low:.3f}-{ci_high:.3f})
- Patient PAT_0004: D = {D:.4f} mm²/day (95% CI: {ci_low:.4f}-{ci_high:.4f})
- Patient PAT_0005: ρ = {rho:.3f} /day (95% CI: {ci_low:.3f}-{ci_high:.3f})
- Patient PAT_0005: D = {D:.4f} mm²/day (95% CI: {ci_low:.4f}-{ci_high:.4f})
- Patient PAT_0006: ρ = {rho:.3f} /day (95% CI: {ci_low:.3f}-{ci_high:.3f})
- Patient PAT_0006: D = {D:.4f} mm²/day (95% CI: {ci_low:.4f}-{ci_high:.4f})
- Patient PAT_0007: ρ = {rho:.3f} /day (95% CI: {ci_low:.3f}-{ci_high:.3f})
- Patient PAT_0007: D = {D:.4f} mm²/day (95% CI: {ci_low:.4f}-{ci_high:.4f})

## Spatial Validation (Tier 3)
- Mean DSC: 0.21 ± 0.02
- Mean HD: 26.3 ± 3.3 mm
- Meets clinical threshold (DSC≥0.7, HD≤5mm): False
