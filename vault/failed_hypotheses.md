# Failed and open hypotheses log (masterplan §195)

Created 2026-10-04. Never delete a row. Add rows when a result comes in. Numbers: [[PAPER_FINDINGS_LEDGER]].

| Hypothesis | Result | Status | Why / evidence |
|---|---|---|---|
| DTI improves forecast | −0.001 Dice, p 0.97 | Rejected (development only) | `forecast_labels_core/results.json`; external not run |
| Growth-only PDE beats no-change on shrinking cores | −0.039 Dice, n=50 | Rejected | same file |
| Growth-only model beats no-change on LUMIERE volumes | Worse (n=40 p 0.0003; n=50 p 0.004) | Rejected | `lumiere_volume_forecast.json` |
| Age is a robust external predictor | HR 1.009, Holm p 0.082 (CGGA pooled) vs 1.031 TCGA | Weakened | `cgga_validation.json` |
| MGMT predicts TTP | HR 1.11, p 0.60 | Not shown; underpowered (needs ~2,807 events) | `limit_checks.json` |
| C-GAT beats scVI at zone classification | Random split only; patient-ID lookup 82.8% | Unreliable | script 86; script 88 pending |
| RL > simple rule | RL 0.2% of oracle vs rule | Rejected | script 68 |
| Pacing gain on real parameters | +2.1 d | Mostly explained (shared assumed kill rates) | script 87 |
| Resistance-driven adaptive gain on real patients | 0 of 6 cells (script 83) | Rejected | `plastic_resistance/results.json` |
| Dual-agent "rescue" is adaptivity | MTD + same drug 358.9 d vs dual-adaptive 305.6 d | Rejected | `dual_mtd_baseline.json` |
| Fractal D_f tracks anisotropy | wrong direction | Rejected | `fractal_aniso_vs_iso.json` |
| Track A inflammation score informs Track B | score = 1.0 for all 61 | Rejected (null link) | `adaptive_cohort_summary.json` |
| rho > 0.02 /day threshold | 10/10 vs 0/51, p 5e-5 | Exploratory (post hoc, n=10) | `rho_threshold_check.json` |
| Anisotropy elongates simulated tumour | +0.178, p 4e-19 | Supported (simulation only) | `fractal_aniso_vs_iso.json` |
| Adaptive keeps resistant fraction lower | 54/61; robust at ratios 500–2000 | Supported (simulation, assumed kill scale) | `emax_sweep_summary.json` |
| Molecular programs transfer to PDE parameters | not tested | Open | Phase 6 |
| Selector (mechanistic vs persistence) helps | not tested | Open | Phase 4 |
| Sequential updating beats fixed fit | not tested | Open | Phase 4 |
| Robust MPC ≥ RL under uncertainty | not tested | Open | Phase 5 |
| H1 PDE beats persistence (broadly positive) | GBM-only +0.0130 (CI -0.0043 to +0.0302); all-eligible +0.0088 (-0.0061 to +0.0235) | Not supported | `output/pde_manifest/results.json` |
| PDE beats the no-PDE geometric rule | GBM-only -0.0019 (-0.0176 to +0.0141) | Not supported | `output/pde_manifest/results.json` |
| H3 tract anisotropy improves the forecast | old solver: aniso r=1 vs iso_same +0.0002 (-0.0002 to +0.0005); exact solver: -0.0004 (-0.0009 to -0.0001), r=10 -0.0030 (-0.0050 to -0.0009) | Rejected: the small old edge was a clamp artefact; with atlas tensors, direction slightly hurts (ledger 2ak) | `output/pde_manifest/results.json`, `output/pde_manifest_monotone/results.json` |
