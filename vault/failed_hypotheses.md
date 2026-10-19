# Failed and open hypotheses log (masterplan §195)

> **HISTORICAL BACKGROUND (corrected 2026-10-19):** Dated record of the older GBM platform, created 2026-10-04. It is not current project status. The current goal is probabilistic RANO (`START_HERE.md`). Rows with phase labels ("Phase 4-6") describe the old Track plan as of 2026-10-04, not current milestones. All treatment results are simulations with assumed parameters unless stated.

Created 2026-10-04. Never delete a row. Add rows when a result comes in. Numbers: [[PAPER_FINDINGS_LEDGER]].

| Hypothesis | Result | Status | Why / evidence |
|---|---|---|---|
| DTI improves forecast | −0.001 Dice, p 0.97 | Rejected (development only) | `forecast_labels_core/results.json`; external not run |
| Growth-only PDE beats no-change on shrinking cores | −0.039 Dice, n=50 (descriptive subgroup defined by the observed outcome: `descriptive_shrank_or_same`, CI −0.061 to −0.019) | Rejected | same file |
| Growth-only model beats no-change on LUMIERE volumes | Worse in volume error, not Dice (HD-GLIO-AUTO n=40 p 0.0003; DeepBraTumIA n=50 p 0.004) | Rejected | `lumiere_volume_forecast.json` |
| Age is a robust external predictor | HR 1.009, Holm p 0.082 (CGGA pooled) vs 1.031 TCGA; adjusted for IDH + treatment 1.011 (1.002-1.021, p 0.018; ledger 2ar; `cgga_adjusted_age_mgmt.json` `Y5_age_adjusted.terms.age.p` = 0.01787 checked, HR/CI taken from ledger 2ar) | Weakened: real but about 1/3 of the TCGA size | `cgga_validation.json` |
| MGMT predicts TTP | HR 1.11, p 0.60; CGGA OS (overall survival, not TTP) within TMZ 0.80 (0.60-1.05), MGMT x TMZ 0.60 (0.33-1.10) (ledger 2ar; these HR/CI were not re-read from JSON; the adjusted interaction p is 0.0995) | Not shown; direction matches literature; underpowered | `limit_checks.json` |
| C-GAT beats scVI at zone classification | Random split only; patient-ID lookup 82.8%. Script 88 is complete (ledger 2o): leak-free C-GAT 63.8% over five patient-level folds (`cgat_leak_free.json`) | Unreliable | script 86; script 88 |
| RL > simple rule | The one-line rule matches RL and costs 0.2% of oracle computation (ledger C1; script-68 output key not re-read) | Rejected | script 68 |
| Pacing gain with real fitted growth rates and common assumed kill rates | +2.1 d (simulated model result) | Mostly explained (shared assumed kill rates) | script 87 |
| Resistance-driven adaptive gain, simulation on the `real_test` fitted-growth cohort (assumed treatment parameters; no patient received these policies) | 0 of 6 tested simulation cells (script 83) | Rejected | `plastic_resistance/results.json` |
| Dual-agent "rescue" is adaptivity | 8-patient simulation, TTP capped at 360 d, assumed drug model: MTD + same drug 358.9 d vs dual-adaptive 305.6 d | Rejected | `dual_mtd_baseline.json` |
| Fractal D_f tracks anisotropy | wrong direction (simulation; primary fractal-dimension difference -0.018, `primary_fd_150d`) | Rejected | `fractal_aniso_vs_iso.json` |
| Track A inflammation score informs Track B | score = 1.0 for all 61 was a code default (0 ID overlap; zone files likely synthetic, unconfirmed) | Not tested (ledger 2ap) | `adaptive_cohort_summary.json`; src/44_adaptive_therapy.py:749-764 |
| rho > 0.02 /day threshold | 10/10 vs 0/51, p 5e-5 | Exploratory (post hoc, n=10) | `rho_threshold_check.json` |
| Anisotropy elongates simulated tumour | +0.178, p 4e-19 | Supported (simulation only) | `fractal_aniso_vs_iso.json` |
| Adaptive keeps resistant fraction lower | 54/61 at ratio 1000, 7 ties, 0 reverse; only ratios 500, 1000, 2000 tested (MTD more resistant 57/54/51) | Supported (simulation, assumed kill scale) | `emax_sweep_summary.json` |
| Molecular programs transfer to PDE parameters | not tested | Open | Phase 6 |
| Selector (mechanistic vs persistence) helps | not tested | Open | Phase 4 |
| Sequential updating beats fixed fit | not tested | Open | Phase 4 |
| Robust MPC ≥ RL under uncertainty | not tested | Open | Phase 5 |
| H1 PDE beats persistence (broadly positive) | GBM-only +0.0130 (CI -0.0043 to +0.0302); all-eligible +0.0088 (-0.0061 to +0.0235; arm/key not mapped to `results.json` in the 2026-10-19 check, where the all-pairs mean read +0.0106; unverified) | Not supported | `output/pde_manifest/results.json` |
| PDE beats the no-PDE geometric rule | GBM-only `aniso\|r=1` vs geometric_train_rate -0.0019 (-0.0176 to +0.0141; mean checked, CI not re-read) | Not supported | `output/pde_manifest/results.json` |
| H3 tract anisotropy improves the forecast | old solver: aniso r=1 vs iso_same +0.0002 (-0.0002 to +0.0005); exact solver: -0.0004 (-0.0009 to -0.0001), r=10 -0.0030 (-0.0050 to -0.0009) | Rejected: the small old edge was a clamp artefact; with atlas tensors r=1 is near zero and r=10 is negative (ledger 2ak) | `output/pde_manifest/results.json`, `output/pde_manifest_monotone/results.json` |
