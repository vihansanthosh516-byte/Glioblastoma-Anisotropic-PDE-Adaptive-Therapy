# Claim ledger (masterplan §95)

One block per claim. Numbers come from [[PAPER_FINDINGS_LEDGER]]; this file adds the wording rules.
Levels: 1 Demonstrated, 2 Supported, 3 Hypothesised, 4 Future use (never in Results).
Created 2026-10-04. Update when a result changes.

## CL-B2 Forecast beats no-change (MU core)
- Evidence: Dice 0.259 vs 0.233, diff +0.026 (CI 0.008–0.044), Holm p 0.0015, n=133. `output/forecast_labels_core/results.json`
- Dataset / unit: MU-Glioma-Post / patient. Level 1 in development cohort only.
- External: not tested on masks. LUMIERE volumes show the opposite (growth model worse; ledger 2l).
- Known limit: median diff +0.005; loses on shrinking cores (−0.039, n=50); target chosen after a volume probe; absolute Dice low.
- Allowed: "Small mean improvement in next-scan overlap in the development cohort."
- Forbidden: "Clinically improves tumour forecasting." "Validated."

## CL-B3 Anisotropy elongates simulated tumours
- Evidence: 1.21 vs 1.03, 98/103, p 4.1e-19. `output/fractal_aniso_vs_iso.json`
- Level 1 for simulation vs simulation only.
- Limit: not tested against real tumour shape; tract alignment is lower under aniso (0.45 vs 0.50).
- Allowed: "In simulation, anisotropic diffusion elongates the tumour."
- Forbidden: "Anisotropy explains real tumour shape."

## CL-B1 Adaptive therapy selects less resistance than MTD (simulation)
- Evidence: MTD ends with more resistant cells in 54/61, 0 reverse; adaptive ends with higher final mass (VR 2.39). `output/adaptive_cohort_summary.json`, `output/emax_sweep_summary.json`
- Level 1 for the simulator; Level 3 for biology.
- Limit: kill scale and resistance are assumed (E_MAX = rho×1000); robust only for ratio 500–2000; Track A→B link is null.
- Allowed: "In this model, adaptive dosing keeps resistant fraction lower at the cost of larger final tumour mass."
- Forbidden: "Adaptive therapy works in GBM."

## CL-C1 Kill-rate-conditioned RL wins at equal budget
- Evidence: 100% win on day-90 volume (`output/rl_kill_conditioned/evaluate.json`); never wins TTP; a one-line rule matches it.
- Level 1 as simulation; adds nothing over the rule.
- Limit: `real_test` uses the same assumed kill rates for all 21 patients; oracle information.
- Allowed: "With known kill rates, RL matches a simple efficiency rule."
- Forbidden: "RL improves treatment." "Real-drug test."

## CL-C3 Pacing extends TTP
- Evidence: +48.5 / +52.5 / +2.1 / +15.7 d. `output/probe_paced_policies/evaluate.json`
- Limit: real-parameter set only +2.1 d; day-90 volume loses.
- Allowed: "Pacing lengthens simulated time to progression; the gain is small with real growth rates."

## CL-A1 C-GAT zone accuracy (UNRELIABLE)
- Evidence: 78.7% vs scVI 73.0% under random cell split. Patient-ID lookup 82.8%. Patient-level LR 60.3%, RF 62.0%, scVI-LR 66.2% (script 86). Leak-free C-GAT: `output/cgat_leak_free.json` (script 88).
- Level: not demonstrated. Leak-free C-GAT (script 88, 5 folds) is 63.8% (SD 13.3). A cells-per-patient-only classifier scores 64.6% (SD 27.5) under the same folds (script 95). Healthy class = 3 donors with other diagnoses (L33).
- Allowed: "Under patient-level folds, a leak-free graph model reached 63.8% (SD 13.3), no better than non-biological baselines within fold spread."
- Forbidden: "C-GAT classifies tumour zones at 78.7%."

## CL-A3 / CL-CGGA Age and survival
- Evidence: TCGA age HR 1.031 (n=518); CGGA pooled 1.009 (1.001–1.018), Holm p 0.082. `output/survival_stats_summary.json`, `output/cgga_validation.json`
- Level 2 for direction; effect size weaker externally.
- Allowed: "Age direction replicated; effect was weaker and not significant after Holm in each batch."
- Forbidden: "Age is a validated predictor in our pipeline."

## CL-LUM LUMIERE volume forecast
- Evidence: growth model error 3.28 vs no-change 1.89 (n=40, p 0.0003); 1.59 vs 1.05 (n=50, p 0.004). `output/lumiere_volume_forecast.json`
- Level 1 (negative), exploratory (seen before plan).
- Limit: volumes only; noisy (median no-change error 0.8–1.8 log units).
- Allowed: "A growth-only model did not beat no-change on LUMIERE volumes."

## CL-TRK Track A→B link
- Evidence: inflammation score 1.0 for all 61; zone files likely synthetic (L32).
- Allowed: "No link was established."
- Forbidden: "Single-cell data inform the PDE."
