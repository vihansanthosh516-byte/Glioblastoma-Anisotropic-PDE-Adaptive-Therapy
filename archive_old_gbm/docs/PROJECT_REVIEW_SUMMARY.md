> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Historical GBM model review: corrected evidence summary

Revised 2026-10-19. Former generated date 2026-07-25 and “Phase 6 complete”
labels do not describe the later edits or current project. The old framework
explored in-silico policy optimization under assumed treatment effects. It did
not establish superiority to standard care, a validated patient digital twin,
clinical robustness or readiness for regulatory translation. Current goal:
`../../START_HERE.md` (probabilistic RANO).

## Traceable historical simulations

Paths below refer to root `../../output/`. They describe different runs and
must not be pooled or substituted for one another.

| Run / keys | Existing values | Meaning and limits |
|---|---|---|
| `phase5_adaptive_metrics.json`: `rl_final_volume_mm3`, `stupp_final_volume_mm3`, `rl_toxicity_reduction_pct`, `rl_recurrence_delay_days` | 0.2471128 vs 11.0059328 mm3; toxicity reduction 0%; recurrence delay 0 days | Single assumed-model artifact. Neither 1.04 mm3/10.6x clearance nor >90-vs-42-day TTP is supported. Final volume does not establish TTP benefit. |
| `phase6_sensitivity_metrics.json`: `total_scenarios_evaluated`, `parameter_correlations_with_rl_volume.*.pearson_r` | 30 scenarios; rho 0.9547642, D_w 0.1466039, alpha_sens -0.0867225 | Exploratory model sensitivity, not a measured patient biomarker or causal clinical result. |
| Same file: `rl_win_rate_pct`, `mean_rl_volume_mm3`, `mean_stupp_volume_mm3`, `drug_auc_ratio` | 10.0%; 15.083309 vs 14.468607 mm3; 0.9921429 | Does not support 36.7% wins, >80% above a rho cutoff, 13% less drug, or general RL superiority. |
| `biomarker_stability_metrics.json`: `model_threshold.rho_star_per_day`, `bootstrap.interpolated.ci95` | rho*=0.0751/day; CI [0.0746,0.0820] (rounded) | Distinct script-62 model threshold, only four of 64 above it (ledger 2). Obsolete 0.024/day rule is withdrawn. No mapping to Ki-67/PET or clinical policy is established. |

Later equal-budget policy comparisons are indexed by
`../../vault/PAPER_FINDINGS_LEDGER.md` C1 and open-conflict item 4: PPO 12.67
versus Stupp 14.47 mm3 describes its named script-66 simulated run. It does not
replace every phase-5/6 value. The later 25% cohort64 win rate is simulated;
20/21 on fitted real_test growth rates still uses assumed kill rates, measures
day-90 volume and never wins on TTP. A simple rule matches the policy. No
independently measured drug response is demonstrated.

## Numerical model scope

The old sections described registration, latent dynamics, reaction-diffusion,
mechanics, tensor construction and RL environments. Implementation descriptions
are not evidence that every module was integrated or calibrated. The following
specific numerical assertions are withdrawn pending exact source artifacts:

- Phase 1 SyN/GP registration and velocity summaries, NEB energy 5.74, dual
  attractors and 32-dimensional flux require code/run identifiers and saved results.
- Patient-specific multi-omic/DTI calibration and retrospective eight-patient
  calibration require a cohort manifest, input lineage and fitted-parameter output.
- Phase 3 128-cubed grid, 500-step timing, DSC>=0.70 and HD95 recurrence claims
  require masks, comparison protocol, hardware and metric artifacts.
- Phase 4 N=500 virtual-trial, non-inferiority, correlation and p-value claims
  require a named output, assumptions, endpoints and registered margins.
- RL observation/action/reward details, training/evaluation seconds, “ultra-fast”
  runtimes and batch times need exact code/configuration and benchmark logs.
  Conflicting 10/50/595-second anecdotes are not a benchmark.
- Untreated/Stupp/RL cellularity, CV, BED, TTP and exposure tables lack a single
  coherent run source. They are removed rather than mixed with other experiments.

Ledger 2e identifies script-48 legacy parallel/perpendicular diffusivities
0.013/0.0013 versus alternate-regime 0.13/0.013 mm2/day. These are selected model
inputs, not established patient-specific calibration. Ledger B3 records the
fractal comparison moving in the unfavorable direction; it does not validate
anisotropic biology. The former 6.7% clinical-necrosis calibration and tract
correlations lack matching calibration evidence here and are withdrawn.
Ledger 2h (`invasion_summary.json`) records a simulated free speed 26.592883
versus analytic 28.284271 um/hr: 5.98% relative error. This tests selected model
arithmetic, not clinical wave-speed accuracy. Ledger 2q
(`solver_verification.json`) covers synthetic numerical checks, clamp artifacts
and resolution limits; passing them does not validate biology.

## Reproducibility and translation limits

No repository-wide “all seed 42,” complete JSON/PNG coverage, deterministic
provenance certification or end-to-end success is established without a manifest
and command logs. `../../LICENSE` contains MIT text for the code; external inputs
have separate licenses. No clean install or full pipeline was run for this review.

Historical proposals for geometry, heterogeneity, toxicity, large-scale computing
and external clinical validation remain untested work proposals. Patient-specific
calibration, transportable forecasts, independent therapy outcomes and biomarker
measurements are missing. The old “strong/novel” heading, clinical recommendations,
biophysical-fidelity claim and translational conclusion are withdrawn. Literature
novelty remains unverified. The old GBM program is background; current state and
next steps are in START_HERE. No new analysis or spent-set reuse is authorized.
