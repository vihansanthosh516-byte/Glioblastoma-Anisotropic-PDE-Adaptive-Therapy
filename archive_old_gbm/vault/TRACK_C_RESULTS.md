> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Track C — Historical model results

**Status:** Historical completion account (commit e638c4f), corrected2026-10-19; not blanket verification

## Overview
Collection of historical fitting/MPC/RL/sensitivity analyses; integration and cross-run comparability are not established by listing modules. All therapy outcomes use assumed effects, including fitted-growth real_test inputs; no validated patient digital twin.

## Script 51 — Inverse Estimation
- 154 patients refit on all scans
- 124/154 sit at model's lower growth bound (observed growth below min)
- 24 have interior rho estimates, not proven “proper” biological fits: median0.0115561/day, script72 cross-check Spearman0.914783 on24 (`inverse_est_metrics.json:summary`). D hits a bound in147/154; no valid parameter-accuracy claim follows.
- D cannot be identified from volumes alone
- Source: output/inverse_est_metrics.json

## Script 52 — Robust MPC
- 68.9% dose sparing (equal budget)
- Saved variance fields conflict: variance_reduction_pct=-31.7562 versus cost_variance_reduction_pct=+0.27573; robust cost_variance_mean0.631810 vs standard0.633557. These are distinct recorded summaries; no blanket variance-improvement inference.
- Source: output/robust_mpc_benchmark.json

## Script 59 — Drug-Budget Artifact
- 100% → 10% win rate at equal budget
- Historical day35/55-day schedule account needs exact actions/configuration; phase6 JSON supports10% RL wins but is not schedule proof.
- Source: output/phase6_sensitivity_metrics.json

## Script 60 — Original low-diffusivity ablation

Numeric comparisons are selected assumed-model endpoints; later Swanson/spatial grid results differ (Script60-66/78 notes). A small effect on total volume does not establish general DTI/shape irrelevance.
- DTI, mechanics, diffusion all <0.01% effect on outcome
- Model is per-voxel logistic growth with no spatial structure
- PPO 100% win; DAgger reaches oracle ceiling; heuristic 43%
- Source: output/ablation_and_baselines_metrics.json

## Script 62 — Biomarker Stability
- Modeled rho*=0.0751/day crossover, not a validated biomarker (only4/64 above it)
- Early start beats Stupp for 60/64 patients
- Bootstrap 95% CI: [0.0746, 0.0820]
- Source: output/biomarker_stability_metrics.json

## Script 64 — Virtual Cohort
- Unconditioned PPO loses: 25% win rate
- Model dependence on varying assumed kill rates is an interpretation; not measured causal patient response.
- Source: output/phase8_cohort_metrics.json

## Script 66 — RL Equal-Budget
- Legacy LHS30 final-volume PPO100%; separate cohort64 blindPPO25% is synthetic, not “real patients.” Fitted-growth real_test PPO20/21 uses shared assumed kill rates. Exact sets/arms differ (source66/67 notes).
- Source: output/rl_equal_budget/evaluate.json

## Script 67 — Conditioned Policy
- Oracle conditioning changes cohort64 simulated win rate25%→100%, on day90 volume; no measured patient efficacy
- Simple rule matches oracle to 0.2%
- Policies without kill rates fail to infer them
- Source: output/rl_kill_conditioned/evaluate.json

## Script 68 — Time-to-Progression Endpoint
- Day-90 optimal schedule vs Stupp on TTP: 0.0 d (cohort64), **-29.4 d (lhs60)**, -10.7 d (real_test), -7.9 d (synth_test). It never wins (corrected 2026-10-04 from `ttp_equal_budget/evaluate.json`; earlier note quoted only the lhs60 value).
- **`real_test` uses real growth rates but the same assumed kill rates for all 21 patients** (script 87; [[Script-75-Probe-Paced]]).
- TTP oracle gains: +78d (cohort64), +50d (lhs60), +3.9d (real_test), +24.9d (synth)
- CEM2/9/lower-bound assertion unverified without exact check artifact and optimization definition; no rigorous bound is established.
- Source: output/ttp_equal_budget/evaluate.json

## Synthesis
In these simulations, assumed case-specific kill effects change policy outcomes, and a simple efficiency rule matches RL. This does not establish what determines clinical response or what every valid adaptive RL system requires.

## Limitations
- Model lacks clonal competition between sensitive/resistant populations
- D unidentifiable from volume data alone
- 124/154 patients fit at the model's lower bound
- TTP search optimum/bound status is not established by an unlinked CEM anecdote.
- Script 52's robust MPC variance regression not fully diagnosed

## Related Vault Notes
- [[Script-51-Inverse-Estimation]]
- [[Script-59-Sensitivity-RL]]
- [[Script-66-RL-Equal-Budget]]
- [[Script-67-Conditioned-Policy]]
- [[Script-68-TTP-Endpoint]]
- [[Repo]]
Artifact check2026-10-19: all named inverse/MPC/sensitivity/ablation/biomarker/cohort/66–68 JSONs exist in root output, contrary to the audit missing-file assertions. Inverse summary keys confirm154 refits,124 rho lower-bound,24 interior rho and Spearman; MPC and phase6 keys were read directly. Preservation of a result is not independent clinical validation.
