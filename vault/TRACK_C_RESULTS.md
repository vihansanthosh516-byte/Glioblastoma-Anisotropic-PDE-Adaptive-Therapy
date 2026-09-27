# Track C — Verified Results

**Status:** Complete (commit e638c4f)

## Overview
Digital twin pipeline: inverse estimation → MPC → RL adaptive steering → 
sensitivity → endpoint comparison.

## Script 51 — Inverse Estimation
- 154 patients refit on all scans
- 124/154 sit at model's lower growth bound (observed growth below min)
- 24 fit properly: median ρ = 0.0116/day, agrees with script 72 (Spearman 0.91)
- D cannot be identified from volumes alone
- Source: output/inverse_est_metrics.json

## Script 52 — Robust MPC
- 68.9% dose sparing (equal budget)
- Caveat: robust MPC variance reduction is -31.8% (worse than standard)
- Source: output/robust_mpc_benchmark.json

## Script 59 — Drug-Budget Artifact
- 100% → 10% win rate at equal budget
- Greedy heuristic front-loads by day 35, leaves 55 days untreated
- Source: output/phase6_sensitivity_metrics.json

## Script 60 — Ablation
- DTI, mechanics, diffusion all <0.01% effect on outcome
- Model is per-voxel logistic growth with no spatial structure
- PPO 100% win; DAgger reaches oracle ceiling; heuristic 43%
- Source: output/ablation_and_baselines_metrics.json

## Script 62 — Biomarker Stability
- ρ* ≈ 0.075/day crossover
- Early start beats Stupp for 60/64 patients
- Bootstrap 95% CI: [0.0746, 0.0820]
- Source: output/biomarker_stability_metrics.json

## Script 64 — Virtual Cohort
- Unconditioned PPO loses: 25% win rate
- Cause: varying kill rates across patients
- Source: output/phase8_cohort_metrics.json

## Script 66 — RL Equal-Budget
- PPO 100% on synthetic LHS, 25% on real patients
- Source: output/rl_equal_budget/evaluate.json

## Script 67 — Conditioned Policy
- Conditioning on kill rates flips 25% → 100% win
- Simple rule matches oracle to 0.2%
- Policies without kill rates fail to infer them
- Source: output/rl_kill_conditioned/evaluate.json

## Script 68 — Time-to-Progression Endpoint
- Day-90 optimal schedule LOSES on TTP (-29.4 days)
- TTP oracle gains: +78d (cohort64), +50d (lhs60), +3.9d (real_test), +24.9d (synth)
- CEM beat the oracle on 2/9 checks — oracle is a lower bound
- Source: output/ttp_equal_budget/evaluate.json

## Synthesis
Track C's central finding: patient-specific drug response, not policy 
complexity, determines outcome. Drug-budget-matched evaluation and 
patient-conditioned rules are essential for valid RL in adaptive 
therapy.

## Limitations
- Model lacks clonal competition between sensitive/resistant populations
- D unidentifiable from volume data alone
- 124/154 patients fit at the model's lower bound
- TTP oracle is a lower bound, not the true optimum
- Script 52's robust MPC variance regression not fully diagnosed

## Related Vault Notes
- [[Script-51-Inverse-Estimation]]
- [[Script-59-Sensitivity-RL]]
- [[Script-66-RL-Equal-Budget]]
- [[Script-67-Conditioned-Policy]]
- [[Script-68-TTP-Endpoint]]
- [[Repo]]