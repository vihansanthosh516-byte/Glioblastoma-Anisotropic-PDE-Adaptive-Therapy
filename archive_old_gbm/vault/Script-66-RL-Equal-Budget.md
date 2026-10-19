> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Script 66 — RL at Equal Drug Budget

**Status:** Committed (95a73a1, a5cbb6a, 697c470)

## Original legacy-run results (not the Swanson rerun)

Source `output/rl_equal_budget/evaluate.json:sets.lhs30`: stupp,heuristic59,ppo_final_s0/s1/s2,bc_dagger_final; fields mean_final_mm3,win_rate_final_pct,mean_auc. This assumed-treatment LHS30 simulation is distinct from cohort64 and the later Swanson rerun (separate directory).

| Arm | Win % | Mean final volume | Drug AUC |
|---|---|---|---|
| Stupp | — | 14.47 mm³ | 35.0 |
| Old heuristic (script 59) | 10% | 15.08 mm³ | 34.7 |
| PPO (3 seeds) | 100% | 12.67 mm³ | 35.0 |
| DAgger from oracle | 100% | 11.43 mm³ | 35.0 |
| Legacy-reward PPO | 0% | 14.89 mm³ | 34.9 |

Simulated evaluation of21 fitted-growth real_test cases with assumed kill rates: PPO about95.2% on day-90 volume, DAgger100%. No measured drug efficacy or TTP advantage follows.

## Four Findings

1. **Original LHS30 simulation** — listed final-volume PPO arm wins100%; blind cohort64 only25% in the separate conditioned-policy evaluation
2. **Winning policies are fixed schedules, not adaptive** — PPO converges to "combo from day 1", DAgger to "combo days 1-35"
3. **The model has no adaptive advantage** — single homogeneous population, no clonal competition
4. **Endpoint tradeoff in this simulation** — the historical55-day untreated interpretation requires exact schedule/action output; final-volume optimization does not establish TTP benefit.

## Model issues found

- `alpha_sens` is stored but never used in `pde_step`
- Real per-patient D is derived as ρ/10, not independently fit

## Literature

Adaptive therapy requires clonal competition (Gatenby 2009, Leder 2014). 
This selected model lacks that competition structure; it does not establish that all adaptive policies lack benefit in other models or patients.

## Files

- src/66_rl_equal_budget_study.py
- output/rl_equal_budget/ (JSON + PNG)