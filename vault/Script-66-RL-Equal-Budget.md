# Script 66 — RL at Equal Drug Budget

**Status:** Committed (95a73a1, a5cbb6a, 697c470)

## Key Results

| Arm | Win % | Mean final volume | Drug AUC |
|---|---|---|---|
| Stupp | — | 14.47 mm³ | 35.0 |
| Old heuristic (script 59) | 10% | 15.08 mm³ | 34.7 |
| PPO (3 seeds) | 100% | 12.67 mm³ | 35.0 |
| DAgger from oracle | 100% | 11.43 mm³ | 35.0 |
| Legacy-reward PPO | 0% | 14.89 mm³ | 34.9 |

Held-out real patients (n=21): PPO wins 95%, DAgger wins 100%.

## Four Findings

1. **RL beats Stupp at equal budget** — 100% win rate
2. **Winning policies are fixed schedules, not adaptive** — PPO converges to "combo from day 1", DAgger to "combo days 1-35"
3. **The model has no adaptive advantage** — single homogeneous population, no clonal competition
4. **Day-90 volume is a poor metric** — the "optimal" schedule leaves tumor untreated for 55 days

## Model issues found

- `alpha_sens` is stored but never used in `pde_step`
- Real per-patient D is derived as ρ/10, not independently fit

## Literature

Adaptive therapy requires clonal competition (Gatenby 2009, Leder 2014). 
Our model lacks this structure, so adaptive control has nothing to exploit.

## Files

- src/66_rl_equal_budget_study.py
- output/rl_equal_budget/ (JSON + PNG)