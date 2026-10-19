> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Scripts 67-68 — Kill-conditioned policy and the TTP endpoint

**Sources:** `output/rl_kill_conditioned/evaluate.json` (script 67), `output/ttp_equal_budget/evaluate.json` (script 68). Values read 2026-10-04.

All rates below are simulations with assumed kill effects, including real_test fitted-growth cases. No clinical treatment response is measured.

## Script 67 — day-90 win rate vs Stupp, equal budget, full 64^3 PDE
| Set | blind `ppo66` | `ppo_cond_ft` | `dagger_cond` | `efficiency_rule` | `oracle` |
|---|---|---|---|---|---|
| cohort64 (n=20) | **25.0%** | 100% | 100% | 100% | 100% |
| lhs60 | 100% | 100% | 100% | 100% | 100% |
| real_test (n=21) | 95.2% | 100% | 100% | 100% | 100% |
| synth_test | 47.5% | 100% | 97.5% | 100% | 100% |

- On cohort64, the conditioned policy gets mean log(final/oracle) of 0.0013. The simple `efficiency_rule` gets0.0016. Exact key `sets.cohort64.arms.{ppo_cond_ft,efficiency_rule}.mean_log_ratio_final_vs_oracle`:0.00127370/0.00158714; the audit incorrectly treated these fields as absent. **A one-line rule matches the RL policy.** The paper must say RL adds nothing over the rule here.
- Blind PPO is already 95-100% on lhs60 and real_test. The 25% failure is cohort64-specific (that set varies kill rates across patients).

## The caveat the paper must carry
- These wins are on **day-90 volume**. Script 68 shows the day-90-optimal schedule never beats Stupp on time to progression (TTP): 0.0 d (cohort64), **-29.4 d (lhs60)**, -10.7 d (real_test), -7.9 d (synth_test). Read from `ttp_equal_budget/evaluate.json` 2026-10-04.
- So "100% win" is a statement about one endpoint that the vault itself calls a poor metric. The optimal schedule leaves the tumour untreated for about 55 days.

## Script 68 — TTP oracle gains (from the vault)
+78 d (cohort64), +50 d (lhs60), +3.9 d (real_test), +24.9 d (synth). The historical CEM2-of9/lower-bound assertion is unverified without the nine-check artifact and exact search definition; no rigorous bound is established.
