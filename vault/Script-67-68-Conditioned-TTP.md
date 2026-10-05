# Scripts 67-68 — Kill-conditioned policy and the TTP endpoint

**Sources:** `output/rl_kill_conditioned/evaluate.json` (script 67), `output/ttp_equal_budget/evaluate.json` (script 68). Values read 2026-10-04.

## Script 67 — day-90 win rate vs Stupp, equal budget, full 64^3 PDE
| Set | blind `ppo66` | `ppo_cond_ft` | `dagger_cond` | `efficiency_rule` | `oracle` |
|---|---|---|---|---|---|
| cohort64 (n=20) | **25.0%** | 100% | 100% | 100% | 100% |
| lhs60 | 100% | 100% | 100% | 100% | 100% |
| real_test (n=21) | 95.2% | 100% | 100% | 100% | 100% |
| synth_test | 47.5% | 100% | 97.5% | 100% | 100% |

- On cohort64, the conditioned policy gets mean log(final/oracle) of 0.0013. The simple `efficiency_rule` gets 0.0016. **A one-line rule matches the RL policy.** The paper must say RL adds nothing over the rule here.
- Blind PPO is already 95-100% on lhs60 and real_test. The 25% failure is cohort64-specific (that set varies kill rates across patients).

## The caveat the paper must carry
- These wins are on **day-90 volume**. Script 68 shows the day-90-optimal schedule **loses on time to progression (TTP): -29.4 days** (`TRACK_C_RESULTS.md`; per-set TTP numbers in `ttp_equal_budget/evaluate.json` not re-read).
- So "100% win" is a statement about one endpoint that the vault itself calls a poor metric. The optimal schedule leaves the tumour untreated for about 55 days.

## Script 68 — TTP oracle gains (from the vault)
+78 d (cohort64), +50 d (lhs60), +3.9 d (real_test), +24.9 d (synth). CEM beat the oracle on 2 of 9 checks, so the oracle is a lower bound.
