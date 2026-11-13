# Track C: the five open decisions, with evidence (2026-11-12)

The user asked for the most ambitious choice that the science supports for each decision in `paper/trackC_unblock_plan_2026-11-12.md`. This record gives the choice, the evidence, and the cost. Nothing here changes a frozen rule; the 14-day target, the success rule, and the `observed_two_scan` stop rule stay. These choices become binding only when written into a dated amendment before any policy outcome.

## Summary

| # | Decision | Earlier recommendation | **Evidence-based choice** | Main reason |
|---|---|---|---|---|
| D1 | Controller model grid | 81 models | **135 models (5 q levels)** | With 81, the true signal quality is off the controller grid in half of all worlds |
| D2 | Forecast length for 112-day scans | 8 cycles | **8 cycles (224 days)** | The stop rule needs two consecutive scans; 6 cycles see only one |
| D3 | Patients per cohort | 6 | **Set by an outcome-free development pilot that minimizes total compute** | Required worlds scale with the variance of world differences, and cohort size controls part of that variance |
| D4 | Duke aggregate | No | **No** | Unresolved licence; spent cohort; a defined fallback already exists |
| D5 | Interval review | Machine review | **Machine review + reference-code cross-check + a wider outcome-free screen that also tests Gaffke's interval** | Adds an independent implementation check and a candidate that may need fewer worlds |

## D1: 135 models

Evidence on disk:
- The controller grid in code is `Q_GRID = (0.5, 0.7, 0.9)` (`src/212_trackc_signal_fisher_grid.py` line 30), giving 3 x 27 = 81 models.
- The simulator draws each world's true signal quality uniformly from **0.6, 0.7, 0.8, 0.9** (`configs/trackABC_integrated_v2.json`, `signal_model.primary_q_signal_sampling_levels`).
- So in worlds with q = 0.6 or 0.8 (half of all worlds by design), neither controller has the true q in its model set. Both arms are misspecified there, which blurs exactly the contrast Track C tests.
- The config's `synthetic_q_signal_grid` is 0.5, 0.6, 0.7, 0.8, 0.9 (5 x 27 = 135 models; `candidate_model_set_maximum` 135). It contains every true q level plus the 0.5 no-information level.

Choice: **135**. Cost: controller work grows about 135/81 = 1.67 times (estimate; to be timed in plan Step 6). Weakest point: a larger model set also makes minimax more conservative, which may shrink any gain; that is a real outcome, not a reason to keep the mismatch.

## D2: 8 forecast cycles for the 112-day sensitivity

- Cycle = 28 days; primary forecast = 6 cycles = 168 days (config `controller_comparison`).
- A 112-day scan interval is 4 cycles, so a 168-day forecast contains only one future scan.
- The frozen stop rule needs **two consecutive** visible readings at or above the threshold (config `stop_rule`). A forecast that sees one scan cannot evaluate when stopping would occur.
- 8 cycles (224 days) is the smallest whole-cycle forecast that contains two 112-day scans. `PolicyBenchmarkSettings` currently rejects `scan_every=4` with 6 cycles (`src/211_trackc_policy_benchmark.py` line 86), which is consistent with this.

Choice: **8 cycles**, for the 112-day sensitivity only; the primary (84-day) analysis keeps 6. Cost: about 8/6 = 1.33 times forecast work in that sensitivity arm only.

## D3: patients per cohort, chosen by a development pilot

Evidence:
- The GE-betting screen shows the needed world count scales with the square of the world-difference SD: about 350 worlds per family at SD 30 days, about 1,400 at SD 60, about 3,900 at SD 100 (`paper/trackc_ge_betting_coverage_audit_result_2026-11-09.md`, section 3; extrapolated by 1/sqrt(n), not tested cells).
- Each world's difference is an average over m virtual patients. Its variance is approximately `tau^2 + sigma^2 / m`, where tau^2 is between-world variance and sigma^2 is within-world patient variance (standard two-level variance decomposition).
- Runtime per world grows about linearly with m. Total compute is about `n(m) x m`. If sigma^2 dominates, larger m lowers n a lot; if tau^2 dominates, larger m wastes compute.
- Neither tau nor sigma is known. They can be estimated only from policy utilities.

Choice: run a **development pilot** on a separate, registered development seed block (the config already distinguishes development from evaluation worlds; `not_yet_locked_before_evaluation` item 2). About 20 development worlds with m = 12 patients, both arms. Estimate tau^2 and sigma^2 by a one-way random-effects decomposition. Then pick m from {2, 4, 6, 8, 12, 16} that minimizes predicted `n(m) x m x time(m)` subject to the 14-day criterion, and fix n. Development results are never used as evaluation data and are labelled development in every report.

Cost: about 20 world pairs at m = 12 (estimate 1 to 2 days on a DGX Spark at the timings below). Weakest point: 20 worlds give a noisy tau estimate; the rule therefore uses the upper 80% confidence bound of the variance, which errs toward more worlds.

## D4: no Duke aggregate

- The Duke licence text was never saved (gap B6, U1 in `paper/gap_status_2026-11-09.md`).
- Duke is an exposed, spent cohort (gap B1).
- The config already defines the fallback MRI baseline as a standard normal score (`input_rules.track_b.fallback_mri_baseline`).

Choice: **No.** Closes config item 7 with no loss: the MRI-like signal in Track C is synthetic in either case.

## D5: interval method review, cross-check, and a wider screen

Evidence:
- GE-betting is the "hedged capital" interval of Waudby-Smith and Ramdas, *Estimating means of bounded random variables by betting*, JRSSB 2024 ([arXiv 2010.09686](https://arxiv.org/abs/2010.09686)). The authors' Python implementation is the ConfSeq package, module `src/confseq/betting.py` ([GitHub](https://github.com/gostevehoward/confseq)).
- The project's GE-betting passed 30/30 synthetic coverage cells; its margin was thin in the maximum-variance cell (lowest Wilson lower bound 0.9499) (`paper/trackc_ge_betting_coverage_audit_result_2026-11-09.md`).
- Ming, Ramdas, Shen, Wang and Waudby-Smith (2026, arXiv 2607.18661) report that Gaffke's interval for bounded means is inadmissible but asymptotically efficient. Not read in full here; any use requires reading the paper and its validity conditions first.

Choice:
1. Cross-check `src/trackc_ge_betting.py` against ConfSeq on fixed synthetic inputs to a stated tolerance.
2. An independent fresh-context review of the code against the paper, written up.
3. Re-run the outcome-free coverage screen with the real paired-arm coupling (after C-8 is fixed) and add Gaffke's interval as a candidate if its validity for this estimand is confirmed in the paper.
4. Adopt, before any outcome, the valid method with the smallest predicted world count. Method choice uses synthetic data only.

Weakest point: no human statistician has reviewed the method. The report will say so.

## Cost at the chosen settings (estimates)

Base timing: 0.793 hours per world pair at 2 patients per cohort, one observation (`paper/gap_status_2026-11-09.md` C6). Assumed DGX Spark speed-up across 20 CPU cores: about 15 times (not measured). The Track C code is CPU-only NumPy.

| Setting | Factor | Hours per world pair on the Spark |
|---|---|---|
| m = 2, 81 models (base) | 1 | about 0.05 |
| m = 6, 135 models | 3 x 1.67 = 5 | about 0.26 |

At m = 6 and 135 models: 350 worlds per family is about 90 hours; 1,400 is about 15 days. The development pilot fixes the true numbers. A total above the runtime cap means the plan reports "not feasible at this budget" before any evaluation run, which the config's non-pass rule already allows.

## What these choices do and do not change

They remove scientific mismatches (D1, D2), replace a guess with a measured design choice (D3), close a licence risk (D4), and strengthen the interval method (D5). They do not guarantee a positive result. After these choices and the code steps in the unblock plan, compute and a declared runtime cap are the only remaining blockers.
