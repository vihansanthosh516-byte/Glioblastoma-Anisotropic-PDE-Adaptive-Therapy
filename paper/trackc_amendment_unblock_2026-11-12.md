# Track C amendment: unblock decisions and PC-side gate work (2026-11-12)

Dated before any v2 evaluation world exists. The user asked to "do everything you can on my computer" and to run the DGX Spark later. That instruction is the approval of the five decisions in `paper/trackC_five_decisions_evidence_2026-11-12.md`, with the corrections below. Unchanged: the 14-day half-width target, the 95% level, the success rule, the `observed_two_scan` stop rule, the utility, and the two response families.

## Decisions adopted

| # | Decision | Where it is in code and config |
|---|---|---|
| D1 | Controller model set = 5 q levels x 27 = **135 models** | `src/trackc_research_grid.py`; config `signal_model.controller_model_set_decision`. Script 212 keeps its 81-model grid as the historical record. |
| D2 | 112-day scan sensitivity uses **8 forecast cycles** (224 days); the other intervals keep 6 | config `controller_comparison.forecast_cycles_by_scan_interval_days`; `cell_spec("scan112")` |
| D3 | Cohort size and world count chosen by a **development pilot** with the rule fixed in `src/trackc_dev_pilot_design.py` | The pilot runs on the Spark; see the correction below |
| D4 | **No Duke aggregate** | config `input_rules.track_b.duke_aggregate_decision`; gate removed |
| D5 | Interval method: independent review + cross-check + outcome-blind screen | See "Interval method" below |

C-8 (paired-arm coupling): **keep** the existing coupling. All arms share scenarios, pre-drawn scan noise and the evaluation seed; later particle-filter draws may desynchronize between arms. This keeps each arm's marginal law and only affects the variance of the world difference (`paper/trackc_population_mean_feasibility_audit_2026-11-08.md`, "Policy and observation streams"). The interval's unit is the world, and worlds are independent, so validity holds. The pilot measures the variance that includes this coupling.

## Corrections found while doing the work

1. **Interval source.** The project GE-betting follows Martinez-Taboada and Ramdas (2026, arXiv:2608.21694), not the Waudby-Smith and Ramdas (2024) hedged-capital interval named in the D5 evidence note. ConfSeq implements the 2024 method. So the ConfSeq check is a comparison with a second valid method, not a numerical match.
2. **Cohort-size rule.** If time per world grows in exact proportion to cohort size m, the cheapest choice is always the smallest m. The reason: the needed worlds n grow with the variance tau^2 + sigma^2/m, so n x m grows with m. The best m depends on the fixed time per world. The rule is now: time per world = a + b x m, fitted from pilot runs at two cohort sizes (m = 2 and m = 12). This is the standard optimal-cluster-size result, m* = sqrt(sigma^2 a / (tau^2 b)).
3. **World count was over-estimated.** The earlier figures (about 350 worlds at SD 30, 1,400 at SD 60, 3,900 at SD 100) used a 1/sqrt(n) extrapolation. The direct screen below shows much lower counts.

## Interval method (C-7)

- Independent review: **PASS WITH NOTES** (`paper/trackc_ge_betting_independent_review_2026-11-12.md`). Every checked element matches the paper. Notes and answers are in that file: world differences must be fed in world-index order; n is fixed before evaluation; an empty interval (`None`) stops inference.
- ConfSeq cross-check (`output/trackc_interval_crosscheck/result.json`; 12 cells, 400 repeats each, synthetic only): all four checks pass. GE coverage was 1.0 in 11 cells and 0.9975 in one (nominal 0.975 per contrast). No empty intervals. GE median half-width was about half of ConfSeq's (for example 10.9 vs 20.5 days at n = 350, SD 30).
- Weakest point: coverage near 1.0 means the interval is conservative here. The wide support [-790, 790] drives this. No human statistician has reviewed the method.
- Decision: **adopt GE-betting** as the Track C v2 interval, delta 0.025 per family contrast (Bonferroni over 2 families). Gaffke's interval is not added: it was not needed to reach a feasible world count, and adding a method now would add a choice point.

## Outcome-blind precision screen (C-6)

`output/trackc_precision_screen/result.json`. Rule: smallest n per family at which the 80th percentile of the GE half-width is at most 14 days (200 synthetic repeats per point; normal and skewed shapes; worst shape shown).

| SD of world differences (days) | 15 | 20 | 30 | 40 | 60 | 80 | 100 | 120 |
|---|---|---|---|---|---|---|---|---|
| Worlds per family | 276 | 276 | 276 | 284 | 312 | 364 | 440 | 544 |

Below SD 30 the count stays at 276: the support width sets a floor. The real SD comes from the pilot.

## Seeds (C-9)

Reserved in `src/trackc_seed_registry.py`, checked free of every earlier seed:
- development pilot: scenario 30,000,000-30,099,999, policy 31,000,000-31,099,999 (`trackc_v2_dev_pilot`);
- evaluation: scenario 40,000,000-40,999,999, policy 41,000,000-41,999,999 (`trackc_v2_evaluation`).
The driver derives each world's 11 seeds from its index (`build_schedule`); worlds are nested, so world i keeps its seeds whatever n is.

## Driver (C-3, C-4, C-11, C-12, C-16)

`src/trackc_research_driver.py`: reads settings from the config; passes the stop rule and the signal f0 range explicitly; runs the seven cells (primary, scan28, scan56, scan112, miss10, penalty0, penalty10); runs worlds in parallel; writes one atomic, hash-checked checkpoint per world; resumes; `--max-hours` stops cleanly; `--verify-only` checks every planned world; refuses the evaluation phase while any gate is open; refuses the max-T method. Tests: `tests/test_trackc_research_driver.py`, `tests/test_trackc_dev_pilot_design.py`.

## Test status (scratch copy, Python 3.14.4)

First pass: 555 passed, 2 skipped, 1 xfailed, 7 failed. Six failed only because the scratch copy lacked `paper/`. One (`test_trackc_policy_context.py`) also failed on unchanged HEAD code; fixed in the second pass below. Final full-suite count: see `paper/trackc_spark_runbook_2026-11-12.md`, step 1.

## Pre-Spark fixes (same day, second pass)

- `test_trackc_policy_context.py`: the test's fake simulator still returned the retired `"either"` stop rule and no `observed_stop_day`, which the benchmark correctly refuses. The test was fixed; the code was not changed.
- Frozen analysis written before any evaluation world: `src/trackc_v2_analysis.py` (tests `tests/test_trackc_v2_analysis.py`). It applies the config `positive_rule` to the primary cell only, feeds worlds in schedule order, and returns "inconclusive" on an incomplete run, an empty interval, or an out-of-support value.
- GE-betting adopted in the config (`outcomes_and_inference.interval_method_adopted`, `inference_status`).
- Driver: refuses to resume or verify over checkpoints made with different code or config (each checkpoint stores file hashes); one math-library thread per worker so 18 workers do not oversubscribe the Spark's cores.
- Seed tables written: `output/trackc_v2_seed_tables/dev_pilot_50.json` and `evaluation_2000.json` (nested; covers any n up to 2,000 per family).

## Gates closed and still open

Closed 2026-11-12 (config `closed_gates_2026_11_12`): Duke aggregate, q grid, interval method, seed table, research driver.

Still open, all needing the Spark: cohort size and world count (pilot design), runtime cap, clean-install verification. No v2 evaluation world has been run.

**Freeze rule.** After the pilot design, the only allowed config change is writing `choice_m`, `choice_n_per_family` and the runtime cap, and removing those gates. After that, no code or config change: the driver refuses to mix checkpoints from different code.

## Cost

Measured local timing: see `paper/trackc_local_timing_2026-11-12.md` (one development world at m = 2, 135 models).
