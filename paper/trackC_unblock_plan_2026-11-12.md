# Track C unblock plan: make compute the only blocker (2026-11-12)

Goal set by the user: remove every Track C blocker except compute, so that the confirmatory evaluation can run on a DGX Spark.

Unchanged rules: the 14-day half-width target, the success rule, and the stop rule (`observed_two_scan`) stay as frozen. No policy outcome is generated until every gate below is closed and the run is authorized. Sources: `configs/trackABC_integrated_v2.json` (`not_yet_locked_before_evaluation`), `paper/gap_status_2026-11-09.md` (C3 to C13), `paper/trackC_integration_followup_2026-11-10.md`, `paper/status_map_2026-11-10.md`.

## One honest point first

Compute can make the precision target **reachable**. It cannot promise a **positive** result. With enough worlds, the result will be one of: positive, negative, or (if the target is still not met in some cells) inconclusive. All three are valid outcomes. The plan does not tune anything toward a positive result.

## All open items (from the files)

| # | Open item | Source | Type | Fix |
|---|---|---|---|---|
| C-1 | Controller model grid: 3 q levels (81 models) in code vs 5 levels (135) in config | gap C9; config `synthetic_q_signal_grid` | **User decision** | D1 below; then align code and config; remove strict xfail test |
| C-2 | 112-day scan sensitivity has no forecast length (`scan_every=4` rejected with 6 cycles) | gap C10 | **User decision** | D2 below; then fix `PolicyBenchmarkSettings`; remove xfail |
| C-3 | `RobustBeliefMPC` default stop rule is legacy `"either"` | gap C11 | Code | Driver passes `observed_two_scan` explicitly; test asserts it |
| C-4 | No module reads the config at run time (values duplicated in code) | gap C8 | Code | Runtime config loader; driver reads every setting from config; parity tests stay |
| C-5 | Virtual-patient count per selection, calibration, and evaluation cohort | config list item 1 | **User decision** | D3 below (sets runtime) |
| C-6 | Number of worlds by an outcome-blind precision and timing rule | config item 2; gap C4 | Code + compute | Synthetic precision screen (Step 5) picks n before any policy outcome |
| C-7 | Interval method not adopted (GE-betting passed 30/30 synthetic cells, no independent review) | config item 3; gap C3 | Review + amendment | Step 3 |
| C-8 | Paired-arm RNG coupling: particle-filter draws can desynchronize | `paper/trackc_population_mean_feasibility_audit_2026-11-08.md` | Code + test | Separate per-arm filter streams; include coupling in the interval screen |
| C-9 | Seed table not materialized; registry not reserved | config item 4; gap C7 | Code | Step 4 |
| C-10 | Runtime cap and feasibility decision | config item 5; gap C6 | Timing + decision | Step 6 on the Spark |
| C-11 | Research-scale driver (runs all worlds, verifies outputs) | config item 6; gap C7 | Code | Step 2 |
| C-12 | Missingness (0.0/0.1), course-penalty sweep (0/10), bootstrap resamples not orchestrated | gap C12 | Code | Inside the driver |
| C-13 | Duke aggregate permission and checksum | config item 7 | **User decision** | D4: recommend **not using it** (standard-normal fallback already defined) |
| C-14 | Clean-install environment verification | config item 8 | Code + test | Step 7 (lock file; fresh install on the Spark, ARM Linux) |
| C-15 | q_signal grid and primary sampling rule | config item 9 | Same as D1 | Fixed with C-1 |
| C-16 | Rejected max-T helper still callable | gap C13 | Code | Driver refuses it; test |

Known weaknesses that are **not** blockers (they stay in every report): the particle filter undercovers (joint 0.30 to 0.46 versus a nominal 0.90; `paper/trackc_prior_recovery_diagnostic_result_2026-11-10.md`); signal quality and dependence are assumed; the simulator builds in dependence, so the `d = 0` arm is misspecified by design; results are simulation only.

## Decisions only the user can make

- **D1 (C-1, C-15):** 81 models (3 q levels, as the code has now) or 135 models (5 levels, as the config has). Recommendation: **81**. It is the tested code path and costs less compute. 135 adds a finer grid and about 1.7 times the controller cost.
- **D2 (C-2):** forecast length for the 112-day scan sensitivity. Recommendation: **8 cycles (224 days)**, the smallest value that holds two 112-day scans.
- **D3 (C-5):** patients per cohort. Recommendation: keep the declared design value if one exists in the config; otherwise **6 per cohort**. Runtime scales about linearly with it.
- **D4 (C-13):** use the Duke aggregate? Recommendation: **no**. It avoids the licence question and keeps Track C independent of a spent cohort.
- **D5 (C-7):** reviewer route for GE-betting. Options: (a) a person outside this project (statistician, teacher); (b) a documented machine review: an independent fresh-context reviewer plus a cross-check against the authors' reference implementation (the `confseq` package by Waudby-Smith and Ramdas). Recommendation: (b) now, and (a) if a person becomes available. The report must say which route was used.

## Steps, in order

Each step: write the code, add synthetic tests, run them in a scratch copy, commit and push. No policy outcome is produced before Step 8.

1. **Record the decisions** (D1 to D5) in a dated amendment `paper/trackc_amendment_unblock_<date>.md`. Update the config to match. Close C-1, C-2, C-13, C-15 in the config list.
2. **Build the driver** `src/trackc_research_driver.py`:
   - reads every setting from the config (C-4);
   - passes the stop rule explicitly (C-3) and refuses the max-T helper (C-16);
   - runs worlds in parallel with `multiprocessing` (worker count set by flag);
   - writes one checkpoint file per world, so a crash or restart resumes without rerunning finished worlds;
   - runs the missingness, penalty, and bootstrap settings (C-12);
   - verifies at the end that every planned world has a complete, hash-recorded output (C-11);
   - keeps the existing `mode="evaluation"` refusal until the config gate list is empty.
   Tests: stub benchmark end to end; resume after a killed worker; output verification catches a missing world; stop-rule and max-T refusals.
3. **Interval method review and adoption** (C-7):
   - cross-check `src/trackc_ge_betting.py` against `confseq` on fixed synthetic samples (same inputs, same interval, to a stated tolerance);
   - an independent fresh-context review of the math against the paper, written to `paper/trackc_ge_betting_independent_review_<date>.md`;
   - if both agree, adopt GE-betting by amendment. If not, stop and report.
4. **Seed table** (C-9): materialize every world, cohort, signal, and policy seed with the declared `PCG64DXSM` derivation; reserve the block in `src/trackc_seed_registry.py`; give separate filter streams per arm (C-8). Tests: no reuse across worlds and arms; paired arms share exogenous scan draws; registry refuses overlaps.
5. **Outcome-blind precision screen** (C-6): with GE-betting and the real coupling, simulate world-level differences from the declared synthetic families (no policy outcomes) over a grid of world counts. Pick the smallest n per family at which the predicted half-width is at most 14 days in **every** declared cell. Record n and the cells where it is impossible at any n under a stated cap. This sets the compute bill.
6. **Timing on the DGX Spark** (C-10): run 3 world pairs with the real benchmark at the D3 cohort size, single worker and all workers. Record hours per world pair. Compute total hours = n x time. Write the runtime cap and the go/no-go into the amendment **before** Step 8.
7. **Clean environment** (C-14): pinned `requirements-trackc.lock`; fresh virtual environment on the Spark (ARM Linux); run the full Track C test suite there; record the result.
8. **Empty the gate list and run.** When C-1 to C-16 are closed, the config gate list is empty and `mode="evaluation"` unlocks. Run once with the fixed n. No reruns with changed settings.
9. **Report** with all arms, the null-information controls, the success rule as frozen, and every weakness listed above.

## What will be true after Steps 1 to 7

Only compute stands between the project and the Track C run. The remaining risks are honest outcomes, not blockers: the run may show no gain, a loss, or an inconclusive result in cells where 14 days cannot be reached.

## Rough effort

Steps 1 to 5 and 7 are code and review on this PC: about 2 to 4 working sessions. Step 6 needs about 1 hour on the Spark. Step 8 is the compute: estimated 11 to 19 hours for 207 to 350 world pairs at 2 patients per cohort, multiplied by the cohort-size factor from D3 (estimate only; see `paper/gap_status_2026-11-09.md` C6).
