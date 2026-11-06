# Track C paired-policy benchmark test record (2026-11-06)

## Scope

This run checked synthetic implementation mechanics only. It did not compare policy outcomes for research, read patient data, or authorize a sealed evaluation.

Code commit: `871bd3d` (`feat(trackc): add paired shared-endpoint policy benchmark`).

## Implemented and checked

`src/211_trackc_policy_benchmark.py` runs six policies on one virtual evaluation cohort and one shared baseline/RANO-like treatment path:

1. Leave-one-out dependence-aware signal MPC.
2. Matched leave-one-out conditional-independence MPC (`d = 0`).
3. No-signal robust MPC.
4. Shared-endpoint Q-learning trained on selection scenarios.
5. Fixed schedule selected on selection scenarios.
6. Standard six-course schedule.

Selection and evaluation scenarios use separate inputs. The two leave-one-out arms use unlabeled baseline signals and exclude each patient's signal from that patient's peer-weight fit. Evaluation outcomes do not fit the signal weights.

The paired score is the earlier of the two progression times, restricted at `tau_days`, minus the course penalty for courses whose cycle starts before `tau_days`. The simulator now returns `courses_through_cutoff` when `course_cutoff_days` is supplied. This matches the shared-endpoint Q-learning reward at the restriction boundary. Total courses remain available separately.

## Test evidence

- The first full-harness smoke check exposed a budget-grid mismatch: the default fixed-schedule catalog included schedules over a reduced six-course test budget. The runner now filters the default catalog to the configured budget and still refuses over-budget custom catalogs.
- Scratch suite: 121 passed in 608.87 seconds across the new benchmark, simulator, Q-learning, robust MPC, integrated-signal conditioning, paired-signal, and leave-one-out suites.
- After tightening the course-cutoff boundary assertion, the new benchmark suite passed again: 5 passed in 10.63 seconds.
- SHA-256 comparison found 17 source/test files identical between the scratch copy and the working tree.
- The scratch copy was `C:\Users\vihan\AppData\Local\Temp\trackc-policy-benchmark-7e59e425143d4a56be726f89a5dc1c8a`.

## Limits and next work

The six-arm harness is not the full baseline panel in `trackc_candidate_methods_design_2026-11-04.md`. It does not yet run the scan-rule and state-estimation policies, equal-course fixed-schedule frontier, dose-dense proxy, or six-versus-twelve-course comparison. It also has no multi-world aggregation, precision gate, or frozen outcome protocol. These must be resolved before any policy outcome comparison.

All tests use assumed virtual patients. Passing tests establish only that the implemented paths follow their contracts. They do not show that any policy improves treatment utility, describes real biology, is unique, or is better than published methods.
