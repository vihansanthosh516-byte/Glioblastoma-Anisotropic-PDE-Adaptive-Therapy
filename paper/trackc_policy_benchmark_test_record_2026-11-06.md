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

### Follow-up audit: signal-prior alignment

A code audit found that the LOO signal controllers used the signal model's broad reference `f0` bounds by default, while the simulator drew `f0` from the narrower log-uniform range `[0.001, 0.1]`. Because the joint signal likelihood integrates over `f0`, this gave the controllers a nuisance prior that did not match the virtual-patient generator. The benchmark now defaults `signal_f0_bounds` to the simulator range, refuses selection or evaluation scenarios outside the declared range, and passes the same bounds to both dependence-aware and conditional-independence controllers. A regression test captures both constructor calls and verifies the bounds. This changes no outcome result; none has been run.

The updated six-suite scratch run passed **80 tests in 25.40 seconds**. SHA-256 matched all seven relevant source/test files between the tested copy and working tree. This is implementation evidence only; it does not show policy benefit.

### Follow-up audit: common observation schedule

A second audit found that `_evaluate_policy` did not attach `scan_every` to the evaluation scenarios. The simulator treats a missing field as a scan every 28-day cycle, while the shared Q-learning evaluator explicitly used the configured interval (normally 84 days). Thus the RANO-like endpoint and observation schedule could differ across arms. The benchmark now copies each evaluation cohort and explicitly sets the same configured scan interval before scoring every arm. The regression test records each simulator call during a complete small harness run and verifies that fixed-schedule selection, Q-learning training, and all six evaluation arms receive the same interval. The test failed before the fix and passed after it.

After this correction, the six focused scratch suites passed **80 tests in 25.41 seconds**. SHA-256 matched all seven relevant source/test files. No policy outcomes were run; earlier outcome comparisons were not made with this harness.

## Limits and next work

The six-arm harness is not the full baseline panel in `trackc_candidate_methods_design_2026-11-04.md`. It does not yet run the scan-rule and state-estimation policies, equal-course fixed-schedule frontier, dose-dense proxy, or six-versus-twelve-course comparison. It also has no multi-world aggregation, precision gate, or frozen outcome protocol. These must be resolved before any policy outcome comparison.

All tests use assumed virtual patients. Passing tests establish only that the implemented paths follow their contracts. They do not show that any policy improves treatment utility, describes real biology, is unique, or is better than published methods.
