# Track C paired-policy benchmark test record (2026-11-06)

**Current-status correction (2026-11-08):** The integrated world runner selects the finite-set minimax versus matched `d = 0` minimax contrast as primary. Older sections below document prior LOO-primary mechanics and are superseded for the integrated v2 path; LOO remains a separately reported contrast.

## Scope

This run checked synthetic implementation mechanics only. It did not compare policy outcomes for research, read patient data, or authorize a sealed evaluation.

Code commit: `871bd3d` (`feat(trackc): add paired shared-endpoint policy benchmark`).

## Interpretation correction (2026-11-06)

The phrase “shared path” below refers to scoring the baseline and RANO-like endpoints on the same path **within one policy**. Different policies produce different action-dependent paths. The paired benchmark gives policies the same initial virtual-patient inputs, precomputed observation-noise stream, scan schedule, and global treatment budget; it does not force their realized disease trajectories to match. The paired utility is calculated per patient as the earlier capped endpoint time minus that policy's courses through the cutoff, then averaged within each world. It is not the smaller of the two cohort-average endpoint times.

The current worktree further adds ordered synthetic-patient IDs, disjoint selection/evaluation validation, ID fingerprints, caller-declared scenario seeds, and retained endpoint-cap flags. This follow-up is verified in the current mechanics record below; it does not retroactively change the historical test evidence below.

### Follow-up audit: latent simulator inputs at policy reset (2026-11-06)

`simulate_two_endpoints` now defaults to `policy_context="observable"` and passes policies only a copy of the scan-precision field `sigma` during `reset`. Full scenario context requires explicit `policy_context="full"` and is reserved for oracle diagnostics. The older single-endpoint `simulate` function is unchanged. The active benchmark and fixed-schedule selection explicitly request observable context; the result record names the exposed field and records that latent scenario fields were not exposed. `FixedSchedulePolicy` and `HolidayContainment` initialize from sigma only. `CourseCappedPolicy` strips extra fields before it resets its inner policy, even if a caller passes a full scenario mapping directly.

The scratch copy passed **12 focused tests in 2.69 seconds**: six context tests, paired-endpoint simulator checks, and fixed-schedule/course-cap checks. SHA-256 matched six source/test files between scratch and worktree. An earlier broader attempt was stopped before completion; the scratch copy initially lacked `trackc_oracle.py`, so it produced no valid aggregate suite result. No policy outcomes, patient data, or sealed evaluations were used. This verifies the reset boundary and compatibility only; it does not show treatment benefit.

## Implemented and checked

`src/211_trackc_policy_benchmark.py` runs six policies on one virtual evaluation cohort. Each policy has its own action-dependent treatment path; the baseline and RANO-like endpoints are measured on that policy's path:

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

A code audit found that the LOO signal controllers used the signal model's broad reference `f0` bounds by default, while the simulator drew `f0` from the narrower log-uniform range `[0.001, 0.1]`. Because the joint signal likelihood integrates over `f0`, this gave the controllers a nuisance prior that did not match the virtual-patient generator. The benchmark now defaults `signal_f0_bounds` to the simulator range, refuses selection or evaluation scenarios outside the declared range, and passes the same bounds to both dependence-aware and conditional-independence controllers. A regression test captures both constructor calls and verifies the bounds. No research-scale or fresh-seed policy comparison used this harness; the tiny contract test only checked plumbing.

The updated six-suite scratch run passed **80 tests in 25.40 seconds**. SHA-256 matched all seven relevant source/test files between the tested copy and working tree. This is implementation evidence only; it does not show policy benefit.

### Follow-up audit: common observation schedule

A second audit found that `_evaluate_policy` did not attach `scan_every` to the evaluation scenarios. The simulator treats a missing field as a scan every 28-day cycle, while the shared Q-learning evaluator explicitly used the configured interval (normally 84 days). Thus the RANO-like endpoint and observation schedule could differ across arms. The benchmark now copies each evaluation cohort and explicitly sets the same configured scan interval before scoring every arm. The regression test records each simulator call during a complete small harness run and verifies that fixed-schedule selection, Q-learning training, and all six evaluation arms receive the same interval. The test failed before the fix and passed after it.

After this correction, the six focused scratch suites passed **80 tests in 25.41 seconds**. SHA-256 matched all seven relevant source/test files. No research-scale outcome comparison was run; the n=2 mechanics test invokes the paths but does not assess performance direction. Earlier research-scale comparisons were not made with this harness.

### Follow-up audit: broader comparator coverage

The benchmark now defines 9 scan-rule comparators, 4 state-estimation comparators, a 24-cell raw-threshold grid, and a 128-cell filtered-threshold grid taken from existing Track C code. At a full benchmark run, all 165 comparator settings are scored on selection patients using mean robust utility. The scan-rule, state-estimation, raw-threshold, and filtered-threshold family winners are recorded; only the selected cells from the two grids are run on evaluation patients. The other declared baselines are evaluated individually. The fixed-schedule catalog already covers the front-loaded and paced schedule references.

Every comparator uses the configured scan interval, the same endpoint definitions within its own path, the global treatment budget, and its explicit arm-specific cap. The historical filtered-threshold setting keeps its old fixed configuration and provenance. The policy selection uses only the supplied selection cohort; the function receives no evaluation cohort for selection.

The integration test checks the family winners, saved settings and provenance, arm-specific course cap, shared scan interval, and evaluation of one selected cell from each threshold grid. It checks mechanics on a two-patient synthetic cohort only. It does not compare performance direction or use patient data.

Before adding the raw-threshold grid, 10 focused scratch suites passed 153 tests in 647.09 seconds. The benchmark suite passed 8 tests in 27.31 seconds after adding both threshold grids. Final full-suite verification is recorded below.

#### Final mechanics verification after both threshold grids

- The scratch copy ran the 10 focused suites listed below. Result: **153 passed in 642.21 seconds**.
- SHA-256 comparison matched all **439 Python files** under `src/` and `tests/`: no missing, extra, or different files.
- `graphify update .` completed after the code changes. It rebuilt the graph with 19,292 nodes, 31,654 edges, and 1,571 communities.
- The integration test uses a two-patient synthetic cohort and reduces each threshold grid to two settings. It verifies selection/evaluation separation, family-winner reporting, configuration uniqueness in the full catalog, shared scan schedule, and course caps. It does **not** execute all 165 settings on a full selection cohort or compare policy outcomes.
- No patient data, research-scale outcome comparison, or sealed evaluation was used.

Focused suites: `test_trackc_policy_benchmark.py`, `test_trackc_sim.py`, `test_trackc_qlearning_baseline.py`, `test_trackc_robust_belief_mpc.py`, `test_trackc_integrated_signal_conditioning.py`, `test_trackc_paired_signal_benchmark.py`, `test_trackc_loo_model_average.py`, `test_trackc_policies.py`, `test_trackc_regime.py`, and `test_trackc_scan_endpoint.py`.

## Limits and next work

### Current mechanics verification after cohort-identity and utility corrections

- The ten focused scratch suites listed above passed **159 tests in 640.91 seconds** after the current code and test files were copied to the scratch tree.
- SHA-256 comparison matched all **439 Python files** under `src/` and `tests/` before and after the run.
- New checks cover patientwise worst-endpoint utility (including a counterexample where it differs from the minimum of endpoint means), overlapping selection/evaluation IDs rejected before fitting, unique/canonical IDs, ID and scenario-seed provenance fields, and preservation of endpoint-cap flags.
- The run used assumed virtual patients only. It did not run a research-scale policy comparison, patient analysis, or sealed evaluation.
- `graphify update .` completed after the source change.

At this checkpoint, the harness contained scan rules, state-estimation policies, and selection across the 24-setting raw-threshold and 128-setting filtered-threshold grids. The fixed-course references and equal-course frontier were added in later follow-ups below. Multi-world aggregation, a precision target, exact held-out world anchors, and a frozen interval method remain open. A dose-dense comparator remains infeasible at the current 28-day time resolution. Resolve the open items before research-scale outcome comparison. Cohorts from one assumed world support development mechanics only.

### Follow-up: explicit fixed-course references (2026-11-06)

The shared benchmark now evaluates three explicit schedule references when the configured budget permits them: six courses at the first six 28-day cycle starts, twelve courses at the first twelve cycle starts, and twelve courses spaced every other cycle (56-day intervals). The twelve-course arms are omitted when the configured global budget is below twelve. Result metadata records planned course counts, cycle indices, intervals, and the non-clinical role of these synthetic schedules. The primary paired contrast is unchanged.

The 56-day schedule is a pacing control, **not** a dose-dense regimen. The simulator's fixed 28-day cycle cannot represent a treatment interval shorter than 28 days; a genuine dose-dense comparator would need a separately reviewed time-resolution change. That follow-up checked schedule construction, budget handling, and runner inclusion. No policy outcome comparison was run. The equal-course frontier prototype is documented in the follow-up below. Multi-world aggregation, a precision target, exact held-out world anchors, and a frozen interval method remain open.

### Follow-up: equal-course fixed-schedule frontier (2026-11-06)

The benchmark now builds an upper fixed-schedule frontier from selection patients, estimates each feedback policy's mean courses through the horizon on an independent calibration cohort, and selects a randomized mixture of adjacent frontier schedules at that expected dose. It evaluates the mixture only on the disjoint evaluation cohort. Targets outside the observed fixed-schedule dose range are marked unmatchable; the code does not extrapolate. Evaluation records actual course-use mismatch and paired differences in worst-endpoint restricted days and robust utility. The mixture does not match each patient's realized dose.

The result fingerprints selection, calibration, and evaluation scenarios, IDs, and calibration/evaluation signals. It records distinct scenario-generation, signal-generation, policy-simulation, evaluation, and mixture-master seeds. Each schedule assignment uses a stable seed derived from the mixture master seed and policy name. No evaluation outcome enters frontier selection, dose calibration, or schedule assignment.

The focused scratch run passed **22 tests in 58.57 seconds** across `test_trackc_policy_benchmark.py` and `test_trackc_equal_course_frontier.py`. It used two synthetic patients in each cohort and reduced threshold grids. This verifies code paths only; it does not show that any policy performs better. The successful graph update used `--code-only` because the full update required an unavailable language-model key.

Remaining limits: the project has not run the benchmark at research scale or across held-out virtual worlds. It still needs exact world-anchor enforcement, a precision target, a frozen interval method, and registered development seeds before outcome comparison. The broad project title remains unchanged; this mechanics work does not establish global uniqueness, biological validity, clinical benefit, or superiority.

All tests use assumed virtual patients. Passing tests establish only that the implemented paths follow their contracts. They do not show that any policy improves treatment utility, describes real biology, is unique, or is better than published methods.

### Follow-up: primary-contrast alignment (2026-11-08)

The integrated v2 configuration names finite-set minimax versus its matched `d = 0` control as the primary contrast. The benchmark previously always wrote the LOO contrast into `primary_paired_contrast`, while the world runner also ran the finite-set arms and labeled them exploratory. That output label conflicted with the integrated protocol even though no result had been scored for a research comparison.

The runner now accepts an explicit primary contrast. The integrated world runner selects `minimax`; the benchmark records both minimax and LOO contrasts separately and writes the chosen one to `primary_paired_contrast`. World summaries infer the recorded primary contrast when callers omit `contrast`, reject mixed primary labels, and still allow explicit secondary LOO summaries. The minimax summary checks recorded paired differences against the underlying candidate and control arrays.

The three focused scratch suites passed **34 tests in 48.87 seconds**. Before and after tests, SHA-256 matched all **351 source files** and **116 test files** between the worktree and scratch copy. The first rerun exposed a missing provenance fixture in the synthetic adapter test; the fixture was corrected, then the suites passed. No policy outcomes were generated or inspected, and the integrated evaluation remains unfrozen.

### Full Track C suite after manifest-test fixture repair (2026-11-06)

- The full scratch-copy run covered all 35 `test_trackc_*.py` files: **377 passed, 1 skipped in 850.89 seconds**.
- The skip was `test_trackc_seed_registry.py:109`, because the scratch copy intentionally omitted `output/`.
- I then copied only the referenced frozen-protocol evidence file into the scratch tree and reran that suite: **13 passed, no skips**.
- The first full run exposed two stale unit-test fixtures that hard-coded seed ranges no longer reserved for the test purpose. The manifest implementation was unchanged; the unit tests now use small synthetic seed values and a mocked registry, while the separate registry suite continues to exercise the real registry.
- SHA-256 comparison matched all **340 source files and 103 test files** between the scratch copy and worktree, excluding `__pycache__`.
- These remain mechanics tests only. No research-scale policy outcomes, patient data, or sealed evaluation were run.

### Follow-up: family-level max-T interval proposal (2026-11-08)

`src/trackc_world_summary.py` now provides a development helper for the two response-family contrasts named in the integrated v2 configuration. It accepts one equal-weight world-level mean difference per independent world and family, resamples worlds separately within each family, holds each observed family standard error fixed, and uses the higher empirical quantile of the maximum absolute centered bootstrap statistic. The reported intervals are symmetric around each observed family mean.

This is an **observed-SE max-T bootstrap proposal**, not a bootstrap-t method. It assumes independent world samples across response families; fingerprints cannot prove that assumption. It has not received independent statistical review and is not frozen. The helper is not wired into a research-scale driver or result writer. No policy outcomes were run or inspected.

Coverage remains unvalidated at the planned world counts. Westphal and Zapf's simulation study found small-sample error inflation for maxT in diagnostic-accuracy co-primary endpoint settings; those endpoints differ from Track C's world-level utility means, so the finding is a caution rather than direct evidence about this helper ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11025299/)).

Focused scratch tests cover deterministic family ordering, the two-family positive rule, a non-passing family, a missed precision target, and invalid inputs. The suite passed **22 tests in 0.13 seconds** in `C:\Users\vihan\AppData\Local\Temp\trackc-max-t-8lsunbxt`. SHA-256 matched all **467 files** under `src/` and `tests/` before and after the test run. No policy outcomes were run.

\n
