# Timing-only follow-up amendment: calibration-world `f0` envelope

**Status:** frozen for a new synthetic timing-only pilot. This does not amend any policy-performance protocol.

## Reason for amendment

The first pilot stopped before policy simulation because `run_shared_policy_world` rejected generated scenario `f0` values outside the default `PolicyBenchmarkSettings.signal_f0_bounds` of 0.001–0.1. The calibration generator shifts the `f0` range across worlds. The failed run and its seed use remain recorded in `trackc_robust_policy_timing_failure_2026-11-06.md` and its manifest. Its seeds will not be reused.

## Frozen correction

The follow-up will set `signal_f0_bounds` to the common calibration-world envelope `F0_REFERENCE_BOUNDS = (0.001 × exp(-0.7), 0.9)`. This is the minimum lower endpoint and maximum capped upper endpoint allowed by the declared `log_f0_shift` ranges in both the `single` and `mixture` calibration families. The same envelope will be used in both worlds; no world-specific bounds will be supplied to the policy. `n_f0_nodes` and every other `PolicyBenchmarkSettings` field will remain at its code default.

This corrects the input contract. It is not selected from policy outcomes. The policy models, 81-point full signal grid, 27-point matched `d=0` grid, signal parameters, cohort size, and timing-only objective will stay unchanged.

## Fresh seed allocation

| Response family | World seed | Scenario and signal seeds | Policy seeds |
| --- | ---: | --- | --- |
| `single` | 20261952 | 20261953–20261957 | 20261958–20261962 |
| `mixture` | 20261963 | 20261964–20261968 | 20261969–20261973 |

The seeds will be assigned in the same order as the original timing protocol: selection, calibration, evaluation scenarios, evaluation signal generation, calibration signal generation; then schedule selection, Q-learning training, calibration simulation, evaluation simulation, and schedule-mixture assignment. They are distinct from the failed attempt and from the reserved Stage 3 block.

## Stop and reporting rules

The two worlds will run sequentially from a scratch copy on an otherwise idle PC. If either run fails, the error will be recorded and the follow-up will stop; the run will not be retried with changed settings or reused seeds. The script will discard each returned benchmark object without reading policy outcome fields. Its output will save only environment, code/protocol hashes, settings, seed assignments, elapsed seconds, and run status. No utilities, endpoint times, course counts, policy ranks, or other outcome scores will be written or inspected.

The result will estimate runtime only. It cannot establish that a controller improves treatment outcomes, is best, or is novel. Any later outcome comparison will need its own frozen worlds, estimand, anchors, precision target, and analysis method.
