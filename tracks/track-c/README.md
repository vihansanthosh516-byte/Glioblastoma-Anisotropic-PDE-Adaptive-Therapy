# Track C — mathematical simulation and treatment policy

## Current role

Track C develops a two-compartment tumor-growth and treatment-suppression simulator and compares adaptive policies in synthetic worlds. Several prototypes and development outcomes already exist. The independent-world integrated benchmark is still a draft; no new sealed evaluation is authorized by this project map.

## Primary records

- [Independent-world benchmark draft](../../paper/trackc_world_benchmark_protocol_draft_2026-11-07.md)
- [World-seed capacity audit](../../paper/trackc_world_seed_capacity_2026-11-07.md)
- [Integrated-control prior-art update](../../paper/trackc_integrated_control_prior_art_update_2026-11-06.md)
- [Finite-set signal-ambiguity prototype](../../paper/trackc_robust_signal_ambiguity_prototype_2026-11-06.md)
- [Full prior-art review](../../paper/trackc_prior_art_full_2026-11-04.md)
- [Integrated plan](../../paper/trackABC_integrated_research_plan_master_2026-11-07.md)

## Main implementation and checks

- `src/trackc_sim.py` — simulator and particle filter.
- `src/trackc_integrated_signal_model.py` — molecular-like and MRI-like signal generation and joint likelihood.
- `src/216_trackc_robust_signal_ambiguity.py` — finite-set dependence-aware minimax controller.
- `src/208_trackc_robust_belief_mpc.py`, `src/211_trackc_policy_benchmark.py` — MPC and shared policy benchmark.
- `src/trackc_world_schedule.py`, `src/trackc_world_runner.py`, `src/trackc_world_summary.py`, `src/trackc_seed_registry.py` — development orchestration, summaries, and seed bookkeeping.
- `tests/test_trackc_integrated_signal_conditioning.py`, `tests/test_trackc_robust_signal_ambiguity.py`, `tests/test_trackc_policy_benchmark.py`, `tests/test_trackc_world_schedule.py`, `tests/test_trackc_world_runner.py`, `tests/test_trackc_world_summary.py` — focused checks.
- `output/trackc_stage2_*/` and `output/trackc_*sensitivity*/` — exposed development records; never treat them as fresh evaluation data.

## Confirmatory evaluation stopped (2026-11-09)

The 14-day precision target is unchanged and cannot be guaranteed at current seed capacity: the [feasibility audit](../../paper/trackc_population_mean_feasibility_audit_2026-11-08.md) needs at least 207 worlds per family and the seed block allows 7. The GE-betting interval matched its published equations and passed a predeclared synthetic coverage screen (30/30 cells), but it is not adopted, has no independent review, and reached 18.3 days at 207 worlds per family in its most favourable cell; no method met 14 days in any cell ([result](../../paper/trackc_ge_betting_coverage_audit_result_2026-11-09.md)). `tests/test_trackc_config_parity_extended.py` pins further settings against the config and marks three known drifts as expected failures. Confirmatory evaluation stays stopped until the user authorizes more worlds with fresh reserved seeds and a runtime cap, or explicitly changes the target. See the [open-gap record](../../paper/gap_status_2026-11-09.md).

## Gate before a research-scale evaluation

The driver, exact cohort sizes and world count, precision-based planning rule, fresh seed allocation, final inference settings, runtime limit, and dated protocol seal must be completed before evaluation outcomes are generated. If the predeclared precision or resource gate fails, report feasibility or an inconclusive result; do not change settings after seeing evaluation outcomes.
