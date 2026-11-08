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

## Gate before a research-scale evaluation

The driver, exact cohort sizes and world count, precision-based planning rule, fresh seed allocation, final inference settings, runtime limit, and dated protocol seal must be completed before evaluation outcomes are generated. If the predeclared precision or resource gate fails, report feasibility or an inconclusive result; do not change settings after seeing evaluation outcomes.
