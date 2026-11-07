# Observation-based treatment stop rule proposal (2026-11-06)

**Status: proposed; not yet wired into the policy benchmark.** This proposal changes how simulated treatment stops. It does not change the project title, patient-data analyses, or Track A/B results.

## Why change it

The current policy benchmark uses `stop_endpoint="either"`. The simulator stops later treatment when hidden true progression reaches either latent endpoint. This keeps the hidden event from entering the policy, but gives the simulated treatment executor perfect knowledge of progression. That is an unrealistic source of information and can change course use and policy utility.

The simulator now has an optional `stop_endpoint="observed_two_scan"` mode. It uses only the noisy MRI-like log-burden readings that occur on scheduled scan visits:

1. The first finite scheduled scan sets the observed baseline and nadir.
2. A follow-up scan is positive when its observed log-burden is at least `log(1.4) + max(observed_nadir, observed_baseline + log(0.10))`.
3. Treatment stops before a course is given at the second consecutive positive scheduled scan.
4. A scheduled scan with a missing reading resets the confirmation count. Non-scan cycles do not count as scans and do not reset it.
5. A negative scan resets the count and may lower the observed nadir.

The baseline and RANO-like progression outcomes remain latent simulator outcomes measured on the policy-specific trajectory. They do not control treatment in this mode. The rule is a deliberately simple noisy-volume trigger. It is **not** formal RANO, a clinical decision rule, or a claim that RANO always requires confirmation. RANO 2.0 uses multiple clinical and imaging criteria and limits confirmation requirements to specified settings ([Wen et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10860967/)).

## Evidence and code status

- Added to `src/trackc_sim.py` as an opt-in mode. Existing `baseline`, `rano`, and `either` modes remain available.
- Added three tests for two-scan confirmation timing, missing-scan handling, and independence from the hidden progression threshold.
- Scratch verification: **62 tests passed in 599.02 seconds** across the new stop-rule tests, `test_trackc_sim.py`, and `test_trackc_scan_endpoint.py`. SHA-256 matched the four copied simulator source modules. No outcome benchmark, patient-data analysis, or sealed evaluation ran.
- After a docstring-only clarification, the three new tests passed again in **2.23 seconds**; SHA-256 matched the simulator source and all three copied test files.
- Benchmark-compatibility check: **16 tests passed in 71.36 seconds**; SHA-256 matched ten copied Track C source modules. The active benchmark still exercised its existing `either` mode.
- The active policy benchmark still uses `stop_endpoint="either"`. Its controller planning, Q-learning training, schedule selection, evaluator validation, and output metadata have not yet been changed to use the new rule.

## Required before any policy outcome comparison

1. Decide whether this observable trigger is the primary execution rule. Keep latent-oracle stopping as a labeled sensitivity analysis only.
2. Pass the same stop rule through policy training, candidate-action simulation, comparator selection, fixed-schedule selection, and held-out evaluation. Refuse a run if one arm uses a different rule.
3. Test that policy resets and action selection never receive latent progression times or hidden parameters. Test that changing only hidden endpoint thresholds cannot change actions or treatment stops under the observable rule.
4. Freeze the trigger, scan interval, missing-scan behavior, endpoint estimands, utility weight, world-level precision target, uncertainty intervals, multiplicity rule, and fresh seed blocks before outcome scoring.
5. Report results across every predeclared signal anchor, including harm. Do not change settings to seek a positive result.

## Research-position check

The umbrella combination remains established prior art. M4RL combines a multiscale GBM biology model with reinforcement-learning drug scheduling ([Liu et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC13155510/)); a 2026 preprint combines reaction-diffusion, learned correction, adaptive updating, and treatment-schedule MPC ([Liu et al., 2026](https://arxiv.org/abs/2607.13877)); and a 2026 preprint studies robust treatment-policy learning under model and measurement uncertainty in prostate-cancer simulation ([Xie et al., 2026](https://arxiv.org/abs/2610.22336)). A 2026 GBM immunotherapy systems paper is a conceptual review, not a validated treatment algorithm ([Chen et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13349897/)).

The narrower dependence-aware molecular-like/MRI-like policy question remains a possible distinction, not a verified novelty claim. Its signal-to-resistance relationship is simulated, and its prior predictive gate failed. This stop-rule improvement addresses simulator observability and fairness; it does not establish novelty, superiority, treatment benefit, or a positive result.
