# Observation-based treatment stop rule proposal (2026-11-06)

**Status: proposed and wired as an opt-in development setting; not adopted as the primary rule.** The benchmark default remains `stop_endpoint="either"`. This proposal changes only how simulated treatment stops. It does not change the project title, patient-data analyses, or Track A/B results.

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
- Added a shared stop-state helper used by both the simulator and belief-MPC forecasts, so execution and planning use the same threshold and missing-scan rule.
- Added the setting to `PolicyBenchmarkSettings`. Fixed-schedule selection, baseline selection, Q-learning training, belief-MPC candidate paths, calibration, evaluation, result validation, and saved settings now receive the same selected stop rule. A mismatch is rejected.
- Q-learning no longer truncates observed-stop training at latent endpoint times. It assigns endpoint reward at the observed stop or the restricted horizon and applies course costs through that point.
- Belief MPC forecasts noisy future scans with common noise draws across candidate treatment patterns, updates the same two-scan state, and stops each candidate path from that observed state. Forecasting assumes future scans are present; missing-scan prediction is not modeled by MPC.
- Added tests for helper behavior, hidden-event independence in MPC stopping, observed-stop Q-learning termination, and consistent end-to-end benchmark propagation.
- Scratch verification: **62 tests passed in 599.02 seconds** across the new stop-rule tests, `test_trackc_sim.py`, and `test_trackc_scan_endpoint.py`. SHA-256 matched the four copied simulator source modules. No outcome benchmark, patient-data analysis, or sealed evaluation ran.
- After a docstring-only clarification, the three new tests passed again in **2.23 seconds**; SHA-256 matched the simulator source and all three copied test files.
- Latest scratch verification: **71 tests passed in 75.57 seconds** across `test_trackc_observed_stop.py`, `test_trackc_qlearning_baseline.py`, `test_trackc_robust_belief_mpc.py`, `test_trackc_robust_signal_ambiguity.py`, and `test_trackc_policy_benchmark.py`. The end-to-end benchmark test used two patients and tiny settings to test wiring; it is a mechanics check, not a research result. Hash comparison against the scratch copy remains to be recorded.
- Broader default-path regression check: **80 tests passed in 1,331.23 seconds** across `test_trackc_sim.py`, `test_trackc_scan_endpoint.py`, `test_trackc_policy_context.py`, `test_trackc_loo_model_average.py`, and `test_trackc_integrated_policy_utils.py` in the same scratch worktree. SHA-256 matched all nine copied source/test files.
- No research-scale multi-world outcome evaluation, held-out analysis, patient-data analysis, or sealed run was performed. Existing default behavior remains `either` unless the caller explicitly selects `observed_two_scan`.

## Required before a research-scale outcome comparison

1. Decide whether this observable trigger is the primary execution rule. Keep latent-oracle stopping as a labeled sensitivity analysis only.
2. Freeze the trigger, scan interval, missing-scan behavior, endpoint estimands, utility weight, world-level precision target, uncertainty intervals, multiplicity rule, and fresh seed blocks before outcome scoring.
3. The MPC now has an opt-in iid future-scan-missingness assumption, shared across candidate schedules. Before any outcome comparison, choose its values using development worlds only and test mismatch against state-linked missingness. The default remains complete future follow-up.
4. Add the registered independent-world precision and held-out scoring path. Report results across every predeclared signal anchor, including harm.
5. Do not change settings after viewing outcomes to seek a positive result. A positive or negative outcome must come from the frozen comparison.

## Research-position check

The umbrella combination remains established prior art. M4RL combines a multiscale GBM biology model with reinforcement-learning drug scheduling ([Liu et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC13155510/)); a 2026 preprint combines reaction-diffusion, learned correction, adaptive updating, and treatment-schedule MPC ([Liu et al., 2026](https://arxiv.org/abs/2607.13877)); and a 2026 preprint studies robust treatment-policy learning under model and measurement uncertainty in prostate-cancer simulation ([Xie et al., 2026](https://arxiv.org/abs/2610.22336)). A 2026 GBM immunotherapy systems paper is a conceptual review, not a validated treatment algorithm ([Chen et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13349897/)).

The narrower dependence-aware molecular-like/MRI-like policy question remains a possible distinction, not a verified novelty claim. Its signal-to-resistance relationship is simulated, and its prior predictive gate failed. This stop-rule improvement addresses simulator observability and fairness; it does not establish novelty, superiority, treatment benefit, or a positive result.

## Proposed MPC missing-scan sensitivity (2026-11-06)

`RobustBeliefMPC` accepts `assumed_future_scan_miss_probability`; the paired benchmark exposes it as `mpc_assumed_future_scan_miss_probability`. For each replan, the MPC samples future scan availability under this fixed iid probability and reuses the same availability draws across candidate treatment patterns. A missing simulated scan resets two-scan confirmation. The setting is valid only with `stop_endpoint="observed_two_scan"`.

The default is `0.0`, which assumes complete future follow-up and preserves the seeded treatment plan. The setting does not read the scenario's hidden `miss_p` or `miss_p_hi`. A value of `1.0` makes all future forecast scans unavailable. This is a planner assumption, not a learned patient-specific missingness model. It does not represent state-dependent or informative missed scans inside the forecast. It remains opt-in and proposed; no policy outcome comparison has used it.
