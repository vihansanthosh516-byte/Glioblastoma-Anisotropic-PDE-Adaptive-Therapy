# Track C candidate-method and baseline protocol (draft, 2026-11-04)

**Status note (2026-11-06):** This is still a draft; no outcome evaluation has run. The endpoint-specific treatment paths and course counts described in the original version are superseded by the integrated plan's shared-path rule: treatment stops after the first simulated event under either endpoint, and both times use one course count. The event uses latent simulator state, so this is an idealized treatment-stop rule, not a clinical stopping rule. See `paper/trackABC_integrated_research_plan_submission_short_2026-11-05.md` for the current proposed wording.

## Why these candidates

Track C remains a simulation on assumed parameter ranges. Test H found value in knowing the simulated responder class in the mixture family, but it gave the class exactly, used one endpoint, and did not model an assay. The median class share of the oracle gap was 0.50, 0.21, and 0.52 at assumed course costs of 0, 5, and 10 days in `mixture`; for `mixture_mix` it was 0.63, 0.15, and 0.32. These exploratory values come from 11 and 20 controlled worlds (`output/trackc_stage2_H_value_of_information_mixture/results.json`, `output/trackc_stage2_H_value_of_information_mixture_mix36/results.json`; ledger C24 and C30). Test K found that no controlled world in either mixture configuration passed the MGMT strata and marker-prevalence checks together (`output/trackc_stage2_K_stratified_anchor_mixture/results.json`, `output/trackc_stage2_K_stratified_anchor_mixture_mix36/results.json`; ledger C14 and C19). The responder label is not validated as MGMT.

The first candidate is a hypothetical imperfect pre-treatment marker that selects a fixed schedule by marker result. The second is a belief-state MPC that uses noisy scans and selects a short treatment plan against both progression definitions. A generic hybrid or RL claim is not treated as novel because prior work already covers adaptive treatment, belief-state control, MPC, and Q-learning (`paper/trackc_prior_art_full_2026-11-04.md`). These candidates test narrower combinations; neither is assumed to be novel or to work.

## Frozen evaluation scope

- Use development controlled worlds only: `mixture` and `mixture_mix` for the marker candidate; include single-family configurations as negative controls if the run budget allows. Controlled-world files: `output/trackc_stage2_D_negative_controls_mixture/controlled_worlds.csv`, `output/trackc_stage2_D_negative_controls_mixture_mix36/controlled_worlds.csv`, `output/trackc_stage2_D_negative_controls/controlled_worlds.csv`, and `output/trackc_stage2_D_negative_controls_mix36/controlled_worlds.csv`.
- Use separate selection and evaluation patients within each world. Freeze each selected policy before evaluating it. Give every policy the same evaluation patients and the same scan-noise stream. The marker-error stream is independent of scan noise.
- Report RANO-like progression and baseline-referenced progression separately. The RANO-like rule is a simplified simulated rule, not clinical RANO: true burden at a scan is compared with 1.4 times the smallest earlier scanned burden, subject to the assumed 0.1-baseline floor. The baseline rule is continuous true burden at 1.2 times baseline.
- Report restricted mean time to progression at 365 and 730 days, mean courses, cost-adjusted utility at 0, 5, and 10 assumed days per course, the arm's gap above the full fixed-schedule frontier at its mean course use, and the paired worst-decile difference against the standard arm. Keep patient pairs together for all arm comparisons. Report intervals over worlds and patient-level paired bootstrap intervals; show each world because mixture controlled sets are small.
- Report sensitivity and specificity separately over a predeclared grid from 0.50 to 1.00. Keep simulated marker prevalence visible. Call the marker hypothetical; do not call it MGMT or imply that it measures a real response class.
- Use development seeds starting at 20261880 and noise keys starting after 1892400, as recorded in `src/trackc_seed_registry.py`. No sealed seed or old H-2/H-3 seed may be used.

## Candidate 1: imperfect marker-conditioned schedules

For each selection patient, derive the latent class only to simulate an assay result. The policy receives the noisy binary marker result, never the latent class or any true tumor parameter. Draw marker results with the declared sensitivity and specificity. In each marker group, select one schedule from the existing 56-schedule library by maximizing the mean robust training utility

`min(RMST_RANO, RMST_baseline) - lambda * courses`,

with `lambda` in 0, 5, and 10 days per course. Use the global selected schedule when a marker group has no selection patients. Freeze the two schedules before evaluating. Evaluate each endpoint and horizon separately, and compare the candidate with the global best fixed schedule and the full fixed frontier at the candidate's mean course count.

Report a sensitivity/specificity surface, marker prevalence, per-world gains, selected schedules, course use, and worst-decile paired outcomes. The exact-class policy may appear only as a non-buildable upper reference. The parameter oracle remains an upper reference only where the existing design defines it.

## Candidate 2: endpoint-robust belief MPC

At visible scans, update the existing particle filter from noisy burden readings and known past actions. At each decision, compare the available short treatment plans using posterior-predictive utility from both progression times on one shared forecast trajectory and one expected course count. The forecast will stop later treatment after the first latent event under either endpoint, finish any current 28-day cycle, and apply no treatment in later cycles. The policy will not receive the hidden event time. This is an idealized simulation rule, not a clinical stopping rule. Apply the first action and replan at the next visible scan. The policy must not access the simulator's true parameters or latent class. Record particle count, prior ranges, horizon, terminal approximation, and all assumed costs in the run protocol. If the implementation uses a surrogate terminal value rather than full posterior rollouts, label that approximation and keep its score separate from observed evaluation outcomes.

If a scheduled scan has no usable reading, do not make a likelihood update or choose a new plan for that patient. Propagate the belief under known actions and hold treatment until the next usable scan; record the missed scan and the resulting no-treatment interval. Use `evaluate_robust_belief_mpc` so the simulator receives the policy's configured scan interval and course budget. Give every comparator the same missing-reading pattern.

## Baseline set

Run all applicable policies on the same evaluation patients and score them under both endpoints:

- Clinical schedule arms: six-course standard TMZ, dose-dense exposure proxy using the existing assumed multiplier, and 12-versus-6 course comparison at the cycle-six landmark. The dose-dense multiplier is not a verified pharmacologic equivalence.
- Open-loop schedules: front-loaded and paced schedules, the 56-schedule selection arm, and the fixed frontier at equal mean courses.
- Scan rules: Zhang/Gatenby threshold rule, raw threshold rule, nadir containment, and response-guided stopping. Label the Zhang/Gatenby rule as an adaptation of the published PSA policy to simulated burden readings.
- State-estimation rules: current filtered threshold and risk-containment policies from `src/trackc_sim.py`.
- Learning baseline: adapted Q-learning inspired by Yauney and Shah's lower-grade-glioma paper; label it an adaptation, not an exact replication. Train on selection patients only. Do not use true patient parameters. The endpoint-specific legacy trainer remains available, but the integrated comparison will use `train_q_learning_shared_endpoints`: one simulated treatment path, treatment stopping at the first simulated event under either endpoint, and a terminal reward based on the earlier restricted progression time minus course cost through that event or the restriction horizon. Its evaluator will score both endpoint times and the common course count from that same path. The shared-endpoint path has contract tests but is not yet wired into a benchmark or outcome-tested.
- Oracle references: exact-class and known-parameter policies, if included, are not buildable methods and cannot establish a deployable win.

Response criteria papers define outcome measurements, not dosing algorithms. Digital-twin and reaction-diffusion forecasting papers are not directly comparable policy arms unless they specify a treatment action rule that can be implemented without inventing one.

## Win rule and limits

A candidate counts as a development signal only if it beats the full fixed-schedule frontier at its own mean course use under both progression definitions at 730 days, has a paired world-level 95% interval above zero, and does not worsen the paired worst-decile outcome against the standard arm. Report 365-day results and all simple-policy comparisons even if this rule fails. Do not choose sensitivity, specificity, priors, or cost after seeing evaluation results.

This is a development rule, not a clinical threshold. A pass on assumed worlds would not show patient benefit, validate the simulator, establish MGMT as the marker, or prove novelty against every method in the field. The literature review is not exhaustive, and the policy set covers methods that can be represented in this simulator without changing their scientific meaning.

## Implementation update (2026-11-06; no outcome evaluation)

`src/207_trackc_qlearning_baseline.py` now includes a separate shared-endpoint trainer and evaluator. The trainer receives the caller's selection scenarios and uses an objective that matches the project's single-path scoring rule. It does not change the endpoint-specific legacy functions. The focused Q-learning, MPC, signal-conditioning, and paired-signal contract suites passed 65 tests in a scratch copy after this addition, including a small real-simulator mechanics check. No treatment-policy outcomes were compared. Next, wire fixed schedules, Q-learning, no-signal MPC, conditional-independence MPC, and dependence-aware signal MPC into one patient-paired development harness; freeze its estimand, split, training settings, uncertainty summaries, and stop rules; then register fresh development seeds before any outcome run.
