# Track C amendment: particle-filter prior support (revision 3, 2026-11-10)

**Status: ADOPTED for the future integrated benchmark only, 2026-11-10,** under the user's standing delegation to Codex's decisions (relayed by the user in chat on 2026-11-10): V3 for every particle-filter arm in the integrated benchmark; the section 4 diagnostic runs report-only; legacy defaults and historical results are preserved; confirmatory policy evaluation stays stopped. Revision 3 changes from revision 2: world-cluster bootstrap; exact narrow-prior sensitivity; accurate pairing statement; baseline filter arms receive V3 (D3 decided). Bounds and settings below were fixed before the diagnostic ran and are not chosen from its results. Track C confirmatory evaluation stays stopped. This change does **not** clear the 14-day precision blocker (gap C4). No endpoint, threshold, seed, cohort, success rule, or the 14-day target changes. No integrated-benchmark policy outcome exists, so the proposal is made before any outcome.

## 1. Problem

The integrated generator can produce some patient values of r and B0 that the primary particle-filter prior excludes. A filter cannot place posterior mass outside its prior. This may distort parameter estimates, and through them the controllers' decisions, for the affected patients. It does **not** prove a bias in the primary policy contrast, and its direction and size on utility are unknown. Both arms of the primary contrast use the same prior (section 3), so any effect may partly cancel in the paired difference; that is also unverified.

Observed counts, one synthetic batch only (100 worlds per family, 12 patients each, world seeds 9,200,000+i, cohort seeds 9,300,000+i, i = 0..99; `tests/test_trackc_prior_support.py::test_fraction_outside_prior_is_recorded`): 246 of 2,400 patients outside the r prior and 324 of 2,400 outside the B0 prior. These are batch counts, not exact population rates. The batch compares against the config prior, which equals the integrated benchmark prior for r and B0.

## 2. Support table (verified against source)

Filter parameters are `PARAMS = ("r", "m", "f0", "c", "eps", "B0")` (`src/trackc_sim.py:39`). Log-scale parameters are reflected at the bounds (`src/trackc_sim.py:190-196`). Linear parameters pass through a bounded logit with z clipped to ±20 (`src/trackc_sim.py:197-199`), so they also cannot leave the bounds. (Revision 1 wrongly said all parameters are reflected.)

Current integrated prior = `INTEGRATED_BELIEF_RANGES` (`src/208_trackc_robust_belief_mpc.py:35-38`): `tc.RANGES` with m replaced by (0.1, 20) and f0 replaced by `F0_REFERENCE_BOUNDS` (`src/trackc_integrated_signal_model.py:18`). The config block `particle_filter_prior` (`configs/trackABC_integrated_v2.json:228-252`) states the same values.

| Param | Scale | Generator support (both families) | Derivation (`src/trackc_integrated_sampler.py`) | Current prior | Proposed V3 prior | Covered now? |
|---|---|---|---|---|---|---|
| r (1/day) | log | 0.000497 – 0.1230 | base 0.001–0.05 (`trackc_sim.py:42`) × e^log_r_shift, shift in [−0.7, 0.9] (lines 28, 34, 86, 106) | 0.001 – 0.05 | 0.00049 – 0.124 | No |
| m (k = m·r) | log | 0.12 – 20 | single: 6–20 × e^log_m_scale, scale in [ln 0.05, 0] → 0.3–20; mixture nonresponders: 6–20 × [0.02, 0.30] → 0.12–6; mixture responders: 6–20 (lines 29, 35, 96-101, 124-133) | 0.1 – 20 | 0.1 – 20 (unchanged) | Yes |
| f0 | log | 0.000497 – 0.9 | 0.001–0.1 × e^log_f0_shift, shift in [−0.7, 2.3], capped at 0.9 (lines 30, 36, 88-89) | 0.000497 – 0.9 | unchanged | Yes |
| B0 | linear | 0.02 – 0.65 | B0_center in [0.12, 0.55] ± 0.10 (lines 25, 31, 37, 92-93) | 0.05 – 0.5 | 0.02 – 0.65 | No |
| c | linear | 0 – 0.25 | `tc.RANGES["c"]`, no world shift (line 110) | 0 – 0.25 | unchanged | Yes |
| eps (epsilon) | linear | 0 – 0.2 | `tc.RANGES["eps"]`, no world shift (line 111) | 0 – 0.2 | unchanged | Yes |

Notes. The V3 r bounds round the analytic support outward. For linear parameters at an exact bound the logit maps to ±infinity and is clipped (`trackc_sim.py:183`, q clipped to [1e-9, 1−1e-9]); that is an existing edge behavior, not changed here.

World-level variables that are **not** particle-filter parameters: family, log_r_shift, log_m_scale / log_m_scale_non, log_f0_shift, B0_center, pi_response, mu, nu, q_signal (`trackc_integrated_sampler.py:26-40, 55-68`); per-patient sigma (assigned, treated as known: `ParticleFilter(scen["sigma"], ...)`, `208_...:172-173`); a_sr and a_rs (fixed at 1, sampler lines 137-138); the mixture responder flag (not returned to policies, sampler docstring lines 73-78).

## 3. Which arms get V3

Traced in `src/211_trackc_policy_benchmark.py`:

- `common_mpc["ranges"] = mpc.INTEGRATED_BELIEF_RANGES` (line 993) is passed to every MPC arm built in `make_signal_policies` (lines 995-1026):
  - `finite_set_minimax_mpc` and `finite_set_d0_minimax_mpc`, both `SignalAmbiguityBeliefMPC` (`216_...:112`) built with `**common_mpc` (lines 1014-1025). **These two are the planned primary estimand** (config `primary_comparator` line 276, `primary_contrast` line 305). They must receive V3.
  - `loo_dependence_aware_mpc`, `loo_conditional_independence_mpc` (`214_...`; lines 996-1005) and `no_signal_mpc` (`RobustBeliefMPC`, line 1006). Secondary or control arms. They should also receive V3 so that every MPC arm in the same world uses the same prior.
- `SignalConditionedBeliefMPC` (`208_...:432`) defaults to `INTEGRATED_BELIEF_RANGES` (line 471) but is **not** instantiated in the integrated benchmark path traced above. It is not part of the planned estimand.
- `FilteredThreshold`, `RiskContainment`, and `HolidayContainment` baseline arms (lines 455-535 region) are built by `factory()` with no `ranges` argument, so they use `tc.RANGES` (r 0.001–0.05, m 6–20, f0 0.001–0.1, B0 0.05–0.5; `trackc_sim.py:347, 383`). They are established comparators, not the primary estimand. **D3 decided (2026-11-10):** in the integrated benchmark they receive V3 explicitly. Their historical configuration (legacy `tc.RANGES`, or `RANGES_WIDE_M` for `wide=True`) stays the default everywhere else.
- `q_learning` and `FixedSchedulePolicy` arms use no particle filter.

Historical defaults are preserved: `tc.RANGES`, `BELIEF_RANGES`, `INTEGRATED_BELIEF_RANGES`, and the `SignalAmbiguityBeliefMPC` default (`216_...:143-148`) are not edited. V3 would be a new named constant passed explicitly. The narrow-prior sensitivity `PF_PRIOR_NARROW_V2` is defined exactly as `INTEGRATED_BELIEF_RANGES`: m (0.1, 20, log), f0 (0.000497, 0.9, log), c (0, 0.25, lin), eps (0, 0.2, lin) identical to V3, with only r narrowed to (0.001, 0.05, log) and B0 to (0.05, 0.5, lin). It is **not** legacy `tc.RANGES` (which also has m 6-20 and f0 0.001-0.1). It is never the primary analysis.

V3 exactly: r (0.00049, 0.124, log), m (0.1, 20, log), f0 (0.001·e^-0.7, 0.9, log), c (0, 0.25, lin), eps (0, 0.2, lin), B0 (0.02, 0.65, lin).

## 4. Frozen outcome-blind parameter-recovery diagnostic (frozen 2026-11-10, before running)

Purpose: describe how the filter's parameter recovery changes between the narrow-prior sensitivity (`PF_PRIOR_NARROW_V2`) and V3. It runs no controller, computes no utility, and compares no policies. Script: `src/221_trackc_prior_recovery_diagnostic.py` (first written as `219_...`; renamed after the run to avoid a number collision with `219_trackc_minimax_action_stability.py`; content unchanged). Output: `output/trackc_prior_recovery_diagnostic_20261110/diagnostic.json`.

- **Generator:** `trackc_integrated_sampler.draw_world_parameters` and `draw_patient_cohort`, both families.
- **Size:** 50 worlds per family, 4 patients per world, 400 patients total. Fixed; not changed after results.
- **Seed allocation (recorded before the run; `check_free` returned no conflict on 2026-11-10):**
  - world seeds 20263000 + i, i = 0..99 (single i = 0..49, mixture i = 50..99);
  - cohort seeds 20263100 + i;
  - filter seeds (policy stream) 20263200 + i, one per world, shared by both priors;
  - observation-noise keys 2630000 + i, one per world;
  - world-cluster bootstrap key 2630100.
  - A `dev_used` registry entry is added after the run with the output as evidence.
- **Sigma:** levels (0.05, 0.10, 0.20) assigned by sorted patient index within each world, cycling in that order (4 patients: 0.05, 0.10, 0.20, 0.05).
- **Fixed treatment and scans (no controller):** 28-day cycles; a 5-day course in cycles 0, 2, 4, 6, 8, 10 (6 courses); filter updates at cycles 0, 3, 6, 9, 12 (days 0, 84, 168, 252, 336); stop after the cycle-12 update. Observation y = log(true burden) + sigma·N(0,1), matching the filter likelihood (`trackc_sim.py:236-237`). No scan missingness. Imaging burden only; no molecular or MRI-like signal channel.
- **Pairing (accurate statement):** each patient's true parameters, treatment schedule, and observation values are identical under both priors. Each world's filter starts from the same seed under both priors, so the initial uniform draws are the same, but they map to different parameter values because the bounds differ. After the first resampling the two filters consume random numbers in different amounts and orders (resampling happens only for filters whose ESS drops below P/2), so later random draws are **not** guaranteed identical. The pairing is therefore common truth and common data, not common filter randomness.
- **Particles:** 500 (`PolicyBenchmarkSettings.mpc_particles` default), h = 0.90, step 0.25 days (filter defaults).
- **Interval:** weighted 5% and 95% quantiles of each marginal particle distribution after the cycle-12 update (equal-tailed rule as in `ParticleFilter.fr_interval`).
- **Summaries:** for each of the six parameters, per family and per prior: marginal 90% coverage of the truth; joint coverage (all six covered); median interval width (log width for log parameters, raw width for linear). Each summary also split by whether the truth lies inside or outside the narrow prior's r and B0 bounds. Filter health: ESS 1/Σw² recorded at each update before resampling (minimum and median per patient), resampling count, wall-clock time per prior.
- **Uncertainty:** paired differences V3 minus narrow, with a **world-cluster bootstrap**: within each family, resample the 50 worlds with replacement, keeping each sampled world's 4 patients together; 2,000 resamples; key 2630100; percentile 95% intervals.
- **Reporting rule: report-only.** No defensible pass threshold exists: a 500-particle filter with five scans is not expected to reach nominal 90% coverage under either prior, and no threshold was declared before now. The diagnostic cannot change the bounds, which follow from the generator support in section 2. All results are reported, including any loss of precision under V3.
- **Not allowed:** choosing bounds, particle counts, or schedules from the results; running policies or utilities.

## 5. Implementation checks

- The two strict expected failures in `tests/test_trackc_prior_support.py` flip to passes; remove their `xfail` marks in the same commit.
- Update the config `particle_filter_prior` with a `version` field; move the old values under `sensitivity_priors`. Update `tests/test_trackc_config_parity_extended.py`.
- A test that every MPC arm in `make_signal_policies` receives V3.

## 6. Unchanged

Endpoint, utility, stop rules, world sampling law, seeds, controllers, success rule, 14-day target, GE-betting candidate status, and the stop on confirmatory evaluation.

## 7. Decision record

- Revision 1 drafted 2026-11-10 at user request. Revision 2 after Codex review, 2026-11-10.
- Revision 3 and adoption for the future integrated benchmark: 2026-11-10, user standing delegation to Codex decisions (D1 V3 all filter arms; D2 diagnostic report-only; D3 baseline filter arms get V3).
