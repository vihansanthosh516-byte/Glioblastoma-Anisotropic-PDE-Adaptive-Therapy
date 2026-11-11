# Track C amendment draft: particle-filter prior support (2026-11-10)

**Status: DRAFT, NOT ADOPTED.** The user asked for this draft on 2026-11-10. It takes effect only after the user approves it in writing in a later dated note. No policy outcome exists for the integrated benchmark, so this change is made before any outcome is seen. No endpoint, threshold, seed, cohort, success rule, or the 14-day target changes.

## 1. Problem (outcome-blind evidence)

- The integrated generator (`src/trackc_integrated_sampler.py`) scales the base ranges by world shifts. Analytic support: r 0.000497 to 0.12298 per day (0.001·e^-0.7 to 0.05·e^0.9); B0 0.02 to 0.65 (centre 0.12 to 0.55, ±0.10).
- The primary particle-filter prior (`configs/trackABC_integrated_v2.json` `particle_filter_prior`; `src/trackc_sim.py` `RANGES`) is r 0.001 to 0.05 and B0 0.05 to 0.5.
- Synthetic draw, 2400 patients: 246 (10.3%) outside the r prior, 324 (13.5%) outside the B0 prior (`paper/status_map_2026-11-10.md`; `tests/test_trackc_prior_support.py`).
- The filter reflects particles at the prior bounds (`src/trackc_sim.py` line 194). It cannot reach true values outside them. That biases the primary contrast by construction.

## 2. Proposed change

1. New primary prior, named `PF_PRIOR_INTEGRATED_V3`:
   - r: log-uniform [0.00049, 0.124] per day (generator support rounded outward).
   - B0: uniform [0.02, 0.65].
   - m: log-uniform [0.1, 20] (unchanged; already covers 0.12 to 20).
   - f0: log-uniform [0.001·e^-0.7, 0.9] (unchanged).
   - c: uniform [0, 0.25]; epsilon: uniform [0, 0.2] (unchanged).
2. Keep `trackc_sim.RANGES` unchanged. It is also the generator base range and is used by historical sealed scripts (150-171, 208, 211). Editing it would change both generator and old results. The new prior is passed explicitly through the existing `ranges=` arguments.
3. The old prior (r 0.001-0.05, B0 0.05-0.5) stays available only as a separately named misspecification sensitivity, `PF_PRIOR_NARROW_V2`. It is never the primary analysis.
4. The config key `particle_filter_prior` gets the new values and a `version` field; the old values move under `sensitivity_priors`.

## 3. Checks required before adoption counts as implemented

- `tests/test_trackc_prior_support.py`: the two strict expected failures must flip to passes. Remove their `xfail` marks in the same commit.
- Config parity test (`tests/test_trackc_config_parity_extended.py`) updated to the new prior.
- Outcome-blind filter check on synthetic patients only: coverage of the true r and B0 by the filter's 90% interval, old vs new prior. No policy, no utility. Report also the runtime change, if any.
- Known risk: a wider prior spreads 500 particles over a larger space. Filter accuracy per particle may drop. The check above must report it; the particle count is not changed by this amendment.

## 4. What stays unchanged

Endpoint, utility, stop rules, world sampling law, seeds, controllers, success rule, 14-day half-width target, GE-betting candidate status, and the stop on confirmatory evaluation (gap C4).

## 5. Decision record

- Drafted: 2026-11-10, at user request.
- Approved: _pending user_.
