# Track B compute audit (2026-11-02)

Label for every projection below: **projection, not a bound; this PC (8 CPU threads, Windows, Python 3.14.4); synthetic full-grid cases only.** No patient data was read. This is not a Duke estimate: Duke grid, gaps, eligible N and central-gap subset are unknown.

Sources: `src/fk_lumiere_protocol.py`, `src/fk_lumiere_runner.py`, `src/fk_lumiere_pde_provider.py`, `src/fk_hybrid_forecast.py`, `src/fk_lumiere_stats.py`; timing records in `output/track_b_timing_pilot_20261031/` (`lumiere_timing_record.json`, `lumiere_history_gap_stress.json`, `film_unet_timing_record.json`, `timed_fixed_growth_timing_record.json`) and `output/track_b_pde_only_grid_timing_20261103.json`; new run `output/track_b_feasibility_pilot_20261102/pilot.json` (script `src/fk_track_b_feasibility_pilot.py`, run once in a scratch copy, 694 s).

## A. Call-count formulas (from the code)

Notation: N eligible patients, N_c patients in the central-gap subset (forecast gap 12-14 week labels), K outer folds (5), I inner calibration folds (5), P prior scenarios (3), R central-gap refits (1), E epochs (300), S U-Net seeds (3), B bootstrap resamples (10,000), G fixed-growth rate-grid size (20).

Per analysis on m patients (one `run_lumiere_cv`); the sum of outer-training sizes over folds is (K-1)m:

| Item | Count | Code basis |
|---|---|---|
| Full forecasts (one per held-out patient) | m | `pde_provider.forecast`, once per test patient |
| Inner calibration-only forecasts | (K-1) m | `fit_fold` -> `run_cross_fitted_boundary_calibration`: each outer-training patient is forecast once, held out of its inner fold. Independent of I. Upper bound: only training patients with three non-empty, non-full masks are used. |
| PDE solver-posterior calls total | K m | sum of the two rows; each call = coarse grid n_coarse^2 + refine grid + up to n_draws unique solver fields |
| Ell-candidate scoring | 0 extra solver calls | scored from stored pre-boundary probabilities |
| Volume-nuisance fits | K (outer) + K I (inner) | cheap, not timed |
| FiLM U-Net optimiser steps per training set | (K-1) m E S | `train_film_models`: every training pair, every epoch, every seed |
| FiLM models per training set | K S | |
| Fixed-growth forecasts | (K-1) m G + m | `select_rate` per fold plus test forecasts |
| Bootstrap resamples | B | `max_t_pass_test`, once per analysis |

Whole protocol (`run_lumiere_protocol`):

- PDE calls = K (P N + R N_c). The three prior scenarios use different bounds, so no PDE reuse; the refit re-solves everything on its subset.
- FiLM steps = (K-1) E S (N + R N_c) with the existing shared cache across the three prior runs (the fingerprint does not include the provider); without the cache, replace N by P N. The central-gap refit has its own cache, so it trains again.
- Fixed-growth forecasts = (K-1) G (P N + R N_c) + P N + R N_c; bootstrap resamples = (P+R) B.
- Example N=28, N_c=22 (old LUMIERE scale): 15(28)+5(22)=530 PDE calls (full 3(28)+22=106, inner calibration-only 424). The old projection used 4 x 140 = 560 because it counted the central-gap run at 28. Consistent with the plan's note.

What the formulas do not include: a refit of U-Net after any change, ablation arms beyond the listed ones (the matched-volume ablation and PDE-only arm share the one PDE call), PDE failures (they remove full forecasts, not inner calls), count-model fitting and file I/O.

## B. Measured pilot (`pilot.json`)

Full grid 91 x 109 x 91 at 2 mm. A hard 25-minute budget with per-item child-process timeouts was enforced; nothing was censored.

| Item | Measured | Note |
|---|---|---|
| FiLM U-Net training step (padded 112x112x96) | **2.99 s** | difference of 3-pair and 1-pair epoch-1 trainings; the record `film_unet_timing_record.json` measured 4.14 s mean step, and the plan notes ~2x session variation |
| One FiLM prediction (single model) | 1.0 s | ~3 s for three seeds |
| Fixed-growth rate selection, 22 training patients | 36.0 s (1.64 s per training patient) | older record 19-24 s: the same ~2x variation |
| Bootstrap B=10,000, 40 patients | 0.034 s | negligible |
| PDE call, **reduced** settings (n_coarse 5, n_refine 3, n_draws 20, n_anchors 64), history 6 wk / forecast 13 wk | full 99.2 s; calibration-only 93.8 s | 38 solver fields counted |
| Same, history 23 wk / forecast 13 wk (extreme) | full 189.1 s; calibration-only 211.6 s | 38 solver fields |

Scaling the reduced PDE call to the frozen full settings multiplies by 762/38 = 20.05 (762 = 21^2 + 11^2 + 200, an upper-ish field count). The scaled numbers are therefore an **upper-leaning scenario**, not a central estimate. Implied full-settings call: about 1.99 ks (representative, full forecast) and 3.79 ks (extreme). The recorded full-settings calls were 1.25 ks (full) and 1.54 ks (calibration-only) at a 30-day history gap (`lumiere_timing_record.json`), and 3.18 ks / 2.52 ks at 23 weeks (`lumiere_history_gap_stress.json`). The scaled values are 1.2x-1.7x above the recorded ones (they overestimate), so the scaling is order-of-magnitude only; the recorded-time column is the **lower scenario**. Both are projections, not bounds.

## Projection table (hours, serial, this PC)

Assumptions: K=5, I=5, P=3, R=1, E=300, S=3, B=10,000, N_c = 0.8 N, FiLM reused across the three prior runs. PDE hours use all calls at one gap. "Pilot" = scaled pilot numbers (FiLM 2.99 s/step). "Recorded" = old records (PDE 1250.4 s full / 1544.3 s calibration-only; FiLM 4.14 s/step).

| N | PDE calls | FiLM steps | Pilot, representative gap | Pilot, extreme (23 wk) gap | Recorded-time version |
|---|---|---|---|---|---|
| 20 | 380 | 129,600 | 308 h | 546 h | 306 h |
| 30 | 570 | 194,400 | 463 h | 819 h | 459 h |
| 40 | 760 | 259,200 | 617 h | 1,092 h | 612 h |
| 61 | 1,160 | 396,000 | 942 h | 1,667 h | 935 h |

Components at N=61 (pilot, representative): PDE 613 h, FiLM 328 h, fixed growth 0.4 h, bootstrap ~0. At N=20: PDE 201 h, FiLM 107 h. Compute is dominated by the PDE posterior calls, then FiLM training. With the recorded step (4.14 s) FiLM is 456 h at N=61.

The old LUMIERE four-analysis scenario of 557-931 h (N=28) lies inside this range, as it should.

## Weakest points

1. PDE timings use reduced posterior settings scaled by solver-field count. Cost per field depends on D (CFL step count) and horizon, so the 20x factor is approximate; the reduced grid includes the endpoints of the D range (including D=10), which may over-weight expensive fields. Cross-check against recorded full-settings calls: scaled values run 1.2x-1.7x higher.
2. One synthetic tumor shape per gap, two gaps. Real tumors, real history gaps (the stress record shows 155-week gaps costing 16-19.5 ks per call, 5-8x the 23-week case; the 155-week gap was not timed in the new pilot and is not in the projection table) and real forecast gaps are not covered; a long-gap patient dominates cost. No Duke gap distribution is assumed.
3. A single run on this PC. FiLM step and fixed-growth timings varied about 2x versus earlier sessions. These projections can be off by that factor in either direction.
4. Serial execution assumed; the PDE batch runs and torch use 8 threads and compete if run together. No parallel or cached path is validated.
5. The call-count formulas give an upper bound on inner calibration calls (patients with unusable masks reduce them) and ignore reruns after crashes.
6. N_c = 0.8 N is a placeholder from the LUMIERE 22/28 ratio.
7. Not timed: count model, PDE-only cutoff selection, FiLM prediction (about 1 s per model), I/O, input hashing. Each is small relative to the totals.
8. This audit does not decide feasibility; the gate-4 budget is a human decision. If the approved budget is below these figures, scope must be reduced by a dated amendment before outcomes, without weakening only the strongest baseline.
