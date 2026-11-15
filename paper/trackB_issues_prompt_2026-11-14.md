# Prompt: help fix Track B (glioblastoma MRI tumor-volume forecasting)

You are reviewing a research track that has failed every locked test. Find the real cause and propose a fix. Do not change the question. Do not suggest tuning on a test set. Every number below cites a file in the repo.

## The locked question (cannot change)
From a post-treatment MRI scan (plus scan dates and treatment dates), predict the tumor-core volume at the NEXT scan.
- Target: y = log(core voxels + 1), core = mask labels 1 (necrosis) + 3 (enhancing).
- Score: mean absolute error (MAE) of y, averaged per patient first.
- Must beat ALL standard methods: persistence (last size), population growth rate, two-point exponential, empirical-Bayes random slope, Gompertz, gradient boosting; Fisher-KPP PDE named (an earlier MU mask test did not beat persistence).
- The project rule: the final method must be a NEW HYBRID method and must win on outside data.

## Data
- Development: MU-Glioma-Post, 334 pairs, 135 patients, 5 patient folds. Automatic (computer-made) tumor masks. Reused cohort (used before in the project).
- Outside test 1: UCSD-PTGBM, 46 pairs, 33 patients. Used once before for a different target. Now used twice (v8, then the hybrid finalist).
- Outside test 2: UCSF-ALPTDG, 298 patients, exactly 2 scans each, EXPERT manual masks. 59 patients were used before (H-3/H-4); 137 untouched patients were eligible for the main test.
- Rejected or spent: LUMIERE (volumes used in an earlier H-3v test; registration no-go), RHUH (13 usable pairs; minimum is 25), ACRIN 6677 (failed area forecast; 15 patients for sRCBV), Duke (exposed), others in `paper/trackB_fresh_cohort_screen_2026-11-11.md`.

## What was tried (all locked before scoring)
Table: `paper/trackB_mu_method_tries_log_2026-11-14.md`.
- 8 single methods (v1-v8): treatment-phase rate shrinkage, noise-aware latent size, median regression, median boosting, input-mask composition, edge shape, MRI brightness, two-part collapse model.
- 6 hybrids by parallel agents (`paper/trackB_hybrid_search_protocol_2026-11-14.md`): PK/PD Gompertz + residual learner, regime-switching state-space, gated mixture of experts, 3-compartment kinetic ODE + residual, segmentation-noise deconvolution, physics-informed stacking.
- A multiple-try correction (CI level 1 - 0.05/k) was applied.

## Results
- MU dev: best baseline pop_rate MAE 1.528. Best method v8 MAE 1.449 (diff -0.079, 99.375% CI -0.287 to +0.121; not passed). All 14 methods landed between 1.45 and 1.64. (`output/track_b_mu_v8_20261114/result.json`, `output/track_b_hybrid_*/result.json`)
- UCSD (33 patients): v8 1.219 vs pop_rate 1.441; diff -0.222, 95% CI -0.700 to +0.066; not passed. (`output/track_b_v8_ucsd_external_20261114/result.json`)
- UCSF (137 untouched patients): persistence 0.832, pop_rate 0.851, v8 1.072. v8 minus pop_rate +0.221, 95% CI +0.102 to +0.340: v8 is significantly WORSE. (`output/track_b_v8_ucsf_alptdg_20261114/result.json`)

## The core problems
1. **Noise floor.** On MU, the median absolute log change between scans is about the size of the measurement noise. The straight two-scan trend (MAE 2.72) is the worst method, so past growth predicts noise, not biology. A two-tool comparison on LUMIERE shows the tools disagree about change by more than the change itself (`vault/PAPER_FINDINGS_LEDGER.md`, entry 2ab).
2. **Mask-source shift.** MU masks are automatic; UCSF masks are expert. Baseline error is 1.53 on MU but 0.83 on UCSF. Corrections learned on MU noise hurt on clean UCSF masks. This is the most likely reason v8 flips from best (MU, UCSD) to worst (UCSF).
3. **Outliers dominate the score.** On MU, 32 pairs with log error > 3 hold 36% of total error. 11 pairs fall to zero volume (15% of error). These are likely surgeries or segmentation drop-out, and later surgery dates are not in the data.
4. **Too few outside patients where methods win.** UCSD has 33 patients, so the CI is too wide.
5. **No prior scan in UCSF.** Each UCSF patient has only 2 scans, so history-based parts of any method switch off.
6. **Missing and assumed treatment data.** UCSD and UCSF lack exact RT end and progression-before-scan dates; phase windows were assumed (RT = 42 days, adjuvant = 180 days).
7. **Reuse.** MU is reused; UCSD is now used twice; 59 UCSF patients are used. No large clean cohort with automatic masks is left.
8. **Selection pressure.** 17 locked tries now; the correction makes a pass very hard on MU.

## What I need from you
1. Is the noise floor the true limit, or is there a method that removes mask-source shift (for example, re-segmenting ALL cohorts with one fixed public tool so masks are comparable) while keeping the question unchanged?
2. Propose ONE new hybrid method that should work on both automatic and expert masks, with a locked plan: development data, outside test data (must be untouched), comparator, and success rule.
3. Say plainly if no fair positive result is possible with this data, and what data would make it possible.
4. Check every claim above against the cited files before you agree or disagree.
