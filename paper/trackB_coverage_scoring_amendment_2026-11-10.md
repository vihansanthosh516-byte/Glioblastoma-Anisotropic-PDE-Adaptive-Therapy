# Track B amendment: coverage-aware scoring (frozen 2026-11-10, before code)

**Scope.** Future Track B scoring code only. Synthetic tests only. The completed Duke 10-mm result, its private records, and the historical scorers (`src/track_b_transition_estimand.py::summarize_window`, `src/fk_mu_scoring.py::dice_with_empty_rule`) are unchanged and are not rerun. No patient image or mask is read. Fixes gap B4 (outside-field voxels counted as non-converted) for any future, different cohort.

## Rules

1. **Support source.** Observable support comes only from scan coverage metadata: the grid shapes and the exact output-to-input transform used to resample each visit (`geometric_coverage`). Tumor labels, intensities, or predictions are never used to build support. Code enforces this by accepting only a `ScanCoverage` object built from geometry.
2. **Multi-visit transitions.** A location is scored only if it is covered at every visit the quantity uses (previous, middle, and target for the transition estimand; every history visit and the target for a forecast). Common support = logical AND of the per-visit coverage masks.
3. **Same support for every method.** All forecast methods in a window are scored on one common support. A method's prediction outside support is ignored, never counted as a false positive.
4. **Outside support is unobserved.** It is neither negative nor positive. For the transition estimand, risk voxels outside support are removed before band eligibility and counting.
5. **Fail closed.** Scoring raises if any required visit lacks a coverage mask, masks are not boolean, shapes differ, the common support is empty, or (transition) no band remains eligible after restriction is reported as ineligible (`None`), not as zero.
6. **Empty-mask rule** on support is unchanged: Dice 1 if prediction and target are both empty on support, 0 if exactly one is.
7. **Reporting.** Every scored window reports the support voxel count and the fraction of the target grid it covers.

## Tests required

Partial coverage; lesions outside coverage; one common support across methods; labels or predictions changed outside support cannot change any score; missing or invalid coverage fails closed.

## Implementation note (2026-11-10, no rule change)

Integrated into the future forecast scoring path `src/fk_lumiere_runner.py::run_lumiere_cv` (required comparators copy-last, time-aware growth, PDE-only, FiLM U-Net, candidate, ablation). Coverage is checked before any target read.

**Correction (2026-11-10).** An earlier version of this note said targets were "restricted" to support. The code did that by setting unobserved target voxels to zero, and then used them as negatives in training and selection. That broke rule 4. It is fixed:

- Targets are no longer zeroed. Each target carries a geometry-derived validity mask (the common support).
- Growth-rate selection scores Dice only on valid voxels (forecast and target both masked).
- FiLM training uses a masked Dice + BCE loss. Invalid voxels have weight 0 and add nothing. The mask is part of the reuse fingerprint. Under full coverage the frozen training path is used unchanged; a test on a padded grid shows identical trained weights. Under partial coverage, padding voxels keep weight 1 as frozen background; only unobserved voxels get weight 0. When partial coverage leaves arms unscored, the result status is "incomplete: partial coverage ..." and `comparison_complete` is false; the PDE provider's fit and forecast are never called (tested).
- Final Dice is scored on support for every method (predictions and target both masked).
- The PDE nuisance fit, PDE boundary calibration, PDE-only cutoff and count arm have no defined estimand for partly observed targets. When any eligible patient has partial coverage, these arms are not fitted or scored (fail closed, no fallback). See `paper/trackB_partial_coverage_pde_count_amendment_draft_2026-11-10.md` (draft, not adopted). Under full coverage they run unchanged.
- Forecast domains and history inputs are unchanged.

Synthetic tests: `tests/test_fk_lumiere_coverage_scoring.py` and `tests/test_track_b_validity_masks.py`. Changing labels or predictions only outside support (supported labels fixed) does not change the masked loss, FiLM weights, selection scores, selected rate, fold records, count results or final Dice. The MU/Duke runner and the Duke result are not changed.
