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
