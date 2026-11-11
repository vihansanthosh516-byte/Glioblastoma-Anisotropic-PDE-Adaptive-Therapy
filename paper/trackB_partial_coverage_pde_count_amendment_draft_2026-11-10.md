# Track B amendment draft: PDE and count arms under partial scan coverage (2026-11-10)

**Status: DRAFT, NOT ADOPTED.** Must be decided before any patient data is used with partial coverage. Until then the runner refuses to fit or score these arms when coverage is partial (fail closed).

## Problem

Under the coverage amendment (`paper/trackB_coverage_scoring_amendment_2026-11-10.md`), target voxels outside the common scan coverage are unobserved. Three arms use target masks in ways that do not have a defined meaning for a partly observed target:

1. **PDE nuisance** (`src/fk_lumiere_pde_provider.py::fit_fold`): `fit_mu_training_volume_nuisance` fits on whole target volumes; a truncated lesion has an unknown volume, so the fitted nuisance would change estimand.
2. **PDE boundary calibration** (`run_cross_fitted_boundary_calibration`): residuals are measured on the target boundary; a boundary cut by the coverage edge is not an observed tumor boundary.
3. **Count arm** (`new_label_pieces`, `fit_count_model`, `count_metrics`): it counts 26-connected new pieces. The coverage edge can split, hide, or create pieces, so the count is a different quantity.

The PDE-only, candidate, and matched-volume ablation arms depend on (1)-(2).

## Current implementation (fail closed, no estimand change)

- Full coverage (common support covers the whole forecast domain for every eligible patient): every arm runs as before.
- Partial coverage for any eligible patient: the PDE provider is not fitted; the candidate, PDE-only, and matched-volume ablation arms are reported as not scored (None) with the reason; the count arm is reported as not scored. Copy-last, time-aware growth, and FiLM U-Net run with masked selection and loss. The comparison is then marked incomplete.

## Options to decide (not chosen here)

- A. Exclude patients with partial coverage from all arms (changes the cohort; must be frozen before targets are read).
- B. Keep partial patients for mask arms only and report PDE/count arms on the full-coverage subset (two cohorts; needs a declared primary).
- C. Define censored-volume and censored-count versions of the PDE nuisance and count model (new estimands; needs method work and validation before use).

No option is adopted. No patient data may be scored with partial coverage until one is.
