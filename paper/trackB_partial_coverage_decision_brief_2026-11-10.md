# Decision brief: PDE and count arms under partial scan coverage (2026-11-10)

**Status: brief only. No option adopted. No code changed.** Source: `paper/trackB_partial_coverage_pde_count_amendment_draft_2026-11-10.md`. Current code fails closed: with any partial-coverage patient, the candidate, PDE-only, ablation and count arms are not scored and the result is `incomplete` (tests: `tests/test_track_b_validity_masks.py`).

**Unknown that drives feasibility:** how many LUMIERE patients have partial coverage. This has not been counted.

**Correction (2026-11-10).** An earlier version said the count needs only headers and transforms, with no image intensities read. That was wrong. The coverage mask and the count use geometry only. But the exact target-to-origin transform they need does not exist yet:

- It is produced by `src/fk_longitudinal_prep.py::register_rigid_mri`, which uses visit-2 and visit-3 CT1 MRI intensities (no label input).
- No LUMIERE transform is recorded. `output/track_b_lumiere_inventory_20261031/header_check.json` stores no affine matrices (keys: source, read, files_per_visit, n_eligible, n_excluded, grid_counts, problems, patients). `src/fk_nifti_header_audit.py` reads affines to check agreement but does not save them.
- `output/cfb_registration_qc_20261028/registration_qc.json` is a different cohort (CFB-GBM) and stores no transform matrices.

So the transform must first come from the approved MRI-only registration step. Tumor-label files stay closed throughout.

## Option A: exclude partial-coverage patients from all arms

- **Cohort:** one cohort. Only patients whose common support covers the whole forecast domain. Smaller than now.
- **Common comparisons:** every arm is scored on the same patients. Candidate vs. every comparator stays paired.
- **Estimand:** unchanged (frozen PDE nuisance, boundary calibration, count arm, Dice). Population changes: "patients with full scan coverage".
- **Code/validation:** small. Add a geometry-only eligibility filter that runs before any target read, plus tests that it reads no labels. Existing full-coverage path is already tested.
- **Feasibility:** high, if enough patients remain. The frozen stop rule (`minimum` eligible) still applies; if too few remain, the run stops with no scoring.
- **Weakest point:** possible selection bias. Full coverage may relate to scan protocol or lesion location. Results do not transfer to partially covered scans.

## Option B: two cohorts

- **Cohort:** all patients for mask arms (copy-last, time-aware growth, FiLM); full-coverage subset for PDE, candidate, ablation and count arms.
- **Common comparisons:** **different patient sets.** Candidate (subset) cannot be compared directly with mask-arm scores on the full cohort. A fair candidate comparison needs mask arms rescored on the subset as well, so the subset comparison is in effect Option A plus an extra report.
- **Estimand:** unchanged per arm, but two populations. Needs a declared primary cohort before targets are read.
- **Code/validation:** moderate. Two scoring passes, two sets of statistics and reports, tests that cohorts never mix in one comparison.
- **Feasibility:** medium. Most work is reporting. Easy to misread; a reviewer may quote the wrong cohort.
- **Weakest point:** risk of unpaired comparisons and of choosing the more favourable cohort after the fact. Must be frozen first.

## Option C: censored estimands

- **Cohort:** all patients, all arms.
- **Common comparisons:** one cohort, all paired.
- **Estimand:** **new.** Censored-volume PDE nuisance, censored boundary calibration, censored count of new pieces. These are not the frozen quantities.
- **Code/validation:** large. New likelihoods for truncated lesions, new boundary handling at the coverage edge, a new count model for pieces split or hidden by the edge. Each needs synthetic validation (known truth, recovery checks) before any patient use. Prior-art check not done.
- **Feasibility:** low in the near term. Weeks of method work; the result would be a new method, not the frozen protocol.
- **Weakest point:** many new modelling choices made after the protocol was frozen. High risk of an unvalidated estimand.

## Comparison

| | A | B | C |
|---|---|---|---|
| Patient sets | one (full coverage) | two (different sets) | one (all) |
| Paired candidate comparison | yes | only on subset | yes |
| Estimand change | none | none | yes, new |
| Code | small | moderate | large |
| Feasibility | high, if n ≥ minimum | medium | low |

## Recommendation

**Option A.** It keeps every frozen estimand, keeps one cohort, and keeps all comparisons paired. It is the smallest change. Its cost (smaller, possibly biased population) is stated, not hidden. Option C can be kept as later, separately frozen method work.

Option A is conditional on meeting the frozen minimum (`MIN_ELIGIBLE_PATIENTS = 25`, `src/fk_lumiere_runner.py:34`).

## Decision sequence for Option A

1. Establish each patient's exact transform with the approved MRI-only registration step (`register_rigid_mri`, CT1 intensities only). Record full matrices and QC. No tumor label is opened.
2. Build coverage masks from those transforms and the grids (geometry only). Count patients whose common support covers the whole forecast domain. No tumor label is opened.
3. Apply the frozen 25-patient minimum. If fewer than 25 remain, the run stops with no scoring, and the option choice is reopened.
4. If 25 or more remain, freeze the cohort (patient list and transforms, hashed) before any target-mask read.

Do not switch options after seeing any target.

## Decision needed from you

Reply with one of:

1. "Adopt Option A, conditional on the frozen 25-patient minimum: exclude patients whose geometry-derived common support does not cover the whole forecast domain, from all arms. Approve steps 1-4 above: MRI-only registration to record transforms, then a geometry-only coverage count with tumor labels closed, then the minimum check, then a cohort freeze before any target-mask read." (recommended)
2. "Adopt Option B with primary cohort = [all patients / full-coverage subset]."

Option C would need a separate method plan; it is not a short decision.
