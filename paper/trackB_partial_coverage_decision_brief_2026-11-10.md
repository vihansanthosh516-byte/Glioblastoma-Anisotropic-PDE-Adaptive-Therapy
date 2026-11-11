# Decision brief: PDE and count arms under partial scan coverage (2026-11-10)

**Status: brief only. No option adopted. No code changed.** Source: `paper/trackB_partial_coverage_pde_count_amendment_draft_2026-11-10.md`. Current code fails closed: with any partial-coverage patient, the candidate, PDE-only, ablation and count arms are not scored and the result is `incomplete` (tests: `tests/test_track_b_validity_masks.py`).

**Unknown that drives feasibility:** how many LUMIERE patients have partial coverage. This has not been counted.

**Correction (2026-11-10).** An earlier version said the count needs only headers and transforms, with no image intensities read. That was wrong. The coverage mask and the count use geometry only. But the exact target-to-origin transform they need does not exist yet:

- Its visit-3 rigid link is produced by `src/fk_longitudinal_prep.py::register_rigid_mri`, which uses visit-2 and visit-3 CT1 MRI intensities (no label input). This link alone is not enough; see the transform-chain audit below for the other links.
- No LUMIERE transform is recorded. `output/track_b_lumiere_inventory_20261031/header_check.json` stores no affine matrices (keys: source, read, files_per_visit, n_eligible, n_excluded, grid_counts, problems, patients). `src/fk_nifti_header_audit.py` reads affines to check agreement but does not save them.
- `output/cfb_registration_qc_20261028/registration_qc.json` is a different cohort (CFB-GBM) and stores no transform matrices.

So the transform must first come from the approved MRI-only registration step. Tumor-label files stay closed throughout.

## Transform-chain audit (2026-11-10; docs, code, headers and zip file names only)

**Result: the chain is incomplete. The coverage count stays BLOCKED.** A single visit-3-to-visit-2 rigid transform is not enough.

### What `coverage_from_geometry` needs

`coverage_from_geometry(visit, moving_shape, output_shape, output_to_input)` (`src/track_b_coverage_scoring.py:59`) needs, for each visit, one map from a voxel of the scoring grid (2 mm atlas grid, 91 x 109 x 91) to a voxel of the grid where that visit was really acquired. Coverage is "where the scanner saw". For LUMIERE this is the **native** scan grid, not the publisher's atlas grid: the atlas images are already resampled onto a full 182 x 218 x 182 grid, so coverage computed on the atlas grid alone would be trivially full and would hide real field-of-view loss.

Full chain per visit: 2 mm scoring voxel -> 1 mm atlas voxel -> (visit 3 only: rigid refinement to visit 2) -> atlas physical space -> inverse of the publisher's native-to-atlas transform -> native physical space -> native voxel (native NIfTI affine).

### Per visit

| Link | Visit 1 (history) | Visit 2 (origin) | Visit 3 (target) | Exists? |
|---|---|---|---|---|
| Native source grid (shape) | `CT1.nii.gz` (native) | same | same | Files present in zip (632 CT1 files). Native headers **not read**; `header_check.json` covers atlas files only. |
| Native affine (full 4x4) | needed | needed | needed | **Missing.** `src/fk_lumiere_headers.py` saves shape, spacing, qform/sform codes only, and reads atlas files only. |
| Publisher native-to-atlas transform | `DeepBraTumIA-segmentation/native/transformation/ct1.tfm` | same | same | Files present (599 per sequence). **Not parsed.** Type (rigid, affine or deformable), direction and coordinate convention (ITK LPS vs. NIfTI RAS) **unverified**. |
| Atlas grid (1 mm) and 2 mm scoring grid | 182 x 218 x 182 at 1 mm; 2 mm = 91 x 109 x 91 | same | same | Shape and spacing verified from headers (84 headers, one grid). Atlas affine **not saved**. Atlas-1 mm to 2 mm mapping: fixed by protocol, not yet written as code. |
| Visit-to-visit MRI registration | **None planned.** Protocol assumes visit 1 is aligned to visit 2 through the atlas only. | Fixed reference | `register_rigid_mri`, visit-3 atlas CT1 onto visit-2 atlas CT1 (`paper/track_b_no_contact_dataset_plan_2026-10-31.md` section 5b) | Visit 3: method frozen, **not run**, no transforms recorded. Visit 1: **no transform and no QC** (assumed identity in atlas space, unverified). |
| Exact `output_to_input` for coverage | composition of the rows above | composition | composition including rigid | **Not built.** No code composes the chain. |

The only saved registration output, `output/cfb_registration_qc_20261028/registration_qc.json`, is a different cohort (CFB-GBM) and stores no matrices.

### Does the proposed MRI-only read approval cover it?

No. The proposed preflight scope (`paper/track_b_readiness_audit_2026-11-01.md`) reads 84 files: visit-2 and visit-3 atlas CT1, plus visit-2 atlas brain mask, for 28 patients. It explicitly excludes visit-1 files and native-space files. So it yields only the visit-3 rigid refinement. It does **not** cover:

1. Native CT1 headers with full affines for visits 1, 2 and 3 (header-only; no voxels needed, but outside the listed scope).
2. The publisher `.tfm` files for visits 1, 2 and 3 (text files, not image voxels; outside the listed scope).
3. Any check that visit 1 is aligned with visit 2. Either an approved visit-1-to-visit-2 MRI-only registration with QC, or a documented decision that the atlas alignment is accepted for visit 1 (with its weakness stated). Neither exists.
4. Code that composes the chain and synthetic tests of it (known transforms, ITK/NIfTI convention, direction).

No tumor-label file (`seg_mask`, native seg masks, HD-GLIO segmentations) is needed for any link.

### What must exist before any count

- Approval of an extended read scope: native CT1 headers (full affines) and `ct1.tfm` for visits 1-3, plus the visit-2/3 atlas CT1 and visit-2 brain mask already proposed. Tumor-label files stay closed.
- A decision on visit 1: MRI-only rigid registration to visit 2 with the same frozen QC, or atlas-only alignment accepted as a stated assumption. Either needs a dated amendment.
- Code to parse `.tfm`, verify its type and convention, and compose the full chain; synthetic tests.
- Then the steps in "Decision sequence for Option A" below.

## Metadata audit of the chain (2026-11-10, after the audit above)

Code: `src/fk_lumiere_transform_audit.py`. Records: `output/track_b_lumiere_transform_audit_20261110/transform_audit.json` (per visit: full qform, sform and chosen affine; .tfm blocks) and `direction_check.json`. Read: 252 files for 28 patients x 3 visits, namely native `CT1.nii.gz` header (348 bytes), atlas `ct1_skull_strip.nii.gz` header (348 bytes), and `native/transformation/ct1.tfm` text. No voxel arrays, intensities or labels. The code asserts that no other file was opened.

| Item | Result | Status |
|---|---|---|
| Files present | 84/84 native CT1, 84/84 atlas CT1, 84/84 ct1.tfm | Complete |
| Atlas grid | 182 x 218 x 182, one affine for all 84 | Verified |
| Native grids | 10 shapes, 8 spacings, incl. thick stacks (224 x 256 x 24 and 448 x 512 x 24 at 6 mm slices) | Recorded. Thin stacks make real partial coverage likely. |
| Native affine | Full 4x4 saved; sform used where set, else qform; 0 qform/sform disagreements | Verified |
| .tfm format and type | ITK `CompositeTransform` with one `MatrixOffsetTransformBase_double_3_3` leaf in all 84; singular values 1.000 +/- 1e-15, so rigid | Verified; no deformable part |
| Coordinate convention | ITK physical space is LPS; converted with diag(-1,-1,1) | Standard ITK rule; code-tested on synthetic data only |
| Direction | **Ambiguous.** ITK files usually store the resampling map (fixed to moving), but which image was fixed in DeepBraTumIA is not documented in files read here. Metadata evidence: mapping the atlas centre through inverse(T) lands nearer the native grid centre in 74/84 visits (median 40.7 mm vs 80.4 mm; outliers up to 695 vs 75 mm). This favours "stored T maps native to atlas". It is evidence, not proof. | **Not frozen** |
| Fits the 4x4 `geometric_coverage` method | Yes, for every link: native affine, rigid .tfm, visit-3 rigid refinement (`rigid_voxel_map` equals `apply_rigid`, tested) and 2 mm downsampling | Representable, once direction is fixed |
| Visit-3 rigid refinement | Not run (needs the approved MRI-only read) | Missing |
| 2 mm block centre offset | `downsample_to_1mm` uses 2i + 0.5 | Protocol choice; not frozen |

Composition per visit (scoring voxel to native voxel): `inv(native_affine) @ atlas_to_native_ras @ atlas_affine [@ rigid for visit 3] @ downsample`. Visit 2 has no rigid factor. Visit 1 depends on the rule below.

Synthetic tests (`tests/test_fk_lumiere_transform_audit.py`, 7 tests, plus `tests/test_fk_longitudinal_prep.py`; 20 passed in the scratch copy): header affines match nibabel for qform and sform; the ITK centre convention; deformable types are refused, not approximated; LPS-to-RAS conversion; rigid map equals `apply_rigid`; known-geometry coverage, where the right direction and the wrong direction give different masks; one exact 4x4 for the full chain.

### Visit-1 alignment choice (recommendation only; no amendment adopted)

- **Option 1: MRI-only rigid registration of visit 1 to visit 2**, with the same frozen method and QC as visit 3.
  - For: the same rule for every non-origin visit. Visit 1 feeds history inputs (FiLM channel, PDE start state). Misalignment there would bias forecasts and coverage. QC also catches bad atlas alignments.
  - Against: it needs 28 more atlas CT1 reads. Some patients may fail QC, which raises the risk of falling below the 25-patient minimum. It changes the history inputs of the frozen protocol, so a dated amendment is needed.
- **Option 2: accept the publisher's atlas alignment for visit 1**, as a stated assumption.
  - For: no new reads, and no change to history inputs.
  - Against: it is unverified, and it treats visit 1 differently from visit 3. A visit-3-style check was judged necessary there, and nothing shows visit 1 is safer.
- **Recommendation: Option 1.** It is the consistent and checkable choice. Its cost (cohort loss) is measured by the QC, not hidden.

### Is a valid geometry-only count possible now?

**No. The count stays blocked.** Missing:

1. A frozen .tfm direction. This needs DeepBraTumIA documentation or source showing the registration's fixed image, or a frozen rule chosen before any count (for example, the metadata evidence above) and stated as an assumption.
2. The visit-3 rigid refinement transforms (the approved MRI-only read of visit-2/3 atlas CT1 and visit-2 brain mask).
3. A frozen visit-1 rule (amendment). Under Option 1, also the visit-1 registrations and QC.
4. The frozen 2 mm block centre convention.
5. A runner that composes the chain per visit and feeds `coverage_from_geometry`, with tests. Not built.

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

1. Establish the full transform chain for visits 1, 2 and 3 (see the transform-chain audit): native affines, publisher `.tfm`, the visit-3 rigid refinement (`register_rigid_mri`, CT1 intensities only), and the approved visit-1 handling. Record full matrices and QC. No tumor label is opened. This needs a wider read approval than the current 84-file proposal.
2. Build coverage masks from those transforms and the grids (geometry only). Count patients whose common support covers the whole forecast domain. No tumor label is opened.
3. Apply the frozen 25-patient minimum. If fewer than 25 remain, the run stops with no scoring, and the option choice is reopened.
4. If 25 or more remain, freeze the cohort (patient list and transforms, hashed) before any target-mask read.

Do not switch options after seeing any target.

## Decision needed from you

Reply with one of:

1. "Adopt Option A, conditional on the frozen 25-patient minimum: exclude patients whose geometry-derived common support does not cover the whole forecast domain, from all arms. Approve steps 1-4 above: a wider read scope (native CT1 headers and `ct1.tfm` for visits 1-3, plus the proposed atlas CT1 and brain-mask files) and a dated visit-1 handling amendment, MRI-only registration and chain code to record transforms, then a geometry-only coverage count with tumor labels closed, then the minimum check, then a cohort freeze before any target-mask read." (recommended)
2. "Adopt Option B with primary cohort = [all patients / full-coverage subset]."

Option C would need a separate method plan; it is not a short decision.
