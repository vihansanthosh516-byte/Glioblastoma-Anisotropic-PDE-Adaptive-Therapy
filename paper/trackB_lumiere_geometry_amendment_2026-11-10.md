# Track B amendment: LUMIERE geometry rules (ADOPTED 2026-11-10, before any MRI intensity read)

Adopted on the coordinator's instruction ("Adopt the 2 mm grid convention ... Adopt MRI-only rigid registration of visit 1 to visit 2"). Written before any CT1 intensity, coverage mask, or eligible-patient count was computed. No rule here was chosen from patient outcomes or cohort size.

## G1. Transform direction (validated, not assumed)

The stored `DeepBraTumIA-segmentation/native/transformation/ct1.tfm` maps native points to atlas points. Coverage uses `inv(T)` as the atlas-to-native map. Evidence: transform-direction validation against paired brain masks under a protocol frozen before reading (`paper/trackB_tfm_direction_check_protocol_2026-11-10.md`; record `output/track_b_lumiere_transform_audit_20261110/tfm_direction_check.json`): 84 scans nested within 28 patients; median Dice 0.9989 (min 0.9808) for `inv(T)` vs. 0.3936 for `T`; `inv(T)` won 84/84. This validates direction only; it is not proof of registration accuracy inside the brain.

## G2. 2 mm scoring grid (supersedes the proposal in the decision brief)

Scoring-grid index i maps to 1 mm atlas index **2i** on each axis (same origin, spacing doubled), matching the existing Track B nearest-neighbour operation (`src/fk_mu_development_data.py::_target_grid` + nibabel `resample_from_to(order=0)`). Grid 91 x 109 x 91. Labels and masks use nearest neighbour, which under this rule is `mask_1mm[::2, ::2, ::2]`. Test: `tests/test_fk_lumiere_transform_audit.py::test_existing_track_b_2mm_grid_maps_scoring_voxel_i_to_1mm_index_2i`.

## G3. Visit-1 registration

Visit 1 is rigidly registered to visit 2 with `src/fk_longitudinal_prep.py::register_rigid_mri`, using the frozen section 5b method and QC of `paper/track_b_no_contact_dataset_plan_2026-10-31.md` unchanged (fixed: visit-2 atlas CT1; fit region: visit-2 atlas brain mask > 0; NCC >= 0.80, max |Euler angle| <= 10 deg, shift <= 10 mm; Powell two-stage; no retry). Moving image: visit-1 atlas CT1. A patient failing visit-1 or visit-3 QC (or the sanity checks) is excluded from all arms before any target read. Visit-1 history labels and the visit-1 coverage chain use the visit-1 rigid transform.

**Code version note.** The section 5b hash (`6588f511...`) predates three committed fixes to the same file: optimiser non-convergence is now a QC failure (57735f2), non-finite parameters are rejected (e027624), and a landmark helper was added (05c1038). Settings and thresholds are unchanged. The registration run uses the current file, SHA-256 `1e6c2ed59f7237b0e1d8be4ff559f72fe34440eee82bd4981d8ac90f5ff460bd`.

## G4. Coverage chain and geometry-only eligibility

For visit k in {1, 2, 3}, the output-to-input voxel map from scoring voxel to native CT1 voxel is

`inv(N_k) @ inv(T_k) @ A @ R_k @ D`

where D = 2i scoring-to-atlas map (G2); `R_k` = the rigid output-to-input voxel map of `apply_rigid` for visit k (identity for visit 2); A = atlas affine; `T_k` = stored ct1.tfm in RAS (G1); `N_k` = native CT1 affine. Coverage = `geometric_coverage(native_shape_k, (91, 109, 91), map)`; common support = AND over the three visits.

**Forecast domain for this count:** the visit-2 atlas brain mask on the 2 mm grid (G2). The frozen support also adds history-label voxels outside the brain mask; these are tumor labels and are not read here, so that part is checked later, before target reads (a remaining gate; a patient whose history-label voxels fall outside common support is excluded then).

**Eligible:** registration QC passed for visits 1 and 3, sanity checks passed, and common support covers every forecast-domain voxel. Frozen minimum: 25 patients; below 25 the run stops with no scoring. The eligible list and transforms are saved and hashed before any target-mask read.
