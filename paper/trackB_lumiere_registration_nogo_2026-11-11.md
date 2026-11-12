# Track B LUMIERE: registration QC no-go under frozen rules (2026-11-11)

**Decision: NO-GO. The LUMIERE eligibility pipeline is stopped.** 7 of 28 patients passed registration QC, below the frozen minimum of 25. No coverage count was run, no cohort was frozen, no target mask was read, and no forecast was scored. QC thresholds and the 25-patient minimum are not relaxed.

## Run

- Code: `src/fk_lumiere_registration_run.py` (commit 13336f3), registration code `src/fk_longitudinal_prep.py` SHA-256 `1e6c2ed59f7237b0e1d8be4ff559f72fe34440eee82bd4981d8ac90f5ff460bd`. Six workers, one run, exit code 0.
- Method: amendment `paper/trackB_lumiere_geometry_amendment_2026-11-10.md` G3 with the frozen section 5b settings. Thresholds read back from code: NCC >= 0.8, max |Euler angle| <= 10.0 deg, shift <= 10.0 mm (unchanged).
- Files read: 112, i.e. 84 atlas `ct1_skull_strip.nii.gz` (visits 1-3) and 28 visit-2 atlas `brain_mask.nii.gz`. Zero `seg_mask` or other tumor-label files. Per-file SHA-256 in each case record.
- Records: `output/track_b_lumiere_registration_20261110/registration_qc.json` and `cases/Patient-*.json` (28/28). All transform parameters finite.

## Result (28 patients)

| | Count |
|---|---|
| Passed both visits | **7** (Patient-023, 032, 042, 043, 058, 062, 068) |
| Failed | 21 |
| Visit-1 registration passed | 9 / 28 |
| Visit-3 registration passed | 19 / 28 |
| Sanity (file/shape/finite/empty-mask) failures | 0 |
| Processing or code errors | 0 |

All 21 failures are QC failures. Reasons, counted per failing check (a patient can have several): visit-1 NCC < 0.8: 17; visit-3 NCC < 0.8: 9; visit-1 rotation > 10 deg: 3; visit-3 rotation > 10 deg: 2; visit-1 shift > 10 mm: 1; visit-3 shift > 10 mm: 1. Per-patient values and reasons are in `registration_qc.json` (`failed` and `cases`).

**The no-go does not depend on the visit-1 amendment.** Visit-3 registration alone (the pre-amendment frozen rule) passed 19 of 28, also below 25.

## Code-defect assessment

No evidence of a code defect: no exceptions, no sanity failures, all parameters finite, and 7 patients pass with NCC 0.81-0.95. Some failures are large (Patient-040 visit 3: 71.8 deg rotation, NCC 0.264; Patient-070 visit 3: 63.8 mm shift; Patient-086 visit 1: NCC 0.124). These may be optimiser local minima or real anatomical change (surgery, resection cavity, follow-up gaps) between scans. Which one is not established. Any investigation must use synthetic data only; patient images are not rerun without a new dated protocol decision.

## Next valid options (none chosen)

1. Close LUMIERE for Track B forecast scoring and report this no-go.
2. A new dated protocol (before any further patient image read) for a different, pre-specified registration method (for example, a multi-start or initialised rigid search). It must be validated on synthetic data first, must be frozen before rerun, and the rerun would be labelled post hoc, since these QC results are now known.
3. Move the Track B forecast question to another cohort, under its own frozen protocol.
