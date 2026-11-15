# Track B: frozen v8 test on UCSF-ALPTDG (time1 -> time2) — locked protocol (2026-11-14)

Written and committed BEFORE any ALPTDG volume forecast or score is computed in this test.

## Same question, same score
Predict the next tumor-core volume. y = log(core voxels + 1), core = labels {1 (NCR), 3 (ET)} from `<id>_time1_seg` / `<id>_time2_seg`. Score: patient MAE of y (1 pair per patient). Bootstrap 5000, seed 20261114.

## Frozen models
v8 and the 6 baselines fit ONCE on all MU primary pairs, exactly as `src/track_b_v8_ucsd_external.py` (`fit_predict`). No refit, no change.

## Cohorts
- PRIMARY: ALPTDG patients NOT in `output/predictgbm/h3_final/ucsf_test_pairs.csv` and NOT in the "UCSF-PDGM Mapping" sheet (those were used in H-3/H-4). Expected about 236-239 patients.
- SECONDARY (labelled "used before"): the mapped patients.
- Pair rules (as MU): dt 14-365 days, time1 core non-empty.
- Minimum 25 primary patients.

## Comparator and success (fixed now)
Comparator: pop_rate (best MU baseline). Success: upper bound of the 95% patient-bootstrap CI of (v8 MAE - pop_rate MAE) < 0 on PRIMARY. Also reported: v8 vs persistence and vs every baseline.

## Inputs from `UCSF_PostopGlioma_Table S1 R1 V5.0_UNBLINDED_FINAL.xlsx`, sheet "Clinical Info"
- day_in = "Days from 1st surgery/DX to 1st scan"; dt = "Days from 1st scan to 2nd scan".
- No earlier scan exists: has_prior = 0 for all (v8's no-prior path).
- RT window: [day_in + "Days from 1st scan to 1st RT start", +42 days].
- Adjuvant window: [day_in - "Days from 1st chemo start to 1st scan", +180 days], outside the RT window.
- post_progression: not derivable (table gives only progression AFTER scan 1) -> never set. Assumption.
- age = "Patient Age" (UCSF-LPTDG sheet, time 1); mgmt = missing; GBM flag = "WHO 2021 Diagnosis" contains "Glioblastoma".
- Composition features (v5) from the time1 seg only.

## Weak points (stated in advance)
- Same hospital as the 59 previously used patients; diffuse gliomas of all grades.
- No prior scan, so v8's latent-size part reduces to y_in.
- Phase windows assumed; post-progression phase unavailable.
- Expert manual masks here vs automatic masks in MU: label shift.
