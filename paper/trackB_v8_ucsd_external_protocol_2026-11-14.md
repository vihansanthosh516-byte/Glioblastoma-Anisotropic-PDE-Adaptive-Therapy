# Track B: frozen v8 external test on UCSD-PTGBM — locked protocol (2026-11-14)

Written and committed BEFORE any UCSD volume, forecast or score is computed in this test.

## Same question, same score
Does v8 predict the next tumor-core volume better than standard methods? Target y = log(core voxels + 1), core = labels {1, 3}. Score: patient-averaged MAE of y. Bootstrap 5000, seed 20261114.

## Frozen model
- v8 (`src/track_b_mu_v8.py`, amendment `paper/trackB_mu_v8_amendment_2026-11-14.md`) and all 6 baselines are fit ONCE on ALL MU primary pairs. No refit, no change on UCSD.
- Comparator fixed in advance: population growth rate (`pop_rate`), the best baseline on MU. All 6 baselines are also reported.
- Success: upper bound of the 95% patient-bootstrap CI of (v8 MAE - pop_rate MAE) < 0. One test, no multiple-try split.
- Minimum: 25 patients with at least one usable pair. Fewer -> no-go, no score reported as a result.

## UCSD pairs
- Sessions: `UCSD-PTGBM-XXXX_NN` with `*_BraTS_tumor_seg.nii.gz` under `data/external/UCSD-PTGBM-v3` (v1 and BraTS-2024 test-set folders; duplicates by session ID read once).
- Scan day (since diagnosis) = -("Days from Acquisition to Date of initial surgery, treatment or diagnosis") from `research_notes/trackB_confirmatory_audit_ucsd/UCSD_PTGBM-clinical-information_v3_2026-12-Mar.xlsx`.
- Pairs: chronologically adjacent sessions of one patient; dt 14 to 365 days; input core non-empty (MU primary rules). Prior-scan rate from the previous adjacent session when it exists.
- All diagnoses kept (as in MU).

## Phase mapping (assumptions; UCSD lacks MU's exact fields)
- RT window: [first-radiation day, first-radiation day + 42]. (RT end date not given; 42 days = standard 6-week course.)
- Adjuvant window: [first-chemo start, first-chemo start + 180], outside the RT window.
- post_progression: "Progression Free survival" known and <= scan day of input. Unknown -> not post-progression.
- Clinical day fields are per session; each session's own row is used.

## Weak points (stated in advance)
- UCSD was exposed before in this project for a different target (annotated core additions, `src/166_trackb_ucsd_external_test.py`). Volume forecasting with v8 has not been run on it. Still not a fully untouched cohort.
- Different hospital, scanners and segmentation pipeline from MU: a real shift test, but label definitions may differ.
- Phase mapping uses assumed windows.
