# Track B: bounded fresh-cohort screen (2026-11-11, metadata only)

**Result: no candidate passes every gate. Recommendation: close patient-level Track B forecasting for now; keep the LUMIERE registration work as a feasibility result** (`paper/trackB_lumiere_registration_nogo_2026-11-11.md`). No code or frozen rule was changed. Nothing was downloaded, no one was contacted, and no image or mask voxel was read.

## Gates (frozen before the search; unknown is not a pass)

G1 repeated scans supporting a forecast interval (at least three labelled visits per patient: two history, one target); G2 patient-level scan and segmentation links; G3 usable tumor labels with documented definitions; G4 enough eligible patients (frozen project minimum 25); G5 access and reuse terms that permit the analysis without a pending request; G6 no overlap with datasets or patients already used here.

## Already screened or used (not repeated)

Used in this project, so G6 fails: MU-Glioma-Post, LUMIERE, RHUH-GBM, TUM-GBM, UCSD-PTGBM, CFB-GBM, BraTS, Duke. PREDICT-GBM is TUM-GBM + LUMIERE + RHUH (G6 fail). Rejected in `paper/track_b_external_cohort_screen_2026-10-28.md`: UTSW-Glioma, UCSF-ALPTDG, SAILOR, MOTUM, UPENN-GBM, BrainMD, CPTAC-GBM, UCSD-BMETS-Longitudinal.

## New candidates this screen

| Candidate | Source | Eligible-patient evidence | G1 | G2 | G3 | G4 | G5 | G6 |
|---|---|---|---|---|---|---|---|---|
| OpenNeuro ds007045 "Real-world Multi-Institutional Glioblastoma MRI Dataset" (v2.0.1) | https://openneuro.org/datasets/ds007045/versions/2.0.1 (README) | 363 patients, 8 centres; README: "underwent preoperative neuroimaging before surgery"; subjects labelled `-000` (one session) | **Fail** (single preoperative timepoint) | Pass (BIDS per subject) | Pass (SegResNet, BraTS labels, neuroradiologist-reviewed) | Fail (0 with three visits) | Unknown (OpenNeuro click-through; licence not checked) | Pass (sites not used here) |
| TCIA Brain-Tumor-Progression | https://www.cancerimagingarchive.net/collection/brain-tumor-progression/ | 20 newly diagnosed GBM patients; "Two MRI exams are included for each patient" | **Fail** (two exams; no two-history window) | Unknown | Unknown | **Fail** (20 < 25) | Unknown | Unknown |
| CFB-GBM v2.0 (augmented, 2026) | https://arxiv.org/html/2608.17884v1 | 264 patients (same centre) | Unknown | Unknown | Unknown | Unknown | Unknown | **Fail** (CFB-GBM already used in a project side study) |

Other hits (MU-Glioma-Post 2025 paper, UCSD-PTGBM 2026 paper, LUMIERE) are already-used cohorts (G6 fail).

## Limits

One bounded search session (three web searches, two page reads); not a systematic review. Controlled-access releases needing a data-use request fail G5 by the project's no-request rule and were not pursued.
