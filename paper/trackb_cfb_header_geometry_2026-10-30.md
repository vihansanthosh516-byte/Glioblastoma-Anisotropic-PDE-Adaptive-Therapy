# CFB candidate MRI and brain-mask header check

2026-10-30. This is a file-geometry check, not a patient forecast or registration-quality assessment.

## Input and procedure

The existing main-worktree crosswalk `research_notes/GBM CFB scoring readiness/CFB_v4_candidate_crosswalk_20261028.csv` supplies 38 distinct candidate patients whose table-derived t1-to-t2 intervals are 76–104 days. These intervals originate from week-resolution data, not verified exact scan dates. Crosswalk SHA-256 begins `05a4ead1`.

For each t0, t1 and t2 visit, the local T1-Gd path was taken from the crosswalk and the corresponding `brain_mask.nii.gz` path was derived by replacing the filename suffix `_t1gd.nii.gz` with `_brain_mask.nii.gz`. Python 3.14 and nibabel 5.4.2 opened NIfTI headers only. No MRI voxel data, brain-mask voxel data, or GTV arrays were read. Each T1-Gd shape and affine were compared with its same-visit brain mask. T1-Gd shapes and affines were also compared across the three visits. Affine equality used absolute tolerance 0.0001 and zero relative tolerance.

## Observed result

- All 114 MRI–brain-mask pairs opened successfully and had matching shapes and affines.
- All 38 patients had matching T1-Gd shapes and affines across their three visits.
- No header read error or within-visit geometry mismatch occurred.

An independent second header-only check reproduced these counts. The crosswalk's source-table provenance fields classify 36 t0 candidate GTVs as human and two as generated; all 38 t1 and all 38 t2 candidate GTVs are generated. Individual review status remains unknown. These provenance counts are metadata, not an inspection or assessment of the GTVs.

These checks establish matching header geometry only. They do not validate brain-mask contents, anatomical registration, tumor-boundary alignment, tumor-label coverage, negative labels, exact dates or scoring eligibility. The original CFB preflight remains unpassed; no patient scoring is authorized by these observations.
