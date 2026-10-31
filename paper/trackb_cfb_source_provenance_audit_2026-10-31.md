# Track B CFB-GBM source and label-provenance audit

**Date:** 2026-10-31

**Scope:** Public-source and file-metadata audit only. No GTV or brain-mask voxel arrays were read. No patient forecast was scored.

**Decision:** CFB-GBM remains a candidate cohort, not cleared for Track B patient scoring. For the selected 38 cases, all t1 and t2 GTV targets are marked `genSegmentation`; future scores would measure agreement with these generated labels unless their inclusion in the expert-reviewed subset can be source-verified.

## What primary sources establish

The [TCIA CFB-GBM collection record](https://www.cancerimagingarchive.net/collection/cfb-gbm/) describes Version 4 as a 264-patient collection. Version 4 corrected RANO TSV values and renumbered some folders; TCIA states that image data were not changed. The image/segmentation package is CC BY 4.0 and requires the dataset citation.

The [CFB-GBM v2.0 paper](https://arxiv.org/html/2608.17884) says that imaging modalities were rigidly registered and resampled to baseline t0 T1-Gd using SimpleITK, a three-step multi-resolution scheme, and Mattes mutual information. This supports the stated cohort-level preprocessing method. The paper does not report case-specific transforms, registration scores or independent per-patient acceptance records for the 38 selected patients.

The paper describes mixed GTV provenance:

- The original visible-tumor GTVs were manually delineated under ESTRO guidance.
- Missing longitudinal GTV annotations were completed with nnU-Net, pre-trained on BraTS 2021 and fine-tuned on CFB-GBM ground-truth contours.
- The authors sampled 44 generated GTV annotations across the available timepoints for review and correction by five radiation oncologists; the model was then applied to generate missing annotations across t0, t1 and t2.

This is useful cohort-level evidence, but the paper alone does not identify which selected t0/t1/t2 files are original manual contours, generated contours, or part of the reviewed subset. Do not treat all released GTV files as uniformly manual ground truth. The metadata-only mapping below resolves the manual-versus-generated category for the 38 candidates, but not membership in the expert-review subset.

### Metadata-only mapping for the 38 selected candidates

I downloaded and read the official [TCIA treatment-imaging availability TSV](https://www.cancerimagingarchive.net/wp-content/uploads/CFB-GBM_treatment_imaging_availability_v03_20260619.tsv) and [MRI-availability TSV](https://www.cancerimagingarchive.net/wp-content/uploads/CFB-GBM_mri_availability_v02_20260129.tsv). This examined availability, `gtv_type` and week-valued timing only; no image or mask voxels were read. The treatment-imaging TSV SHA-256 is `33d5b6e185cc9582af98285d1d199aa48ccd1a97e0c6e83e0defbd7f2dd515d8` and the MRI-availability TSV SHA-256 is `74c40a04a5e0f328232e639e080dda8fbcc0d0297771d11eef3dab5e8abed389`.

Across the 38 candidate patients, all 114 t0/t1/t2 rows mark a GTV file as present. The source `gtv_type` values are:

- t0: 36 `humanSegmentation`, 2 `genSegmentation` (patient IDs 18 and 197).
- t1: 38 `genSegmentation`.
- t2: 38 `genSegmentation`.

The metadata does not identify which generated files, if any, belong to the 44-mask expert-review subset. Do not describe all 76 future-visit targets as expert-reviewed or manual. The per-candidate metadata-only join is preserved in [trackb_cfb_selected_candidate_metadata_2026-10-31.csv](trackb_cfb_selected_candidate_metadata_2026-10-31.csv).

The MRI-availability table records `time_diff_t0 (weeks)`, not exact acquisition dates. For these candidates, t1-to-t2 gaps computed from that field are 77, 84, 91 or 98 days. The table does not establish exact day-level windows.

The paper defines CFB's binary GTV as enhancing tumor plus necrotic tumor, excluding edema. Before any CFB scoring, the protocol must explicitly map that binary clinical GTV to the Track B source-defined target; do not silently assume the two definitions are identical.

## Linked registration code inspection

TCIA links the [author CFB-GBM repository](https://github.com/AurelienCD/CFB-GBM), whose README identifies `02_MRItoMRI.ipynb` as a 3D Slicer MRI-to-MRI registration notebook. Static inspection of the [raw notebook](https://raw.githubusercontent.com/AurelienCD/CFB-GBM/main/02_MRItoMRI.ipynb) found:

- The notebook computes a transform and writes resampled NIfTI images, but the shown cells do not persist the per-case transform or a case-level QC/acceptance record.
- The same variable name `files` is used for the directory listing and as the loop variable in consecutive loops. After the first loop, the later t1/t2 loops iterate over the last filename string rather than the original directory listing. If run verbatim in order, those later loops would not process the t1/t2 file list as intended.
- The notebook's stored execution counts are null. This source-code issue does not show that the released images are unregistered or that this exact notebook produced Version 4. It shows that the linked notebook, as currently published, is not sufficient evidence of the actual per-patient transform lineage.

## Current local audit boundary

The local paginated review aid covers 38 candidate IDs and contains T1-Gd display pages only. Its manifest states that GTV and brain-mask voxel arrays were not opened. Header/grid agreement and sampled MRI displays do not prove anatomical or tumor-boundary registration. No independent reviewer has assessed the candidate set; see [the review packet](trackb_cfb_independent_review_packet_2026-10-31.md) and [blank pair-level review form](trackb_cfb_registration_review_template_2026-10-31.csv). The per-case `gtv_type` and week-derived timing metadata are in [the selected-candidate metadata table](trackb_cfb_selected_candidate_metadata_2026-10-31.csv).

The dataset paper and code do not establish that zero-valued voxels are valid negative tumor labels over the intended PDE evaluation grid. A brain mask alone cannot prove that. Valid-negative support, per-file GTV source, and correspondence to the planned Track B target remain separate unresolved gates.

## Next permitted sequence

1. Search the public source package for a per-file reviewer flag identifying which selected generated GTVs, if any, were in the 44-mask radiation-oncologist review subset. The currently checked public metadata lacks that field; if no public mapping exists, report it as unavailable and do not infer membership.
2. Freeze a defensible registration-review rubric, then obtain an independent review of full-resolution t0-to-t1 and t0-to-t2 T1-Gd volumes. The sampled 5-mm PNG pages are navigation aids only.
3. Decide whether the study will target the CFB binary GTV (enhancing plus necrotic tumor, edema excluded) and how that maps to the current Track B endpoint. Establish from source evidence that included voxels have valid positive/negative meaning. Do not infer valid negatives from zero values or brain masks.
4. Only after these gates pass, freeze eligible cases, target definition, baselines, analysis code and failure rules; then consider patient-mask scoring.

If source metadata cannot map label provenance or define negative support, record the failed gate and keep CFB out of patient-level PDE scoring. Do not substitute a post hoc mask rule to obtain a positive score.
