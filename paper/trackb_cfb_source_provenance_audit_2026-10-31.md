# Track B CFB-GBM source and label-provenance audit

**Date:** 2026-10-31  
**Scope:** Public-source and file-metadata audit only. No GTV or brain-mask voxel arrays were read. No patient forecast was scored.  
**Decision:** CFB-GBM remains a candidate cohort, not cleared for Track B patient scoring.

## What primary sources establish

The [TCIA CFB-GBM collection record](https://www.cancerimagingarchive.net/collection/cfb-gbm/) describes Version 4 as a 264-patient collection. Version 4 corrected RANO TSV values and renumbered some folders; TCIA states that image data were not changed. The image/segmentation package is CC BY 4.0 and requires the dataset citation.

The [CFB-GBM v2.0 paper](https://arxiv.org/html/2608.17884) says that imaging modalities were rigidly registered and resampled to baseline t0 T1-Gd using SimpleITK, a three-step multi-resolution scheme, and Mattes mutual information. This supports the stated cohort-level preprocessing method. The paper does not report case-specific transforms, registration scores or independent per-patient acceptance records for the 38 selected patients.

The paper describes mixed GTV provenance:

- The original visible-tumor GTVs were manually delineated under ESTRO guidance.
- Missing longitudinal GTV annotations were completed with nnU-Net, pre-trained on BraTS 2021 and fine-tuned on CFB-GBM ground-truth contours.
- The authors sampled 44 generated GTV annotations across the available timepoints for review and correction by five radiation oncologists; the model was then applied to generate missing annotations across t0, t1 and t2.

This is useful cohort-level evidence, but it does not identify which of the selected patients' t0/t1/t2 files are original manual contours, generated contours, or part of the reviewed subset. Do not treat all released GTV files as uniformly manual ground truth. The per-file source map remains an eligibility requirement.

## Linked registration code inspection

TCIA links the [author CFB-GBM repository](https://github.com/AurelienCD/CFB-GBM), whose README identifies `02_MRItoMRI.ipynb` as a 3D Slicer MRI-to-MRI registration notebook. Static inspection of the [raw notebook](https://raw.githubusercontent.com/AurelienCD/CFB-GBM/main/02_MRItoMRI.ipynb) found:

- The notebook computes a transform and writes resampled NIfTI images, but the shown cells do not persist the per-case transform or a case-level QC/acceptance record.
- The same variable name `files` is used for the directory listing and as the loop variable in consecutive loops. After the first loop, the later t1/t2 loops iterate over the last filename string rather than the original directory listing. If run verbatim in order, those later loops would not process the t1/t2 file list as intended.
- The notebook's stored execution counts are null. This source-code issue does not show that the released images are unregistered or that this exact notebook produced Version 4. It shows that the linked notebook, as currently published, is not sufficient evidence of the actual per-patient transform lineage.

## Current local audit boundary

The local paginated review aid covers 38 candidate IDs and contains T1-Gd display pages only. Its manifest states that GTV and brain-mask voxel arrays were not opened. Header/grid agreement and sampled MRI displays do not prove anatomical or tumor-boundary registration. No independent reviewer has assessed the candidate set; see [the review packet](trackb_cfb_independent_review_packet_2026-10-31.md) and [blank pair-level review form](trackb_cfb_registration_review_template_2026-10-31.csv).

The dataset paper and code do not establish that zero-valued voxels are valid negative tumor labels over the intended PDE evaluation grid. A brain mask alone cannot prove that. Valid-negative support, per-file GTV source, and correspondence to the planned Track B target remain separate unresolved gates.

## Next permitted sequence

1. Reconcile the 38 selected patients and each t0/t1/t2 GTV filename against source metadata that identifies manual, generated, and expert-reviewed annotations. Keep this metadata-only; do not read mask voxels.
2. Freeze a defensible registration-review rubric, then obtain an independent review of full-resolution t0-to-t1 and t0-to-t2 T1-Gd volumes. The sampled 5-mm PNG pages are navigation aids only.
3. Decide the scoring support and establish from source evidence that included voxels have valid positive/negative meaning. Do not infer valid negatives from zero values or brain masks.
4. Only after these gates pass, freeze eligible cases, target definition, baselines, analysis code and failure rules; then consider patient-mask scoring.

If source metadata cannot map label provenance or define negative support, record the failed gate and keep CFB out of patient-level PDE scoring. Do not substitute a post hoc mask rule to obtain a positive score.
