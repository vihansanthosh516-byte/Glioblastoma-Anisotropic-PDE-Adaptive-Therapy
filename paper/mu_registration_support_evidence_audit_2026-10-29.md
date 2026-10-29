# MU Registration and Target-Support Evidence Audit

Updated 2026-10-29. Scope: determine whether public MU-Glioma-Post documents clear the existing pre-mask scoring gates. This is a source and metadata audit. No tumor-mask voxel arrays were read and no patient scores were calculated.

## Findings

### Longitudinal registration

The dataset paper says images were resampled to 1 mm isotropic resolution and rigidly co-registered to the SRI24 atlas with CaPTk/Greedy. It also says scans were reviewed to exclude significant motion or suboptimal quality. These statements establish a described preprocessing method and general image-quality review. They do not give the selected patients' applied transforms, visit-to-visit registration measurements, pass thresholds, or tumor-boundary alignment checks. Atlas registration of each scan does not by itself prove that each longitudinal triplet is aligned well enough for tumor-growth scoring.

### Segmentation review and valid support

The paper defines labels 1–4 as non-enhancing tumor core, surrounding non-enhancing FLAIR hyperintensity, enhancing tissue, and resection cavity. It reports manual refinement and neuroradiologist consensus for segmentations. This supports annotation review, not registration quality. The paper does not define label 0's spatial meaning or say that every zero-valued voxel is a valid negative throughout the whole NIfTI grid.

The paper also says the released MR images are skull-stripped and extracranial structures are removed. Therefore, treating the entire grid as a valid negative region needs explicit justification. This is a support-definition concern, not proof that any released mask is wrong. The current `entire_image_grid` gate may be stricter than the scientifically needed support. Do not loosen it until the protocol defines a source-backed evaluation support and keeps that support separate from the origin-only PDE domain.

### Public pipeline code

The linked public repository provides a generic wrapper that calls the FeTS preprocessing command. The wrapper does not expose patient-level transform records or registration-QC acceptance records. It cannot establish that no separate records exist.

## Local gate status

The loader requires a hashed evidence JSON tied to the exact manifests, MRI-only report, and selected patient IDs. Registration evidence must cover every selected patient and attest a passed QC status, method, reference, and a hash-bound local evidence file. Target-support evidence must cover every selected patient, define the support, establish valid negatives, and provide a hash-bound local evidence file. The current schema specifically requires `coverage_scope="entire_image_grid"`. See `src/fk_mu_development_data.py:7-28` and `:240-300`.

No such evidence JSON or supporting MU evidence file was found in the project evidence folders. The latest metadata-only preflight stopped at the missing evidence document after checking 72 NIfTI headers. It did not hash mask files, read target-mask arrays, or start scoring. See `output/fk_mu_pipeline_verification_20261029.json`.

The existing MRI-only report defines nonzero T1n intensity as a *support proxy*, explicitly says it is not a validated acquisition mask, and records that no registration or valid-support gate passed. See `src/fk_mu_mri_only_qc.py:211-251` and `output/fk_mu_mri_only_qc_20261028_v2/mri_only_qc.json`. It cannot be promoted to source-backed support by itself.

## Required next evidence

Ask the dataset team whether they can provide or point to:

1. Registration method and reference for each included timepoint, plus any available case-level transform/QC records or an accepted review procedure that covers the selected longitudinal triplets and tumor-relevant anatomy.
2. A definition of label 0 and any invalid, unannotated, or out-of-support voxels in the released segmentation volumes.
3. Any per-timepoint brain/support masks or source-backed instructions that define a valid evaluation region, including how that region relates to skull-stripping and the registration reference.

If records exist, preserve source files, hashes, exact cohort coverage, and human review in the evidence JSON. If not, keep patient-mask arrays closed. Revise the support protocol first; then design a prospective MRI-only registration review with frozen criteria. The current descriptive overlap/correlation screen has no pass threshold and does not clear that gate.

## Message draft for user to send

**Subject:** MU-Glioma-Post longitudinal registration and segmentation support documentation

Hello Dr. Isufi and Dr. Nada,

I am a student researcher at **[school or program]**, working with **[teacher or research mentor]**. My project is an exploratory study of patient-level PDE-based tumor-growth forecasts from serial postoperative MRI. I am using the public MU-Glioma-Post release for development only. I will follow the TCIA data-use terms, cite the dataset and paper, and will not attempt to identify or contact any participant. The work is not for clinical use.

I read that the images were rigidly registered to SRI24 and skull-stripped, and that the tumor segmentations were manually reviewed. Could you point me to any available documentation that answers these questions?

1. Were each patient's timepoints registered independently to SRI24, or was there also longitudinal visit-to-visit registration? Are case-level transforms, registration QC measures, acceptance criteria, or review records available for the released timepoints?
2. What does label 0 mean in the segmentation NIfTI files? Does it mean no tumor/other labeled tissue throughout the full image grid, or are any voxels outside the valid annotation/support region?
3. Are per-timepoint brain/support masks or instructions available that define which voxels should be treated as valid for longitudinal tumor-mask evaluation, especially outside the skull-stripped brain?

I do not need any identifying information. If any records are available, I can match them to the anonymized released timepoint IDs under the dataset's terms. If no such records exist, a brief clarification of the intended registration and support interpretation would still help me set a cautious evaluation protocol.

Thank you for your time,

**[name]**

**[school/program and mentor contact, if required]**

## Sources

- Mahmoud et al., *MU-Glioma Post: A comprehensive dataset of automated MR multi-sequence segmentation and clinical features*, 2025. Methods and QC: sections “MRI acquisition and imaging data,” “Preprocessing pipeline and data deidentification,” and “Technical Validation.” [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12635165/).
- Authors' public [MU_Glioma_Post_segmentation_pipeline](https://github.com/amnmed09/MU_Glioma_Post_segmentation_pipeline), especially `run_fets_preprocessing.py`, which calls the general FeTS preprocessing command.
- Official [TCIA MU-Glioma-Post collection](https://www.cancerimagingarchive.net/collection/mu-glioma-post/).
