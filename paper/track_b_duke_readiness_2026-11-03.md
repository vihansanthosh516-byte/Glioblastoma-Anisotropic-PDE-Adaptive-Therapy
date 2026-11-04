# Track B Duke readiness (updated 2026-11-04; metadata audit completed 2026-11-03)

## Current direction

The proposed Track B primary question is biological and radiographic: within the same patient and physical-distance band, do newly appearing nonenhancing FLAIR-labeled regions convert to a next-visit necrotic/enhancing core label more often than persistent FLAIR-labeled regions? ADC prediction is secondary. Reaction-diffusion forecasting is a separate secondary benchmark. The competition-form draft is `paper/trackB_research_plan_draft_2026-11-03.md`; it is proposed, not adopted or frozen.

The study concerns transitions among automated MRI labels under observed treatment. It cannot establish histologic infiltration, a causal PDE effect, clinical benefit, or independent external confirmation without the relevant evidence.

## Archive and metadata status

- The user reports having Duke download access. The exact accepted written agreement and its output-sharing terms are not saved in the project.
- The local `DATASET_FINAL.zip` is 17,043,569,844 bytes and matched the MD5 `99580fb8ccb3daed186dc6afe5b5bf2f` supplied with the authorized download. The restricted Zenodo API did not expose a file checksum for independent comparison. The local SHA-256 recorded by the inventory is `0b8c2bfa7d881bd0e81e87e52646d676c314e73e5e60290fc8685a5c0ebc867a`.
- The metadata-only inventory is complete. It listed ZIP central-directory entries and matched metadata session keys to archive session folders. It read the approved anonymized patient/session, relative-date, diagnosis/IDH, and treatment metadata. Its patient-level output remains outside the repository.
- No NIfTI member was extracted or decompressed. No MRI or segmentation-mask voxel array and no landmark-coordinate table was opened. No patient score was calculated.
- The metadata audit is not a final eligibility or independence determination. The proposed inclusion rule uses records coded as glioblastoma and IDH `Negative`, excludes `Not Provided`, `indeterminate`, and mixed/ambiguous diagnoses, and does not claim comprehensive IDH1/2 confirmation. Patient overlap with prior project cohorts has not been reconciled. Image quality and label usability have not been checked.

## Source facts and limits

The [Zenodo record](https://zenodo.org/records/20548528) describes 61 patients, 503 sessions, six MRI sequences, and automated labels: 1 necrotic core, 2 nonenhancing T2-FLAIR abnormality, 3 enhancing tumor, and 4 resection cavity. These published totals do not establish the eligible three-visit cohort.

The [author guide](https://github.com/deshanyang/GBM-DIR-QA) says separate sessions are not longitudinally registered and that extracted dates may be off by several days. The archive inventory found the expected six MRI sequence filenames and a segmentation file per session; it found no released skull-stripped brain mask. Therefore, the research-plan draft defines an operational MRI foreground from history T1-c and FLAIR images. That foreground is not a validated anatomical brain mask and must be checked during a separately approved image-only preflight.

The dataset uses automated segmentation labels. The exact trained segmentation checkpoint is not stated in the record, and label 0 is not defined. The study will call its outcomes automated-label transitions, not tumor truth. The terms for derived masks, patient-level predictions, figures, and model weights depend on the exact accepted agreement and remain open.

## Gates before patient analysis

1. Save the written Duke grant and exact accepted data-use terms, including allowed analysis and retention or sharing of outputs.
2. Save written school or competition review, with date, reviewer, and conditions.
3. Complete the allowed patient-overlap check and freeze the eligibility and minimum-cohort rules. Do not call the cohort independent if overlap cannot be resolved.
4. Adopt and freeze a Duke-specific protocol, including the primary transition estimand, labels, registration, operational support, exclusion and fallback rules, folds, all secondary models, uncertainty methods, precision rules, software, and runtime ceiling.
5. Complete synthetic integration and full-pipeline timing. Do not weaken comparators or thresholds to fit the available compute after seeing outcomes.
6. Only after the source, school, and protocol gates pass, approve and run an MRI-only registration/support QC preflight. Keep target segmentation labels closed during that step. Freeze the QC-passing cohort and hashes before outcome scoring.
7. Seal code, protocol, cohort, folds, and allowed input fingerprints. Score once, report negative and incomplete results unchanged, and keep patient-level outputs local unless the accepted agreement permits sharing.

Until these gates pass, patient MRI/mask/landmark reads and scoring remain blocked. Metadata inventory completion alone does not authorize them.
