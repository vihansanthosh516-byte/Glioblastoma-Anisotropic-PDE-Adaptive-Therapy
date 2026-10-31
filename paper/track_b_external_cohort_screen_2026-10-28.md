# Track B PDE external cohort screen (2026-10-28)

**Status (updated 2026-10-31):** metadata-only source review. Duke declined the first access request because planned-use and academic-affiliation details were missing. The user reports resubmitting on 2026-10-29 with the requested research-purpose and license-compliance information; dataset access and school/fair permitted-use approval remain pending. No restricted files, patient-level outcomes, images, masks, or scores were accessed for this screen.

## Strongest current lead: Duke longitudinal GBM MRI

The [Zenodo record](https://zenodo.org/records/20548528), *A Longitudinal Glioblastoma MRI Dataset with Anatomical Landmark Pairs for DIR Validation*, describes 61 patients with 3–18 sessions each (503 total), skull-stripped 1-mm isotropic atlas-aligned NIfTI, and T1, T1-c, T2, FLAIR, ADC, and DWI at each session. Relative scan and treatment dates are available. Every session has an nnU-Net mask with necrotic core, nonenhancing T2/FLAIR abnormality, enhancing tumor, and resection cavity labels. The record was published 2026-06-12 and modified 2026-09-23; its data files are restricted and total 119.3 GB.

The record offers CC BY-NC-ND 4.0 or a custom Duke license. A request asks for academic affiliation, if any, and a research-purpose description; access is not guaranteed. Duke declined the initial request and the user reports resubmitting on 2026-10-29. School/fair acceptance of the selected license and derivative-data rules, and access approval, remain unconfirmed. Do not download restricted files before both approvals are documented.

## Fit and limitations

- Multi-visit masks make this a plausible lead for a **new, separately amended patient-level PDE growth-forecast evaluation**. A candidate design uses two historical masks to forecast the next eligible mask. The visit schedule is irregular, so freeze whether the estimand is next-observed-visit conditional on the realized gap or a prespecified fixed-day window supported by the available dates. Do not describe a realized-gap analysis as a scheduled-at-origin forecast.
- The published record says full cohort-selection details will appear in a future paper. Its recurrence category can be assigned by tumor-volume change or clinical progression, so a recurrence-conditioned analysis can select on tumor change. The authors' [dataset guidance](https://github.com/deshanyang/GBM-DIR-QA) adds that separate sessions have not been registered to one another, and recommends rigid or affine registration before downstream use. It says scan dates were extracted automatically from exam/report dates and may be imperfect or off by a few days; its recurrence workflow can label a scan after >40% enhancing-volume growth from the postoperative nadir. These create registration, timing, and outcome-selection risks for a longitudinal forecast. The masks are automated nnU-Net outputs, not independent expert ground truth. Treat a mask-based endpoint as released-label change unless validated against an independent reference.
- The 61 patients are not a power guarantee. Before outcome access, audit patient overlap and prior exposure against MU, UCSD-PTGBM, CFB-GBM, and other project cohorts; estimate precision under the final patient-level estimand; and freeze patient-grouped splits, baseline comparators, parameter-identifiability and calibration checks, and technical-failure rules.
- UCSD-PTGBM is a separate, partially exposed annotation-change side study with a negative exploratory result; it is not the PDE core. Duke would require a new protocol/amendment and cannot retroactively change that result.
- CFB-GBM is a separate annotation-transition candidate with unresolved review lineage and registration QC. A CFB result cannot complete the PDE parameter-estimation and forecasting program.

## PREDICT-GBM is not a fresh cohort

The new [PREDICT-GBM platform](https://doi.org/10.1038/s41746-026-03194-0) reports 243 patients across TUM-GBM (142), LUMIERE (61), and RHUH-GBM (40), and releases processed data and an evaluation pipeline. Those cohorts already appear in this project’s development or prior analyses: the local PREDICT-GBM work used the 142 TUM patients, and LUMIERE/RHUH were also examined. Treat the platform as useful for reproducibility and baseline review, not as untouched external validation. Its primary endpoint is geometric coverage of later recurrence for radiotherapy planning, which is related to but distinct from a PDE forecast of future tumor-core geometry. See the local exposure record in `README_STATUS.md` and the [PREDICT-GBM data release](https://huggingface.co/datasets/LZimmer/PREDICT-GBM).

## Decision

Duke remains the strongest metadata-level lead for a future multi-visit Track B PDE study, but it is **not ready for download or scoring**. Access approval and school-use approval are pending. Its source says the 61 patients have 3–18 sessions and automated labels, but it does not establish within-patient registration between sessions. If approval arrives, first audit cohort overlap, dates, registration, support, mask provenance, and metadata-only eligibility. Do not inspect target-mask values until those gates and the frozen protocol pass.

## 2026-10-31 Track B gate audit

### Forecast target and protocol status

The proposed primary forecast uses each patient's first two eligible chronological masks as history, with the second visit as origin and the immediate third eligible mask as target. Issue the forecast at origin day 90; include a patient only when the observed second-to-third interval is 76–104 days. Do not substitute a later target. The target-window rule is explicit in `paper/track_b_pde_external_protocol_draft_2026-10-28.md`, but that document remains a candidate draft, not a frozen scoring protocol.

### MU-Glioma-Post gates

MU is publicly available under CC BY 4.0, but this project has already used it for patient-level PDE development and scoring. Keep all new MU results exploratory. The paper describes 1-mm resampling and rigid atlas registration, and describes the released tumor labels and an automated-segmentation workflow with review/refinement. It does not provide case-level longitudinal registration acceptance records or establish tumor-boundary alignment for each selected triplet. It also does not define label 0 as a valid negative throughout the complete NIfTI grid. The images are skull-stripped, so a whole-grid negative interpretation needs source-backed support guidance. Matching headers, image intensity, brain-mask overlap, or visual similarity do not clear these gates.

The read-only MU preflight record (`output/fk_mu_pipeline_verification_20261029.json` at sibling commit `0cf951b`) reports that 72 NIfTI headers were checked before the loader stopped at a missing hashed source-evidence document. It reports no mask-file hashing, no target-mask voxel reads, and no scoring. The source evidence still needed is: (1) case-level registration method/reference and QC or an accepted review procedure; (2) label-0 meaning and any invalid or unannotated regions; and (3) source-backed valid-support masks or instructions. The existing entire-grid support gate remains in force; do not relax it to obtain a score.

### Access and mask provenance

Duke files are restricted and require request approval under CC BY-NC-ND 4.0 or a custom agreement. Approval and school-use terms remain pending; no Duke files were downloaded. The dataset description says each patient has 3–18 sessions and an nnU-Net mask per session, but the authors' guidance says sessions are not registered to each other and scan dates may be off by a few days. Relative dates alone do not pass scan-date or registration QC.

CFB-GBM remains a separate annotation-change side study. A newer sibling checkout contains an MRI-only, paginated display aid for 38 cases: its manifest records 760 pages and 5,668 sampled slices. The rendering script records T1-Gd volumes as the only image inputs and marks GTV and brain-mask arrays as unopened. No independent reviewer has assessed the display set and no registration pass exists. Crosswalk metadata records mixed t0 mask provenance and generated t1/t2 candidate masks, with individual review status unknown. This establishes no Track B PDE result and grants no scoring permission. The image outputs remain in that sibling checkout and are not copied here.

### Fresh external cohort search

The 2026-10-31 search found no fresh, public, no-request cohort that meets the three-mask, timing, label, reuse-rights, and prior-exposure gates.

| Candidate checked | Metadata finding | Gate result |
| --- | --- | --- |
| UTSW-Glioma | Public TCIA CC BY 4.0 release with expert-reviewed subregion labels, but source describes pre-surgical MRI selection, not three longitudinal visits. | Fails longitudinal history/target requirement. |
| UCSF-ALPTDG | Expert voxelwise masks at exactly two consecutive post-treatment visits per patient; access requires agreeing to a data-use agreement. | Fails three-mask and no-request requirements. |
| SAILOR | High-grade glioma series with multiple visits and masks, but EBRAINS lists controlled access under a data-use agreement. | Fails no-request requirement. |
| MOTUM | 67 patients, including 29 high-grade gliomas; source describes one MRI study and segmentation pair per patient. | Fails longitudinal mask requirement. |
| UPENN-GBM | Public release describes pre-operative scans and follow-up scans for some second-resection cases, but does not establish three dated, labeled visits per patient. | Not shown to pass; metadata eligibility remains unresolved. |
| BrainMD | The NeurIPS release describes longitudinal MRI and records, but the full dataset requires a data-use agreement and the documentation does not provide tumor voxel masks. | Fails no-request and mask requirements. |
| CPTAC-GBM | TCIA describes follow-up imaging and preserves relative timepoint metadata, but its MRI package is under NIH controlled access and the collection page does not document tumor-mask coverage. | Fails no-request and label-coverage requirements. |
| UCSD-BMETS-Longitudinal | Public CC BY 4.0, two to four masked timepoints, but cohort is brain metastases. | Fails GBM disease requirement. |

Previously examined MU, LUMIERE, RHUH, UCSD-PTGBM, CFB-GBM, and BraTS material are not untouched confirmation. PREDICT-GBM is composed of TUM-GBM, LUMIERE, and RHUH cohorts already present in this project. No candidate above is recommended for patient scoring.

### Comparator implementation gate

The current checkout's `src/fk_mu_scoring.py` thresholds the hybrid and no-PDE inclusion fields separately at 0.5. A newer sibling checkout contains a scorer change that ranks no-PDE probabilities within support to match the hybrid point mask's voxel count, with stable tie handling. That change is not present in this checkout. Review and verify this protocol-alignment fix before any scoring; do not describe the current point comparator as volume-matched until then.

### Next exact gate and unsent author message

First send the authors the documentation request below, and wait for source-backed answers. In parallel, wait for documented Duke access approval and school/fair permission. Then complete metadata-only overlap, date, and eligibility checks; resolve registration, valid support, and mask-provenance evidence; review the comparator change; and freeze code, cohort, folds, tuning, target rule, metrics, and failure rules. Score only after every gate passes.

**Subject:** MU-Glioma-Post longitudinal registration and segmentation support documentation

Dear MU-Glioma-Post authors,

I am a student researcher studying patient-level reaction-diffusion tumor-growth forecasts from serial MRI. I am using the public MU-Glioma-Post release for exploratory development only and will follow the dataset and TCIA terms. I will not attempt to identify or contact any participant.

Could you point me to documentation that answers these questions?

1. Were each patient's visits registered independently to SRI24, or was there also visit-to-visit longitudinal registration? Are case-level transforms, registration-quality measurements, acceptance criteria, or review records available for released timepoints?
2. What does label 0 mean in the segmentation files? Does it mark a valid negative throughout the image grid, or can it include voxels outside the annotation or acquisition support?
3. Are per-timepoint brain/support masks or source-backed instructions available to define the valid region for longitudinal tumor-mask evaluation, including the relationship to skull-stripping and atlas space?

I do not need identifying information. Aggregate guidance or records that can be matched to anonymized released timepoint IDs would be sufficient. If these records do not exist, a brief clarification of the intended label and support interpretation would still help me keep the evaluation cautious.

Thank you,

[Name]

**Status:** draft prepared; not sent.

## Sources

- Primary dataset record and license: [Zenodo 10.5281/zenodo.20548528](https://zenodo.org/records/20548528).
- Dataset-specific registration and date guidance: [GBM-DIR-QA README](https://github.com/deshanyang/GBM-DIR-QA).
- Current Track B boundaries and CFB preflight: `START_HERE.md` and `paper/track_b_cfb_scoring_protocol_draft_2026-10-27.md`.
- Historical cohort and access audit: `README_STATUS.md` section 5a.1 and ledger entry 2bm.
- MU public access and dataset description: [TCIA MU-Glioma-Post](https://www.cancerimagingarchive.net/collection/mu-glioma-post/) and [dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12635165/).
- New candidate screen: [UTSW-Glioma TCIA](https://www.cancerimagingarchive.net/collection/utsw-glioma/), [UCSF-ALPTDG paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11294954/), [SAILOR EBRAINS record](https://search.kg.ebrains.eu/instances/cae85bcb-8526-442d-b0d8-a866425efff8?noSilentSSO=true), [MOTUM paper](https://doi.org/10.1038/s41597-024-03634-0), [UPENN-GBM TCIA](https://www.cancerimagingarchive.net/collection/upenn-gbm/), [BrainMD data-use terms](https://proceedings.neurips.cc/paper_files/paper/2024/file/b53513b83232116ae25f57a174a7c993-Supplemental-Datasets_and_Benchmarks_Track.pdf), [CPTAC-GBM TCIA](https://www.cancerimagingarchive.net/collection/cptac-gbm/), and [UCSD-BMETS-Longitudinal TCIA](https://www.cancerimagingarchive.net/collection/ucsd-bmets-longitudinal/).
