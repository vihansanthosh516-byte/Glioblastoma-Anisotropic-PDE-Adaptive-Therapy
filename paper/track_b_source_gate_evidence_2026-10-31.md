# Track B source-gate evidence record (2026-10-31)

> **HISTORICAL for the MU/Duke path (updated 2026-11-01).** The gates below apply to MU-Glioma-Post and Duke only. They are not the LUMIERE gates. Current Track B work uses LUMIERE (exploratory only); its gates are in `paper/track_b_no_contact_dataset_plan_2026-10-31.md` §8 and `paper/track_b_readiness_audit_2026-11-01.md`.

**Metadata and published documents only.** No target-mask voxel was read. No restricted file was downloaded. No author or school was contacted. No patient mask was scored. Amendment A-2026-10-31-R was a draft when this record was written; the user adopted it on 2026-11-01 (reuse only; `paper/track_b_amendment_decision_2026-10-31.md`). The 64³ crop stays out of the required U-Net comparison.

## Sources checked

| ID | Source | Accessed |
|---|---|---|
| S1 | Mahmoud et al., "MU-Glioma Post: A comprehensive dataset of automated MR multi-sequence segmentation and clinical features", *Sci Data* 12, 1847 (2025), doi:10.1038/s41597-025-06011-7. Full text from Europe PMC, PMC12635165 (`fullTextXML`). | 2026-10-31 |
| S2 | TCIA collection page, https://www.cancerimagingarchive.net/collection/mu-glioma-post/ (Version 1, updated 2025/03/21). | 2026-10-31 |
| S3 | Duke Zenodo record https://zenodo.org/records/20548528 and the GBM-DIR-QA README (already in `paper/track_b_external_cohort_screen_2026-10-28.md`). Not re-fetched. | 2026-10-28 |
| S4 | Local preflight record `output/fk_mu_pipeline_verification_20261029.json` (as cited in the cohort screen). | — |

## Gate status

| Gate | Evidence found | Status |
|---|---|---|
| **Cavity-label meaning (MU)** | S1 and S2: "Label 4: Resection Cavity (RC)". It covers post-surgical changes, including recent blood products and air foci and chronic CSF-isointense material. S1: automated labels were manually validated and refined by neuroradiologists, and cavity errors were corrected. | **Definition: source-backed (new).** The appendix rule (core = labels 1 or 3; cavity 4 excluded from the core) matches the source definitions. Limit: one released mask per visit, no inter-rater data; the cavity edge is not error-free truth. |
| **MU label-zero meaning** | S1 and S2 define labels 1–4 only. Label 0 / background is not stated. Images are skull-stripped (S1: deep-learning brain extraction). | **Open.** No source says label 0 is a valid negative over the whole grid or inside the brain. |
| **Valid target support (MU)** | No brain or support mask is released (S2). S1 describes skull-stripping only; no extraction mask or QC is listed. | **Open.** Nonzero-MRI support stays withdrawn: no source establishes it. |
| **MRI-only registration evidence (MU)** | S1: each image was resampled to 1 mm³ and "co-registered to … SRI24 via rigid-body transformation using the CapTK/Greedy registration tool". S2 does not say whether timepoints were registered to each other. No per-case transform, QC metric or acceptance record is released. | **Open.** Atlas alignment per timepoint is documented. Visit-to-visit alignment and its QC are not. |
| **Duke access** | S3. No local approval record found (`paper/`, `data/manifests/`, `vault/` checked by file search). The user reports a resubmission on 2026-10-29. | **Open (pending).** |
| **School-use approval** | No written record on file. User report (2026-10-31, contact log below): the school said this use is OK. | **Open** until the written approval (date, approver, conditions) is saved in the project. |

**New fact, not a gate:** the S1 article text is CC BY-NC-ND 4.0. The data on TCIA are CC BY 4.0 (S2). Do not copy article figures or text into project outputs.

**Source inconsistencies seen in S2 (for the record):** the 3T scan count appears as both n=318 and n=261 at 45.08%; the 2D/3D shares sum to 100.10%. These do not affect the gates.

## Exact missing items

1. MU: the meaning of label 0, and whether any region of the grid is unannotated or invalid.
2. MU: a per-timepoint valid-region (brain/support) mask, or a written rule for one.
3. MU: whether visits were registered to each other, or only to SRI24 one by one; and any per-case registration QC or acceptance record.
4. Duke: written access approval and the license that applies.
5. School: written approval for this use of the selected data and license.

The no-scoring rule stays. The runner still refuses to score until the loader records the registration and target-support gates as passed.

## Contact log

- 2026-10-31, user report: the school said this use is OK. No written copy is on file yet, so the school gate stays open until it is saved.
- 2026-10-31, user report: the user contacted the MU authors and Duke about one week earlier. No reply yet. No written evidence has arrived, so no gate changes. The school has not been contacted.

## Pass evidence per gate

A gate passes only on written evidence that meets its rule. Record the date, sender and file of each reply. An informal or partial answer leaves the gate open.

| Gate | Pass evidence | Not enough |
|---|---|---|
| Cavity label meaning | Passed (S1, S2). Limit: one released mask per visit, refined by neuroradiologists, no inter-rater data, so the cavity edge is not error-free truth. | — |
| MU label 0 | A written statement from the dataset authors or TCIA that defines label 0, and says where in the grid it is a reviewed "no tumor" value. | Our own reading of masks; BraTS convention by analogy. |
| Valid target support | A released per-timepoint support/brain mask, or a written support rule from the authors or TCIA that we can apply without target labels. | Nonzero MRI voxels, image-intensity rules, atlas brain masks, or any rule we define ourselves. |
| Registration evidence | A written statement of how the visits of one patient relate (only per-visit rigid to SRI24, or also visit-to-visit), plus per-case QC records or a documented accepted review procedure. | Header match, visual overlap, or a registration quality score we compute ourselves. |
| Duke access | A written approval from Duke that names the license, the allowed use, and the rules for derived data. | A pending request; a Zenodo page alone. |
| School approval | A written approval from the school for this secondary-data use, or a written statement that no review is needed, with any conditions. | Verbal approval. |

## Draft messages (reviewed 2026-10-31; ready for the user to send; NOT sent)

Corrections from the earlier draft:
- MU question 2 no longer offers nonzero MRI voxels as a support rule. That rule is withdrawn and must not be suggested.
- MU now asks for per-case QC records directly, and for a written rule if no mask exists.
- Duke now asks for use and derived-data conditions.
- School now asks what review is needed, not only for a yes or no.

### A. To the MU-Glioma-Post authors (contact via the corresponding author of S1 or the TCIA helpdesk)

> **Subject:** MU-Glioma-Post (TCIA, Version 1): label 0, valid region and longitudinal registration
>
> Dear Dr. Mahmoud and colleagues,
>
> I am a student doing independent computational research with MU-Glioma-Post (TCIA, Version 1, updated 2025/03/21). Thank you for releasing it. Your *Scientific Data* paper (12, 1847, 2025) defines labels 1–4. I want to evaluate tumor-growth forecasts between visits, and I have three questions.
>
> 1. **Label 0.** What does label 0 mean in the segmentation files? Is it a reviewed "no tumor" value everywhere in the image grid, only inside the brain, or can it include regions that were not reviewed?
> 2. **Valid region.** Is a per-timepoint brain or valid-region mask available (for example, the skull-stripping mask)? If not, is there a rule you recommend for the region where the segmentations can be treated as complete?
> 3. **Registration between visits.** The paper says each image was rigidly registered to SRI24 with CaPTk/Greedy. Was each timepoint registered to SRI24 on its own, or were a patient's timepoints also registered to each other? Are per-case registration quality checks or acceptance records available, or is there a documented review procedure?
>
> Short answers, or a pointer to documentation, are enough. I do not need any identifying information.
>
> Thank you,
> [name]

### B. To the Duke dataset contact (via the Zenodo access-request thread for record 20548528)

> **Subject:** Access request status and conditions: Longitudinal Glioblastoma MRI Dataset (Zenodo 20548528)
>
> Hello,
>
> I resubmitted my access request on 2026-10-29 with the research purpose and license information you asked for. Could you please tell me:
>
> 1. The current status of the request.
> 2. Which license applies if it is approved: CC BY-NC-ND 4.0 or the custom Duke license.
> 3. The conditions on use and on derived data. For example: may I report summary statistics and model results computed from the data? May I keep model weights trained on it? May I show derived images or masks in a written report?
>
> I will not download any file before written approval.
>
> Thank you,
> [name]

### C. To the school (research coordinator or the person who approves student research)

> **Subject:** What review is needed for computational research on de-identified MRI datasets?
>
> Hello,
>
> I am doing independent computational research that uses existing, de-identified brain MRI datasets. No new data is collected and no people are contacted. The datasets are:
>
> 1. MU-Glioma-Post from The Cancer Imaging Archive, public under CC BY 4.0.
> 2. If access is granted, a Duke longitudinal glioblastoma MRI dataset under CC BY-NC-ND 4.0 or a Duke license.
>
> I will analyse the data only on my own computer and will not share the data or derived images. Could you tell me:
>
> 1. What review this secondary use of de-identified data needs (for example, an IRB or school review), or whether none is needed.
> 2. Whether you can give written approval, and any conditions on it.
>
> Thank you,
> [name]
