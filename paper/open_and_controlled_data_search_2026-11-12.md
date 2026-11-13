# Open and controlled data search for Tracks A, B, C and integration (2026-11-12)

Metadata only. No image, expression value, outcome, or response value was read. Sources were read 2026-11-12. Small metadata files are in git-ignored `data/external/`.

## Need for each part

- Track A: tumor biology with a measured or credible drug-response link.
- Track B: a fresh longitudinal MRI cohort (at least 25 patients with three labelled post-treatment visits, per `paper/trackB_fresh_cohort_screen_2026-11-11.md`).
- Track C: real data to set the drug-suppression part of the simulator (not to evaluate policies; evaluation stays stopped).
- Integration: biology and MRI in the same patients, with treatment and outcome.

## New open candidates checked today

| Dataset | Verified facts | Open? | Verdict |
|---|---|---|---|
| **BraTS-TCGA-GBM** + TCGA-GBM (GDC) | 135 patients: processed co-registered pre-operative MRI with segmentations and radiomic features, **CC BY 3.0** ([TCIA](https://www.cancerimagingarchive.net/analysis-result/brats-tcga-gbm/)). Same patients have open TCGA RNA expression, MGMT, treatment and survival in GDC ([GDC](https://portal.gdc.cancer.gov/projects/TCGA-GBM)). Original TCGA-GBM DICOMs are controlled. | **Yes** | **Best open integration candidate.** Biology (bulk RNA) and MRI in the same patients, with treatment and survival. Limits: one pre-operative scan (not longitudinal); outcome is survival, not measured drug response; radiogenomics on TCGA is heavily published, so novelty needs a narrow, pre-registered question. |
| **WAW-Glioma** (Bartnik et al. 2026) | 192 adults, grades 2-4 (106 grade 4); 156 pre-op and 161 post-op exams; manual segmentations; molecular markers (IDH, MGMT, EGFR, CDKN2A, chr7/10); progression and survival; **CC BY 4.0, open download** ([Zenodo 22205631](https://doi.org/10.5281/zenodo.22205631)). Never used by this project. Post-op exams per patient (from file names): 1: 56, 2: 12, 3: 19, 4: 2, 5: 3. Grade 4 with 3+ post-op exams: **22**. | **Yes** | Fresh and open. **Fails Track B's 25-patient, three-visit gate narrowly (22 grade 4; 24 all grades).** Usable for a smaller, explicitly exploratory study, or for an A-B marker link (molecular markers, not expression). |
| UPENN-GBM | 630 patients, CC BY 4.0; from the open availability table: 41 patients with two timepoints, 0 with three ([TCIA](https://www.cancerimagingarchive.net/collection/upenn-gbm/)). | Yes | Fails Track B (no three-visit patients). |
| GBM-MR-NER-Outcomes | 45-52 TCGA patients, baseline features, PFS/OS, CC BY 3.0 tables | Yes | Small, baseline only. |
| Spitzer/Nomura et al. 2025, GSE274546 | 59 patients, paired primary and recurrent snRNA-seq with treatment and MGMT ([paper](https://www.nature.com/articles/s41588-025-02168-4)) | Yes (processed) | Track A longitudinal biology under treatment. No MRI, no measured drug response. |
| HGCC (done) | 83 lines, measured TMZ AUC | Yes | Used: the 119-gene panel showed no association (`paper/trackA_hgcc_tmz_102_result_2026-11-12.md`). Lines could serve a different, new pre-registered question. |
| ACRIN 6677 tables (done) | 123 patients, central-reader areas, arm, outcomes | Yes | Used for Track B (failed). Can still give Track C real treatment-response ranges (exploratory calibration, not evaluation). |

## Controlled sets worth requesting (ranked)

| Rank | Dataset | Gives | How to get it | Who can file |
|---|---|---|---|---|
| 1 | **BraTS 2024 post-treatment glioma** (Synapse syn59059776) | Large multi-site post-treatment MRI with labels; may include multiple visits | Synapse account plus data-use terms. **No dbGaP sponsor stated** | Likely the user directly (check terms; under-18 rules may apply) |
| 2 | **UCSF-ALPTDG** | Expert voxel masks at two consecutive post-treatment visits | Data-use agreement ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11294954/)) | Check whether an individual can sign |
| 3 | Burdenko-GBM-Progression | 180 patients, 1-8 follow-ups, MGMT, response labels | dbGaP (phs004225) | Needs a faculty sponsor |
| 4 | Ivy GAP images | Region RNA plus MRI, 10 tumors | dbGaP | Needs a faculty sponsor |
| 5 | CPTAC-GBM images | Proteogenomics plus MRI with follow-up | dbGaP | Needs a faculty sponsor |
| 6 | Moniz Garcia 2026 (6 patients) | Biology, MRI, measured TMZ tests | Email authors (draft exists) | User |

## Recommended next step

Two open studies can start without any request:

1. **Integration (A + B, same patients):** BraTS-TCGA-GBM MRI features plus TCGA RNA plus treatment and survival. Needs a narrow pre-registered question and a prior-art check first, because TCGA radiogenomics is crowded.
2. **Track B (exploratory, fresh):** WAW-Glioma post-operative series (22 grade-4 patients with 3+ exams), labelled exploratory because it is below the 25-patient floor.

In parallel, the user can try rank 1 (BraTS 2024 post-treatment on Synapse), which may not need a sponsor.

None of these gives measured drug response in patients. That remains unavailable in open data.
