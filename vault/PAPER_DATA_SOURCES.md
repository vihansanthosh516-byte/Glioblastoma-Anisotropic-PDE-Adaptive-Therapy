> **HISTORICAL BACKGROUND (corrected 2026-10-19):** Old-track inventory, not current RANO dataset guidance. Workbook `data/tcia/MU-Glioma-Post_ClinicalData-July2025.xlsx` exists (101568 bytes), contrary to the audit absence claim. Stored data presence does not establish provenance or clinical eligibility. External download dates/sizes and unkeyed access/license statements remain unverified without dated inventories. No input was downloaded or analysis rerun.

# Paper: Data Sources (the four real cohorts)

Written 2026-10-04. Paths/counts were checked where available; this is not a complete provenance audit. Recheck input versions and eligibility before reuse.

| # | Cohort | What it is | Size used | Where on disk | Used by |
|---|---|---|---|---|---|
| 1 | **multiomic-gbm** (UCSC Cell Browser) | Single-cell and spatial GBM atlas. Source paper: Nature Communications, DOI 10.1038/s41467-026-69716-2 (source and attribution corrected in `archive_old_gbm/docs/dataset_info.md`). | Full atlas: 223,113 scRNA-seq cells, 25 patients. Script 01 keeps **140,355 cells** (Core 80,828 / Periphery 47,547 / Healthy 11,980). Scripts 04-16 use a **15,000-cell subsample**. | `output/01_filtered_three_class.h5ad` (559,657,115 bytes, about 560 decimal MB; checked 2026-10-19). Untracked external input. Full-atlas counts/size require the source release; filtered class counts are in `output/01_class_counts.tsv` and latent array headers corroborate retained totals. | Track A, scripts 01-34 |
| 2 | **TCGA-GBM** | Clinical data from cBioPortal; expression from Xena HiSeqV2. | 518 clinical patients. 150 patients have both expression and survival. | `data/tcga_gbm_clinical.csv`, `data/tcga_gbm_expression.tsv` | Track A, scripts 35-41 |
| 3 | **MU-Glioma-Post** (TCIA) | Longitudinal post-treatment MRI, clinical data (MGMT, TTP), tumour masks. | 203 raw patients. 154 with fitted growth. 103 in the stromal/shape workflow (exact inclusion/skip provenance missing here); 152 initial forecast pairs, with 133 seedable cellular-core cases and 19 empty-core skips (`output/forecast_labels_core/results.json:n_scored,n_skipped`), a distinct workflow. 130 GBM patients with MGMT for script 82. | `data/tcia/MU-Glioma-Post/`, `data/tcia/MU-Glioma-Post_ClinicalData-July2025.xlsx`, `output/mu_glioma_params_real.csv` | Track B and C |
| 4 | **UCSF-PDGM** (TCIA) | Pre-operative diffusion MRI. Historical grade 2-3 IDH-mutant attribution is unverified without its clinical eligibility manifest. | 62 result rows per mode (42 train / 20 test, seed 42; `output/cv_split.json` and three result CSVs). | `output/results_standard.csv`, `results_atlas.csv`, `results_dki.csv`, `cv_*.csv`. Imaging is not in git. | Tensor-construction study (see [[UCSF-PDGM-Tensor-Study]]) and DTI atlas |

## Facts the paper must state
- **MU-Glioma-Post initial forecast uses one input/target scan pair per patient, not only one MRI per patient.** The forecast scores scan 1 to scan 2. Details in [[NEGATIVES_REVISITED]].
- **UCSF-PDGM has one scan per patient.** It cannot give a forecast target. It supplies tensors and a shape test only.
- **IvyGAP lineage is not established here.** Script 40 describes molecular subtype as a zone proxy; exact redistribution/license and input lineage need the retained manifest. Do not call this measured spatial-zone validation.
- **BraTS 2021 data exists in `data/brats/`.** No result use was identified in reviewed vault notes; this is not an exhaustive source-code lineage check. Do not cite it as a cohort unless a script is found that uses it.
- **Track C RL/PPO "synthetic" sets** (cohort64, lhs60, synth_test) are model-generated. `real_test` (21 patients) uses real MU-Glioma growth rates but the same assumed kill rates for every patient. The paper must label which is which.
- **Track A to Track B link is absent in the recorded analysis (ledger A24).** The recorded score is 1.0 for all 61; exact per-patient artifact provenance remains missing here.
- **External cohorts (downloaded 2026-10-04, in `data/external/`, git-ignored):** CGGA mRNAseq_693 and mRNAseq_325 clinical and expression (WHO IV usable: 237 and 137); LUMIERE (91 patients; pyradiomics volume tables from HD-GLIO-AUTO and DeepBraTumIA, RANO ratings, demographics); RHUH-GBM (40 patients; clinical CSV and 2.9 GB NIfTI with segmentations; no MGMT column). Used by scripts 90 and 91; RHUH historical not-analysed status is superseded by H-2 scoring/H-4 training (ledger 2ax/2bf); both old test sets are spent. Burdenko open clinical CSV was checked (2az); image access is a separate unresolved item. Sources: CGGA (cgga.org.cn), LUMIERE (figshare), RHUH-GBM (TCIA).

## Cohort counts that look alike (do not mix)
| Number | Meaning and source |
|---|---|
| 203 | raw MU cohort, `data/manifests/consort_mu.json:mu`; directory count is not the eligibility rule |
| 154 | retained rows in `output/mu_glioma_params_real.csv` and inverse summary 154; R2 filter provenance not established by row count alone |
| 152 | initial forecast pairs (`output/forecast_labels_core/results.json`, 133 scored and 19 skipped) |
| 130 | MGMT/TTP analysis, `output/mgmt_resistance/results.json:clinical.n` (79+51) |
| 111 | patients across the 4 Track C evaluation sets (20 + 30 + 21 + 40; the 111 in the vault note) |
| 103 | stromal CSV 103 rows and script 74 shape cohort; not the forecast eligibility denominator |
| 61 | adaptive workflow positive-rho subset, `output/adaptive_geometry_metrics.json`, ledger 2j; 42 negative-rho cases excluded |
