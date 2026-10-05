# Paper: Data Sources (the four real cohorts)

Written 2026-10-04. Every row was checked against a file on disk.

| # | Cohort | What it is | Size used | Where on disk | Used by |
|---|---|---|---|---|---|
| 1 | **multiomic-gbm** (UCSC Cell Browser) | Single-cell and spatial GBM atlas. Source paper: Nature Communications, DOI 10.1038/s41467-026-69716-2 (as stated in `docs/dataset_info.md`). | Full atlas: 223,113 scRNA-seq cells, 25 patients. Script 01 keeps **140,355 cells** (Core 80,828 / Periphery 47,547 / Healthy 11,980). Scripts 04-16 use a **15,000-cell subsample**. | `output/01_filtered_three_class.h5ad` (560 MB). Not in git. External dataset, about 1.4 GB. | Track A, scripts 01-34 |
| 2 | **TCGA-GBM** | Clinical data from cBioPortal; expression from Xena HiSeqV2. | 518 clinical patients. 150 patients have both expression and survival. | `data/tcga_gbm_clinical.csv`, `data/tcga_gbm_expression.tsv` | Track A, scripts 35-41 |
| 3 | **MU-Glioma-Post** (TCIA) | Longitudinal post-treatment MRI, clinical data (MGMT, TTP), tumour masks. | 203 raw patients. 154 with fitted growth. 103 with usable masks and scan pairs. 152 scan pairs for the forecast (133 scored). 130 GBM patients with MGMT for script 82. | `data/tcia/MU-Glioma-Post/`, `data/tcia/MU-Glioma-Post_ClinicalData-July2025.xlsx`, `output/mu_glioma_params_real.csv` | Track B and C |
| 4 | **UCSF-PDGM** (TCIA) | Pre-operative diffusion MRI. WHO grade 2-3 IDH-mutant astrocytoma (per README). | 62 patients (42 train / 20 test, seed 42). | `output/results_standard.csv`, `results_atlas.csv`, `results_dki.csv`, `cv_*.csv`. Imaging is not in git. | Tensor-construction study (see [[UCSF-PDGM-Tensor-Study]]) and DTI atlas |

## Facts the paper must state
- **MU-Glioma-Post has only one MRI per patient pair for the forecast.** The forecast scores scan 1 to scan 2. Details in [[NEGATIVES_REVISITED]].
- **UCSF-PDGM has one scan per patient.** It cannot give a forecast target. It supplies tensors and a shape test only.
- **IvyGAP spatial zones are not redistributed.** Script 40 uses molecular subtype as a stand-in for zone. Say so.
- **BraTS 2021 data exists in `data/brats/`.** Nothing in the vault says a result uses it. Do not cite it as a cohort unless a script is found that uses it.
- **Track C RL/PPO "synthetic" sets** (cohort64, lhs60, synth_test) are model-generated. Only `real_test` (21 patients) uses real MU-Glioma parameters. The paper must label which is which.

## Cohort counts that look alike (do not mix)
| Number | Meaning |
|---|---|
| 203 | patients in the raw MU-Glioma directory |
| 154 | patients with a fitted growth rate (R2 > 0) |
| 152 | patient pairs in the forecast study (script 81 scores 133) |
| 130 | GBM patients with MGMT status (script 82) |
| 111 | patients across the 4 Track C evaluation sets (20 + 30 + 21 + 40; the 111 in the vault note) |
| 103 | patients with usable masks (scripts 43, 74) |
| 61 | patients with positive fitted growth rate (script 44) |
