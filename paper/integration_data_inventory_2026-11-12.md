# Integration study: data inventory (2026-11-12)

IDs and column names only. No survival value, expression value, or image intensity was read for this record.

## Linking file

`data/external/brats2021_map/BraTS2021_MappingToTCIA.xlsx` (TCIA, RSNA-ASNR-MICCAI-BraTS-2021 page; git-ignored). It maps each BraTS 2021 case to its source collection and patient ID. BraTS 2021 images and expert tumor masks (enhancing, necrotic core, edema) are already on disk in `data/brats/` (1,266 cases; training cases have masks).

## Patients with pre-operative MRI + tumor masks on disk, by source

| Source | In BraTS 2021 | On disk with masks | Other data | Use |
|---|---|---|---|---|
| TCGA-GBM | 135 | 102 | 91 also in TCGA U133A microarray (539 samples; `data/external/tcga_gbm_u133a/HT_HG-U133A.gz`, UCSC Xena, open). Only 20 are in the older RNA-seq file `data/tcga_gbm_expression.tsv`. Clinical: `data/tcga_gbm_clinical.csv`. | **Development cohort: MRI + RNA + survival** |
| UPENN-GBM | 447 | 403 | Survival, age, sex, IDH1, MGMT, KPS, resection (`data/external/upenn_meta/clin.csv`). No RNA. | **External test of MRI-based parts** |
| CPTAC-GBM | 39 | 33 | RNA and clinical in GDC (CPTAC-3); availability not yet checked. | Possible external MRI + RNA test |
| IvyGAP | 34 | 30 | Regional RNA (`data/external/ivygap/`), survival, MGMT, subtype (`tumor_details.csv`). | Possible external MRI + RNA test (regional RNA, small) |

## Prior exposure in this project

- TCGA-GBM survival was analysed before with clinical covariates (age, sex, subtype; ledger A3, n = 518) and with 4 genes (script 38). Age is a known strong factor there. MRI features were never tested against survival in this project.
- UPENN-GBM survival: no use found in `src/` (search for `upenn_meta/clin` and `Survival_from_surgery`). Treated as fresh.
- CPTAC-GBM and IvyGAP survival: IvyGAP clinical was ingested by historical script 35; treat IvyGAP survival as exposed. CPTAC unused.

## Size limits

The TCGA development cohort is small (91). Any survival model built on it alone will have wide uncertainty. The main strength is UPENN-GBM (403): a large, fresh, independent cohort for any MRI-based component.
