# Track A dataset audit — Core vs Periphery, donor-held-out (2026-10-31)

## Decision: NO-GO. Download nothing.
No public, no-contact cohort has 6 or more eligible paired donors in one coherent assay and cell population. The appendix rule (fewer than 6 eligible donors: feasibility report only) applies. The Track A question stays unchanged. No amendment is proposed here.

## What was done
- Read the appendix Track A section (eligibility rule: section 2, lines 41–49) and the prior audits `paper/trackA_v2_data_check.md` and `paper/trackA_v2_gate_report.md`.
- GEO web pages returned a robot check (CAPTCHA). It was not bypassed. Records were read through the official NCBI E-utilities API (`esummary`, db=gds) and the GEO FTP SOFT files.
- Downloaded 3 small metadata files only (git-ignored, `data/external/trackA_audit_meta/`):

| File | Source | Bytes | SHA-256 |
|---|---|---|---|
| GSE273274_family.soft.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273274/soft/ | 2,896 | 11832a094b8cc639f4453161712fd83f831bf6cb79607a8f57c5096f95f3bf36 |
| GSE162631_family.soft.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162631/soft/ | 2,069 | 4dd8a5c96bc1c739329dd1d3ecd4404c637e4a52ae7bb97f01ac1ecd371a5ae4 |
| GSE318562_family.soft.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE318nnn/GSE318562/soft/ | 4,874 | ab534c8d01ac2f229c62183dfc32be4913b310176a155d292776b928d805db87 |

- No expression file was downloaded. Nothing was scored. No person was contacted.

## Comparison

| Cohort | Assay / population | Eligible paired donors (verified) | Label source | Blocker |
|---|---|---|---|---|
| GSE273274 (PMID 39933531) | 10x scRNA-seq, "Glioma cells" | **≤3, donor identity unverified** | GEO field `tissue: GBM centre` / `GBM surround` | Samples are titled "replicate 1–3" (files R1/R2/R3). GEO does not state that the replicates are 3 patients. Supplement: `GSE273274_RAW.tar`, 706 MB. |
| GSE162631 (PMID 34228647) | 10x scRNA-seq, **CD31+ endothelial-enriched** (Dynabeads) | 4 (R1–R4, "4 glioblastoma patients" in design) | Sample title only | Endothelial cells only. Every sample, including "peripheral tissue", has `Sample_description = tumor core` (metadata conflict). Supplement: `GSE162631_raw_counts_matrix.zip.gz`, 806 MB. |
| GSE318562 (no PMID yet; released 2026-06-12) | 10x Visium V1 spots, mixed OCT and FFPE | 2 by GEO names (WU1221 Core2/Edge1, WU1227 Core5/Edge2) | Sample title only | Design says 22 sections, 9 patients, "edge-core axis", but 18 sections carry no region label on GEO. Spots, not cells. Fixation differs by batch. `GSE318562_RAW.tar`, 674 MB. |
| GSE237183 | Visium spots | 5 | Navigation-guided location | Fails donor and topic gates (`output/trackA_v2_data_check/alternate_cohort_screen_20261030.json`). |
| GBM-Space | snRNA / spatial H5AD | 4 (strict same-section core + leading edge) | Feature names only | Below 6; nucleus-to-region values not verified (`output/trackA_v2_data_check/s14_metadata_check.json`). |
| HRA009021 | scRNA-seq | 10 donors, map unknown | Paper tables by donor code only | Sample-to-region map not public; needs author contact, which the rule forbids (`paper/trackA_v2_gate_report.md`). |
| UCSC multiomic-gbm (historical) | sc multi-omic | 6–7 core+periphery | Source label | Region = library confound; already exposed. Leak-free C-GAT 0.638 vs cells-per-patient baseline 0.646 (`vault/claim_ledger.md` line 43; scripts 88, 95). Not validation. |

## Why the cohorts cannot be pooled
- GSE273274 + GSE162631 gives at most 7, and only if the GSE273274 replicates are distinct patients (unverified).
- The two studies profile different cells (glioma cells vs CD31+ endothelium). Region would be nested in study, and study is fully confounded with cell population and lab. A pooled classifier would learn the study, not the region.
- Visium spots (GSE318562, GSE237183) are a different unit and modality from dissociated cells. Mixing them changes the question.

## Go / no-go against appendix thresholds
- Fewer than 6 eligible paired donors in every coherent cohort → **stop as infeasible for the planned analysis.**
- No descriptive nested evaluation is run. No superiority claim is possible.

## Gaps (not checked)
- Not searched in this pass: CELLxGENE Census, HCA Data Portal, BioStudies/ArrayExpress single-cell, Zenodo/Figshare, HTAN, CPTAC, pathology image sets.
- GSE273274 and GSE318562 papers' supplementary tables were not read. GSE318562 has no linked paper yet.
- GEO data are public; reuse follows NCBI GEO policy and citation of the source paper. No file-level license was found. Not legal advice.

## What would change the decision (no contact needed)
1. A published GSE318562 paper whose supplement maps all 22 sections to core or edge per patient. That would allow a Visium-based section-level amendment (9 patients, still below 10: descriptive only).
2. A new public single-cell cohort with ≥6 paired core/periphery donors in one assay, found in the unsearched repositories above.
