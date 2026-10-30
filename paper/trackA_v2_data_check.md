# Track A v2 data and label audit

**Status:** metadata screening only. No patient expression values were downloaded, opened, or scored. The GSA workbook confirms donor/sample/run links. The approved paper supplement gives HBI/LBI/Tumor labels for the same 10 donor codes, but not for each GSA sample name. Evidence is in `output/trackA_v2_data_check/hra009021_metadata_audit.json` and `output/trackA_v2_data_check/hra009021_supplement_audit.json`.

## HRA009021: conditional main-test cohort

The article reports 10 patients, 30 single-cell tissues, and 209,762 cells across LBI, HBI, and tumor core (`output/trackA_v2_data_check/dataset_candidate_inventory.json`). HBI/LBI regions were selected using preoperative MRI blood-flow maps, reviewed by radiologists and a neurosurgeon, sampled with image-guided surgery, and checked with histology. These labels do not come from gene expression.

The user-supplied `HRA009021.xlsx` is 17,181 bytes (SHA-256 `0aabd48105df82ad705ac4ebfc5482e9233602b43616834c9bceb7a900bb7f03`). It contains 10 individual records, 30 sample records, 30 experiment records, and 30 run records. All 30 sample-to-donor-to-run links are complete. It has no region-label columns and no HBI/LBI/core label values (`output/trackA_v2_data_check/hra009021_metadata_audit.json`). Thus the region-to-sample map remains unverified.

The approved Supplementary Information PDF is 3,127,913 bytes. Supplementary Table 4 lists 10 HBI, 10 LBI, and 10 Tumor sample labels. Its 10 donor codes match all 10 GSA donor codes. The PDF contains none of the GSA sample names `GBMxx01/02/03`, so it does not show which sample belongs to which region (`output/trackA_v2_data_check/hra009021_supplement_audit.json`).

The workbook's `Download Read file2` link points to the first-read file for all 30 runs. Do not use those second-read links. The separate raw archive inventory lists distinct first- and second-read URLs (`output/trackA_v2_data_check/hra009021_file_inventory.json`).

The archive lists 60 raw FASTQ files. Their measured total is 764,601,130,806 bytes. The exact name, URL, and size for each FASTQ and the checksum file are in `output/trackA_v2_data_check/hra009021_file_inventory.json`. No raw reads or expression matrix have been downloaded.

The main weak point is that the study reports exactly the minimum donor count, and no source links GSA sample IDs to regions. The full raw FASTQ set is large. One possible next file is the article Source Data workbook (`41467_2024_55243_MOESM4_ESM.xlsx`, 50,969,169 bytes; URL in `output/trackA_v2_data_check/hra009021_supplement_manifest.json`). It may not contain the crosswalk and may include expression-derived values. It has not been downloaded. Do not download any file without the user's explicit OK.

## GBM-Space section audit

The user's `GBM_space_snRNA.h5ad.gz` has 12 donor IDs and 57 donor/site combinations in its cell metadata (`output/trackA_v2_data_check/metadata_check.json`). Table S14 lists 55 single-nuclei donor/site blocks (`output/trackA_v2_data_check/s14_metadata_check.json`). Neither count establishes core/edge eligibility.

Table S14 ties histopathology tags to serial Visium/H&E sections. The nuclei were pooled from thick sections. There is no verified crosswalk from those tags to each individual nucleus. A conservative site screen finds one possible donor with both core-only and edge-only sites. This is not an author-validated nucleus label. Weakest point: the possible pair depends on our tag grouping and lacks a cell-level map.

The supplement workbook is 309,634 bytes with SHA-256 `15cf8a2e7d75e11ad9abbabf0ed00af7dd3f584903072116c09b92b2469ef423`; an earlier range check reported 136,557 bytes (`output/trackA_v2_data_check/s14_metadata_check.json`). The mismatch remains unresolved. Only Index, Table S1, and Table S14 were read. No expression sheet or matrix was read.

## Other screened datasets

- GSE286419 reports 6 snRNA-seq donors, with 4 core/margin assignments concordant with the authors' spatial check. It is below the donor minimum; weakest point: only 4 donors remain after the concordance check (`output/trackA_v2_data_check/dataset_candidate_inventory.json`).
- GSE316917 reports 9 patients and spatial samples; its dissociated single-cell reference does not have independently sourced region labels for the main test. Weakest point: the spatial samples are not independent donors (`output/trackA_v2_data_check/dataset_candidate_inventory.json`).
- E-MTAB-17183 has 12 donors, but no verified per-nucleus region mapping. Weakest point: pooled nuclei lack a cell-to-region crosswalk (`output/trackA_v2_data_check/dataset_candidate_inventory.json`).

Counts, source links, label descriptions, terms, and file sizes are in `output/trackA_v2_data_check/dataset_candidate_inventory.json`. The other shortlisted options fail the paired-region gate. HRA009021 reaches the numerical minimum only conditionally; weakest point: its crosswalk is unchecked.

## Prior art and next gate

PubMed, Europe PMC, and Google Scholar screens found prior work on HBI/LBI cell composition and on hybrid GBM cell states. The proposed paired hybrid-state comparison is only a candidate, not a novelty claim. The three closest papers and differences are in `output/trackA_v2_data_check/prior_art_review.json` and `paper/trackA_v2_proposal.md`.

Next gate: resolve mapping from GSA sample names to HBI/LBI/Tumor labels. Before scoring, freeze the full protocol, gene lists, settings, seeds, and failure rules; hash and commit them. Use `src/donor_cv.py` so each donor stays wholly in one fold. Do not join Track A donors with Track B or Track C.
