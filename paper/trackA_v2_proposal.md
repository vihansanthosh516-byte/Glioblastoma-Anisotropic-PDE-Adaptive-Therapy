# Track A v2 proposal

**Status: proposal only.** No model has been trained and no patient expression has been scored. HRA009021 is the only screened cohort that reaches the 10-donor minimum on the article's report. Its run-to-donor-to-region mapping is not yet verified. Evidence: `output/trackA_v2_data_check/dataset_candidate_inventory.json`, `output/trackA_v2_data_check/hra009021_file_inventory.json`, and `output/trackA_v2_data_check/proposal_settings.json`.

## Dataset shortlist

| Candidate | Donors and labels | Terms and files | Decision and weakest point |
|---|---|---|---|
| **HRA009021 — recommended, conditional** | The paper reports 10 patients and 30 single-cell tissue samples across LBI, HBI, and tumor core. HBI/LBI labels use preoperative MRI blood-flow maps, radiologist and neurosurgeon review, image-guided surgery, and histology. Labels do not come from gene expression. GSA sample records expose individual IDs, but the full run-to-donor-to-region crosswalk still needs checking. | GSA-Human marks it Open access; no separate file-level reuse license was found. The raw archive has 60 FASTQ files totaling 764,601,130,806 bytes. The required workbook, `HRA009021.xlsx`, has no published transfer size. Exact URLs and sizes are in `output/trackA_v2_data_check/hra009021_file_inventory.json`. | It reaches the minimum only if the crosswalk confirms 10 paired HBI/LBI donors. Weakest point: exactly the minimum donor count and an unverified run map (`output/trackA_v2_data_check/dataset_candidate_inventory.json`). |
| GBM-Space, E-MTAB-17183 | 12 donors. Histology tags are from serial Visium/H&E sections, but pooled nuclei lack a verified nucleus-to-region map. One donor was a possible pair in a conservative site screen (`output/trackA_v2_data_check/dataset_candidate_inventory.json`). | Public record; file reuse terms are not confirmed. `GBM_space_snRNA.h5ad.gz`, 11,267,029,458 bytes, is already on the user's computer (`output/trackA_v2_data_check/dataset_candidate_inventory.json`). | Fails the label gate. Weakest point: no verified nucleus-level region label. |
| Pai et al., GSE286419 | 6 donors reported; 4 had region assignments concordant with the spatial check. Core and margin were sampled using MRI and checked by neuropathology (`output/trackA_v2_data_check/dataset_candidate_inventory.json`). | Public GEO files; no file-specific license was found. `GSE286419_RAW.tar`, 583,915,520 bytes (`output/trackA_v2_data_check/dataset_candidate_inventory.json`). | Fails the donor minimum. Weakest point: only 4 donors pass the conservative region check. |

GSE316917 was also screened. It has 9 patients and spatial samples, but does not provide a 10-donor region-labeled dissociated single-cell cohort (`output/trackA_v2_data_check/dataset_candidate_inventory.json`). No new patient-expression files have been downloaded.

## Bar 1: method

**Base model:** scVI, a published neural generative model for single-cell counts. The scvi-tools package uses the BSD-3-Clause license; the exact links are recorded in `output/trackA_v2_data_check/prior_art_review.json`.

**One change:** add a gradient-reversal donor adversary to scVI's latent space.

1. scVI models single-cell count data with a neural network.
2. The donor adversary tries to remove donor identity from the learned representation while the region head predicts HBI versus LBI.
3. This targets the old donor-leakage failure. The adversary could also remove real biology, so it must beat the original model and an ablation.

Use the same nested donor-held-out folds from `src/donor_cv.py` for original scVI, adversarial scVI, the adversary-off ablation, and every baseline. Keep all samples and cells from a donor together. Fit preprocessing and feature selection inside the training folds only. Aggregate cell predictions to sample scores.

**Metric:** donor-paired concordance. For each held-out donor, compare the HBI score with the LBI score; then give each donor equal weight. This keeps the evaluation unit as the donor.

**Baselines:** prevalence; the old cells-per-donor rule rebuilt from its source; sample cell count alone; cell-type composition; and regularized pseudobulk classification. The old result is motivation only and is not directly comparable to this new endpoint. The original and improved models must use identical donor splits and the same metric.

**Bar 1 passes only if:** at least 10 donors have verified paired HBI/LBI labels; adversarial scVI beats original scVI and every baseline on the same held-out donors; the Holm-adjusted donor-bootstrap 95% interval for every difference is above zero; and disabling the adversary returns performance to within ±0.02 of original scVI. The draft uses 10,000 donor bootstrap resamples. These settings are in `output/trackA_v2_data_check/proposal_settings.json` and are not sealed. CPU runtime is a risk; the article reports 209,762 cells and the raw FASTQs total about 765 GB (`output/trackA_v2_data_check/dataset_candidate_inventory.json`, `output/trackA_v2_data_check/hra009021_file_inventory.json`). Weakest point: only 10 donors are reported, so any uncertainty interval may be wide.

## Bar 2: biology

**Question:** among malignant cells from paired LBI and HBI peritumoral samples, does the fraction of cells co-expressing at least two fixed Neftel cell-state programs differ within donor after accounting for hypoxia, normal-brain content, and donor composition?

The improved scVI representation will support cell-state scoring. Use fixed published Neftel state gene lists and a hybrid-state rule frozen before scoring. Compute one hybrid-cell fraction per sample and compare HBI with LBI within each donor. Repeat with direct fixed-signature scoring as a sensitivity check. Adjust for hypoxia and normal-brain content, and keep donor paired throughout.

### Prior art check

PubMed, Europe PMC, and Google Scholar were searched on 2026-10-30. Search terms and scope are in `output/trackA_v2_data_check/prior_art_review.json`. These searches were not exhaustive.

1. **Wang et al. (2024), “Turning attention to tumor-host interface and focus on the peritumoral heterogeneity of glioblastoma.”** This is the closest region-specific paper. It compares HBI and LBI and reports malignant-cell proportions and expression differences. The proposed hybrid-cell fraction is a different endpoint, but the regional overlap is substantial.
2. **Nomura et al. (2025), “The multilayered transcriptional architecture of glioblastoma ecosystems.”** It studies malignant hybrid states in single-nucleus GBM data, but not paired MRI-defined HBI and LBI samples.
3. **Neftel et al. (2019), “An integrative model of cellular states, plasticity and genetics for glioblastoma.”** It defines the canonical malignant programs and reports co-expression between states. It does not compare hybrid-state frequency between HBI and LBI within donor.

The hybrid-state concept is already published. Only the paired HBI/LBI question may be new. This remains a candidate, not a novelty claim. The supplements still need review.

**Bar 2 passes only if:** at least 10 donors have verified paired HBI/LBI samples; the prespecified within-donor HBI-minus-LBI difference in hybrid-cell fraction has a donor-bootstrap 95% interval excluding zero after the confound checks; and the direction stays the same with a second malignant-cell definition and in leave-one-donor-out analyses. These draft settings are in `output/trackA_v2_data_check/proposal_settings.json`. If any check fails, report the result as negative or inconclusive without changing it. Weakest point: prior studies already describe hybrid states and HBI/LBI cell-state differences; the new endpoint may still be too close to existing work.

## Downloads and next gate

The exact 60-file FASTQ list, URLs, byte sizes, and checksum file are in `output/trackA_v2_data_check/hra009021_file_inventory.json`. The full raw set is about 765 GB. The GSA page also lists `HRA009021.xlsx`, needed to identify which runs belong to which donor and region. The repository streams that workbook without reporting its size. No FASTQ or expression matrix has been downloaded.

Do not download the raw reads yet. First verify the workbook's donor/run/region map, choose only the paired HBI/LBI samples, then return with the exact subset list and total size for approval. Do not seal or score the protocol until the label map and usable donor count pass. Do not join Track A data to Track B or Track C.
