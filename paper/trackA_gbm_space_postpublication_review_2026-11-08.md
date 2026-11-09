# Track A GBM-Space post-publication metadata review

**Date:** 2026-11-08
**Scope:** public source and repository metadata only. No expression values, patient images, or new archives were downloaded or read. No score or figure was changed.

## Decision

GBM-Space should not be treated as an external test of the frozen Core-versus-Periphery Track A endpoint. The new Nature paper and official portal establish a substantial, multi-region, IDH-wildtype cohort and malignant-cell annotations. They do not, in the public sample-level metadata checked here, establish the exact nucleus-to-Core/Periphery labels and joins required by that endpoint. The source paper also directly studies spatial organization of malignant cell states, so a broad regional-biology rediscovery would not be novel.

The earlier note that the entire GBM-Space dataset fails because EGA requires an application is too broad. The raw sequencing layer is controlled, but processed single-cell and spatial data are described as downloadable through the GBM-Space portal and deposited in ArrayExpress. This changes the access classification, not the Core/Periphery eligibility decision.

## Verified metadata

- The 2026 Nature article reports **12 primary IDH-wildtype glioblastomas**, **57 tumour sampling sites**, more than **1,025,000 nuclei**, and paired single-nucleus multiome and Visium profiling across sites. The portal reports 28 malignant cell states and makes processed atlas data available for exploration and download.
- The paper distinguishes malignant cells using inferred chromosome 7 gain and chromosome 10 loss together with expression programs; it then assigns malignant-state annotations. This is a documented author-generated annotation, not an independently validated drug-response label.
- Metadata for ArrayExpress **E-MTAB-17183** was inspected through the BioStudies study API. Its sample-characteristics fields include donor (`Individual`), broad `Sampling site`, modality, and data-accession identifiers. The sample table has 24 donor-by-modality rows. The listed site values are broad lobes such as frontal, temporal, parietal, and occipital neocortex; the inspected table does not give a nucleus-level Core/Periphery field or a direct nucleus-to-anatomical-region crosswalk.
- This metadata check agrees with the earlier local audit: the available single-nucleus records did not establish six paired donors with exact, validated Core/Periphery labels. The Visium and histopathology results in the paper do not by themselves transfer binary regional labels to every single nucleus.
- Raw sequencing records are controlled through EGA. The current Nature data-availability statement names EGAD00001015526 (10x multiome), EGAD00001015527 (10x Visium), EGAD00001016155 (LCM-WGS), and EGAD00001016126 (Parse scRNA-seq). Processed data are available through the official portal and ArrayExpress. The checked ArrayExpress metadata did not expose an explicit record-level license field. EMBL-EBI says its Terms of Use do not add restrictions beyond those set by the data owner and may apply to records without an individual license. Therefore, access to processed files is not evidence that every derived output may be redistributed; retain attribution and do not redistribute source matrices or derived cell-level data without a clear basis.

## Prior-art and scientific fit

The article directly reports spatially patterned malignant-cell states, including developmental-like states near infiltrated/leading-edge tissue and gliosis/hypoxia states nearer the tumour core, across this same 12-tumour atlas. That preempts a broad claim that Track A newly discovers malignant-cell spatial organization in GBM. A fixed TMZ in-vitro-response-associated gene score could support a separate exploratory association with the published state axis, but this would not validate resistance: GBM-Space does not provide matched TMZ dose-response or treatment-response outcomes for these nuclei.

A targeted post-publication search found a closer treatment-response precedent: Suter et al. (2026) introduced scFOCAL, which integrates GBM single-cell states with drug perturbation signatures and reports in-vitro, ex-vivo, and in-vivo validation. The proposed 119/121 score is narrower and uses a different source (a published culture-viability association panel), but the broad idea of connecting GBM states to drug-response signatures is already established. Any state-score result here is therefore an exploratory transfer check, not a new drug-response method or a claim of superiority.

The state labels and the response-associated score both derive from expression data. A difference between groups could reflect shared transcriptomic structure or annotation construction and cannot establish functional TMZ sensitivity. The study's state annotations are useful for an association analysis, but they are not independent drug-response labels.

## Best next Track A route

Keep the 119/121 score and A5 result intact. A separate state-axis protocol is at `paper/trackA_gbm_space_tmz_state_protocol_2026-11-08.md`. The metadata-only preflight at `output/trackA_gbm_space_tmz_state/metadata_preflight.json` found 12 donors with at least 1,215 malignant nuclei in each primary group and 115 exact shared panel genes (47 positive-viability and 68 negative-viability directions). The local H5AD size and SHA-256 are frozen in `paper/trackA_gbm_space_tmz_state_input_fingerprint_2026-11-08.json`; this identifies the local bytes but does not authenticate them against an authoritative publisher checksum. Thirteen synthetic tests passed in a scratch copy, including backed-AnnData access, stored raw-count guards, explicit-zero handling, and exclusion of all measured source-panel symbols from control backgrounds. The code, protocol, source files, and input fingerprint must be sealed before scoring. The actual analysis has not run. Every raw-count check must pass during scoring. Treat donors, not nuclei, as the replication units. Report any result as an association between a transcriptomic proxy and author-annotated malignant states, not as drug resistance, treatment response, external confirmation of A5, new spatial biology, or a novel drug-response method.

If the single-cell state annotation cannot be linked reproducibly to donor and modality, stop. Do not substitute broad lobe labels for Core/Periphery or use the paper's Visium map as though each nucleus had a measured spatial coordinate.

## Sources

- de Jong, G. et al. (2026). *A spatiotemporal cancer cell trajectory underlies glioblastoma heterogeneity*. Nature. https://doi.org/10.1038/s41586-026-11234-8
- Suter, R. K. et al. (2026). *Drug and single-cell gene expression integration identifies sensitive and resistant glioblastoma cell populations*. Nature Communications, 17, 99. https://doi.org/10.1038/s41467-025-67783-5
- Official GBM-Space portal: https://www.gbmspace.org/
- BioStudies/ArrayExpress records: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-17183 and https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-17164
- Current EGA records named by the article: https://ega-archive.org/datasets/EGAD00001015526, https://ega-archive.org/datasets/EGAD00001015527, https://ega-archive.org/datasets/EGAD00001016155, and https://ega-archive.org/datasets/EGAD00001016126.
- BioStudies metadata API responses inspected without retrieving assay files: https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-17183 and https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-17164.
- EMBL-EBI licensing policy: https://www.ebi.ac.uk/licencing/
- Prior local metadata audit: `paper/trackA_v2_data_check.md`, `output/trackA_v2_data_check/region_provenance.json`, and `output/trackA_v2_data_check/s14_metadata_check.json`.
