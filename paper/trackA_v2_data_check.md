# Track A v2: GBM-Space data check

Metadata inspection only; no patient expression values were read or scored. Evidence: `output/trackA_v2_data_check/metadata_check.json`. The downloaded archive and expanded file remain in Downloads.

- The archive matches the repository's stated byte length and passed gzip decompression/CRC. A local SHA-256 is recorded in the output. Weakest point: no publisher checksum was located for independent comparison.
- The cell metadata contains 12 nonblank donor IDs and 57 donor/site combinations. Weakest point: donor IDs are provided by the dataset; eligibility for paired core and edge is unknown. See the output JSON.
- No explicit core or edge column appears in the cell metadata. `site_id` identifies sampling sites, not a verified anatomical label. No independent mapping has been established; the minimum paired-donor gate is therefore unproved. Do not infer labels from gene expression.
- The H5AD has an `X` matrix and `raw/X`, but no layers. Their presence does not establish that either contains unnormalized integer counts. Count suitability for scVI remains unverified; do not inspect expression or score until the applicable protocol and seal are committed and pushed.
- [BioStudies terms](https://www.ebi.ac.uk/biostudies/submit) say new datasets are released under CC0. [BioStudies help](https://www.ebi.ac.uk/biostudies/RH3R/help) says that where no license is stated, EBI terms permit reuse and redistribution. The exact deposited asset's license field and any donor-specific restrictions have not been independently confirmed.
- The prior-art search remains incomplete: a broad Europe PMC query is recorded in `output/trackA_v2_data_check/europepmc_search.json`, but only its first results were reviewed. Google Scholar and full-text review are still needed. No biology novelty or model advantage is established.

Next: source a non-expression site-to-region mapping and verify at least ten eligible paired donors without downloading new data; then resolve raw-count provenance, terms, and prior art before writing and sealing the protocol. Additional files require the user's explicit approval.

## Site-to-region documentation check

The [GBM-Space preprint](https://doi.org/10.1101/2025.05.13.653495) reports pathologist-drawn IvyGAP regions on H&E images, independently validated by other pathologists. These labels are on Visium tissue sections; the nuclei were collected from pooled thicker serial sections of each block. Distinct anatomical regions can coexist within a block, so a matching site ID does not establish a core or edge label for an individual nucleus. A label also cannot be derived from gene-expression niches without making the target circular.

The preprint reports 39 histology-annotated Visium sections; see `output/trackA_v2_data_check/region_provenance.json`. Weakest point: we did not find a usable site-to-region crosswalk, review supplement S14, or establish the number of eligible paired donors. This is an unresolved gate, not proof that no suitable mapping exists. The spatial annotations may support a different design, but would need an appropriate independent-label and donor audit. No additional patient data were downloaded or scored.

## Additional screening

The [preprint methods](https://www.researchgate.net/publication/391763862_A_spatiotemporal_cancer_cell_trajectory_underlies_glioblastoma_heterogeneity) say the authors used scVI to integrate their nuclei, with reaction and donor-related covariates. That published embedding and any donor-wide feature selection cannot serve as leakage-free inputs to our held-out test. Retraining and all feature selection must occur inside the training folds. The methods do not identify which H5AD slot holds counts suitable for a new scVI fit.

The [EGA regional dataset](https://ega-archive.org/datasets/EGAD50000000778) has four patients, not eleven; the latter is the number of tissues. The [Khan spatial study](https://www.biorxiv.org/content/10.1101/2025.04.04.647096v3.full) reports nine patients in its spatial cohort; its separate single-cell reference does not establish a paired eligible cohort. Both are below the ten-donor gate as described. Evidence and limitations: `output/trackA_v2_data_check/candidate_method_checks.json`.

Preliminary PubMed, Europe PMC, and Google Scholar query results are in `output/trackA_v2_data_check/`. These are limited searches; full-text and supplement checks remain incomplete. Related adversarial single-cell models have already been published, including [scCobra](https://pmc.ncbi.nlm.nih.gov/articles/PMC11825689/); do not describe the general adversarial idea as novel.

The [preprint supplement page](https://www.biorxiv.org/content/10.1101/2025.05.13.653495v1.supplementary-material) links `media-3.xlsx`, an accompanying tables workbook that may contain the sectioning plan, Table S14. A headers-only request reports 136,557 bytes (see `output/trackA_v2_data_check/supplement_manifest.json`). Approval is required before downloading. Examine sample provenance only; do not inspect expression results before the protocol is sealed. A sectioning plan alone may not prove a nucleus-level anatomical label.
