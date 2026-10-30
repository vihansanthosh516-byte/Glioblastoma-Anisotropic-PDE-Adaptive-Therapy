# Track A v2: GBM-Space data check

Metadata inspection only; no patient expression values were read or scored. Evidence: `output/trackA_v2_data_check/metadata_check.json`. The downloaded archive and expanded file remain in Downloads.

- The archive matches the repository's stated byte length and passed gzip decompression/CRC. A local SHA-256 is recorded in the output. Weakest point: no publisher checksum was located for independent comparison.
- The cell metadata contains 12 nonblank donor IDs and 57 donor/site combinations. Weakest point: donor IDs are provided by the dataset; eligibility for paired core and edge is unknown. See the output JSON.
- No explicit core or edge column appears in the cell metadata. `site_id` identifies sampling sites, not a verified anatomical label. No independent mapping has been established; the minimum paired-donor gate is therefore unproved. Do not infer labels from gene expression.
- The H5AD has an `X` matrix and `raw/X`, but no layers. Their presence does not establish that either contains unnormalized integer counts. Count suitability for scVI remains unverified; do not inspect expression or score until the applicable protocol and seal are committed and pushed.
- [BioStudies terms](https://www.ebi.ac.uk/biostudies/submit) say new datasets are released under CC0. [BioStudies help](https://www.ebi.ac.uk/biostudies/RH3R/help) says that where no license is stated, EBI terms permit reuse and redistribution. The exact deposited asset's license field and any donor-specific restrictions have not been independently confirmed.
- The prior-art search remains incomplete: Google Scholar and Europe PMC checks and full-text review are still needed. No biology novelty or model advantage is established.

Next: source a non-expression site-to-region mapping and verify at least ten eligible paired donors without downloading new data; then resolve raw-count provenance, terms, and prior art before writing and sealing the protocol. Additional files require the user's explicit approval.
