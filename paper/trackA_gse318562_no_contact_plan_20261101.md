# Track A no-contact cohort plan: GSE318562

**Status:** Candidate plan only, not adopted (2026-11-01). No person was contacted. No dataset file was downloaded, and no expression values were read.

## Recommendation

Use GSE318562 as the first public candidate for an optional, descriptive spatial replication of a signature already frozen from the existing Track A training data. It is not a replacement source for a malignant-cell-only external test, and it cannot support a new-biology or novelty claim.

GEO describes 22 Visium V1 spatial samples from nine primary IDH-wildtype GBM patients. The corresponding preprint reports core-to-edge sampling and manual neuropathology annotation of tissue regions. Public files are listed at [GSE318562](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318562), with additional processed files at [Zenodo record 20709378](https://zenodo.org/records/20709378). GEO lists the processed `GSE318562_RAW.tar` at about 674 MB. The full Zenodo record is about 51.9 GB, so a future run should use only the smallest sufficient GEO files if the user chooses to proceed.

The study is a close topic match: its preprint already analyzes cell states and multicellular niches along the GBM core-to-edge axis. Any use here must be described as a fixed-signature replication on a separate cohort, with no claim that the core-edge question or its biology is new. See the [preprint](https://www.biorxiv.org/content/10.1101/2025.04.04.647096v3.full).

## Eligibility limits

- **Donor and region map:** GEO identifies nine patients and 22 samples, but the sample-title list alone does not link every section to patient and region. Before expression access, verify the public supplement's sample-to-patient and sample-to-histology map, then count donors with usable Core and Edge observations. Stop if the map is incomplete; do not contact the authors.
- **Cell specificity:** Visium measures mixed spots. Histology and spot-level CNV estimates can support a tumor-enriched spatial analysis, but they do not establish a malignant-cell-only expression test. Report this as spatial, not malignant-cell, evidence.
- **Data terms:** GEO says it places no restrictions on GEO data use, while warning that submitters may assert rights. The Zenodo record's license field appears blank in this screen. Recheck the exact files and any stated rights before download. If permission or terms remain unclear, stop without contacting anyone.
- **Prior exposure and overlap:** Compare the mapped patient/sample IDs against the existing Track A training and prior-use records before analysis. Do not count any overlapping or previously scored samples as independent replication.
- **Prior art:** The source study itself studies core-to-edge organization. A result here can only test whether a pre-frozen Track A gene set reproduces spatially in these samples; it cannot establish novelty.

## If the user later chooses to proceed

1. Complete the metadata-only checks above. Keep the assay, donor, section, region, tumor-content, and prior-use map in a small audit record. Do not read expression values during this step.
2. Confirm that the frozen signature's genes and the planned comparator sets are represented in the public assay files. Check the exact GEO/Zenodo terms and attribution requirements. No author contact.
3. Before downloading or reading expression, write and commit a dated Track A amendment that freezes the spatial endpoint, the patient-level analysis, controls, code, and pass rule. Preserve the existing Yu 2020 result and its separate protocol.
4. If the map and terms pass, score the frozen signature per donor in histology-labeled Core and Edge spots. Treat the donor—not the spot—as the analysis unit. Compare it with normal-brain and cell-composition controls, a prior Core/Periphery signature, and size-matched random gene sets. Account for estimated tumor content and report all donor values. Use no p-values or cell-level pseudo-replication.
5. Keep the claim descriptive. Because there are nine patients and mixed spots, report uncertainty and limitations; do not call it malignant-cell transfer or a new discovery.

## Current decision

GSE318562 is worth keeping as an optional spatial candidate, pending the public donor-to-region map and terms check. It does **not** pass the existing malignant-cell external-validation gate. The closer single-nucleus option already screened in `paper/trackA_v2_proposal.md`, GSE286419, has only four region-concordant Core/Margin donors and is also too small for the intended primary confirmation. If the project requires malignant-cell-level confirmation, the current decision stays **no-go without another eligible public cohort**.

## Sources checked

- [NCBI GEO GSE318562](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318562): public status, study design, sample records, and processed-file listing.
- [Khan et al. preprint](https://www.biorxiv.org/content/10.1101/2025.04.04.647096v3.full): core-edge study design and histology-based annotation; preprint, not treated as peer-reviewed evidence here.
- [Zenodo record 20709378](https://zenodo.org/records/20709378): public processed spatial and inferred-CNV file listing; no license displayed in the record view checked.
- [NCBI GEO disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html): NCBI places no restrictions on GEO data use but notes that submitters may claim rights.
