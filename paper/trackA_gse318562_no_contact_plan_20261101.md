# Track A no-contact cohort plan: GSE318562

**Status:** Metadata-only candidate screen done (2026-11-01); not selected for malignant-cell confirmation. Public NCBI BioSample records map all 22 Visium sections to nine donor IDs and sample-region labels. No person was contacted, no dataset file was downloaded, and no expression values were read. The separate alternative-dataset plan is in `trackA_alternate_dataset_plan_20261101.md`.

## Recommendation

Do not select GSE318562 for the primary malignant-cell external test. It can remain an optional, descriptive spatial replication candidate. The public sample labels give only **five unambiguous donor-level Core/Edge pairs**. It cannot support a new-biology or novelty claim.

GEO describes 22 Visium V1 spatial samples from nine primary IDH-wildtype GBM patients. The corresponding preprint reports core-to-edge sampling and manual neuropathology annotation of tissue regions. Public files are listed at [GSE318562](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318562), with additional processed files at [Zenodo record 20709378](https://zenodo.org/records/20709378). GEO lists the processed `GSE318562_RAW.tar` at about 674 MB. The full Zenodo record is about 51.9 GB, so a future run should use only the smallest sufficient GEO files if the user chooses to proceed.

### NCBI sample-to-donor map

NCBI BioSample `isolate` and `sample type` attributes give this map. Exact Core/Edge pairs occur for B186, B189, WU1220, WU1221, and WU1227. B185 has a Core sample plus `Transition/Core` and `Transition/Edge` samples; it is not counted as an exact pair. B172, B176, and B178 have transition-labeled samples only. The 22 sample labels comprise 7 Core, 6 Edge, and 9 transition or mixed-region samples.

| NCBI BioSample | Sample | Donor | Sample region |
|---|---|---|---|
| [SAMN54992638](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992638) | B172-01 | B172 | Transition |
| [SAMN54992639](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992639) | B172-02 | B172 | Transition |
| [SAMN54992640](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992640) | B176-01-02 | B176 | Transition, mostly Core |
| [SAMN54992641](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992641) | B176-02-01 | B176 | Transition |
| [SAMN54992642](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992642) | B176-02-04 | B176 | Transition |
| [SAMN54992643](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992643) | B178-02-01 | B178 | Transition, mostly Core |
| [SAMN54992644](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992644) | B185-01 | B185 | Core |
| [SAMN54992645](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992645) | B185-02 | B185 | Transition/Core |
| [SAMN54992646](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992646) | B185-03 | B185 | Transition/Edge |
| [SAMN54992647](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992647) | B186-01 | B186 | Edge |
| [SAMN54992648](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992648) | B186-02 | B186 | Edge |
| [SAMN54992649](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992649) | B186-03 | B186 | Core |
| [SAMN54992650](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992650) | B186-04 | B186 | Core |
| [SAMN54992651](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992651) | ST073021-02 | B189 | Edge |
| [SAMN54992652](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992652) | ST073021-04 | B189 | Core |
| [SAMN54992653](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992653) | WU1220-01 | WU1220 | Transition |
| [SAMN54992654](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992654) | WU1220-02 | WU1220 | Edge |
| [SAMN54992655](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992655) | WU1220-03 | WU1220 | Core |
| [SAMN54992656](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992656) | WU1221-Core2 | WU1221 | Core |
| [SAMN54992657](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992657) | WU1221-Edge1 | WU1221 | Edge |
| [SAMN54992658](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992658) | WU1227-Core5 | WU1227 | Core |
| [SAMN54992659](https://www.ncbi.nlm.nih.gov/biosample/SAMN54992659) | WU1227-Edge2 | WU1227 | Edge |

The study is a close topic match: its preprint already analyzes cell states and multicellular niches along the GBM core-to-edge axis. Any use here must be described as a fixed-signature replication on a separate cohort, with no claim that the core-edge question or its biology is new. See the [preprint](https://www.biorxiv.org/content/10.1101/2025.04.04.647096v3.full).

## Eligibility limits

- **Donor and region map:** NCBI BioSample metadata link all 22 Visium sample records to donor and sample-region labels. Only five donors have exact sample-level Core and Edge pairs. The preprint says H&E sections were annotated for high- and low-tumor-content regions, but availability of a public per-spot histology map has not been verified. Do not count mixed-region labels as exact Core/Edge pairs. If the public spot map is absent, keep the cohort at five paired donors.
- **Cell specificity:** Visium measures mixed spots. Histology and spot-level CNV estimates can support a tumor-enriched spatial analysis, but they do not establish a malignant-cell-only expression test. Report this as spatial, not malignant-cell, evidence.
- **Data terms:** GEO says it places no restrictions on GEO data use, while warning that submitters may assert rights. The Zenodo record's license field appears blank in this screen. Recheck the exact files and any stated rights before download. If permission or terms remain unclear, stop without contacting anyone.
- **Prior exposure and overlap:** Compare the mapped patient/sample IDs against the existing Track A training and prior-use records before analysis. Do not count any overlapping or previously scored samples as independent replication.
- **Prior art:** The source study itself studies core-to-edge organization. A result here can only test whether a pre-frozen Track A gene set reproduces spatially in these samples; it cannot establish novelty.

## If the user later chooses to proceed

1. Confirm whether public files include the per-spot histology map and whether it can be joined to the listed Visium samples without expression values. Keep the assay, donor, section, region, tumor-content, and prior-use map in a small audit record. Stop if that map is not public; do not contact anyone.
2. Confirm that the frozen signature's genes and the planned comparator sets are represented in the public assay files. Check the exact GEO/Zenodo terms and attribution requirements. No author contact.
3. Before downloading or reading expression, write and commit a dated Track A amendment that freezes the spatial endpoint, the patient-level analysis, controls, code, and pass rule. Preserve the existing Yu 2020 result and its separate protocol.
4. If the map and terms pass, score the frozen signature per donor in histology-labeled Core and Edge spots. Treat the donor—not the spot—as the analysis unit. Compare it with normal-brain and cell-composition controls, a prior Core/Periphery signature, and size-matched random gene sets. Account for estimated tumor content and report all donor values. Use no p-values or cell-level pseudo-replication.
5. Keep the claim descriptive. Because there are nine patients and mixed spots, report uncertainty and limitations; do not call it malignant-cell transfer or a new discovery.

## Current decision

GSE318562 is not selected for the primary test. The metadata map is complete at the sample level, but the strict paired count is five, and public per-spot histology labels remain unverified. It does **not** pass the malignant-cell external-validation gate. GSE286419 has only four region-concordant Core/Margin donors. E-MTAB-8743 is already used in the Track A TAM exploration and cannot serve as independent confirmation. GSE317289 is a new public spatial lead with five IDH-wildtype GBM patients, below the six-donor floor; spatial spots also cannot validate a malignant-cell-only result. See `trackA_alternate_dataset_plan_20261101.md`; the primary malignant-cell confirmation remains **no-go**. No outside dataset was added to the current figure.

## Sources checked

- [NCBI GEO GSE318562](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318562): public status, study design, sample records, and processed-file listing.
- [NCBI BioSample PRJNA1372229 records](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1372229): public donor (`isolate`) and sample-region (`sample type`) attributes used for the 22-row map above.
- [Khan et al. preprint](https://www.biorxiv.org/content/10.1101/2025.04.04.647096v3.full): core-edge study design and histology-based annotation; preprint, not treated as peer-reviewed evidence here.
- [Zenodo record 20709378](https://zenodo.org/records/20709378): public processed spatial and inferred-CNV file listing; no license displayed in the record view checked.
- [NCBI GEO disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html): NCBI places no restrictions on GEO data use but notes that submitters may claim rights.
