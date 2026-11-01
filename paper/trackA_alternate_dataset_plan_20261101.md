# Track A no-contact alternate dataset plan

**Status:** Plan only (2026-11-01). No one was contacted. No expression data were opened or downloaded. No figure was changed. This plan does not authorize data access or scoring.

## Aim

The primary Track A question remains whether the already frozen malignant-cell Core/Periphery signature transfers to independent patient samples. The current best no-contact lead, E-MTAB-8743, can support only a separate bulk-tissue regional replication; it cannot answer the malignant-cell-specific question. Treat any future use as replication only. Do not claim a new core-edge discovery or add outside data to the current figure.

## Candidate screen

| Candidate | Publicly reported structure | Access screen | Decision |
|---|---|---|---|
| **E-MTAB-8743 (2023)** | Ten newly diagnosed, treatment-naive IDH-wildtype GBM patients; the paper reports bulk RNA-seq from spatially distinct Core, Rim, and Invasive-margin tissue, plus FACS-sorted 5ALA+/- cells from the margin. | The paper deposits raw RNA-seq at the public BioStudies/ArrayExpress record. BioStudies documents public access without login and anonymous FTP. Its help says records with no stated license use EMBL-EBI Terms of Use; the exact record-level license and file manifest have not been surfaced here. The paper reports at least 70 million paired 75-bp reads per sample; exact download size is unverified. | **Best current no-contact option for a separate regional RNA replication.** The ten-patient cohort meets the donor floor on paper, but bulk tissue cannot establish malignant-cell-specific transfer. Verify complete Core/Invasive-margin sample pairing, file manifest, prior-use overlap, and record terms before any data download. |
| **SNUH-GliomaAtlas (2026)** | The paper reports six matched fresh Core/Periphery single-cell donors and 223,113 quality-filtered cells. Metadata from the linked preprint supplement identifies five IDH-wildtype GBM pairs (SNU27, SNU32, SNU35, SNU46, SNU51); the sixth paired donor, SNU40, is IDH-mutant grade 4 astrocytoma. The paper reports CNV checks for neoplastic-cell classification. | The linked GitHub README says raw study data are available upon request, while the atlas portal lists some downloadable objects. Exact sample-level files and terms for the five eligible GBM pairs are unverified. | **Not a primary candidate at present.** The strict IDH-wildtype GBM count is five in the preprint metadata, below the six-donor floor; final-paper sample mapping remains unverified. No contact or request. |
| **GSE182109 (2022)** | The paper reports 18 glioma patients, including 11 newly diagnosed GBM, 5 recurrent GBM, and 2 lower-grade glioma; ten patients had multiregion sampling. The paper's regional classes include enhancing, necrotic, infiltrative, and margin/non-enhancing tissue. | GEO lists 44 processed 10x samples; the paper states raw single-cell data are public without restrictions. The GEO sample rows do not expose the needed donor-to-region crosswalk. | **Unverified (U).** Potential cell-level lead, but six eligible IDH-wildtype GBM Core/Periphery donor pairs are not established from public metadata. No download until the map is verified. |
| **HRA005318** | Open 10x single-cell dataset. Public sample records show three exact tumor/peritumor pairs: G0153, G0162, and G0171; other records are unpaired. | GSA-Human lists open access. The public record links raw sequencing files; a processed expression matrix and reuse license were not verified. | Only three exact pairs. Too small for the planned primary test. Possible descriptive fallback only if the user later accepts a smaller cohort. |
| **GSE84465 (Darmanis 2017)** | 3,589 cells from four IDH-wildtype GBM patients, with Core and Periphery samples in each patient. GEO lists a 19.6 MB processed file. | Direct public GEO download. GEO's general data-use notice still needs to be cited; no separate file license was found in this screen. | Four donors, below the project's six-donor minimum. Exploratory fallback only; no confirmatory or novelty claim. |
| **GSE318562** | 22 Visium sections from nine patients, with five exact sample-level Core/Edge pairs. Spots mix cell types. | Public GEO files; spot-level histology join and exact reuse terms remain unverified. | Not suitable for malignant-cell validation. Keep only as an optional spatial description. See `trackA_gse318562_no_contact_plan_20261101.md`. |
| **GSE286419** | Four region-concordant Core/Margin donors in the existing audit. | Public GEO record. | Below the six-donor minimum; do not use as the primary confirmation. |

The E-MTAB-8743 paper already analyzes expression differences across Core, Rim, and Invasive-margin regions, and the SNUH study analyzes Core-versus-Periphery cell states and composition. Any reuse can test a frozen signature's reproducibility only; neither supports a novelty claim.

## No-contact gates for any candidate

1. Read public metadata only. Verify donor/sample/region links, diagnosis and IDH status, complete paired regions, and the exact eligible donor count. Do not infer pair counts from total atlas size or a study-level statement alone.
2. Verify the exact processed/raw files and the record-level terms. If an account approval, access application, or person-to-person request is required, stop. If use terms remain unclear, stop without contact.
3. Compare patient and sample identifiers with Track A's training and prior-use records. Stop if there is overlap or prior scoring.
4. If those gates pass, prepare and commit a dated amendment before opening expression values. Freeze the endpoint, controls, code, and interpretation limits. Obtain the user's go-ahead before any expression-data download or scoring.
5. Keep the evidence type explicit. E-MTAB-8743 can support only a donor-level bulk-tissue Core-versus-Invasive-margin replication; it cannot validate a malignant-cell-only result. A cell-level validation requires a separate cohort with a defensible malignant-cell label. Treat each patient as one analysis unit and report every eligible paired value.
6. Keep any future result descriptive. Do not add outside data to the current figure; any display would need a separate, clearly sourced figure.

## Current decision

**Download nothing.** E-MTAB-8743 is the strongest currently screened no-contact candidate for a separate regional bulk-RNA replication: the paper reports ten IDH-wildtype GBM patients and Core/Rim/Invasive-margin RNA-seq, but the exact accession manifest, complete-pair map, prior-use overlap, and record-specific terms still need metadata-only verification. It cannot validate a malignant-cell-only signature. SNUH-GliomaAtlas has six matched pairs overall but only five strict IDH-wildtype GBM pairs in the preprint supplement; its raw data are described as available upon request, so it is not a current no-contact primary candidate. HRA005318 and GSE84465 have only three and four exact pairs. No screened alternative currently passes every gate for primary malignant-cell confirmation. No one was contacted, no request was made, no expression was opened, no score was run, and no figure was changed.

## Search coverage and limits

This targeted metadata screen checked NCBI GEO, GSA-Human, EGA, Dryad, BioStudies/ArrayExpress (including E-MTAB-8743), the SNUH-GliomaAtlas site and linked GitHub, plus linked source papers. The follow-up also screened GSE182109, GSE317289, and an eight-patient Dryad imaging-mass-cytometry cohort. No expression or sequence files were downloaded or read. It did not exhaustively search HCA, HTAN, CPTAC, OSF, or every Figshare and Zenodo record. Do not describe the search as exhaustive.

## Sources

- [SNUH-GliomaAtlas paper](https://doi.org/10.1038/s41467-026-69716-2): six matched Core/Periphery patients, cell counts, CNV-supported neoplastic assignment, and study findings.
- [SNUH-GliomaAtlas data/code page](https://github.com/nameetas/TSKGA): public atlas links and statement that raw study data are available upon request.
- [SNUH preprint supplement, patient/sample table](https://assets-eu.researchsquare.com/files/rs-4468724/v1/5044141e7974eafd0cf3cce4.pdf): older preprint-era diagnosis and Core/Periphery sample metadata used for the provisional five-IDH-wildtype-pair count; this does not replace the final paper's supplement.
- [E-MTAB-8743 paper](https://doi.org/10.1186/s13073-023-01207-1): ten IDH-wildtype GBM patients and bulk RNA-seq from Core, Rim, Invasive margin, and FACS-sorted margin cells.
- [E-MTAB-8743 BioStudies record](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8743): public accession for the deposited raw RNA-seq files.
- [BioStudies access and terms](https://www.ebi.ac.uk/biostudies/RH3R/help): public records are accessible without login; anonymous FTP is supported; records with no stated license use EMBL-EBI Terms of Use.
- [GSE182109 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182109) and [source paper](https://doi.org/10.1038/s41467-022-28372-y): 44 single-cell samples from a mixed-grade cohort; ten-patient multiregion design; sample-to-region crosswalk not confirmed in GEO rows.
- [GSA-Human HRA005318](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA005318): open-access single-cell record and download listing.
- [GSA-Human BioProject PRJCA019090](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA019090): public sample names and Core/Periphery-related metadata.
- [NCBI GEO GSE84465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465): four-patient single-cell cohort and processed file listing.
- [NCBI GEO GSE286419](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286419): six-tumor study with only four region-concordant donors in the existing audit.
- [NCBI GEO GSE318562](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318562): Visium cohort screened in the companion no-contact plan.
