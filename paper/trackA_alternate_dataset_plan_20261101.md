# Track A no-contact alternate dataset plan

**Status:** Plan only (2026-11-01). No one was contacted. No expression data were opened or downloaded. No figure was changed. This plan does not authorize data access or scoring.

## Aim

Test whether the already frozen Track A malignant-cell Core/Periphery signature reproduces in a separate patient cohort. Treat this as replication only. Do not claim a new core-edge discovery or add outside data to the current figure.

## Candidate screen

| Candidate | Publicly reported structure | Access screen | Decision |
|---|---|---|---|
| **SNUH-GliomaAtlas (2026)** | Six biologically independent patients have matched fresh Core and Periphery single-cell samples; 223,113 cells passed the study's quality filter. Samples were processed on the same day. The study checked its neoplastic-cell assignment with CNV analysis. | The paper describes a public atlas, but its linked GitHub README says raw study data are available upon request. The exact public download for the six paired expression samples and file-level reuse terms are not verified. | **Best scientific fit, access unverified (U).** Do not download or contact anyone. Use only if a no-contact public file and its terms can be verified. |
| **HRA005318** | Open 10x single-cell dataset. Public sample records show three exact tumor/peritumor pairs: G0153, G0162, and G0171; other records are unpaired. | GSA-Human lists open access. The public record links raw sequencing files; a processed expression matrix and reuse license were not verified. | Only three exact pairs. Too small for the planned primary test. Possible descriptive fallback only if the user later accepts a smaller cohort. |
| **GSE84465 (Darmanis 2017)** | 3,589 cells from four IDH-wildtype GBM patients, with Core and Periphery samples in each patient. GEO lists a 19.6 MB processed file. | Direct public GEO download. GEO's general data-use notice still needs to be cited; no separate file license was found in this screen. | Four donors, below the project's six-donor minimum. Exploratory fallback only; no confirmatory or novelty claim. |
| **GSE318562** | 22 Visium sections from nine patients, with five exact sample-level Core/Edge pairs. Spots mix cell types. | Public GEO files; spot-level histology join and exact reuse terms remain unverified. | Not suitable for malignant-cell validation. Keep only as an optional spatial description. See `trackA_gse318562_no_contact_plan_20261101.md`. |
| **GSE286419** | Four region-concordant Core/Margin donors in the existing audit. | Public GEO record. | Below the six-donor minimum; do not use as the primary confirmation. |

The SNUH study itself analyzes Core-versus-Periphery cell states and composition. Any reuse can test a frozen signature's reproducibility only; it cannot support a novelty claim.

## No-contact gates for SNUH-GliomaAtlas

1. Read public metadata only. Confirm the six patient-to-sample pairs, diagnoses, IDH status, region labels, and which records belong to the newly generated cohort. Do not treat the full 28-patient atlas or its pooled reference as six paired samples.
2. Confirm that the expression object for those six pairs is downloadable without an account approval, application, or person-to-person request. Check the exact file terms. If a request or unclear terms are required, stop and do not contact anyone.
3. Compare patient and sample identifiers with Track A's existing training and prior-use records. Stop if there is overlap or prior scoring.
4. If these gates pass, prepare and commit a dated amendment before opening expression values. Freeze the signature, malignant-cell rule, patient-level endpoint, controls, code, and interpretation limits. Obtain the user's go-ahead before any expression-data download or scoring.
5. Score only CNV-supported neoplastic cells. Aggregate by patient and region; each patient is one analysis unit. Show all six paired patient values. Compare with the normal-brain contamination score, Neftel cell-state scores, a prior Core/Periphery signature, and size-matched random gene sets. Do not treat cells as independent replicates.
6. Limit any result to descriptive replication. Six pairs do not support a strong superiority or novelty claim. Keep external data out of the current figure unless the user later requests a separate, clearly sourced figure.

## Current decision

**Download nothing.** SNUH-GliomaAtlas is the strongest alternative by scientific fit, but it does not pass the no-contact access and terms gates yet. The open HRA005318 and GSE84465 datasets have only three and four exact pairs, respectively. No screened alternative currently passes every gate for primary malignant-cell confirmation. No author contact, data request, expression read, score, or figure edit is part of this plan.

## Search coverage and limits

This targeted metadata screen checked NCBI GEO, GSA-Human, EGA, Dryad, the SNUH-GliomaAtlas site and linked GitHub, plus linked source papers. It did not exhaustively search HCA, BioStudies/ArrayExpress, HTAN, CPTAC, OSF, or every Figshare and Zenodo record. Do not describe the search as exhaustive.

## Sources

- [SNUH-GliomaAtlas paper](https://doi.org/10.1038/s41467-026-69716-2): six matched Core/Periphery patients, cell counts, CNV-supported neoplastic assignment, and study findings.
- [SNUH-GliomaAtlas data/code page](https://github.com/nameetas/TSKGA): public atlas links and statement that raw study data are available upon request.
- [GSA-Human HRA005318](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA005318): open-access single-cell record and download listing.
- [GSA-Human BioProject PRJCA019090](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA019090): public sample names and Core/Periphery-related metadata.
- [NCBI GEO GSE84465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465): four-patient single-cell cohort and processed file listing.
- [NCBI GEO GSE286419](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286419): six-tumor study with only four region-concordant donors in the existing audit.
- [NCBI GEO GSE318562](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318562): Visium cohort screened in the companion no-contact plan.
