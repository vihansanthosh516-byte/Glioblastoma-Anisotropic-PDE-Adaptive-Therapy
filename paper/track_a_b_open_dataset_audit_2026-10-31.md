# Track A/B open-dataset decision

**Date:** 2026-10-31  
**Scope:** Tracks A and B only. Track C is deferred at the user's direction.  
**Rule:** No author contact, access request, download, or patient scoring in this audit.

## Decision

- **Track A:** GBM-Space fails the current six-donor descriptive feasibility gate for the strict same-section `Cellular tumor` plus `Leading edge` feature screen. Best no-request fallback is a possible seven-donor descriptive analysis combining GEO GSE273274 (3 paired core/periphery donors) and GSE162631 (4 paired core/periphery donors). This needs a protocol amendment and metadata feasibility check. One cohort contains broad single-cell samples; the other enriches CD31-positive endothelial cells. Seven donors cannot support the plan's inferential superiority claim, which requires ten.
- **Track B:** No newly identified public GBM MRI cohort clears the current two-history-to-next-mask PDE gates. Do not download Duke or another replacement. CFB-GBM is already local; registration review and target-mask support/provenance gates remain open. Do not score patient masks.

## Track A

### GBM-Space complete metadata inventory

The official E-MTAB-17164 SDRF lists 97 Visium H5AD files from 12 donors. The completed HDF5 metadata-only inventory is [tracka_gbmspace_h5ad_feature_inventory_20261031_v2.json](../output/tracka_gbmspace_h5ad_feature_inventory_20261031_v2.json). It found:

- 97/97 files complete; zero errors.
- 38 sections with at least one histopathology feature name.
- Six sections with both `Cellular tumor` and `Leading edge` feature names, across four donors: AT10, AT15, AT3, and AT4.
- Eighteen sections with both `Cellular tumor` and `Infiltrating tumor` feature names, across five donors. Infiltrating tumor is not the same label as leading edge.
- No expression values, annotation values, or images were read. Feature names do not prove nonzero or usable labels. The strict screen stays below the six-donor feasibility gate.

The [GBM-Space portal](https://www.gbmspace.org/) lists a 3.9 GB processed download and spatial histopathology annotations. The public [E-MTAB-17164 record](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-17164) has no explicit accession-level license field. BioStudies directs users to EMBL-EBI Terms of Use when a record does not state a license; owner terms and attribution still apply. This cohort has already been screened in this project, so it is not untouched external validation.

### Best direct-public paired core/periphery fallback

| Dataset | Source-backed structure | Exact file | Fit and limit |
|---|---|---|---|
| [GSE273274](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273274) | Three patients; three paired tumor-core and periphery samples | `GSE273274_RAW.tar`, 706.1 MB | Broad single-cell samples; direct public GEO access. |
| [GSE162631](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162631) | Four patients; four paired core and periphery samples | `GSE162631_raw_counts_matrix.zip.gz`, 806.1 MB | CD31-positive endothelial-cell enriched; not a broad all-cell cohort. |

The combined maximum is seven paired donors, subject to checking donor overlap, sample mapping, and usable cell types. NCBI GEO says anyone can download public GEO data without an account and describes GEO as unrestricted access; cite the source papers and accessions. The [GEO FAQ](https://www.ncbi.nlm.nih.gov/geo/info/faq.html) supports no case-by-case access request. Because the assays differ, do not pool all cells without a frozen cell-type strategy. This option is descriptive only and cannot establish superiority under the ten-donor rule.

Other candidates do not clear the current design: [GSE318562](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318562) has five represented core/edge groups in the audited spot file and lacks spot-level histology labels there; [GSE237183](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237183) reports 13 GBM spatial samples, not 13 established independent donors, and uses surgeon-region annotations; its accession is already mentioned in the project archive. HRA009021 remains unusable for the current question until its sample-to-region map is source-verified; no access request is needed for its public files.

## Track B

- **CFB-GBM v4:** already downloaded locally. Existing target-blind MRI review aid covers 38 candidate patients and 760 pages, but no independent reviewer has assessed it. MRI rendering is not a registration pass. Valid-negative support and mask provenance remain unresolved. No target-mask values were read and no scoring is authorized.
- **MU-Glioma-Post:** public, but previously used for development; selected-case registration and full-grid target-support evidence remain unresolved. It is not fresh external validation.
- **RHUH-GBM:** TCIA lists 40 patients and public 2.9 GB brain-extracted images/segmentations under CC BY 4.0. Its preoperative, early postoperative, and recurrence stages include surgery between the first two scans; only one later post-surgery interval remains. It does not provide two clean post-surgery history masks before a future target for the current growth-only PDE.
- **UCSD-PTGBM:** public, but its aggregate 178-patient/243-timepoint count does not establish enough eligible three-visit sequences or the required registration/support evidence. The annotation-transition use is already exposed and remains separate from the PDE core.
- **Duke:** requires individual access approval; excluded under the user's no-request rule. **Burdenko:** image access requires a restricted-license request. **LUMIERE/PREDICT-GBM:** previously examined and not fresh confirmation.

### Next step

Do not download Track B data. For Track A, the two GEO files above are the smallest direct-public feasibility pair identified that may reach six donors. No download has started. Before scoring, verify the sample-to-donor map, confirm a comparable cell population across both studies, freeze the split/comparator/failure rules, and amend the protocol. Do not promise a positive result.

## Search coverage

Targeted source checks covered GEO/BioStudies, GBM-Space, TCIA, Zenodo, NGDC, and linked papers. Parallel read-only audits also checked imaging, spatial/molecular, clinical, pathology, proteomics, and cell-line repositories, including GDC, Synapse, PDC/PRIDE, and DepMap. This is a broad targeted search, not proof that every dataset has been found.
