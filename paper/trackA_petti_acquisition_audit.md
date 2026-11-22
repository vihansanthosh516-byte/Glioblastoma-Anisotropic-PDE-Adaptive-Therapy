# Track A Petti acquisition and eligibility audit

Acquired the histology overview, author specimen table, Xenium plotting code, and release and code manifests in the Git-ignored data/external/trackA_petti_15014838 folder. The PDF MD5 matches the [latest Zenodo manifest](https://zenodo.org/records/20709378): b2b335e3b1a56cd9c35e33db7cbcc187. The annotated Xenium object is 4485875917 bytes with published MD5 c6b5d15167b0f82328fc0e29a5cc41e6. It must pass that checksum before use.

## Regional eligibility

The [author specimen table](https://github.com/Petti-Lab/GBM_Spatial_paper_figures_code/blob/444e2757f06608698cc1bf9b6c8b94988f10d2ac/Figure_2/spatial_data_metadata.txt) lists exact core and edge specimen labels under six shared code prefixes: B185, B186, ST073021, WU1220, WU1221, WU1227. This count is a literal table screen. It does not establish approved patient pair eligibility. The [author README](https://github.com/Petti-Lab/GBM_Spatial_paper_figures_code/blob/444e2757f06608698cc1bf9b6c8b94988f10d2ac/README.md) identifies this table as Figure 2 metadata and describes the Xenium ROI boxes as figure annotations. The table and code inspection do not verify an independent pathology-to-Xenium-cell registration. Do not use expression-derived cluster labels as pathology ground truth. Require a released mask or cell-to-ROI crosswalk, coordinate units and transform provenance, and a deduplicated count of untouched eligible patients. Density must be a baseline because it partly defines the pathology regions.

## Measured TMZ eligibility

The [organoid study](https://ega-archive.org/studies/EGAS50000001221) has controlled-access sequencing. Its public catalog does not establish a complete join among untreated expression, measured TMZ viability and independently derived sister cultures by patient. Record the coded crosswalk and request authorized access before evaluation. HGCC remains spent.

## Completion status

The original Artzi benchmark and the frozen simulations ran. Artzi, Yu and HGCC remain spent. This audit closes the public metadata screen for the inspected Petti files. A new real-data region or TMZ test has not been run. No new patient expression values or viability outputs were opened. Positive real-data validation remains open. [Artzi result](trackA_artzi_result_2026-11-17.md) | [Joint error simulation](trackA_joint_error_result.md) | [Recovery plan](trackA_recovery_plan.md)

## Xenium cohort 3 specimen join

The downloaded [GEO family index](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335275&targ=all&form=text&view=full) lists the following joins to the author table above. Normalize leading zeros in terminal specimen numbers only. GEO explicitly supplies the patient code for each specimen. Labels apply at specimen level; a cell-level pathology join remains unverified.

| GEO patient | Edge accession | Core accession | Author specimen labels |
|---|---|---|---|
| B189 | GSM9809293 | GSM9809294 | ST073021-2 / ST073021-4 |
| WU1220 | GSM9809295 | GSM9809296 | WU1220-1 / WU1220-3 |
| WU1221 | GSM9809297 | GSM9809298 | WU1221-Edge1 / WU1221-Core2 |

These are three candidate patients from this series. Even if every paired comparison favored a candidate with no ties, the smallest possible one-sided exact sign-test p-value would be 1 / 2^3 = 0.125 before correcting for multiple comparators. This is a design calculation. The series cannot independently meet the previous 0.05 superiority rule. Subdividing specimens into blocks or counting Visium duplicates as new patients would not fix this.

## Specific data request for closure

For the region endpoint, request a coded specimen-to-patient crosswalk covering the other Xenium cohorts and an independently annotated core/edge ROI export in Xenium coordinates (or masks plus registration transforms), with unit, section and assay identifiers. Request whether samples overlap the above series or prior Track A cohorts. For TMZ, request the deidentified baseline expression file or authorized access, assay AUC and endpoint definition, plus origin-to-culture-to-assay mapping and which cultures arose independently. No emails were sent.
