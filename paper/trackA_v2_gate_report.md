# Track A v2 gate report — HRA009021 HBI/LBI (2026-10-31; updated 2026-11-03)

## Current decision: STOP before data access

**Blocking item:** a public, explicit map from each GSA sample name/accession to donor code and HBI, LBI or Tumor. The existing metadata-only audit links 10 donor codes to the paper's region labels, but not individual GSA samples to regions. Do not infer region from suffix or expression. The user has ruled out author contact; proceed only if the map becomes public. Do not download FASTQs, seal a protocol, or score this cohort.

## Work done in this session
- Checked the stated project facts against the repo files (all match; see table).
- Read the full-text Data availability section (Europe PMC, PMC11685534).
- Read Supplementary Table 3 (pages 15 of `41467_2024_55243_MOESM1_ESM.pdf`, already on disk).
- Checked the NGDC BioProject PRJCA031420 page.
- Tried the GSA-Human page for HRA009021 (timed out twice from this computer).

## Work NOT done
- No FASTQ, count matrix or expression file was downloaded or read.
- No protocol was sealed. No model was trained. Nothing was scored.
- No author contact was made. Do not contact authors.
- No replacement cohort search beyond the 2026-10-30 screen.

## Facts checked

| Claim | Status | Evidence |
|---|---|---|
| 10 donors, 30 samples, 3 per donor | Verified | `output/trackA_v2_data_check/hra009021_metadata_audit.json` (`samples_per_donor`) |
| IDH-wildtype status for scRNA-seq subset | Supported at study-cohort level | The paper reports 26 IDH-wildtype GBM patients and a 10-person scRNA-seq subset from the same study; exact accession-to-paper case mapping remains incomplete. |
| GSA workbook has no region labels | Verified | same file, `region_labels_present: false` |
| Supplement Table 4 labels HBI/LBI/Tumor by donor code only | Verified | `output/trackA_v2_data_check/hra009021_supplement_audit.json` (`direct_gsa_sample_to_region_map_found: false`) |
| Raw FASTQ total 764,601,130,806 bytes (about 765 GB), 60 files | As recorded | `output/trackA_v2_data_check/hra009021_file_inventory.json` |
| GBM-Space nuclei lack a nucleus-to-region map | As recorded | `output/trackA_v2_data_check/s14_metadata_check.json`, `region_provenance.json` |
| SenMayo, OXPHOS and TAM sealed tests failed adjusted criteria | Verified | `output/trackA_senmayo/result.json`, `output/trackA_oxphos/result.json`, commit 2874524 (TAM). Preserved, not rerun. |

## New evidence on the crosswalk (all negative)
1. **Data availability statement** (PMC11685534): raw scRNA-seq at GSA-Human HRA009021; "All other data included in this study are provided in Supplementary data and Source data." No code repository, no processed matrix, no sample metadata table named.
2. **Supplementary Table 3** (per-biopsy CBF and ROI): rows numbered 1–26 with LBI/HBI/Tumor columns. No GSA sample names, no donor codes. Its 26 rows do not match the 10 sequenced donors, so it cannot serve as a map.
3. **BioProject PRJCA031420**: no linked sample records or region attributes shown.
4. **Earlier audits**: the Source Data workbook and GSA workbook contain no sample-to-region link.

Conclusion: no public source maps GSA sample names to regions. Suffix order (01/02/03) must not be assumed. Inferring the region from expression is not allowed, because the label must not come from gene expression.

## Access and reuse terms
- **Access:** the official [GSA-Human record](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA009021) lists HRA009021 as Open access. The [NGDC policy](https://ngdc.cncb.ac.cn/gsa-human/document/Principle_of_Accessing_Human_Genetic_Resource_Data_in_NGDC_V1.pdf) says open-access data may be used directly for research.
- **Use limits:** the policy limits use to the stated research scope and research group/collaborators, bars redistribution and participant identification, and allows submitters to set specific restrictions. No accession-specific license or restriction is displayed in the public record. The article’s CC BY-NC-ND license covers article text, not dataset rights.
- **Status: CONDITIONAL.** No-request access and repository-level research-use terms are documented. Confirm school/fair suitability and any accession-specific restriction before data use.

## Exact counts
- Donor records: 10 (GBM23, 35, 36, 37, 41, 44, 45, 51, 52, 53). Countable paired donors under Track A's exact sample-region rule: 0, because sample accessions lack region labels.
- Samples: 30 (10 HBI + 10 LBI + 10 Tumor by the paper). Per-sample identity: unknown.
- Paired HBI/LBI subset size: cannot be computed until the map exists. Upper bound if all 20 HBI/LBI runs are used: about 2/3 of the 765 GB total (roughly 510 GB, estimate only, assumes even file sizes).

## Contact status
The user has ruled out author contact. The earlier unsent message draft is removed. Continue only if required mapping becomes publicly available.

## Unresolved questions
1. The sample-to-region map (blocking).
2. Cell-level malignant-label fields, provenance, and validation.
3. Patient-level overlap with prior Track A use; exact public identifiers do not establish cross-study identity.
4. Any accession-specific use restriction and school/fair suitability. Repository-level terms are documented above.
5. Storage and runtime, only if all metadata gates later pass and a dated scope amendment authorizes data access.

## If a public map becomes available
1. Check the map against donor, sample, experiment and run accessions in `hra009021_metadata_audit.json`.
2. Recheck the six-donor, IDH, malignant-label, reuse, and independence gates. A public map alone does not qualify the cohort.
3. Report evidence for review. Do not download or score data unless a dated scope amendment and explicit user authorization cover the selected cohort and analysis.

## If the map never comes
Stop this cohort. Broad scouting stays paused. Reopen metadata screening only for a specific named accession that could plausibly pass the current Track A New gates. Keep the six-donor minimum and independence requirement. Record any change of dataset or endpoint in a dated scope amendment.
