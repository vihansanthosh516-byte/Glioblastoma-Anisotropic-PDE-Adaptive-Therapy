# Track A v2 gate report — HRA009021 HBI/LBI (2026-10-31)

## Decision: WAIT for one named item

**Item:** an explicit map from each GSA sample name (`GBMxx01/02/03`, accessions HRS…) to HBI, LBI or Tumor, from the study authors.
Until it arrives: do not download FASTQs, do not seal, do not score. If no map arrives, stop this cohort (see "If the map never comes").

## Work done in this session
- Checked the stated project facts against the repo files (all match; see table).
- Read the full-text Data availability section (Europe PMC, PMC11685534).
- Read Supplementary Table 3 (pages 15 of `41467_2024_55243_MOESM1_ESM.pdf`, already on disk).
- Checked the NGDC BioProject PRJCA031420 page.
- Tried the GSA-Human page for HRA009021 (timed out twice from this computer).

## Work NOT done
- No FASTQ, count matrix or expression file was downloaded or read.
- No protocol was sealed. No model was trained. Nothing was scored.
- No message was sent to the authors. The user must send it.
- No replacement cohort search beyond the 2026-10-30 screen.

## Facts checked

| Claim | Status | Evidence |
|---|---|---|
| 10 donors, 30 samples, 3 per donor | Verified | `output/trackA_v2_data_check/hra009021_metadata_audit.json` (`samples_per_donor`) |
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
- **Download:** earlier audit records GSA-Human status "Open access" for HRA009021. Not re-checked today (page timed out). Re-check before any download.
- **Reuse:** no file-level license found. The article is CC BY-NC-ND (Europe PMC record). That covers the article text, not the data. NGDC/GSA-Human data-use policy still needs reading, plus any school or fair rule on human genomic data.
- **Status: UNRESOLVED.** Public download does not prove reuse rights.

## Exact counts
- Donors: 10 (GBM23, 35, 36, 37, 41, 44, 45, 51, 52, 53). Exactly the minimum.
- Samples: 30 (10 HBI + 10 LBI + 10 Tumor by the paper). Per-sample identity: unknown.
- Paired HBI/LBI subset size: cannot be computed until the map exists. Upper bound if all 20 HBI/LBI runs are used: about 2/3 of the 765 GB total (roughly 510 GB, estimate only, assumes even file sizes).

## Named item: request to the authors (user sends)
Corresponding authors (from PMC11685534): Ying Hu (huying@hit.edu.cn), Hang Ji (jihang0723@163.com), Shaoshan Hu (shaoshanhu421@163.com).

Draft message:
> Dear Dr. Hu and Dr. Ji,
> I am a student researcher reading your 2024 Nature Communications paper "Turning attention to tumor–host interface and focus on the peritumoral heterogeneity of glioblastoma." I would like to reanalyze the public data at GSA-Human HRA009021. The sample names there are GBMxx01, GBMxx02 and GBMxx03. Could you tell me which sample is HBI, which is LBI and which is Tumor for each patient? A table with sample accession (HRS…) and region would be enough. I will cite your paper and follow the GSA-Human data-use terms.
> Thank you,
> [name]

## Unresolved questions
1. The sample-to-region map (blocking).
2. The GSA-Human / NGDC reuse terms, and any school rule (blocking before download).
3. Storage: about 510–765 GB of FASTQ plus alignment space. The current free disk space is not checked.
4. CPU runtime for Cell Ranger alignment of 20–30 runs on this laptop is not estimated.

## If the map arrives
1. Check the map against donor, sample, experiment and run accessions in `hra009021_metadata_audit.json`.
2. List the exact paired HBI/LBI runs and their byte sizes from `hra009021_file_inventory.json`. Ask the user to approve that subset.
3. Write a protocol amendment and a frozen protocol for mentor review before scoring: donor-held-out folds (all cells of a donor in one fold, preprocessing and tuning on training donors only); adversarial scVI vs ordinary scVI vs adversary-off ablation vs listed baselines on identical folds; the hybrid-state biology test kept separate from the model test.

## If the map never comes
Stop this cohort. Search for a replacement only under the same gates: 10 or more donors, single-cell or single-nucleus data, region labels not from expression, a verified sample-to-label map, and clear reuse terms. Record any change of dataset or endpoint as a protocol amendment.
