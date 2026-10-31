# Track A dataset audit — Core vs Periphery, donor-held-out (2026-10-31)

## Decision (pass 3): GO for a descriptive-only analysis of Yu 2020 (GSE117891), with 6 strict / 8 lenient donors. The protocol and pre-expression amendment A1 are frozen in `paper/trackA_yu2020/PROTOCOL.md`. The 39 MB matrix has not been downloaded; no expression has been scored.

## Decision (pass 2, superseded): CONDITIONAL GO for one cohort, Yu 2020 (GSE117891). Download nothing yet.
Pass 1 (below) found no GEO cohort with 6 paired donors. Pass 2 searched CELLxGENE and found that the CC BY 4.0 GBmap atlas carries Core/Periphery `sector` labels for 17 donors in 3 studies. One study, Yu 2020, has **9 paired donors in one assay** (STRT-seq). That puts it in the appendix's 6–9 band: **descriptive nested evaluation only, no superiority claim.** One item must be confirmed before download (see "Remaining blocker").

## Pass 2 — CELLxGENE / GBmap (read remotely; obs metadata only, no expression)
Source: CELLxGENE collection 999f2a15-3d7e-440b-96ae-2c806799c08c, "Core GBmap" (bioRxiv 10.1101/2022.08.27.505439), H5AD 8,127,057,404 bytes, license CC BY 4.0 (CELLxGENE standard). Only the `obs` table was read through HTTP range requests.

`sector` values: unknown 320,127 cells; Core 10,591; Periphery 7,846. Donors with both Core and Periphery:

| Study (assay) | Donor | Core cells (neoplastic) | Periphery cells (neoplastic) |
|---|---|---|---|
| Yu2020 (STRT-seq) | S1 | 253 (231) | 124 (46) |
| | S2 | 196 (121) | 62 (3) |
| | S3 | 201 (127) | 55 (39) |
| | S5 | 291 (260) | 35 (2) |
| | S6 | 242 (192) | 113 (83) |
| | S11 | 766 (723) | 244 (17) |
| | S12 | 483 (158) | 150 (39) |
| | S13 | 866 (757) | 112 (2) |
| | S15 | 270 (0) | 230 (3) |
| Wu2020 (BD Rhapsody) | LWZ, PXC, YCM, ZBY | 85–3,391 | 100–2,762 |
| Darmanis2017 (Smart-seq2) | BT_S1, S2, S4, S6 | 203–957 | 55–550 |

"Neoplastic" = GBmap `annotation_level_1`, an expression-derived annotation; used here only to show composition.

### Yu 2020 provenance map
- Paper: Yu K. et al., PMID 34692159. GEO GSE117891: 137 libraries, 13 glioma + 1 metastasis patients, 73 "surgery points".
- GEO cell table `GSE117891_Sample_barcode_cell_information.xlsx` (231,046 bytes, SHA-256 2669c32aa678d941e60e19dd786cb581c9212bb9e088967cb2bf2af52aa016da, downloaded) maps every cell → patient → point (e.g. S1P1_SC1 → S1 → S1P1). **It has no core/periphery column.**
- GBmap code (`github.com/ccruizm/GBmap`, notebook `01_data_import_and_preprocessing.ipynb`) assigns the region from an author file `clinical_info.xlsx`, sheet "Point and Cell information": "Tumoral" → Core, "Peritumoral" → Periphery. It excludes S4 (gliosarcoma), S7 (IDH-mutant), S8, S9 (OA), S14 (O). So the label is author-supplied per surgical point, not derived from expression.
- Each point is one library. Region is a within-donor contrast, so library is nested in region inside each donor. That is unavoidable and must be reported.
- Prior exposure in this project: none. GSE117891 is not used in any `src/`, `output/` or `paper/` file. Darmanis 2017 is only cited (in `paper/trackA_tam_5ala/PROTOCOL.md`).

### Pass 3 — label table found at the primary source (blocker resolved)
- Source: Yu K. et al. 2020, *Natl Sci Rev* 7:1306, doi 10.1093/nsr/nwaa099, PMC8289159, license **CC BY 4.0** (Europe PMC record). Retrieved via the Europe PMC supplementaryFiles API.
- Files (git-ignored, `data/external/yu2020_supp/`):
  - `nwaa099_supplemental_files.zip`, 14,797,073 bytes, SHA-256 99428d35804417e9576ab28080804109f13ff9f4b8513157597efb8a9ef43db2
  - `Table S1. ClinicalInformation.xlsx`, 15,125 bytes (inside the zip), SHA-256 0676607ff4c57348fc043ab4bace1c3985617d0d90507742344ce19642cb772a
- Sheet "Point and Cell information" gives every sampling point a `Sample type` (Tumoral / Peritumoral / Normal) and an MRI `Contrast enhancement` flag. The paper text says the points were marked on MRI ("Red and yellow dots mark locations of tumoral and peritumoral sampling points"). Point numbers are **not** in region order (e.g. GS1: P1, P4–P7 tumoral; P2, P3 peritumoral), so suffix order must never be used.
- GBmap's sector labels match this table for all 9 donors.

### Donor eligibility from Table S1 (sheet "Clinical_Info")
| Patient | Histology | IDH1 R132 | Tumoral points | Peritumoral points | Eligible |
|---|---|---|---|---|---|
| GS1 | GBM (O included) | − | P1,P4,P5,P6,P7 | P2,P3 | Yes |
| GS2 | GBM | − | P1,P2,P4 | P3 | Yes |
| GS3 | GBM | − | P2–P7 | P1,P8 | Yes |
| GS5 | GBM | − | P1–P4 | P5 | Yes |
| GS6 | GBM (O included) | − | P1–P4 | P5 | Yes |
| GS13 | GBM (O included) | − | P4,P5,P6 | P1,P2 (P3: 0 cells) | Yes |
| GS11 | GBM | not reported | P4,P5,P6 | P2,P3 | Conditional: IDH status missing |
| GS12 | GBM | not reported | P1–P4 | P5 | Conditional: IDH status missing |
| GS15 | **not reported**, grade not reported | not reported | P3,P4,P5 | P1,P2 | **No.** GEO says the cohort includes 1 brain metastasis; GS15 is the only patient with no histology. GBmap includes it; this audit does not. |
| GS4 | Gliosarcoma | n/a | — | — | No |
| GS7, GS8, GS9, GS14 | IDH-mutant or grade II–III | — | — | — | No |

- **Strict count (IDH-wildtype confirmed): 6 donors.** This is exactly the appendix minimum. Leave-one-donor-out leaves 5 training donors and 4 in every inner split, so the inner-fold rule (≥4) holds.
- **Lenient count (GBM histology, IDH not reported): 8 donors.**
- Either way: 6–9 band → **descriptive nested evaluation only.**
- Label and MRI enhancement are nearly collinear (Tumoral mostly "+", Peritumoral "−"). The region label is MRI/surgical, not expression-derived.

### Protocol and next data gate
- **Complete:** eligibility rule is frozen (six strict donors primary, eight lenient donors sensitivity), as are the cell-type-only composition diagnostic and failure rules. Pre-expression Amendment A1 clarifies the exact metadata join and interpretation. The donor table, protocol, and analysis-plan appendix are hash-pinned in `paper/trackA_yu2020/SHA256SUMS`.
- **Not done:** the 39 MB expression matrix has not been downloaded or scored. A later run must first verify the protocol hashes, exact cell-barcode joins, and metadata provenance. Results remain descriptive only; no superiority or tumor-cell biology claim is allowed with six to eight donors.

### Scientific risk (must be stated with any result)
- Composition confound: Core points are mostly neoplastic and Periphery points mostly non-neoplastic (e.g. S13: 757/866 vs 2/112). A classifier can score well by learning cell type, not region. Required control: report accuracy within cell-type strata and against a cell-type-only baseline.
- Few cells per point (35–866). STRT-seq is a low-throughput plate assay.
- Pooling with Wu2020 and Darmanis2017 (17 donors) would mix 3 assays, with study fully confounded with platform. That needs a separate draft amendment and is exploratory only.

### Smallest next download (after the blocker clears and the protocol is frozen)
| File | Source | Size |
|---|---|---|
| GSE117891_all_6148.umi.count.matrix.tsv.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117891/suppl/ | about 39 MB (FTP listing) |
| Point → region table | primary source to be found | small |

Not the 8.1 GB GBmap H5AD.

## Pass 1 decision (superseded): NO-GO
No GEO cohort checked in pass 1 had 6 or more eligible paired donors in one coherent assay and cell population. The Track A question stays unchanged. No amendment is proposed here.

## What was done
- Read the appendix Track A section (eligibility rule: section 2, lines 41–49) and the prior audits `paper/trackA_v2_data_check.md` and `paper/trackA_v2_gate_report.md`.
- GEO web pages returned a robot check (CAPTCHA). It was not bypassed. Records were read through the official NCBI E-utilities API (`esummary`, db=gds) and the GEO FTP SOFT files.
- Downloaded 3 small metadata files only (git-ignored, `data/external/trackA_audit_meta/`):

| File | Source | Bytes | SHA-256 |
|---|---|---|---|
| GSE273274_family.soft.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273274/soft/ | 2,896 | 11832a094b8cc639f4453161712fd83f831bf6cb79607a8f57c5096f95f3bf36 |
| GSE162631_family.soft.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162631/soft/ | 2,069 | 4dd8a5c96bc1c739329dd1d3ecd4404c637e4a52ae7bb97f01ac1ecd371a5ae4 |
| GSE318562_family.soft.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE318nnn/GSE318562/soft/ | 4,874 | ab534c8d01ac2f229c62183dfc32be4913b310176a155d292776b928d805db87 |

- No expression file was downloaded. Nothing was scored. No person was contacted.

## Comparison

| Cohort | Assay / population | Eligible paired donors (verified) | Label source | Blocker |
|---|---|---|---|---|
| GSE273274 (PMID 39933531) | 10x scRNA-seq, "Glioma cells" | **≤3, donor identity unverified** | GEO field `tissue: GBM centre` / `GBM surround` | Samples are titled "replicate 1–3" (files R1/R2/R3). GEO does not state that the replicates are 3 patients. Supplement: `GSE273274_RAW.tar`, 706 MB. |
| GSE162631 (PMID 34228647) | 10x scRNA-seq, **CD31+ endothelial-enriched** (Dynabeads) | 4 (R1–R4, "4 glioblastoma patients" in design) | Sample title only | Endothelial cells only. Every sample, including "peripheral tissue", has `Sample_description = tumor core` (metadata conflict). Supplement: `GSE162631_raw_counts_matrix.zip.gz`, 806 MB. |
| GSE318562 (no PMID yet; released 2026-06-12) | 10x Visium V1 spots, mixed OCT and FFPE | 2 by GEO names (WU1221 Core2/Edge1, WU1227 Core5/Edge2) | Sample title only | Design says 22 sections, 9 patients, "edge-core axis", but 18 sections carry no region label on GEO. Spots, not cells. Fixation differs by batch. `GSE318562_RAW.tar`, 674 MB. |
| GSE237183 | Visium spots | 5 | Navigation-guided location | Fails donor and topic gates (`output/trackA_v2_data_check/alternate_cohort_screen_20261030.json`). |
| GBM-Space | snRNA / spatial H5AD | 4 (strict same-section core + leading edge) | Feature names only | Below 6; nucleus-to-region values not verified (`output/trackA_v2_data_check/s14_metadata_check.json`). |
| HRA009021 | scRNA-seq | 10 donors, map unknown | Paper tables by donor code only | Sample-to-region map not public; needs author contact, which the rule forbids (`paper/trackA_v2_gate_report.md`). |
| UCSC multiomic-gbm (historical) | sc multi-omic | 6–7 core+periphery | Source label | Region = library confound; already exposed. Leak-free C-GAT 0.638 vs cells-per-patient baseline 0.646 (`vault/claim_ledger.md` line 43; scripts 88, 95). Not validation. |

## Why the cohorts cannot be pooled
- GSE273274 + GSE162631 gives at most 7, and only if the GSE273274 replicates are distinct patients (unverified).
- The two studies profile different cells (glioma cells vs CD31+ endothelium). Region would be nested in study, and study is fully confounded with cell population and lab. A pooled classifier would learn the study, not the region.
- Visium spots (GSE318562, GSE237183) are a different unit and modality from dissociated cells. Mixing them changes the question.

## Go / no-go against appendix thresholds
- Fewer than 6 eligible paired donors in every coherent cohort → **stop as infeasible for the planned analysis.**
- No descriptive nested evaluation is run. No superiority claim is possible.

## Gaps (not checked)
- Pass 2 searched CELLxGENE collections (name/description match "glioma"/"glioblastoma"). Also seen: SL040/SL057 (2 donors, 5-ALA multisampling), Darmanis (4), a 60-donor glioma atlas (no region field checked). Still not searched: HCA Data Portal, BioStudies/ArrayExpress single-cell, Zenodo/Figshare, HTAN, CPTAC, pathology image sets.
- GSE273274 and GSE318562 papers' supplementary tables were not read. GSE318562 has no linked paper yet.
- GEO data are public; reuse follows NCBI GEO policy and citation of the source paper. No file-level license was found. Not legal advice.

## What would change the decision (no contact needed)
1. A published GSE318562 paper whose supplement maps all 22 sections to core or edge per patient. That would allow a Visium-based section-level amendment (9 patients, still below 10: descriptive only).
2. A new public single-cell cohort with ≥6 paired core/periphery donors in one assay, found in the unsearched repositories above.
