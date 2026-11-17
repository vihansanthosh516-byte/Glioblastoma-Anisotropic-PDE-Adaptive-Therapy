# CARE external-replication eligibility gate (metadata only)

Date and access date: 2026-11-16. Author: Claude. Follows `paper/recurrence_reliability_care_metadata_audit_2026-11-16.md` and `paper/recurrence_reliability_postrun_review_2026-11-16.md`. The locked Wang study (protocol, config, cohort, run `run_20261116_locked_a1`) is unchanged.

Not opened in this work: any CARE expression file (GEO `*_umi_counts.RDS.gz`, `GSE274546_RAW.tar`, the Smart-seq2 matrix), any CARE per-cell annotation file, state or score file, or any survival or outcome field. Machine-readable record: `output/recurrence_reliability_care_metadata/care_eligibility_manifest.json` (built by `src/recurrence_reliability/care_metadata.py`; source hashes inside).

## Decision

**Update after the blinded read (same day): 33 CARE pairs pass every gate that metadata can settle. Two checks remain and need the count files: the barcode match to GEO matrix columns and the malignant recount after the Wang filters. No CARE score has been computed.** Original decision (before the read): sample-level gates pass for 42 patients; confirmed eligible pairs 0. CARE can plausibly support an external replication, but only after the per-cell annotation file is shown to be safe to read and the per-specimen malignant support is counted.

| Gate | State | Count or evidence |
|---|---|---|
| C1 Specimen identity (Table S1 ↔ GEO) | Pass | 121 of 121 specimens agree on patient, time point, site code, stage, and GEO file name; 0 issues |
| C2 Primary + first recurrence present | Pass for 52 | 49 patients with P+R1, 3 with P+R1+R2, 4 with P+R2 only, 3 with R1+R2 only |
| C3 IDH-wildtype GBM | Pass for all | Table S1: "Glioblastoma, IDH-wildtype", IDH "wildtype" for 121 of 121 (author-reported) |
| C4 Treatment context (Wang rule) | Pass for 42 | Untreated primary and first recurrence after radiotherapy + TMZ only; 10 of 52 fail |
| C5 Access and terms per component | Pass for the needed public components, with caveats | See "Access and terms" |
| C6 Barcode ↔ specimen ↔ malignant-label join | **Specimen level: pass. Barcode ↔ GEO matrix columns: open** | Blinded read (below): 121 Sample codes, each mapped to exactly one P#T# ID present in Table S1/GEO; all 429,305 CellIDs unique and of the form `<Sample>-<16-base barcode>`. The match to GEO matrix column names needs the count files and is checked in the QC step after protocol lock |
| C7 Malignant support ≥ 50 nuclei per specimen | **Pass for 33 pairs at author-QC level**; Wang-filter recount open | 9 of 42 pairs have a specimen with < 50 author-labelled malignant nuclei; the Wang filters (≥ 200 genes, MT ≤ 2.5%) can only remove nuclei, so the final count is ≤ 33 |
| C8 Patient overlap with Wang and prior exposure | Pass at provenance level | No shared institution; no identity crosswalk exists; no prior project use of CARE expression |
| C9 Scoring-reference choice | Resolved on estimand grounds (proposal) | Refit as primary; transport as secondary (see below) |

## C2 correction to the earlier audit

The earlier audit described the 59 CARE patients as having "primary and first-recurrence labels". That is not true for all of them. GEO source labels and Table S1 agree: 52 patients have both a primary and a first-recurrence specimen; 4 have a primary and a second recurrence; 3 have only first and second recurrences. Spitzer et al. state the same ("three pairs that consisted of two recurrent samples and four pairs that consisted of a primary and a second recurrence"). For the 3 patients with three time points (P20/TK01, P32/TK13, P58/NL02; all three pass C4), the primary and first-recurrence specimens define the pair and the second recurrence is not used. This rule follows the locked Wang question (primary vs first recurrence) and was fixed before any CARE cell data were seen.

## C4 treatment mapping

Wang et al. describe their cohort as untreated primaries with first recurrences after standard of care only (TMZ, radiotherapy, surgery). Wang's treatment status is a cohort-level author statement. CARE Table S1 gives per-specimen fields "Radiation before the surgery" and "Alkylate agents before the surgery". The transported rule: primary with both fields "No"; first recurrence with radiation "Yes" and alkylating agent exactly "Yes (TMZ)". 42 pairs pass. Failures: 4 with no radiotherapy or alkylator before the first recurrence surgery; 2 with radiotherapy but no alkylator; 1 with TMZ but no radiotherapy; 3 with lomustine (CCNU) or nimustine (ACNU) added. All 52 primaries were treatment-naive by these fields. Table S1 columns read (and only these): Sample ID, Tumor ID, Patient ID, Time point, Primary or recurrent, Tissue source, Lab processed, 10x, Smart-seq2, Diagnosis, IDH mutation status, Radiation before the surgery, Alkylate agents before the surgery. Not read: sex, age, MGMT, tumor location, site of recurrence, steroids, number of TMZ cycles, surgical interval, overall survival, outcome.

The 42 eligible pairs come from 7 sites (TK 13, SN 7, FR 6, SM 5, NL 5, DU 4, MD 2). Processing: Protocol 1 for 30 pairs, Protocol 2 for 12 (GEO "protocol number"; Protocol 2 equals Lab3 in Table S1). No patient has specimens from two protocols, so the within-patient change is not confounded by protocol.

## C5 access and terms

| Component | Route | Terms found | Needed? |
|---|---|---|---|
| Processed 10x counts | GEO GSE274546, one `GSM…_P#T#_umi_counts.RDS.gz` per specimen (per-GSM FTP folders), public, no login | GEO series: "raw data are not available due to patient privacy concerns". The GEO disclaimer page returned a CAPTCHA on 2026-11-16 and was not re-read; the earlier Wang note records GEO's general statement that NCBI places no use restrictions but submitters may hold rights. | Yes: 84 files, 2,975,201,246 bytes listed for the 42 pairs (RAW.tar 4,272,261,120 bytes) |
| Per-cell cell-type labels | GitHub `dravishays/GBM-CARE-WT` (commit tree `383b0a3`), `data/celltype_meta_data_2025_01_08.RDS` (2,344,381 bytes); same repository archived at Zenodo 10.5281/zenodo.14966015 (`GBM-CARE-WT.zip`, 122,875,674 bytes, md5 eb811bed…) | GitHub: MIT License (written for software). Zenodo record: CC BY 4.0, open access. | Yes |
| Table S1 | Nature Genetics supplementary file `41588_2025_2168_MOESM2_ESM.xlsx` (26,702 bytes, SHA-256 `6bfb81f4…45a0`). A copy is also in the repository (`data/spitzer_supptable1.xlsx`, 26,432 bytes, different hash; content not compared) | Article license CC BY-NC-ND 4.0. Only 13 identity/treatment columns were extracted; the workbook stays local and uncommitted. | Yes (selected fields) |
| Smart-seq2 | GEO GSE274548, one aggregate GSM (14,401 nuclei) | Same as GEO | No (no patient/time point linkage in GEO; overlaps the 10x tumors) |
| Raw reads | DUOS-000475 to DUOS-000480 (controlled) | Consent-limited | No; not requested |
| DNA tables | Synapse `syn26464346` ("CARE IDH-wildtype (GBM)") | Synapse API lists no access requirement on the project entity; download needs a Synapse account and the Synapse Terms of Use | No |

The terms gate passes for local, non-redistributed analysis of the public GEO counts and for reuse of the CC BY 4.0 label file with attribution. Neither the GEO count files nor the label file will be committed.

## C6 and C7: why the cell-level gates stay open

GEO holds only per-specimen count matrices; no per-cell label is deposited there. The author code (`R/GBM-CARE-WT_Preprocess.R`, read as code only) builds cell IDs as `<Sample>-<10x barcode>` from Cell Ranger v3.1.0 filtered matrices, where `Sample` is a site-code-plus-surgery string (the code refers to `NL022` for patient P58 time point 3). The released label file `celltype_meta_data_2025_01_08.RDS` is used in the figure scripts with columns `CellID`, `Sample`, `ID` and `CellType` (value `Malignant` for CNA-classified malignant cells). The preprocessing code also joins clinical fields into its working metadata (`OSmonth`, `PatientStatus`, `MGMTstatus`, `SurgicalIntervalmonth`). It is therefore not established that the released label file is metadata-only: it may carry survival columns. A second file, `malignant_meta_data_2025_01_08.RDS`, carries malignant cell-state assignments (`State`, `isCC`) and is an outcome-adjacent file; it must not be used.

Not verified: (a) that GEO RDS column names equal the label-file `CellID`s; (b) per-specimen counts of malignant nuclei after the Wang filters (≥ 200 genes, mitochondrial fraction ≤ 2.5%). The author code applied its own QC (200 to 10,000 genes, mitochondrial ≤ 3%) and a confidence filter before labelling, so cells absent from the label file are expected and must be handled by a frozen rule (proposed protocol, section 3). One author code comment states that P58T3 has no malignant cells; P58T3 is a second recurrence and is not in the pair set.

Safe ways to close C6/C7, in order of preference (needs your decision):
1. A blinded, column-restricted extraction from `celltype_meta_data_2025_01_08.RDS`: load it, print column names only, write `CellID`, `Sample`, `ID`, `CellType` for the 84 eligible specimens, discard every other column without printing, and report per-specimen malignant counts. Exposure: composition counts (the same kind recorded in the Wang Phase B). Any survival column would be held in memory but never displayed or written.
2. After that, read only the column names (barcodes) of the 84 GEO RDS files to verify the join. This touches expression files and should happen only after the external protocol is locked.

## C8 overlap and independence

CARE institutions (Spitzer et al. Methods; Table S1 "Tissue source"): Tokyo University 30 specimens, Pitié-Salpêtrière 28, Seoul National University 20, St. Michael's 14, NorLux 13, Duke 12, MD Anderson 4. Wang specimens are all UCSF (`SF…`). The CARE authors themselves compare "CARE dataset" against "Wang et al dataset" as separate cohorts (`R/spitzer_edf2b_f.R`). No public patient identity crosswalk exists, so independence is supported by provenance, not proven at patient level. The SNU cohort's DNA data were previously released in GLASS; this concerns DNA, not these snRNA counts. Project exposure: no CARE expression, label, or score has been used in this project (Git search of tracked files on 2026-11-16 found only planning notes). Literature knowledge disclosed: Spitzer et al. report a lower malignant-cell fraction at recurrence; Alexopoulos and Karampelas (2026) analysed CARE AC-like changes. Classification for this study: **fresh** for project outcomes.

## C9 scoring reference: transport versus refit

Only procedure R (Tirosh/Neftel control-gene score) depends on a reference; procedure A (UCell) is reference-free and identical under both options.

- **Refit (same algorithm, CARE reference).** Rebuild aggregate-expression bins, control genes (same seed rule) and gene centres from the CARE eligible malignant pool. Estimand: discordance between the two scoring procedures when each is run as specified inside a new cohort. This is how the Wang estimate was produced (in-cohort pool), so it is the closer replication of the Wang endpoint. Consequence: patient indicators share the CARE pool (the dependence issue raised in the post-run review).
- **Transport (Wang-fitted reference).** Reuse the Wang bins, control-gene sets and centres; drop genes absent from CARE. Estimand: discordance when a reference fitted in one cohort is deployed in another. It tests a fixed instrument across platforms and labs, not the same quantity. The Wang run saved only control-gene counts, so a transported reference would first have to be regenerated deterministically from the Wang inputs (no scores computed) and checked against the saved counts (4,027 and 4,041).

Proposal, made on estimand grounds before any CARE data: refit is primary (replication); transport is a declared secondary analysis (deployment). Neither is chosen by results. Details: `paper/recurrence_reliability_care_protocol_PROPOSED.md`.

## Precision at the sample-level ceiling

If all 42 pairs pass C6/C7 (an upper bound, not a forecast): with 0 discordant patients the two-sided exact 95% upper bound is 0.084, still above the 5% tolerance; the lower bound exceeds 5% from 6 discordant patients. CARE alone therefore cannot confirm reliability below 5%; it can show non-trivial discordance or remain inconclusive.

## Blinded read of the author cell-type table (user-authorised, 2026-11-16)

Script: `src/recurrence_reliability/care_labels.py` with a standard-library RDS reader (`src/recurrence_reliability/rds_reader.py`; installing an R-file package was blocked in this environment). File: `celltype_meta_data_2025_01_08.RDS`, 2,344,381 bytes, SHA-256 `dc65ebe5…6bec` (GitHub tree `383b0a3`; Zenodo CC BY 4.0 archive). Output: `output/recurrence_reliability_care_metadata/care_label_extract_summary.json`; the 4-column extract stays local (`data/recurrence_reliability/raw/care/care_celltype_4col.tsv.gz`, git-ignored).

- **Schema: a tibble with exactly 4 columns, `CellID` (character), `Sample` (character), `ID` (character), `CellType` (factor).** The file has no survival, MGMT, state, score or other clinical column, so it is metadata-only. Nothing was discarded unread because nothing else was there.
- 429,305 rows (equal to the 10x cell count in Spitzer et al.); CellIDs unique; CellType labels: Astrocyte, Endothel, Excitatory neuron, Inhibitory neuron, Lymphocyte, Malignant, OPC, Oligodendrocyte, Other, Pericyte, TAM.
- `Sample` is the processing lab's internal code, not the Table S1 Tumor ID. The first run wrongly compared the two and listed false mismatches; the check was fixed and rerun. Prefix by source and lab: TKU (Tokyo, Lab1, 30 specimens), N and P (Pitié-Salpêtrière: Lab1 12, Lab2 16), SN (Seoul, Lab3, 20), SMMG (St. Michael's, Lab1, 14), NL (NorLux, Lab3, 13), Duke (Duke, Lab1, 4), COLUMBIA (Duke tissue processed by Lab2, 8), MDA (MD Anderson, Lab1, 4). `COLUMBIA` names a processing lab for Duke tissue; it is not a new tissue source and does not change the Wang overlap assessment.
- All 84 specimens of the 42 sample-level pairs are present. Pairs with a specimen below 50 malignant nuclei (9): P5, P23, P31, P32, P38, P58, P69 (first recurrence low) and P66, P72 (primary low). **33 pairs pass:** P2, P4, P6, P8, P9, P20, P21, P22, P24, P25, P26, P27, P28, P29, P30, P34, P36, P39, P40, P42, P45, P47, P48, P50, P52, P62, P63, P65, P67, P68, P70, P75, P78.
- Incidental exposure: per-specimen counts of author-labelled malignant nuclei (composition), the same kind of exposure as the Wang Phase B. It relates to the published CARE finding of a lower malignant fraction at recurrence. No MES or other program value was seen.
- Selection note: 7 of the 9 exclusions are recurrences with few malignant nuclei, so the analysed cohort over-represents recurrences that still contain many tumor cells. Report this as a limitation.
- Precision at 33 pairs (ceiling; may fall after the Wang filters): zero discordant patients give a two-sided upper bound of 0.106; the lower bound exceeds 5% from 5 discordant patients.
