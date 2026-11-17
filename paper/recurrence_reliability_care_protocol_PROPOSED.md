# CARE external replication: PROPOSED protocol (not locked)

Date: 2026-11-16. Status: **proposed; ready to lock on the user's go-ahead.** Metadata gates are settled (33 pairs after the blinded label read, `paper/recurrence_reliability_care_eligibility_gate_2026-11-16.md`). The barcode match to GEO matrix columns and the Wang-filter recount are QC steps in section 2 that run after locking and before scoring. No CARE expression, cell label, score, or outcome has been opened. This protocol does not change the locked Wang protocol, config, cohort, or result. If C6/C7 change a rule below, the change is made and dated before any CARE score is computed.

## 1. Question and type

Replication: in an independent cohort (CARE, GSE274546), how often do the two fixed scoring procedures of the Wang study give opposite directions for a patient's primary-to-first-recurrence change in the Neftel MES-like program?

This is a replication of the Wang endpoint with a refitted reference (section 4). It is not a test of recurrence biology, not a scorer ranking, and not pooled with Wang for the primary reading. Anything not listed in sections 6–7 is exploratory and must be labelled so.

## 2. Cohort rules (frozen from metadata)

1. Start: the 33 pairs that pass the sample-level rules (`care_eligibility_manifest.json`: untreated primary; first recurrence after radiotherapy + TMZ only; IDH-wildtype; Table S1 ↔ GEO identity) and have ≥ 50 author-labelled malignant nuclei per specimen (`care_label_extract_summary.json`): P2, P4, P6, P8, P9, P20, P21, P22, P24, P25, P26, P27, P28, P29, P30, P34, P36, P39, P40, P42, P45, P47, P48, P50, P52, P62, P63, P65, P67, P68, P70, P75, P78. Second-recurrence specimens are never used.
2. Specimen label join: the author label file `celltype_meta_data_2025_01_08.RDS` (SHA-256 `dc65ebe5…6bec`; 4 columns, verified metadata-only; Zenodo 10.5281/zenodo.14966015, CC BY 4.0) joined to GEO matrix columns by `CellID` = `<Sample>-<16-base barcode>`; the label `ID` (P#T#) selects the GEO file `GSM…_P#T#_umi_counts.RDS.gz`. If the GEO column names use another barcode form, only a GEM-suffix difference (`-1`) may be normalised; any other difference stops the study. A specimen passes if at least 95% of its label-file cells are found among its GEO matrix columns and every match is unique. Otherwise the pair is excluded.
3. Analysis nuclei: nuclei present in the author label file (author QC and confidence filter) that also pass the Wang filters (≥ 200 detected genes, mitochondrial fraction ≤ 2.5%). Matrix barcodes absent from the label file are outside the author-QC population: they are not labelled, not analysed, and do not trigger the Wang amendment-1 stop rule (that rule concerned a labelled-by-design population). Their number is reported.
4. Malignant nuclei: `CellType == "Malignant"`. No other label (for example cell states) is read.
5. Support: each specimen needs ≥ 50 malignant analysis nuclei; otherwise its pair is excluded.
6. Freeze the final cohort and all file hashes in a CARE QC report before scoring. Smart-seq2 data are not used.

## 3. Program and gene mapping

Neftel MESlike1 and MESlike2 (`data/external/neftel2019/IDHwt.GBM.MetaModules.tsv`, SHA-256 `046905b6…ce25`). Exact symbol match to the CARE features; genes with zero or several matches are dropped; a module below 80% mapped stops the study. MES-like per nucleus = max of the two module scores.

## 4. Procedures and reference (chosen on estimand grounds)

- **R, primary (refit):** Tirosh/Neftel control-gene score exactly as in the Wang protocol, with bins, control genes (100 per gene, 30 bins, seed 20261116) and gene centres re-estimated from the CARE final-cohort malignant pool. Estimand: discordance when each procedure is run as specified inside a new cohort, matching how the Wang estimate was made.
- **A:** UCell, max rank 1500 (no reference).
- **Transport (declared secondary, S-T):** R with the Wang-fitted reference (bins, control genes, centres), regenerated deterministically from the Wang inputs without computing any score and checked against the saved counts (4,027 and 4,041 control genes). Genes absent from CARE are dropped from module and control sets. Estimand: deployment of a fixed reference across cohorts.

## 5. Direction call

Same as Wang: delta = recurrence mean − primary mean over malignant nuclei; SE from 500 nucleus-bootstrap draws per specimen (seed 20261117); determinate if |delta| > 1.96·SE.

## 6. Primary endpoint and interpretation

Opposite-direction discordance between R (refit) and A over all patients in the final CARE cohort; exact two-sided 95% Clopper-Pearson interval; the Wang reading rule against the 5% tolerance (lower bound > 5%: unreliable beyond tolerance; upper bound < 5%: reliable within tolerance; otherwise inconclusive). Reported next to the Wang result (3/21; 3.0%–36.3%) without pooling. Precision at the 33-pair ceiling: zero events give an upper bound of 10.6%, so a "reliable within tolerance" reading is impossible; the lower bound exceeds 5% from 5 events.

The binomial interval assumes independent patient indicators, but the refit reference is shared. Section 7 S-D reports the size of this dependence; the primary interval is not adjusted.

## 7. Declared secondary analyses (descriptive)

- S-D dependence check: leave-one-patient-out reference (bins, controls, centres refit without patient i), directions recomputed, discordance reported.
- S-T transport (section 4).
- S-P protocol strata: discordance within Protocol 1 and Protocol 2 pairs (no patient mixes protocols).
- S-C calibration: 20 within-specimen split halves, determinate-call rate (expected about 5%).
- S-0 continuous: per-patient deltas and SEs for R and A; patient bootstrap of the mean (10,000 draws, seed 20261118).

No other scorer, program, threshold, or clinical variable is added. Survival, MGMT, surgical interval, and cell-state labels are not read.

## 8. Run protection

Separate config `configs/recurrence_reliability_care.json` and output directory; the Wang config and code are not edited (any shared helper is imported, not changed). Adapters (RDS reading, ID join) are tested on synthetic inputs in a scratch copy. Hashes of all CARE input files are frozen. The scoring run executes once; a bug found later needs a dated correction and disclosure.

## 9. Open implementation points before locking

1. Done: the blinded read showed the label file has 4 columns only (no survival or state columns); 33 pairs pass the author-level malignant floor.
2. Exact `CellID` string form versus GEO RDS column names.
3. A Python reader for the GEO `.RDS.gz` count objects: `src/recurrence_reliability/rds_reader.py` handles data frames; sparse `dgCMatrix` (S4) count objects need an extension, tested on synthetic files.
4. Final cohort size and precision recomputed from the frozen QC report.
