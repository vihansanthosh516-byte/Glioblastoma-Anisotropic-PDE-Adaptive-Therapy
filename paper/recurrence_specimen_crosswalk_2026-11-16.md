# GSE174554 Supplementary Table 1 crosswalk — 2026-11-16

## Result
Table 1 supports **33 metadata-defined pairs with one primary and one recurrent snRNA-seq specimen both labeled IDH wildtype**. Of those, **27 pairs** have both SF specimen IDs matched to exactly one GEO GSM title and each GEO progression label agrees with Table 1. These 27 are metadata candidates, not a final eligibility count: Table 1 pairing labels still need source confirmation, and the sheet does not identify which recurrent specimen is the first recurrence for the extra multi-recurrence group.

## Counts from selected Table 1 fields

- Table 1 contains 111 sample rows; 86 are marked `snRNA-seq = Y`.
- Among the 86 snRNA rows: 83 are `IDH wildtype` (38 primary, 45 recurrent) and 3 are `IDH mutant` (2 primary, 1 recurrent); all list diagnosis `Glioblastoma`.
- Across all IDH labels, 35 Table 1 Pair# groups have exactly one primary and one recurrent sample. Of these, 33 are IDH wildtype at both sample rows, one has discordant IDH labels, and one is IDH mutant at both rows.
- Table 1 Pair# 36 has one primary and two recurrent snRNA rows; Pair# 37 is recurrent-only. The 12 snRNA rows with Pair# `NA` are not paired by this field.
- Pair# 16 has a primary row marked IDH mutant and a recurrent row marked IDH wildtype. Pair# 24 has both rows marked IDH mutant. The 33 fully IDH-wildtype pairs exclude both groups.
- Exact SF-ID matching links 78 of the 86 Table 1 snRNA rows to one GEO sample title; 8 have no exact title-ID match. Two exact matches have conflicting progression labels: SF9798 and SF9494.
- GEO Pair# values are largely renumbered relative to Table 1. Comparing SF-ID membership sets rather than label numbers yields 36 identical group sets among matched records; two Table 1 group sets are unmatched and one GEO group set is unmatched. Numeric label differences alone are not treated as biological pairing conflicts.
- Table 1 assay-presence marks: 86 snRNA-seq, 11 snATAC-seq, 35 RNA-exome, 10 spatial RNA, 6 spatial protein. These are specimen-row marks; modalities may overlap and they are not independent sample counts.

## Interpretation of 27 versus 36/40

The source paper reports 36 analyzed primary–first-recurrence patients; official GEO series summary advertises 40 patient-matched primary/recurrent pairs. Table 1 supplies 35 simple primary/recurrent Pair# groups plus Pair# 36 with two recurrent specimens. Table 1 does not include recurrence order for the two Pair# 36 samples. A literal SF-ID + progression crosswalk leaves 27 of the 33 both-IDH-wildtype simple groups with a unique GEO title match and no stage mismatch. Therefore 27 is a conservative cross-source metadata candidate set, not proof that only 27 patients are usable. The 36 and 40 descriptions remain different scopes; this table does not establish the first-recurrence selection needed to certify a final count.

## Provenance and scope

Read only Table 1 columns `ID`, `Stage`, `Pair#`, `Diagnosis`, `IDH`, and assay-presence columns (`snRNA-seq`, `snATAC-seq`, `RNA-exome`, `Spatial RNA`, `Spatial protein`). Did not inspect `Age`, `Sex`, `Tumor site`, `Elapsed time to recurrence`, `Overall survival`, or any other supplementary tables. No expression matrices, outcomes, or numeric biomarkers were opened.

- Source: Wang et al., [Nature Cancer paper](https://www.nature.com/articles/s43018-022-00475-x), [PMC full text and supplementary workbook](https://pmc.ncbi.nlm.nih.gov/articles/PMC9767870/).
- Workbook: `43018_2022_475_MOESM2_ESM.xlsx`, sheet `Table 1`; SHA-256 `58e42239d4e105cfbbda6441bd640afbd64e0c55f87e07b1e0aaab0d1eff1c42`.
- GEO linkage source: official [GSE174554](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174554) sample manifest preserved in the project metadata folder.
- Machine-readable selected-field mapping and summary: [recurrence_table1_crosswalk.json](../output/paper_strengthening_metadata/recurrence_table1_crosswalk.json).

**Gate:** conditional. Before treating any count as eligible, confirm the Table 1 Pair# interpretation and Pair# 36 recurrence selection against the authors’ analysis/sample manifest; resolve the eight missing exact title-ID joins and two stage-label disagreements.
