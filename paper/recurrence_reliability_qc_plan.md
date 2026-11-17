# Recurrence reliability: QC-only plan (Phase C), frozen before matrices are opened

Date: 2026-11-16. Frozen in a commit before any matrix or feature file was downloaded. QC computes technical quantities only. It does not compute, plot, rank, or test any gene-program score or primary-versus-recurrence difference. Code: `src/recurrence_reliability/qc.py`; output: `output/recurrence_reliability/qc_report.json`.

## Files

Download only, from the per-GSM GEO route (`https://ftp.ncbi.nlm.nih.gov/geo/samples/<GSMnnn>/<GSM>/suppl/`):

- `matrix.mtx.gz` and `features.tsv.gz` of the primary library of the 46 eligible specimens;
- the same two files for the three independent-capture libraries (SF12243, SF10433, SF4209 batch2), for the library-aware secondary analysis only.

98 files, 615,079,991 listed compressed bytes. Same-capture batch2 libraries are not downloaded. Not downloaded: `GSE174554_RAW.tar`, raw reads, spatial files. Files stay in `data/recurrence_reliability/raw/gsm/` (git-ignored). SHA-256 of every file goes into the QC report; any later hash change stops scoring.

## Checks (all technical)

1. Matrix Market header is `coordinate integer` (or `real` with every value integral); dimensions equal (number of features, number of barcodes); every stored value is a positive integer; no duplicate (row, column) entries.
2. Features file: three columns; the third is `Gene Expression` for every row; Ensembl IDs unique; record the number of duplicated symbols. The feature list must be identical (same hash of the ID column) in every library; otherwise record the difference and map by Ensembl ID.
3. Barcodes unique within each library; cell IDs are `<library_id>:<barcode>`, so libraries never collide.
4. Author annotation re-joined as in Phase B; join rate recorded.
5. Fixed technical nucleus filters, taken from Wang et al. Methods (not tuned): keep nuclei with at least 200 detected genes and mitochondrial (`MT-` symbol) fraction at most 2.5%. The GEO matrices are expected to be post-filter already; QC records how many nuclei each filter would remove.
6. After filters, each specimen must keep at least 50 author-labelled tumor nuclei (Phase B floor). Otherwise its pair is excluded.
7. Gene-ID mapping for the frozen program: count how many Neftel MESlike1 and MESlike2 genes (`data/external/neftel2019/IDHwt.GBM.MetaModules.tsv`, SHA-256 `046905b6...ce25`) and Verhaak mesenchymal genes (`data/external/neftel2019/TCGA.Verhaak.GBM.txt`, secondary construct) have exactly one feature with that symbol. Rule: genes with zero or several matching features are dropped; a module with fewer than 80% of its genes mapped stops the study.
8. Library-depth summary per specimen (median UMIs and genes per nucleus, tumor and all nuclei). These are technical depth values for the matched-depth plan, not outcomes.

## What QC must not do

No module score, no mean expression of program genes, no comparison of primary and recurrent specimens on any expression quantity, no survival or clinical field. If a check fails because of a code bug, fix the code and rerun QC; record the fix.

## After QC

Freeze the final eligible cohort in `qc_report.json`, then write and commit the protocol (Phase D) before any scoring.

## QC result (2026-11-16, appended after the run; plan text above unchanged)

Run: `py -3 src/recurrence_reliability/qc.py` (3 min 43 s). Report: `output/recurrence_reliability/qc_report.json`. No program score was computed.

- All 98 downloaded files passed `gzip -t`; sizes match the GEO listing; SHA-256 of all 147 used files is in the report.
- Every matrix header is `coordinate integer general`; dimensions match; all values are positive integers; no duplicate entries; barcodes unique.
- **Deviation from check 2 (file format, not outcome-driven):** the GEO feature files have one column (gene symbol) only, with no Ensembl IDs or feature-type column. The symbol is the only identifier. The feature list (33,694 symbols, 0 duplicated) is identical in every library. Mapping is by symbol, as planned for the modules. Code change: `tenx_io.load_library` accepts the symbol-only format; `qc.py` no longer requires Ensembl uniqueness or a `Gene Expression` type.
- **The matrices are not uniformly post-filter.** The fixed Wang filters remove 3,432 nuclei for < 200 genes and 2,904 for mitochondrial fraction > 2.5% across the 49 libraries. Large removals: SF6118 (565 of 568 nuclei < 200 genes), SF9715 (1,616 of 1,834), SF10432 (541 of 994), SF1199 (512 of 634).
- Malignant support after filters: **SF6118 (0) and SF9715 (4) fall below 50, so pairs 19 and 20 are excluded.** GEO also has separate `SF6118v2` and `SF9715v2` records (independent captures by barcode overlap). Using them would change the frozen primary-library rule after QC, so they are not used.
- Module mapping: MESlike1 50/50, MESlike2 49/50 (ERO1L absent; current symbol ERO1A not substituted, per the exact-symbol rule), Verhaak Mesenchymal 100%. No stop.
- **Final cohort: 21 patients** (Table 1 pairs 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 21, 26, 29, 30, 31, 32, 33, 34). Smallest malignant counts after filters: SF6098 55, SF12115 59, SF7025 77.
- Independent-capture libraries for S3 after filters: SF12243 batch2 1,301 malignant nuclei, SF10433 batch2 132, SF4209 batch2 4 (too few to be informative; reported as is).
