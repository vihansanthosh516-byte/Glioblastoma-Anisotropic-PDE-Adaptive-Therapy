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
