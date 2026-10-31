# Track A frozen protocol — Yu 2020 (GSE117891) Core vs Periphery, donor-held-out

Status: frozen before any expression file is downloaded. Base freeze commit: `f4e592d`. Amendment A1 below further specifies the confound diagnostics before expression data are opened. The current protocol, donor list, and pinned analysis-plan appendix are hashed in `SHA256SUMS`; its paths are relative to this folder. Verify them before every run. Any later change requires a dated, labeled amendment committed before expression data are opened.

## 1. Question (unchanged from the appendix, section 2.1)
Will a sparse, donor-balanced expression classifier improve donor-level Core-versus-Periphery prediction over ordinary logistic regression when every test donor is excluded from training?
With 6–8 eligible donors, the study will be **descriptive only**. It will not make an inferential superiority claim.

## 2. Data
- Expression: `GSE117891_all_6148.umi.count.matrix.tsv.gz` from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117891/suppl/ (about 39 MB, listed size). Its byte size and SHA-256 will be recorded on download.
- Cell → point map: `GSE117891_Sample_barcode_cell_information.xlsx`, sheet `Cell_information` (231,046 bytes, SHA-256 2669c32aa678d941e60e19dd786cb581c9212bb9e088967cb2bf2af52aa016da).
- Point → label map: Yu et al. 2020, *Natl Sci Rev*, doi 10.1093/nsr/nwaa099, Table S1, sheet `Point and Cell information` (15,125 bytes, SHA-256 0676607ff4c57348fc043ab4bace1c3985617d0d90507742344ce19642cb772a). License CC BY 4.0.
- Raw files will stay local under the git-ignored `data/external/`. Nothing will be uploaded to a cloud service.

## 3. Labels and eligibility (frozen in `donors.tsv`)
- Core = `Tumoral`. Periphery = `Peritumoral`. `Normal` points will be excluded.
- Label provenance: MRI-marked surgical sampling points (paper text, Fig. 1). Labels are not derived from expression.
- Point numbers will never be used to infer a label. Only the Table S1 mapping will be used.
- GEO patient `S<n>` will map to Table S1 patient `GS<n>`. The cell table's `point` column (e.g. `S13P4`) will map to Table S1 `GS13P4`. A cell whose point is absent from Table S1 will be excluded and counted.
- **Primary cohort (strict, 6 donors):** GS1, GS2, GS3, GS5, GS6, GS13. GBM, IDH1 R132 negative.
- **Sensitivity cohort (lenient, 8 donors):** strict plus GS11 and GS12 (GBM, IDH not reported).
- Excluded: GS15 (no histology; GEO lists one metastasis), GS4 (gliosarcoma), GS7 (IDH-mutant), GS8, GS9, GS14 (grade II–III or no peritumoral point).
- The run will stop if any primary donor has fewer than 20 QC-passed cells in either class after filtering.

## 4. Preprocessing (fitted on training donors only)
- Cell QC: keep cells with ≥1,000 UMIs and ≥500 detected genes. These thresholds do not use labels.
- Normalize each cell to 10,000 counts, then log1p.
- Gene filter: genes detected in ≥5% of training cells. Select the 2,000 most variable genes on training cells.
- Scale with the training-cell means and SDs. Apply the same transform to test cells.

## 5. Models
- **Candidate:** the appendix section 2.2 model (sigmoid gene gate, L1 on the gate, ridge on W, smooth worst-donor loss with τ = 0.1, training-only 15-NN graph smoothness, test-time neighbor blend). Its hyperparameter grid, optimizer, epoch cap and early-stopping rule will be as in the appendix.
- **Primary comparator:** ordinary L2 multinomial logistic regression, with C in {0.001, 0.01, 0.1, 1, 10}.
- **Secondary baselines (appendix 2.3):** 15-NN vote, linear SVM, gradient-boosted trees. The three ablations will also run.
- **Leakage diagnostics (not competitors):**
  1. Cell-type-only model: use only GBmap Core's `annotation_level_2` from the CELLxGENE collection `999f2a15-3d7e-440b-96ae-2c806799c08c`. Before the count-matrix download, save the metadata-only `obs` snapshot and record its release/version and SHA-256 in the run manifest; do not download the 8.1 GB expression object for this purpose. Join by exact cell barcode to the frozen GEO cell table; do not normalize, fuzzy-match, or infer barcodes. Fit an L2 logistic model on one-hot cell-type labels with `C=1` and balanced class weights, on outer-training donors only. Score it on cells with a unique cell-type match, and score each expression model on that same mapped subset for this comparison. Report unmapped and ambiguous cells by donor and region. This expression-derived annotation is a composition diagnostic only; it is never a feature or target for the expression models.
  2. Point-size model: use only `log1p` of the number of cells in the cell's point, with an L2 logistic model (`C=1`, balanced class weights) fitted on outer-training donors. This is a diagnostic for sampling-depth leakage, not a biological predictor.
- Every model will get the same outer folds and the same fixed inner tuning budget (full grid, no extra search).

## 6. Splits
- Outer loop: leave one donor out (6 folds for the primary cohort, 8 for sensitivity).
- Inner loop: leave one training donor out. Each inner training split will keep ≥4 donors (primary: 4).
- All cells of a donor will stay in one fold. No test-to-test graph edges.
- Seed: 20261031 for all random steps.

## 7. Endpoints
- **Primary:** class-balanced AURC (appendix 2.3), averaged equally over donors. Lower is better.
- **Secondary:** full-coverage balanced accuracy per donor; Brier score.
- Report every per-donor value for every model. Report the mean and the min–max range over donors. No p-values for superiority.

## 8. Required confound checks
- Balanced accuracy of the cell-type-only diagnostic per donor.
- Expression-model balanced accuracy **within** each `annotation_level_2` stratum that has ≥10 cells of both classes in ≥4 donors.
- Point (library) is nested in region inside each donor. This will be stated as a limit. It cannot be removed.

## 9. Pre-stated reading of results
- "Candidate better than logistic regression (descriptive)" only if the candidate's donor-mean AURC is lower **and** it is lower in at least 5 of 6 primary donors.
- If the cell-type-only diagnostic reaches balanced accuracy within 0.05 of the best expression model on the same mapped cells, report that cell-type composition alone gives comparable prediction. Do not say it proves the expression result is explained by composition.
- If no stratum passes the size rule in section 8, the within-cell-type check will be reported as not possible, and no region-biology claim will be made.
- The lenient cohort cannot change the primary conclusion.
- A negative or null result will be reported unchanged. Settings will not be retuned after scoring.

## 10. Failure handling
- Hash mismatch in this folder or in an input file: stop.
- Download size different from the recorded size: download again before any analysis.
- Model fails to converge within the epoch cap: keep the last inner-best weights and report it.
- CPU runtime above 6 hours in a 1-fold pilot: stop and record an amendment before reducing the grid.

## 11. Outputs
- `output/trackA_yu2020/result.json`: donors, cell counts per class, per-donor metrics for every model and diagnostic, the seed, the protocol commit and package versions.
- Before reporting, the file will be checked to exist, to be non-empty and to list 6 primary donors.

## 12. Limits known now
- 6 donors: wide uncertainty; descriptive only.
- STRT-seq plate data, 35–866 cells per donor-class.
- Core is mostly neoplastic and Periphery mostly non-neoplastic (GBmap annotation), so composition is the main expected explanation.
- The labels come from one surgical team's MRI-guided sampling. Core/periphery boundaries are not histology-confirmed per cell.
- The method is not claimed as novel, and no new-biology claim will be made from a favorable score without a literature comparison.

## Amendment A1 — composition and sampling-depth diagnostics (2026-10-31, pre-expression)

This amendment makes these pre-expression changes:

- Pins the cell-type diagnostic to GBmap Core `annotation_level_2` from CELLxGENE collection `999f2a15-3d7e-440b-96ae-2c806799c08c`; requires an exact barcode join and a hashed metadata-only snapshot before the count-matrix download.
- Fixes the cell-type-only diagnostic as L2 logistic regression (`C=1`, balanced class weights), fitted on outer-training donors and scored on uniquely mapped cells. Expression models are also scored on that same mapped subset for this comparison.
- Fixes the sampling-depth diagnostic as L2 logistic regression (`C=1`, balanced class weights) on `log1p` point cell count.
- Changes the interpretation of a similar cell-type-only score to “comparable prediction from measured cell-type composition”; it does not prove the expression result is explained by composition.

It does not change donor eligibility, labels, the expression candidate or comparator, splits, endpoints, or the descriptive-only rule. No expression count matrix has been downloaded or scored.
