# Track A frozen protocol — Yu 2020 (GSE117891) Core vs Periphery, donor-held-out

Status: base protocol frozen before expression download. Base freeze commit: `f4e592d`. Amendment A1 below further specifies the confound diagnostics before expression data are opened. The count matrix was downloaded on 2026-10-31 and its SHA-256 matches the source record; no count values have been read or scored. Amendment A2 below freezes implementation details before the first expression read. The protocol, donor list, pinned analysis-plan appendix, runner, and four source inputs are hashed in `SHA256SUMS`; its paths are relative to this folder. Verify every entry before each run. Any later change requires a dated, labeled amendment committed before expression data are opened.

## 1. Question (unchanged from the appendix, section 2.1)
Will a sparse, donor-balanced expression classifier improve donor-level Core-versus-Periphery prediction over ordinary logistic regression when every test donor is excluded from training?
With 6–8 eligible donors, the study will be **descriptive only**. It will not make an inferential superiority claim.

## 2. Data
- Expression: `GSE117891_all_6148.umi.count.matrix.tsv.gz` from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117891/suppl/. It was downloaded on 2026-10-31 (40,989,643 bytes); the SHA-256 is recorded in `SHA256SUMS`. The file passed gzip integrity checking.
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

It does not change donor eligibility, labels, the expression candidate or comparator, splits, endpoints, or the descriptive-only rule. At the A1 freeze time, no expression count matrix had been downloaded or scored.

## Amendment A2 — frozen implementation choices (2026-10-31, before expression values are read)

This amendment resolves implementation details left open by the appendix and Amendment A1. It does not change the target, donor eligibility, primary endpoint, candidate model, required comparisons, or descriptive-only interpretation.

### A. Input and cell handling

- Before loading expression values, verify every entry in this folder's `SHA256SUMS`, including `PROTOCOL.md`, `donors.tsv`, the pinned appendix, the runner, the count matrix, the GEO cell table, Table S1, and the GBmap metadata snapshot. Any missing file or hash mismatch stops the run.
- Join matrix columns to the GEO cell table by exact cell barcode only. Join GEO points to Table S1 by exact patient and point identifiers only. Map `Tumoral` to Core and `Peritumoral` to Periphery. Exclude `Normal`, other diagnoses, and cells without a unique exact join; count each exclusion reason. Duplicate cell barcodes, duplicate gene identifiers, conflicting point labels, or a mismatch between Table S1 and `donors.tsv` stops the run.
- Validate raw counts as finite, nonnegative integers. Apply the frozen QC rule on raw counts: at least 1,000 UMIs and at least 500 detected genes. Select the primary and sensitivity cohorts only from `donors.tsv`. Require at least 20 QC-passed cells in each class for every primary donor. If a lenient donor misses this threshold, omit the entire eight-donor sensitivity analysis and record why; do not alter the six-donor primary cohort.
- In each inner and outer training split, filter genes detected in at least `ceil(0.05 × number of training cells)` cells, rank by variance on training cells only, and break variance ties by gene identifier. Keep the 2,000 highest-variance genes; stop if fewer than 2,000 remain. Normalize each cell to 10,000 counts, apply `log1p`, then fit feature means and standard deviations on training cells only. Replace a zero standard deviation with 1. Apply the same feature order and transform to validation or test cells.

### B. Nested selection and candidate refit

- Keep donor order from `donors.tsv`; use leave-one-donor-out outer folds and leave-one-training-donor-out inner folds. Use seed `20261031`. Assign each cell a deterministic tie rank by sorting SHA-256 digests of `20261031`, a null byte, and its exact cell ID. Reuse the same ranks for every model and fold, so ties do not favor one model or change across cohorts.
- In each inner fit, early-stop the candidate and ablations on the lowest mean donor validation cross-entropy of the unblended classifier output. Use the appendix limit of 500 epochs and 25 epochs without improvement; on an exact tie, keep the earliest epoch. Select hyperparameters by mean inner-donor class-balanced AURC. Among inner-fold AURC values within `1e-12` of the minimum, select the first setting in the fixed ascending grid order for every tuned model.
- Select candidate regularization and blend settings jointly from the 75 regularization combinations and 20 `(T, gamma)` combinations, using only inner-fold predictions. Blend uses probabilities from training cells only. For each outer refit, train the selected setting on all outer-training donors for the median best epoch of that setting across inner folds. Do not use the outer donor to select an epoch or any setting.
- Tune each ablation independently on its remaining grid: ablation 1 fixes `lambda_1=0`; ablation 2 fixes `lambda_2=0` and `gamma=0`; ablation 3 uses mean donor cross-entropy instead of the smooth maximum. Keep `lambda_W` and any remaining blend settings on their appendix grids.

### C. Comparator details

- Logistic regression uses the appendix C grid, L2 penalty, `lbfgs`, no class weighting, and a fixed 1,000-iteration cap.
- The 15-nearest-neighbor vote uses cosine distance, uniform voting, and the shared deterministic tie order.
- Linear SVM uses `LinearSVC` with the appendix C grid, squared-hinge loss, no class weighting, and a fixed 10,000-iteration cap. For probability metrics, fit `LogisticRegression(C=1, solver='lbfgs', max_iter=1,000, class_weight=None)` to donor-held-out margins generated only within the current training donors. Give each calibration donor equal total weight. In inner folds, calibration is nested inside the inner-training donors. Use calibrated probabilities for Brier score and confidence ranking. A calibration fit failure stops the run; do not substitute a fixed sigmoid.
- Gradient boosting uses `GradientBoostingClassifier` with `log_loss`, 100 estimators, depth-2 trees, `min_samples_leaf=10`, `subsample=1`, and learning rate in `{0.05, 0.1, 0.2}`. Disable internal early stopping and use seed `20261031`. No other tree settings are tuned.
- The cell-type and point-size diagnostics remain fixed L2 logistic models with `C=1` and balanced class weights. Cell type joins remain exact and unique; score cell-type comparisons on the same uniquely mapped cells for all expression models. Point size is `log1p` of QC-passed cells in that source point. These diagnostics are not competitors.
- For each eligible `annotation_level_2` stratum (at least 10 Core and 10 Periphery cells in at least four primary donors), report each expression model's held-out balanced accuracy, class-balanced AURC, Brier score, and class counts by donor. Also report donor-macro balanced accuracy across donors meeting that stratum rule. If no stratum meets it, report that the check was not possible.

### D. Endpoint and runtime

- For each held-out donor, sort cells by `max(p_core, p_periphery)` using the shared tie order. Compute class-specific error among retained cells. If a retained set has no cells from one class, assign error 1 to that class. Compute error as the mean of Core and Periphery errors.
- Use every integer retained-cell count from `ceil(N/2)` through `N`. Record actual coverage `k/N`, and integrate the risk curve with the trapezoid rule over those actual coverages. Lower AURC is better. Balanced accuracy uses the 0.5 probability threshold. Brier score is the per-cell binary squared-probability error, averaged within donor. Report all donor values and the donor mean and range; do not calculate p-values or claim inferential superiority.
- The timing pilot runs the first outer fold of the eight-donor sensitivity cohort when all eight donors pass the cell-count gate; otherwise it runs the first outer fold of the six-donor primary cohort. It includes inner tuning and final fits on outer-training donors, but never predicts, scores, or saves results for the held-out donor. It records timings only. Extrapolate the measured fold time to 14 outer folds when sensitivity is eligible, or 6 when sensitivity is omitted, and add the measured input-load time. If the estimate exceeds six hours on CPU, stop before full scoring and write a new amendment before changing any grid. If it passes, run the six-donor primary cohort first, then the eligible eight-donor sensitivity cohort.
- Write the complete `result.json` only after every required fold finishes. Include all implementation choices, input and code hashes, package versions, cell-flow counts, per-donor metrics, runtime, and limitations. Refuse to overwrite an existing result. The pilot timing file contains no scores.

Commit this amendment, the runner, and their checksum manifest before the first expression-value read. Keep the downloaded matrix local; never upload it.
