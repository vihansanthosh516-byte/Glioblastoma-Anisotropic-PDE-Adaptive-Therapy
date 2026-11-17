# Track A Artzi region benchmark: pre-analysis plan

Status: frozen before opening `large_sparse_matrices2.h5` or calculating any expression score.  
Analysis date: 2026-11-16.  
Review: [full scoped prior-art inventory](trackA_prior_art_full_2026-11-16.md).

## Question and interpretation

Among paired primary GBM samples from one patient, can expression distinguish the pathology-labeled primary core (P-CT) from the primary infiltrative edge (P-IT) in a patient not used to train that fold?

The target is **whole-tissue regional ecosystem classification**. It is not a tumor-cell-intrinsic claim, a clinical response prediction, or a treatment recommendation. The public dataset includes pathology-guided region labels confirmed by a neuropathologist; the exact same-resection linkage is inferred from the study design because the public metadata has no specimen-level surgery ID. The user-supplied article and the [metadata-only audit](trackA_artzi_zenodo_metadata_2026-11-16.md) record this limitation.

Artzi was selected after earlier negative Track A results. Therefore the evaluation is exploratory, even though no Artzi expression matrix or score has been opened. It can estimate within-study generalization to held-out patients, not independent external validation. Yu and HGCC remain spent and are not used as test data or for further tuning. The prior negative A5/TMZ results and the failed previous simulation remain in the [failure audit](trackA_failure_audit_2026-11-13.md), [TMZ records](trackA_tmz_response_119of121_amendment_2026-11-08.md), and [simulation result](trackA_margin_simulation_result_2026-11-13.md).

## Data and eligibility

Use only the public Artzi Xenium release named in the metadata audit: `Final_meta.csv` for cell and sample metadata and `large_sparse_matrices2.h5` for the count matrix. Its license and metadata-only file checks are recorded in the [Artzi audit](trackA_artzi_zenodo_metadata_2026-11-16.md). No expression matrix, gene-level count, or score has been accessed for this plan.

Use primary P-CT and primary P-IT only. Pair by the study’s `Patient` code and `Setting`. Aggregate raw counts over all cells and all fields of view for each patient-setting into one regional pseudobulk. Include all cells; do not use study-derived malignancy, subtype, state, ecosystem, or pooled cluster calls to select or weight cells. This prevents those pooled expression-derived annotations from serving as independent test inputs. A matrix row must be joined to metadata through an exact released cell identifier. If the matrix cannot be joined unambiguously, stop. If a paired regional pseudobulk is empty or non-finite, record the metadata and stop before model scoring; do not drop the patient after seeing predictions.

Keep only patient codes that have one or more primary P-CT and one or more primary P-IT records and whose counts join unambiguously. Aggregate replicate records by patient and setting. Do not use recurrent samples, treatment history, inferred cell states, or outcome labels. The pre-access audit found ten shared primary patient codes; the complete-matrix eligibility count remains unknown until the post-lock matrix audit. [Source: metadata-only audit](trackA_artzi_zenodo_metadata_2026-11-16.md).

Minimum eligible set: nine complete patient pairs. This minimum is fixed because the planned upper bound of twenty-five comparator arms requires at least nine discordant pairs for the best possible Holm-adjusted exact sign-test result to be no greater than 0.05: `25 / 2^9 = 0.0488`; with eight, the best possible value is `25 / 2^8 = 0.0977`. The calculation is implemented by `minimum_possible_adjusted_p()` in [the shared evaluator](../src/trackA_margin_evaluation.py). If fewer than nine pairs are eligible, report a data-gate failure and do not calculate model scores.

## Candidate method: maximin patient-paired rank-pair vote (MP-RPV)

For each training patient `i`, specimen `r`, and gene pair `(g,h)`, calculate `v[i,r,g,h] = sign(x[i,r,g] - x[i,r,h])`. The sign is `-1`, `0`, or `+1`; ties are neutral. Define the within-patient region contrast as `d[i,g,h] = v[i,P-IT,g,h] - v[i,P-CT,g,h]`.

For each pair, choose orientation `s` from `{-1,+1}` to maximize the worst-patient signed contrast `m[g,h] = max_s min_i (s*d[i,g,h])`. Break equal `m` by larger mean oriented contrast, then lexicographic gene-symbol order. Keep the first five pairs with `m > 0`; use equal vote weights and do not allow a gene pair to be duplicated. If fewer than five pairs have positive `m`, use only those pairs. If none have positive `m`, give both specimens the same score and therefore predict a tie. No fallback selector is allowed.

For each held-out specimen, its score is the arithmetic mean of the oriented votes from the selected training pairs. A patient-level prediction is concordant when the held-out P-IT score is greater than the P-CT score. This maximin objective is a narrow method candidate, not an established unique method. TSP, k-TSP, AUCTSP, pGLMM pair selection, and a GBM rank-pair ensemble are close prior art; the [prior-art review](trackA_prior_art_full_2026-11-16.md) explicitly leaves novelty unresolved. The previously proposed weakest-state selector is not substituted for MP-RPV; it failed the prespecified simulation and is reported separately in the [simulation result](trackA_margin_simulation_result_2026-11-13.md).

## Frozen patient split and preprocessing

Use leave-one-patient-out outer folds. In each fold, withhold that patient’s P-CT and P-IT pseudobulk together; use all other eligible patients for training. Each patient receives exactly one out-of-fold score per region. Never split cells, fields of view, or specimens from the same patient across training and test.

For non-rank methods, transform each pseudobulk to `log1p(count / total_count * 10,000)`. This is a specimen-local transform; do not estimate normalization, centering, scaling, gene filtering, or feature weights from the full held-out cohort. Any feature filtering or model fitting is performed inside each outer training fold. Do not impute missing genes from held-out samples. Gene-symbol mapping and signature coverage use the frozen assay feature list, after the plan and comparator definitions are committed.

Use the same eligible patients and the same region labels for every arm. Methods may use their predeclared feature representation, but every method must produce one scalar score for each held-out region and be evaluated by the same patient-level paired concordance. Fit no method using test-patient cell annotations, pathology labels beyond the paired outcome used for final scoring, or test-patient expression distributions. No batch correction, cross-sample normalization, or post hoc orientation is allowed.

## Comparator arms

The candidate must be compared against every portable predictor below. This is an upper bound of twenty-four method/signature arms plus an always-tie reference, for twenty-five comparator arms. Methods that use held-out labels to estimate effects, need unavailable modalities, or have no fixed score or prediction rule are documented as prior art but are not falsely described as out-of-sample classifiers; see the [prior-art review](trackA_prior_art_full_2026-11-16.md).

| ID | Comparator | Frozen application to this target |
|---|---|---|
| C01 | TSP | Select one gene pair by the largest training-only class ordering-proportion difference; use its source-defined vote and training-only direction. |
| C02 | k-TSP | Greedily select five disjoint pairs by the published ordering-proportion procedure; equal-weight pair vote. |
| C03 | AUCTSP | Select the single pair with maximum training-only AUC of its within-specimen ordering indicator. |
| C04 | Top-scoring N/triplet ordering | Apply the source method with the gene prefilter and maximum feature count frozen from the source implementation; any computational cap must be set before matrix access and used inside training folds. |
| C05 | Rank-pair tree ensemble | Reproduce the public ranktreeEnsemble method at a frozen revision, including training-only pair filtering and fixed seed. |
| C06 | GBM rank-pair/Robencla | Reproduce the source procedure from its public code and frozen revision; if its required inputs cannot be represented without held-out data, record it as structurally ineligible before scoring. |
| C07 | Penalized pair-order GLMM | Reproduce the Rashid et al. model using patient-grouped training folds; prediction for a held-out patient uses population-level fixed effects only. Freeze code, optimizer and failure rule before the matrix is opened; if the model cannot fit under its published specification, record a pre-score implementation failure, not a result-based exclusion. |
| C08 | Conditional logistic regression | Fit a patient-stratified conditional logistic model on training-patient core/margin pairs; compare fixed-effect linear predictors for the held-out pair. |
| C09 | Stability-selected conditional logistic regression | Resample training patients only, apply the published stable-CLR procedure with a frozen resample count and selection-frequency threshold, then score with the resulting fixed-effect model. |
| C10 | Paired centroid | Fit the average P-IT minus P-CT log-expression vector on training patients and score each held-out sample by its dot product with that vector. |
| C11 | Ridge logistic regression | Fit on the two regional samples per training patient; fixed `C=1`, balanced classes, and training-only feature standardization. |
| C12 | Linear SVM | Fit on the two regional samples per training patient; fixed `C=1`, balanced classes, and training-only feature standardization. |
| C13 | Random forest | Fixed 500 trees, square-root feature subsampling, balanced class weights, and fixed random seed; no test-driven tuning. |
| C14 | Paired differential-expression score | Rank genes by absolute mean within-patient log-expression difference in the training fold; score the top 50 with equal signed weights from the training direction. |
| C15 | Ivy GAP leading-edge versus cellular-tumor score | Apply the exact published gene list and direction as a fixed signed score. |
| C16 | Darmanis infiltration score | Apply the exact published infiltrating-tumor gene list and direction as a fixed signed score. |
| C17 | Kim regional score | Apply the exact published core/edge gene list and direction as a fixed signed score. |
| C18 | 5-ALA invasive-cell score | Apply the exact published invasive-cell gene list and direction as a fixed signed score. |
| C19 | Greenwald outer-layer score | Apply the exact published layer/state gene list and direction as a fixed signed score. |
| C20 | Harwood region modules | Apply the exact published core-versus-infiltrated module gene lists and directions as one equal-weight signed score. |
| C21 | Pai glioma-infiltration score | Apply the exact published infiltration gene list and direction as a fixed signed score. |
| C22 | GBM matrisome region score | Apply the exact published region-associated gene list and direction as a fixed signed score. |
| C23 | Petry iron/invasion score | Apply the exact published regional gene list and direction as a fixed signed score. |
| C24 | Sonpatki spatial cell/niche score | Apply the exact published regional gene list and direction as a fixed signed score only if its source defines a region-directed score; otherwise record it as not a direct region predictor. |
| C25 | Always-tie reference | Give both specimens the same score; its paired concordance is one-half for every patient. |

For C15–C24, freeze source supplement/table, exact genes, sign, and score code before opening the test matrix. Require at least 70% of a source signature’s genes in the assay. A signature below that threshold is structurally ineligible and is not replaced by another signature. If a source does not publish an unambiguous region direction or fixed gene set, record it as a regional biology precedent but not as an executable predictor. Potential overlap between an external signature’s source cohort and Artzi is unresolved for some studies; these signature comparisons and the full analysis remain exploratory.

Before matrix access, create and commit a machine-readable manifest with exact comparator gene lists, source file/table, direction, code revision/hash, package versions, and eligibility rule. Any direct method that is not defined by that manifest blocks expression access. This requirement does not permit choosing signatures based on Artzi expression, cell labels, predictions, or outcomes.

## Primary score, inference, and pass rule

For arm `a` and patient `i`, assign paired concordance `c[i,a] = 1` if the P-IT score exceeds P-CT, `0` if it is lower, and `0.5` if tied. The primary score is the unweighted mean `C[a] = mean_i c[i,a]`. This scores patients, not cells or fields of view.

For every comparator arm `j`, calculate patient differences `q[i,j] = c[i,MP-RPV] - c[i,j]`. Count wins where `q > 0`, losses where `q < 0`, and ties where `q = 0`. Use a one-sided exact paired sign test of whether MP-RPV wins more often than it loses, excluding only paired ties. Adjust all raw p-values with Holm’s procedure across every comparator arm eligible under the rules above. Report for every arm the mean concordance, per-patient concordance, score differences, wins, losses, ties, raw p-value, Holm-adjusted p-value, and any eligibility or implementation failure.

**Pass only if both conditions hold:**

1. MP-RPV has strictly higher mean patient-level concordance than every eligible comparator.
2. The one-sided exact paired sign test rejects against every eligible comparator after Holm correction at `alpha = 0.05`.

No near-pass, trend, selected subgroup, alternate metric, or post hoc threshold will count as a pass. If the candidate fails, report the negative result and its per-patient comparisons. If an arm or patient is missing, do not silently remove it after score inspection. The conservative nine-patient floor and upper-bound comparator count are checked before scoring using [the shared evaluator](../src/trackA_margin_evaluation.py).

## Labels for reporting

Label the benchmark **exploratory** because Artzi was selected after prior Track A results and is evaluated using within-study held-out patients. Label the previous synthetic method check **simulation** and the earlier project-derived regional/TMZ results **post hoc** where they are summarized. Distinguish those categories in every result table and figure. No clinical response or clinical utility statement is allowed from this endpoint.

## Seal state

This plan fixes the question, data eligibility, split, candidate, comparator upper bound, shared score, multiplicity rule, and pass rule before any test-expression access. The parameter values above are prespecified design constants, not observed results. The matrix-to-metadata join, exact external signature gene tables, executable comparator manifest, code hash and software environment must still be verified and committed before opening the count matrix. If an implementation change alters the estimand, predictors, eligible patients, score, or pass rule, amend and commit the plan before proceeding.
