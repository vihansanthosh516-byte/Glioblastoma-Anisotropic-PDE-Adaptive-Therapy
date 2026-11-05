# Track A Research Plan — Science Fair MD Final

**Scope:** This submission covers Track A A5 only: the frozen Yu 2020 Core-versus-Periphery prediction study. Track A New is a separate proposed biology study and is not part of this plan.

**Status note:** A5 is complete and frozen. The sections below state its protocol in research-plan form; they do not imply that the analysis remains unperformed. The completed result and limits are recorded in the [A5 plan and result record](Track%20A5%20Research%20Plan%20Science%20Fair%20MD%20Final.md).

## Research Question

Will a sparse, donor-balanced expression classifier improve donor-level Core-versus-Periphery prediction over ordinary logistic regression when every test donor is excluded from training?

## Hypothesis

The candidate classifier will achieve a lower mean donor-level class-balanced area under the risk-coverage curve (AURC) than ordinary L2 logistic regression. It will also achieve lower AURC in at least five of the six primary donors. This is a prediction hypothesis, not a claim about malignant-cell biology or causation.

## Procedure

1. **Select the frozen cohort and labels.** The analysis will use the public Yu et al. 2020 single-cell RNA-sequencing dataset, GEO GSE117891, and its published sample mappings. The six primary donors will be GS1, GS2, GS3, GS5, GS6, and GS13, as specified in the frozen protocol. Eligibility will use GBM histology and the recorded IDH1 R132 field. That field will not be described as comprehensive IDH1/2-wildtype confirmation. `Tumoral` sampling points will define Core, `Peritumoral` points will define Periphery, and `Normal` points will be excluded. These are MRI-marked surgical sampling-point labels, not histology labels assigned to individual cells.
2. **Map cells and apply quality filters.** Cell barcodes will be joined to sampling points and donor labels only through exact identifiers in the GEO cell table and Table S1. Unmapped or conflicting records will be excluded and counted. Cells will need at least 1,000 UMIs and 500 detected genes. No region label will be inferred from expression or file order.
3. **Preprocess without donor leakage.** Within each training fold, genes detected in at least 5% of training cells will be retained, and the 2,000 most variable genes will be selected. Counts will be normalized to 10,000 per cell and log-transformed with `log1p`. Feature means and standard deviations will be fitted on training cells only, then applied unchanged to validation and test cells.
4. **Hold out donors.** The primary evaluation will use six outer leave-one-donor-out folds. Inner tuning will also hold out one training donor at a time. All cells from a donor will remain in one fold. Preprocessing, feature selection, hyperparameter selection, and model fitting will use training donors only. The random seed will be 20261031. The held-out donor will not connect to other test cells through graph edges.
5. **Fit the candidate and comparisons.** The frozen candidate will use a sigmoid gene-selection gate, L1 gate regularization, ridge classifier weights, a smooth worst-donor loss with temperature 0.1, training-only 15-nearest-neighbor graph smoothing, and a test-time neighbor-probability blend. Ordinary L2 logistic regression will be the primary comparator. Secondary comparisons will include a 15-nearest-neighbor vote, linear support-vector machine, gradient-boosted trees, three candidate ablations, a cell-type-only diagnostic, a point-size diagnostic, and within-cell-type checks when the prespecified cell-count rule is met. Any comparator that returns constant scores will be labeled a constant-score control.
6. **Measure performance and apply the frozen rule.** The primary measure will be class-balanced AURC, calculated for each donor and averaged with equal donor weight; lower values will be better. Balanced accuracy and Brier score will be secondary measures. The candidate will meet its success rule only if mean AURC is lower than logistic regression and candidate AURC is lower in at least five of six donors. Every donor-level result will be reported. No p-value will be used to claim classifier superiority. If cell-type-only balanced accuracy is within 0.05 of the best expression model on the same mapped cells, the result will be described as comparable prediction from measured cell-type composition; this will not be treated as proof that composition explains the expression result.
7. **Limit interpretation.** The study will be described as a within-cohort prediction analysis using mixed tumor and peritumoral cell expression with MRI-marked sampling-point labels. It will not be described as malignant-cell-specific biology, causation, equivalence, clinical usefulness, or external validation. The six-donor result will be treated as descriptive, and any within-cell-type analysis that misses its cell-count rule will be reported as not possible.

## Materials Required

- **Computer:** Dell 16 Plus DB16250, Windows 11, 32 GB RAM; CPU-only, 8-thread setting.
- **Python and libraries:** Python 3.14.4; NumPy 2.5.1, pandas 2.3.3, scikit-learn 1.9.1, and PyTorch 2.14.0+cpu for A5. Track A supporting work will use SciPy, PyTorch Geometric, Scanpy, h5py, Matplotlib, Seaborn, openpyxl, and pypdf.
- **A5 data:** Yu et al. 2020 / GEO GSE117891 count matrix and cell-to-sampling-point metadata; Table S1 for point-region labels; and the GBmap Core / CELLxGENE `annotation_level_2` metadata snapshot for the cell-type diagnostic only. GBmap expression values will not be used as A5 model inputs.
- **Earlier Track A exploratory data:** UCSC Cell Browser multiomic-GBM, Ivy GAP / GSE107559, Angers / E-MTAB-5092, and E-MTAB-8743. These datasets will be listed as prior exploratory materials and will not be inputs to the A5 model or treated as untouched validation.
