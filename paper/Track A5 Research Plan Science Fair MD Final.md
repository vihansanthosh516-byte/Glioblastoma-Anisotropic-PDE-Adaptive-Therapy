# Track A A5 — Frozen Research Plan and Completed Record

**Status:** Completed and frozen. The candidate failed the prespecified 5-of-6 donor success rule. No retuning or figure change is authorized by this record.

## Research Question

Will a sparse, donor-balanced expression classifier improve donor-level Core-versus-Periphery prediction over ordinary logistic regression when every test donor is excluded from training?

## Hypothesis and Result

The candidate was predicted to have lower mean class-balanced area under the risk-coverage curve (AURC) than ordinary L2 logistic regression and lower AURC in at least five of six primary donors. Mean AURC was 0.1479 for the candidate versus 0.1764 for logistic regression. It won in GS1, GS2, GS5, and GS6, and lost in GS3 and GS13. The frozen 5-of-6 rule failed, so superiority was not demonstrated. Candidate Brier score was worse (0.2335 versus 0.1944), and the no-gate ablation had slightly lower AURC (0.1473 versus 0.1479).

## Procedure

1. **Use the frozen cohort and source mappings.** A5 used Yu et al. 2020 single-cell RNA-sequencing data from GEO GSE117891. The six primary donors were GS1, GS2, GS3, GS5, GS6, and GS13. Exact GEO barcode-to-sampling-point and sampling-point-to-Table S1 joins identified 75 labeled points with no unmapped target points. The cohort rule used GBM histology and the recorded IDH1 R132 field; a dash in that field does not establish comprehensive IDH1/2-wildtype status.
2. **Apply the frozen inclusion, label, and quality rules.** Matched tumor and peritumoral cells were retained after requiring at least 1,000 UMIs and 500 detected genes. The analysis included 4,139 quality-controlled matched target cells, of which 2,519 belonged to the six scored donors. Author-supplied Tumoral points defined Core, Peritumoral points defined Periphery, and Normal points were excluded. Labels represent MRI-marked sampling points, not histology labels assigned to individual cells.
3. **Evaluate by held-out donor.** The analysis used six outer leave-one-donor-out folds and nested leave-one-training-donor-out tuning. Preprocessing and feature selection were fitted on training donors only. No test-to-test graph edges were used.
4. **Fit the frozen candidate and comparators.** The candidate used a sigmoid gene-selection gate with L1 regularization, ridge classifier weights, smooth worst-donor cross-entropy with temperature 0.1, training-only 15-nearest-neighbor graph smoothing, and a test-time neighbor-probability blend. It used Adam at learning rate 0.001, up to 500 epochs, and early stopping after 25 epochs without inner-validation improvement. Logistic regression was the primary comparator, with C values 0.001, 0.01, 0.1, 1, and 10. Secondary controls included a 15-nearest-neighbor vote, linear SVM, gradient-boosted trees, three candidate ablations, cell-type-only prediction, point-size diagnostics, and within-cell-type checks when the frozen sample rule was met.
5. **Report all donor-level outcomes and controls.** Mean balanced accuracy was 0.6786 for the candidate and 0.6621 for logistic regression. Cell-type-only balanced accuracy was 0.6546, while the best expression model on the same mapped cells reached 0.6789, a 0.0243 difference. Under the frozen interpretation, cell-type composition gave comparable prediction; this does not prove composition explains the expression result. Only the Myeloid stratum met the within-cell-type size rule, with five eligible donors. The linear SVM was a constant-score control in all six donors. The point-size diagnostic was a constant-score control in GS1, GS2, GS3, and GS5, and varied in GS6 and GS13.
6. **Apply the frozen analysis and interpretation rules.** Class-balanced AURC was the primary endpoint, averaged equally over six donors; lower AURC was better. Balanced class risk was calculated across retained-cell coverage from ceil(N/2) through N and integrated over actual coverage values with the trapezoid rule. Balanced accuracy and Brier score were secondary outcomes. The candidate could be called better only if it had lower mean AURC and lower AURC in at least five of six donors. The project reported donor-level results without p-values for classifier superiority.
7. **State the limits and deferred analyses.** A5 is a within-cohort result on mixed-cell expression and MRI-marked point labels. It does not show malignant-cell-specific biology, causation, equivalence, clinical usefulness, or external validation. The eight-donor sensitivity cohort was not run. Calibration, within-library, scANVI, and attention analyses are **not run; deferred** pending a future dated amendment. No additional A5 analysis is authorized here.

## Materials

- **Computer:** Dell 16 Plus DB16250, Windows 11, 32 GB RAM; CPU-only, 8-thread setting.
- **Software used for A5:** Python 3.14.4; NumPy 2.5.1; pandas 2.3.3; scikit-learn 1.9.1; PyTorch 2.14.0+cpu.
- **Supporting Track A libraries:** SciPy, PyTorch Geometric, Scanpy, h5py, Matplotlib, Seaborn, openpyxl, and pypdf.
- **Data:** Yu et al. 2020 / GEO GSE117891 count matrix and mapping metadata; Supplementary Table S1 for MRI-marked point labels; GBmap Core / CELLxGENE annotation_level_2 metadata snapshot for the cell-type diagnostic only. GBmap expression values were not A5 features.
- **Protocol and results:** [Frozen A5 protocol](trackA_yu2020/PROTOCOL.md), [result JSON](../output/trackA_yu2020/result.json), and [execution audit](trackA_execution_audit_20261031.md).