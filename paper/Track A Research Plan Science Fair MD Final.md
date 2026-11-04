# Track A Research Plan — Science Fair MD Final

**Status:** Consolidated plan and record of completed Track A work, updated November 3, 2026
**Scope:** Track A only

This plan keeps the frozen A5 research question, protocol, result, and current figures unchanged. A5 is complete. Its candidate did not meet the frozen success rule. The higher-potential biology study below is a separate proposal that would require a dated amendment, a qualifying independent cohort, and a new frozen protocol before any new expression values are read.

## Rationale

Glioblastoma (GBM) cells and their surrounding brain environment vary across tumor core, infiltrative margin, and periphery. This matters because infiltrating cells can remain outside the visible tumor mass. Earlier studies already report regional differences, infiltrating malignant-cell programs, and spatially organized tumor states. Yu et al. (2020) profiled multi-sector glioma biopsies; Darmanis et al. (2017) described infiltrating neoplastic cells and a shared infiltration signature; Wang et al. (2024) studied peritumoral heterogeneity; Pai et al. (2025) defined a tumor-margin infiltration signature; and Sonpatki et al. (2026) reported six matched core–periphery single-cell samples with CNV-supported malignant-cell labels. Recent model studies also investigate local and distal invasion programs and dynamic invasion behavior (Chanoch-Myers et al., 2026; Mangukiya et al., 2026).

Track A first explored cell-level classifiers in a mixed-cell atlas. Those early random-cell results were vulnerable to donor leakage and region/library confounding. A5 addressed the prediction question with donor-held-out testing, but used mixed tumor and peritumoral cells and MRI-marked sampling-point labels. Its result is a within-cohort benchmark, not evidence of malignant-cell-specific biology. The strongest future biological question is therefore whether a predeclared regional signal remains within validated malignant cells and transfers across genuinely independent patients after accounting for normal-brain content, cell composition, and established tumor-cell states. A positive result would support a regional association and its transportability; it would not establish causation or a treatment mechanism.

## Research Question

**Frozen A5 question:** Will a sparse, donor-balanced expression classifier improve donor-level Core-versus-Periphery prediction over ordinary logistic regression when every test donor is excluded from training?

This question is unchanged from the frozen A5 protocol. The separate future biology extension does not replace it.

## Hypotheses and Expected Outcomes

### Frozen A5 hypothesis and observed outcome

The frozen hypothesis predicted that the sparse candidate would have lower mean class-balanced area under the risk-coverage curve (AURC) than ordinary L2 logistic regression and would have lower AURC in at least five of six primary donors. A5 met the mean-AURC condition but won in only four of six donors. The hypothesis was not supported under its frozen success rule.

The candidate mean AURC was 0.1479 versus 0.1764 for logistic regression. Candidate Brier score was worse, 0.2335 versus 0.1944. The no-gate ablation had slightly lower mean AURC than the candidate, 0.1473 versus 0.1479. The result does not demonstrate superiority.

### Proposed future biology hypothesis

If a separate study is adopted, it will test whether a predeclared Core-versus-Periphery transcriptomic signal remains detectable within defensibly labeled malignant cells across at least six new, independent IDH-wildtype GBM donors after accounting for normal-brain content, cell-type composition, and established malignant-cell states.

The outcome may support, weaken, or fail to resolve this hypothesis. No malignant-cell signature is currently frozen, no eligible external cohort currently passes the gates, and no positive result is promised. Even a successful result would establish an association or independent transfer result, not causation, clinical utility, or a new mechanism by itself.

## Track A Work Completed

### Early classifier and composition work

The first Track A classifier work used the UCSC Cell Browser multiomic-GBM atlas. A historical random-cell split reported 78.7% C-GAT accuracy, but a patient-ID lookup reached 82.8%. The original graph also included region-derived edge features, and the model-selection record used test performance. These results do not establish generalization to new donors. A later patient-held-out C-GAT run reached 63.8% accuracy, near the 64.6% cells-per-patient baseline. An exploratory gene-expression region classifier reached AUC 0.904 versus 0.713 for cell-type-only prediction, but region was confounded with sequencing library and malignant periphery cells were sparse. These findings motivated donor-level evaluation and composition controls; they are not confirmatory evidence. Historical details are in the [Track A output record](../archive_old_gbm/vault/TRACK_A_RESULTS.md) and [spatial handoff](trackA_spatial_handoff.md).

### Bulk regional gene-program tests

Three hypothesis-informed bulk-tissue tests compared regional gene programs with expression-matched random sets and known tissue gradients. All failed their prespecified adjusted success rules:

| Test | Cohorts and raw result | Prespecified adjustment | Adjusted result |
|---|---|---|---|
| SenMayo senescence/SASP | Ivy GAP (8 donors) plus Angers E-MTAB-5092 (7 patients); raw Core > Edge, p = 0.0035 | Hypoxia | p = 0.37 |
| Hallmark OXPHOS | Same 15 patients; raw Edge > Core, p = 0.0024 | Hypoxia and neuron-marker scores | p = 0.67 |
| TAM activation | E-MTAB-8743 (10 patients); raw Core > Invasive, p = 0.0039 | Hypoxia, neuron, and myeloid-marker scores | p = 0.44 |

These nonsignificant adjusted tests mean the regional effects were not demonstrated under those models. They do not prove that the measured covariates explain the effects or that the effects are absent. The gene lists followed exploratory inspection, so these tests were hypothesis-informed rather than blind. Full results and limits are in [the Track A composition finding](trackA_composition_finding.md).

### Frozen A5 donor-held-out analysis

A5 used Yu et al. 2020 single-cell RNA-sequencing data from GEO GSE117891. The frozen primary donors were GS1, GS2, GS3, GS5, GS6, and GS13. The cohort rule used GBM histology and the recorded IDH1 R132 field; that field is not comprehensive evidence of IDH1/2-wildtype status.

The exact GEO barcode-to-sampling-point and sampling-point-to-Table S1 joins identified 75 labeled points, with no unmapped target points. The analysis included 4,139 matched tumor/peritumoral target cells that passed frozen quality control; 2,519 belonged to the six scored donors. Quality control required at least 1,000 UMIs and 500 detected genes. Core was defined by the author-supplied Tumoral point label, Periphery by Peritumoral, and Normal points were excluded.

The sealed analysis used six outer leave-one-donor-out folds and nested leave-one-training-donor-out tuning. Preprocessing and feature selection were fitted only on training donors. The candidate used a sigmoid gene-selection gate with L1 regularization, ridge classifier weights, smooth worst-donor cross-entropy loss with temperature 0.1, training-only 15-nearest-neighbor graph smoothing, and a test-time neighbor-probability blend. It used Adam at learning rate 0.001, up to 500 epochs, and early stopping after 25 epochs without inner-validation improvement. Ordinary L2 logistic regression was the primary comparator, with C values 0.001, 0.01, 0.1, 1, and 10. Secondary controls included a 15-nearest-neighbor vote, linear SVM, gradient-boosted trees, three candidate ablations, cell-type-only prediction, point-size diagnostics, and within-cell-type checks where the frozen sample rule was met.

The candidate won in GS1, GS2, GS5, and GS6, and lost in GS3 and GS13. Mean balanced accuracy was 0.6786 for the candidate and 0.6621 for logistic regression. The cell-type-only diagnostic had balanced accuracy 0.6546; the best expression model on the same mapped cells reached 0.6789, a 0.0243 difference. Under the frozen interpretation, cell-type composition gave comparable prediction. This does not prove composition explains the expression result. Only the Myeloid stratum met the within-cell-type size rule, with five eligible donors.

The linear SVM was a constant-score control in all six donors. The point-size diagnostic was a constant-score control in GS1, GS2, GS3, and GS5; it varied in GS6 and GS13. Neither is a substantive model competitor.

A5 is a within-cohort result on mixed-cell expression and MRI-marked point labels. Labels are not histology labels assigned to individual cells. The study does not show malignant-cell-specific biology, causation, equivalence, clinical usefulness, or independent external validation. The eight-donor sensitivity cohort was not run. Calibration, within-library, scANVI, and attention analyses are **not run; deferred**. A future dated amendment must decide their disposition.

The result and protocol audit are available in [the A5 result JSON](../output/trackA_yu2020/result.json), [the execution audit](trackA_execution_audit_20261031.md), and [the frozen protocol](trackA_yu2020/PROTOCOL.md).

## Proposed Higher-Potential Biology Extension

### Question and potential contribution

The proposed question is: **Does a predeclared Core-versus-Periphery signal within validated malignant GBM cells transfer across at least six new, independent IDH-wildtype donors after adjustment for normal-brain content, cell-type composition, and established malignant-cell states?**

Core-versus-margin biology and malignant-cell regional signatures have already been reported. The possible contribution would be an independently tested, donor-level estimate of how much regional signal remains after these controls. A useful biological result would distinguish a change in the proportions of established malignant-cell states from an expression shift within the same state. It would not be a claim that core and periphery differences were discovered for the first time. Current literature also includes cell-state and invasion-program analyses in experimental models, so a broad claim about a novel invasion mechanism is not supported by this computational plan.

No analysis-ready malignant-cell signature is frozen. Before any validation data are opened, a dated amendment must name the discovery source, the untouched validation accession, the exact signature or model, the primary donor-level endpoint, the effect direction, the success and precision criteria, and the permitted data workflow. A5’s failed candidate is not a malignant-cell signature and cannot be relabeled as one. If the only adopted signature comes from published prior work, the study will be framed as a transfer or replication test, not signature discovery.

### Required cohort gates

A named public accession must pass every gate before download or expression review:

1. Public, no-request access and clear data-specific access and reuse terms. The article license and dataset-record license must be checked separately.
2. At least six eligible IDH-wildtype GBM patients, with Core and Periphery or a prespecified equivalent mapped to the same patient in public metadata. More donors are preferred; six donors do not guarantee adequate precision.
3. Exact sample-to-donor and sample-to-region joins, with donor-level IDH1/2 status established for each included patient.
4. Defensible cell-level malignant labels with documented provenance and validation, such as source-supported CNV evidence. Spatial spots alone cannot label malignant cells.
5. Exact sample and patient identifiers cross-checked against the local donor manifest, training sources, and prior-use audits. Independence must be verified, not assumed.
6. A biologically coherent assay and region definition. Do not pool incompatible studies or assays to reach six donors.

The current screen has no cohort that passes all gates. GSE286419 is the closest label-ready lead, with four retained Core/Margin pairs. GSE318564 has three explicit IDH-wildtype Core/Edge pairs and an unresolved annotation join. SNUH-GliomaAtlas has six pairs but overlaps prior project use. HRA009021 lacks a public sample-to-region map and confirmed reuse terms. Harwood Zenodo 22013127, concept record 16957940, remains unconfirmed and ineligible: donor/sample mapping, per-case IDH status, FOV-to-region mapping, annotation fields and validation, and patient-level overlap remain unresolved. Zenodo record metadata reports CC BY 4.0 for the dataset; that record license is separate from the article license. The paper already studies core-versus-infiltrated differences, and CosMx is a targeted spatial assay, so this lead would not establish new biology by itself. The ZIP tail requests fetched some compressed member bytes into memory; none were decompressed, inspected, or saved, and no further archive requests occurred. This was a scope deviation, not a successful listing-only read. Exact retained ranges and the missing Sample C central-directory endpoints are recorded in [the metadata and scope-deviation plan](trackA_alternate_dataset_plan_20261101.md).

Broad dataset scouting remains paused. Reopen metadata screening only for a specific named accession that could plausibly pass all gates. Any missing donor map, IDH status, region label, malignant-cell validation, reuse term, or verified independence is a stop. Do not lower the six-donor minimum.

### Proposed procedure after amendment and all gates pass

1. **Freeze the study design.** Select a permitted discovery source and a separate untouched validation cohort. If a new signature is to be discovered, use only the discovery source for gene selection, model fitting, covariate decisions, and threshold choices. Keep all validation donors out of development. If no suitable discovery source passes its own gates, use a literature-defined signature and describe the work as replication.
2. **Verify records before data access.** Record the accession, source-paper version, licenses, exact file list, donor IDs, sample IDs, regions, per-case IDH status, cell-label method, validation evidence, and local overlap crosswalk. Stop if any item is unknown.
3. **Map cells and regions exactly.** Join each cell to its sample, patient, region, and malignant label through source-supported identifiers. Count exclusions and conflicts. Do not infer regions from file order, expression, suffixes, or cell-state clusters.
4. **Define malignant cells.** Apply the source’s validated cell-level annotation rule. Keep nonmalignant cells out of the malignant-cell analysis. Quantify nonmalignant content separately when the source supports it. Do not use spatial spots as malignant-cell labels.
5. **Construct paired donor measurements.** Aggregate malignant-cell counts by donor and region into pseudobulk measurements, or use another method frozen in the amendment. Cells are subsamples; donors are the independent units. Use the same region definitions in discovery and validation.
6. **Control known alternatives.** Measure normal-brain content and cell-type composition. Compare regional expression within established malignant-cell states and standardize state composition so the analysis can separate state-mixture changes from within-state expression shifts. Compare the candidate with published regional and cell-state signatures and expression-matched gene-set controls. Keep feature selection, scaling, model fitting, and tuning inside discovery donors.
7. **Run the untouched validation once.** Apply the locked gene list or model and normalization to the validation cohort without feature selection, threshold changes, or hyperparameter retuning. Report each donor’s paired effect, the equal-donor summary, uncertainty, and all exclusions. Do not count cells or spots as independent samples.
8. **Interpret the biology.** If a stable adjusted signal remains within malignant cells, describe it as a donor-level regional association. If the signal disappears after state standardization, report that the observed contrast tracked the proportions of established states. If it remains within states, report a residual within-state regional association. A prespecified pathway or regulatory enrichment of the locked residual genes may be reported as exploratory interpretation, with multiple-testing correction; it will not establish a mechanism.
9. **Report both outcomes.** A stable adjusted signal would support transfer across this cohort. An inconsistent or absent signal would limit the hypothesis. Neither outcome would establish causation, therapy response, or clinical utility. Keep any future external result in a separate, accurately labeled figure; do not change the current figure under this plan.

The extension’s primary endpoint and success rule remain open until a dated amendment sets them before validation expression is opened. The amendment must include a donor-level precision or power rationale, a single primary endpoint, a rule for multiple secondary tests, and a stop condition for insufficient malignant cells or failed mapping. Six donors are an eligibility floor, not proof of high power. Do not change the endpoint or success rule after scoring.

## Data Analysis

For completed A5, class-balanced AURC was the primary endpoint and was averaged equally over six donors. Lower AURC was better. Balanced accuracy and Brier score were secondary outcomes. The project reported per-donor results and donor-level summaries, without p-values for classifier superiority. The frozen candidate-success rule required both lower mean AURC and lower AURC in at least five of six donors; it failed.

For the proposed extension, the amendment will specify one donor-level primary measure, such as the paired difference in a locked malignant-cell signature score between Core and Periphery pseudobulk. The analysis will show every donor’s paired value, the equal-donor effect estimate, and uncertainty. A donor-level paired test and a precision criterion will be fixed before expression access. Secondary gene-level results will use a prespecified multiplicity correction. Figures may include paired donor plots and a donor-level effect plot. No cell-level p-value will be treated as independent patient evidence.

## Materials Required

### Computer and software

- **Computer used for A5:** Dell 16 Plus DB16250 running Windows 11, with 32 GB RAM. A5 ran CPU-only with an 8-thread setting; no laboratory equipment was used.
- **A5 software versions:** Python 3.14.4; NumPy 2.5.1; pandas 2.3.3; scikit-learn 1.9.1; and PyTorch 2.14.0+cpu.
- **Other third-party Python libraries used in Track A scripts:** SciPy, PyTorch Geometric, Scanpy, h5py, Matplotlib, Seaborn, openpyxl, and pypdf.
- **Python standard-library modules used in Track A scripts:** argparse, collections, csv, dataclasses, datetime, gc, gzip, hashlib, importlib, itertools, json, math, os, pathlib, platform, random, re, resource, subprocess, sys, time, typing, warnings, zipfile, urllib.parse, and xml.etree.ElementTree.

### Datasets used in completed Track A work

- **Yu et al. 2020 / NCBI GEO GSE117891:** single-cell RNA counts, cell-to-sample and cell-to-point metadata, and Supplementary Table S1 point labels. This was the A5 expression cohort. The project used the source’s Tumoral and Peritumoral point labels as Core and Periphery, respectively.
- **GBmap Core / CELLxGENE:** a metadata-only cell annotation snapshot joined to Yu cells by exact barcode for the A5 composition diagnostic and earlier donor audits. GBmap expression values were not used as A5 features.
- **UCSC Cell Browser multiomic-GBM atlas:** used in the historical exploratory classifier and composition analyses. This cohort is exposed and cannot serve as untouched validation.
- **Ivy GAP / GSE107559 and Angers / E-MTAB-5092:** bulk regional expression data used in the SenMayo and OXPHOS tests.
- **E-MTAB-8743:** bulk regional expression data used in the TAM-activation test.

Metadata-only cohort screens did not count as scored datasets. The full candidate ledger and source evidence are in [the Track A dataset plan](trackA_alternate_dataset_plan_20261101.md). No future validation dataset is currently approved or selected.

## Risk, Privacy, and Data Handling

Track A uses computational analysis of previously collected human tumor data. The project does not recruit people, administer surveys, use animals, culture organisms, or perform wet-lab procedures. It uses no chemicals, biological agents, or hazardous laboratory devices.

The files may originate from human clinical research. The project will use only datasets with permitted public access and clear reuse terms, keep source files local, avoid uploading patient-level data, and follow each repository’s conditions. Any school or institutional review required for future data use will be resolved before that use begins. The results are research observations and are not for diagnosis or treatment.

## Effort and Current Feasibility

A5 is complete; no further A5 scoring or retuning is planned. Its four secondary analyses remain not run; deferred pending a future dated amendment. External confirmation is unavailable today, and broad scouting remains paused.

For a specifically named future accession, a metadata-only gate review is estimated at one to three focused days per candidate. A separate discovery and validation source could require two such reviews. Access approval timing is unknown. Only after all gates pass and a dated protocol amendment is frozen is the analysis estimated at one to three weeks; a discovery-stage signature pipeline could require additional time. No completion date is defensible until a qualifying cohort and frozen signature are identified.

## Bibliography

Chanoch-Myers, R., Hara, T., Greenwald, A. C., Hoefflin, R., Bussema, L., Weisman, H. R., Nieblas-Bedolla, E., Gurramkonda, J., Bekezhankyzy, Z., Jucht, A., Calvo Fernández, E., Spitzer, A., Fung, J. M. H., Chiara, B., Al-Holou, W. N., Camelo-Piragua, S., Behnan, J., Galli, R., Wakimoto, H., . . . Tirosh, I. (2026). A blueprint for local and distal invasion programs in glioblastoma. *Nature Communications, 17*(1), Article 4079. https://doi.org/10.1038/s41467-026-70470-8

Darmanis, S., Sloan, S. A., Croote, D., Mignardi, M., Chernikova, S., Samghababi, P., Zhang, Y., Neff, N., Kowarsky, M., Caneda, C., Li, G., Chang, S. D., Connolly, I. D., Li, Y., Barres, B. A., Gephart, M. H., & Quake, S. R. (2017). Single-cell RNA-seq analysis of infiltrating neoplastic cells at the migrating front of human glioblastoma. *Cell Reports, 21*(5), 1399–1410. https://doi.org/10.1016/j.celrep.2017.10.030

Mangukiya, H. B., Skeppås, M., Lasses Armatowski, A., Uppman, I., Hallin, L., Kruse, F., Nordquist, C., Kundu, S., Berglund, M., Malik, A. A., Krona, C., Jörnsten, R., & Nelander, S. (2026). Reconstructing the single-cell spatiotemporal dynamics of glioblastoma invasion. *Nature Communications, 17*, Article 11245. https://doi.org/10.1038/s41467-026-78908-9

Pai, B., Ramos, S. I., Cheng, W. S., Joshi, T., Özen, E., Kulumani Mahadevan, L. S., Silva-Hurtado, T. J., Price, G. A., Tome-Garcia, J., Nudelman, G., Shroff, S., Beaumont, K. G., Yong, R. L., Sebra, R. P., Zaslavsky, E., & Tsankova, N. M. (2025). Spatial multi-omics defines a shared tumor infiltrative signature at the resection margin in high-grade gliomas. *Cancer Research, 85*(21), 4233–4250. https://doi.org/10.1158/0008-5472.CAN-24-4708

Sonpatki, P., Park, H. J., Xing, Y. L., Han, K. Y., Schroeder, B. A., Yu, H. J., Kim, H. J., Chowdhury, T., Hwang, J. H., Nam, S. M., Byun, Y. H., Kang, H., Lee, J. H., Lee, S. T., Won, J. K., Kim, T. M., Choi, S. H., Ku, J. L., Lee, S., . . . Shah, N. (2026). A spatially resolved human glioblastoma atlas reveals distinct cellular and molecular patterns of anatomical niches. *Nature Communications, 17*(1), Article 2951. https://doi.org/10.1038/s41467-026-69716-2

Wang, F., Dong, J., Xu, Y., Jin, J., Xu, Y., Yan, X., Liu, Z., Zhao, H., Zhang, J., Wang, N., Hu, X., Gao, X., Xu, L., Yang, C., Ma, S., Du, J., Hu, Y., Ji, H., & Hu, S. (2024). Turning attention to tumor–host interface and focus on the peritumoral heterogeneity of glioblastoma. *Nature Communications, 15*(1), Article 10885. https://doi.org/10.1038/s41467-024-55243-5

Yu, K., Hu, Y., Wu, F., Guo, Q., Qian, Z., Hu, W., Chen, J., Wang, K., Fan, X., Wu, X., Rasko, J. E. J., Fan, X., Iavarone, A., Jiang, T., Tang, F., & Su, X.-D. (2020). Surveying brain tumor heterogeneity by single-cell RNA-sequencing of multi-sector biopsies. *National Science Review, 7*(8), 1306–1318. https://doi.org/10.1093/nsr/nwaa099
