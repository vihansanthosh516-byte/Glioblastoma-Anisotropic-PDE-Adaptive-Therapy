# Track A New — Proposed Malignant-Cell Biology Research Plan

**Status:** Proposed study only. No qualifying untouched validation cohort passes the current gates. No malignant-cell signature, primary endpoint, or success rule is frozen. This plan does not amend A5 or authorize data access or analysis.

## Rationale

Glioblastoma varies across core and infiltrative regions. Published studies already describe regional programs, infiltrating malignant cells, tumor margins, spatial niches, and invasion-associated states (Yu et al., 2020; Darmanis et al., 2017; Wang et al., 2024; Pai et al., 2025; Sonpatki et al., 2026; Chanoch-Myers et al., 2026; Mangukiya et al., 2026). The contribution here cannot be to claim that core and periphery differ for the first time. The stronger question is whether a regional expression contrast remains within validated malignant cells after separating changes in malignant-cell-state proportions from within-state expression shifts, then transfers across independent patients.

This plan aims for a donor-level, controlled biological association. It will not establish causation, a new mechanism, treatment response, or clinical utility.

## Research Question

Across at least six independent IDH-wildtype GBM donors with matched Core and infiltrative Periphery samples, does a predeclared malignant-cell expression signal persist within malignant-cell states after accounting for normal-brain admixture and differences in state proportions?

## Hypothesis

After the signature and analysis are fixed in a permitted discovery source, the paired Core-versus-Periphery signal will remain directionally consistent in malignant-cell measurements from at least six untouched validation donors after standardizing malignant-cell-state composition and accounting for normal-brain admixture. A result that disappears after state standardization would instead support a state-mixture explanation. Neither outcome will be interpreted as causal.

## Procedure

1. **Keep A5 separate.** A5 remains a completed mixed-cell prediction study with MRI-marked sampling-point labels. Its research question, cohort, protocol, result, and figure remain unchanged. A5’s failed classifier does not count as biological evidence for this question and will not supply a malignant-cell signature.
2. **Identify a specific cohort only through a named-lead metadata review.** Broad scouting is paused. Before any data download or expression review, the proposed validation cohort must have no-request public access, explicit dataset reuse terms, at least six previously unused IDH-wildtype GBM donors with matched Core and Periphery samples, exact donor/sample/region joins, per-case IDH1/2 evidence, source-supported cell-level malignant labels with validation, coherent assay and region definitions, and verified nonoverlap with the local donor manifest and prior-use audits. Article and dataset licenses will be checked separately. Spatial spots alone will not count as malignant-cell labels. Any unknown gate stops the study.
3. **Use the current screening ledger as the starting point, not as approval.** No screened cohort currently passes all gates. GSE286419 has four retained Core/Margin pairs. GSE318564 has three explicit IDH-wildtype Core/Edge pairs and an unresolved annotation join. SNUH-GliomaAtlas has six pairs but overlaps prior project use. HRA009021 lacks a public sample-to-region map and confirmed reuse terms. Harwood Zenodo 22013127 / concept record 16957940 remains unconfirmed: donor/sample mapping, per-case IDH status, FOV-region mapping, annotation fields and validation, and patient-level overlap are unresolved. The dataset-record license is separate from the article license. The paper already studies core-versus-infiltrated differences and CosMx is targeted; this lead does not establish new biology by itself. The ZIP scope deviation is documented in the [Track A dataset plan](trackA_alternate_dataset_plan_20261101.md). Broad scouting remains paused.
4. **Freeze the design before opening validation expression values.** A dated amendment must name the discovery source, untouched validation accession, exact locked signature or model, region mapping, malignant-label rule, primary donor-level endpoint, effect direction, success and precision criteria, multiplicity plan, and permitted data workflow. The validation set must contain at least six independent donors after all exclusions. If only a published signature is used, describe the work as a transfer or replication test. Do not claim signature discovery without a separate discovery dataset.
5. **Map samples, regions, and malignant cells from source metadata.** Join each cell to its sample, donor, region, and malignant label using documented identifiers. Confirm per-case IDH status and label provenance and validation. Count exclusions and conflicts. Do not infer labels from expression, file order, suffixes, or clusters; do not substitute spatial spots for malignant cells. Stop if any mapping or label gate is unclear.
6. **Build paired donor-level measurements.** Aggregate validated malignant-cell counts by donor, region, and prespecified malignant-cell state into pseudobulk measurements. Donors will be the independent units and cells will be subsamples. Use identical region definitions across discovery and validation. Measure normal-brain admixture and state proportions where source annotations support them.
7. **Separate state composition from within-state expression.** Estimate the total paired regional difference, then standardize malignant-state proportions and estimate within-state regional contrasts. Compare the locked signal with published regional and invasion signatures and expression-matched control gene sets. Fit transformations, feature selection, and any model only in discovery donors. Do not tune on validation donors.
8. **Test once in the untouched cohort and report every donor.** Apply the locked model and preprocessing once. Report each donor’s paired difference, equal-donor estimate, uncertainty, exclusions, state-standardized result, and sensitivity to normal-brain content. Use the prespecified success and precision criteria. Do not treat cells as independent patients or add donors after seeing results.
9. **Interpret both positive and negative outcomes.** A persistent signal supports a replicated within-malignant-cell regional association. A signal that disappears after state standardization indicates that the total contrast tracked state proportions under this design. A failed or imprecise result does not establish equivalence or absence. Pathway enrichment may be shown only as exploratory, with multiple-testing correction. Keep any future result separate from the current figure.

## Data Analysis

The primary endpoint will be a donor-level paired Core-minus-Periphery difference in a locked malignant-cell expression score, computed from pseudobulk and averaged equally across validation donors. The dated amendment must define the score, direction, primary test or interval, minimum precision, exclusions, and multiplicity rule before validation expression is opened. Secondary analyses will partition the total regional contrast into state-composition and within-state components and test sensitivity to normal-brain admixture. No cell-level p-value will be treated as donor-level evidence. Six donors are the eligibility floor, not a guarantee of adequate precision.

## Materials and Access

- **Computer and software:** The existing Track A computer is a Dell 16 Plus DB16250 with Windows 11 and 32 GB RAM. Existing Track A Python libraries include NumPy, pandas, scikit-learn, PyTorch, SciPy, PyTorch Geometric, Scanpy, h5py, Matplotlib, Seaborn, openpyxl, and pypdf. Additional packages, if needed, will be declared in the amendment.
- **Data:** A specific permitted discovery source and a separate untouched public validation accession have not been selected. No dataset is approved for Track A New. No data will be downloaded or scored before all metadata, independence, reuse, label, and school-review gates pass.
- **Data handling:** Use only permitted public data, keep source files local, do not upload patient-level data, and follow repository conditions. Resolve any school or institutional review required for future use before access.

## Stop Rules

Stop before expression access if fewer than six independent matched validation donors remain, any case lacks documented IDH1/2 status, exact donor/sample/region mapping fails, malignant labels lack source-backed validation, reuse terms are unclear, overlap cannot be ruled out, or assay/region definitions cannot be harmonized. Stop the planned confirmation claim if the validation cohort has already been used or if the locked signature, endpoint, or success rule changes after validation is opened.

## Novelty Boundary

Core-versus-margin biology, malignant infiltration signatures, and regional cell states are established topics. The potential contribution is a controlled, donor-level separation of state-mixture effects from within-state malignant-cell expression differences, followed by untouched independent validation. This contribution remains conditional on finding a qualifying cohort and obtaining a result that passes the frozen criteria.

## Bibliography

Chanoch-Myers, R., Hara, T., Greenwald, A. C., Hoefflin, R., Bussema, L., Weisman, H. R., Nieblas-Bedolla, E., Gurramkonda, J., Bekezhankyzy, Z., Jucht, A., Calvo Fernández, E., Spitzer, A., Fung, J. M. H., Chiara, B., Al-Holou, W. N., Camelo-Piragua, S., Behnan, J., Galli, R., Wakimoto, H., . . . Tirosh, I. (2026). A blueprint for local and distal invasion programs in glioblastoma. *Nature Communications, 17*(1), Article 4079. https://doi.org/10.1038/s41467-026-70470-8

Darmanis, S., Sloan, S. A., Croote, D., Mignardi, M., Chernikova, S., Samghababi, P., Zhang, Y., Neff, N., Kowarsky, M., Caneda, C., Li, G., Chang, S. D., Connolly, I. D., Li, Y., Barres, B. A., Gephart, M. H., & Quake, S. R. (2017). Single-cell RNA-seq analysis of infiltrating neoplastic cells at the migrating front of human glioblastoma. *Cell Reports, 21*(5), 1399–1410. https://doi.org/10.1016/j.celrep.2017.10.030

Mangukiya, H. B., Skeppås, M., Lasses Armatowski, A., Uppman, I., Hallin, L., Kruse, F., Nordquist, C., Kundu, S., Berglund, M., Malik, A. A., Krona, C., Jörnsten, R., & Nelander, S. (2026). Reconstructing the single-cell spatiotemporal dynamics of glioblastoma invasion. *Nature Communications, 17*, Article 11245. https://doi.org/10.1038/s41467-026-78908-9

Pai, B., Ramos, S. I., Cheng, W. S., Joshi, T., Özen, E., Kulumani Mahadevan, L. S., Silva-Hurtado, T. J., Price, G. A., Tome-Garcia, J., Nudelman, G., Shroff, S., Beaumont, K. G., Yong, R. L., Sebra, R. P., Zaslavsky, E., & Tsankova, N. M. (2025). Spatial multi-omics defines a shared tumor infiltrative signature at the resection margin in high-grade gliomas. *Cancer Research, 85*(21), 4233–4250. https://doi.org/10.1158/0008-5472.CAN-24-4708

Sonpatki, P., Park, H. J., Xing, Y. L., Han, K. Y., Schroeder, B. A., Yu, H. J., Kim, H. J., Chowdhury, T., Hwang, J. H., Nam, S. M., Byun, Y. H., Kang, H., Lee, J. H., Lee, S. T., Won, J. K., Kim, T. M., Choi, S. H., Ku, J. L., Lee, S., . . . Shah, N. (2026). A spatially resolved human glioblastoma atlas reveals distinct cellular and molecular patterns of anatomical niches. *Nature Communications, 17*(1), Article 2951. https://doi.org/10.1038/s41467-026-69716-2

Wang, F., Dong, J., Xu, Y., Jin, J., Xu, Y., Yan, X., Liu, Z., Zhao, H., Zhang, J., Wang, N., Hu, X., Gao, X., Xu, L., Yang, C., Ma, S., Du, J., Hu, Y., Ji, H., & Hu, S. (2024). Turning attention to tumor–host interface and focus on the peritumoral heterogeneity of glioblastoma. *Nature Communications, 15*(1), Article 10885. https://doi.org/10.1038/s41467-024-55243-5

Yu, K., Hu, Y., Wu, F., Guo, Q., Qian, Z., Hu, W., Chen, J., Wang, K., Fan, X., Wu, X., Rasko, J. E. J., Fan, X., Iavarone, A., Jiang, T., Tang, F., & Su, X.-D. (2020). Surveying brain tumor heterogeneity by single-cell RNA-sequencing of multi-sector biopsies. *National Science Review, 7*(8), 1306–1318. https://doi.org/10.1093/nsr/nwaa099