# GBM Tracks A, B, and C: complete research plan

This plan describes the full GBM Tracks A, B, and C project in one sequence. It states the selected method for each part and places existing evidence and its limits beside that method. Results from earlier analyses are exploratory evidence, not preregistered or untouched confirmation. Before submission, we need to verify every claim, state the student’s own contribution, and follow the fair’s disclosure rules.

## Rationale

Glioblastoma (GBM) varies between patients and within a tumor. Computational models use molecular data, MRI, and treatment information to study tumor location, growth, and therapy (Hormuth et al., 2021). These data sources do not share established patient identifiers, so they cannot support one fused patient-level “digital twin.” We need to evaluate Tracks A, B, and C as separate but related GBM analyses, then develop a mathematical model for postoperative MRI volume that separates latent net growth from measurement error. Existing results give a clear reason to test this: the recorded growth forecasts did not show a reliable gain over simple baselines, while treatment effects in simulations are assumed. Published work now includes a recurrence-planning comparison with a 15 mm margin (Zimmer et al., 2026), so our recurrence result is context and replication, not a novelty claim. Clearer uncertainty may help future studies interpret MRI change, but this retrospective project cannot guide care.

## Research Question/Hypothesis(es)/Engineering Goal(s)

**Primary question:** Does a hierarchical state-space model of latent postoperative log-volume and MRI measurement error improve next-scan prediction for held-out patients over a training-fold population-median log-volume slope (primary) and persistence (secondary)?

**Hypothesis:** The model improves patient-level predictive log score. We need to treat this prediction hypothesis as unsupported if the median difference is not positive or its patient-bootstrap 95% interval includes zero. If parameter recovery fails, we must not claim that measurement error and net growth are separately identified.

**Supporting questions:** How well do molecular labels predict tumor zone in patients excluded from training? Do the recorded MRI growth and recurrence methods improve on simple comparators? What outcomes do adaptive-therapy simulations produce under assumed drug effects?

**Engineering goal:** We need to build an auditable forecasting method that reports predictive uncertainty and makes its measurement-noise assumptions explicit.

**Expected outcome:** Any result, including no gain, poor calibration, or parameters that cannot be separated, is informative. We need to state which values come from patient data, simulations, or assumptions. We must not claim a linked molecular-to-MRI model, treatment benefit, or validated digital twin.

## Materials

We need a computer with 16 CPU cores and 34 GB RAM, Python and the project’s scientific-computing packages, and the existing project files. A Kaggle GPU is optional for neural-network training only if the exact data terms permit that host. No chemicals, organisms, laboratory equipment, or human contact are involved.

The project uses the public UCSC multiomic-gbm atlas (25 donors; 223,113 published cells; Sonpatki et al., 2026); MU-Glioma-Post v1 (203 patients; about 11 GB; CC BY 4.0); UCSF-PDGM (495 subjects; about 142 GB on the collection page; CC BY 4.0; Calabrese et al., 2022); and RHUH-GBM brain-extracted NIfTI derivatives (40 patients; about 2.9 GB; listed CC BY 4.0; raw DICOM is controlled access). LUMIERE (91 patients) is publicly downloadable, but its usage notes say non-commercial use while its Figshare record lists CC0, so we need to verify terms for the exact files. PREDICT-GBM is prior-art context only; its files are not an analysis input. The Abu Khalaf repeatability study (2021) is a sensitivity reference, not project data. Dataset citations and download links are in the appendix. We need to record the exact release and use terms before using or transferring each download.

## Procedures
We need to create a source and patient-visit manifest, then analyze Tracks A, B, and C without joining patients across cohorts. For Track A, we need donor-grouped splits and compare the attention-gated classifier with fixed PCA/logistic regression and a 200-gene panel, reporting macro-F1 and donor-level uncertainty. Existing 24-donor estimates (0.729, 0.765, and 0.710; PCA-minus-attention seed-42 interval −0.009 to +0.064) are exploratory. We need to reconcile the unmatched donor by a fixed source-manifest identifier rule, preserve the existing 16/3/5 donor assignments for resolved donors, and exclude an unresolved donor without reshuffling. For Track B, we need to compare reaction-diffusion forecasts with persistence and population geometric growth on the same held-out patients, and compare the selected H-2 map on LUMIERE/RHUH and H-3 map on UCSF with an equal-volume 15 mm expansion, using preoperative inputs only (enhancing recurrence primary; all recurrence secondary). Reused cohorts make these comparisons exploratory. For Track C, we need to compare adaptive and MTD schedules only in simulation with assumed kill rates and present MGMT survival results as observational. Then we need to check parameter recovery for the state-space model and compare rolling forecasts on MU-Glioma-Post against a training-fold population-median log-volume slope (primary) and persistence (secondary). Prior growth analyses used this cohort, so folds are exploratory; LUMIERE is a sensitivity check. We need to finish with one table of data, estimates, uncertainty, evidence type, and failure modes.
## Risk Assessment

We use de-identified, publicly accessible records and follow the terms for each file. We need to avoid re-identification, guessed patient linkage, redistribution of files without permission, and transfer to Kaggle unless the exact license allows it. Public derivatives and raw DICOM can have different access rules. No participant is contacted and no new human or animal data are collected; we need to ask the fair review committee whether secondary analysis needs prior review.

The main ethical risk is that readers mistake a model for medical advice. We need to label the study retrospective and exploratory, distinguish measured results from simulations and assumptions, and state that no clinical benefit or patient-specific treatment effect was established.

## Data Analysis
We need to split and resample by patient. Track A: report macro-F1 and paired donor-level differences, noting the seed-42 interval omits variation across training seeds. Track B: report paired patient-level Dice and equal-volume recurrence coverage for the 15 mm comparator (enhancing primary; the project's all-recurrence endpoint secondary). Track C: report simulated progression time, final burden, resistance, dose, and censoring together; present MGMT survival and log-rank results only as observational.

For the new model, simulate 200 datasets at each of three fixed noise scales (600 total) using observed visit schedules and first log-volumes; report 90% parameter-interval coverage. For prediction, average each patient's paired rolling one-step-ahead log-score differences against the training-fold median log-volume-slope baseline (primary) and persistence (secondary). Report the median patient effect and a 2,000-resample patient-bootstrap 95% interval. Higher scores favor the model; support requires a positive primary effect and an interval wholly above zero. Repeat with prespecified slope and process-noise prior alternatives; a sign change means the result is prior-sensitive. Poor recovery does not alone prove non-identifiability. Report calibration, exclusions, recovery and patient-data noise sensitivity results separately, and failure rules in tables and plots of observed and predicted volumes and uncertainty intervals.
## Bibliography

Abu Khalaf, N., Desjardins, A., Vredenburgh, J. J., & Barboriak, D. P. (2021). Repeatability of automated image segmentation with BraTumIA in patients with recurrent glioblastoma. *AJNR American Journal of Neuroradiology, 42*(6), 1080–1086. https://doi.org/10.3174/ajnr.A7071

Calabrese, E., Villanueva-Meyer, J., Rudie, J. D., Rauschecker, A. M., Baid, U., Bakas, S., Cha, S., Mongan, J., & Hess, C. P. (2022). The University of California San Francisco Preoperative Diffuse Glioma MRI (UCSF-PDGM) dataset. *Radiology: Artificial Intelligence, 4*(6), e220058. https://doi.org/10.1148/ryai.220058

Hormuth, D. A., II, Al Feghali, K. A., Elliott, A. M., Yankeelov, T. E., & Chung, C. (2021). Image-based personalization of computational models for predicting response of high-grade glioma to chemoradiation. *Scientific Reports, 11*, Article 8520. https://doi.org/10.1038/s41598-021-87887-4

Sonpatki, P., Park, H. J., Xing, Y. L., Han, K. Y., Schroeder, B. A., Yu, H. J., Kim, H. J., Chowdhury, T., Hwang, J. H., Nam, S. M., Byun, Y. H., Kang, H., Lee, J. H., Lee, S. T., Won, J. K., Kim, T. M., Choi, S. H., Ku, J.-L., Lee, S., … Shah, N. (2026). A spatially resolved human glioblastoma atlas reveals distinct cellular and molecular patterns of anatomical niches. *Nature Communications, 17*, Article 2951. https://doi.org/10.1038/s41467-026-69716-2

Zimmer, L., Weidner, J., Balcerak, M., Kofler, F., Krupa, M., Ezhov, I., Cepeda, S., Zhang, R. Z., Lowengrub, J., Menze, B., & Wiestler, B. (2026). PREDICT-GBM: A multicenter platform advancing personalized glioblastoma radiotherapy planning. *NPJ Digital Medicine, 9*, Article 686. https://doi.org/10.1038/s41746-026-03194-0

Dataset-specific citations and links are in the reproduction appendix.
