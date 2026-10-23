# GBM Tracks A, B, and C: research plan

Tracks A, B, and C are equally important parts of this GBM project. Each track uses its own data and endpoint. Earlier exploratory work is disclosed, and previously used test groups are not presented as fresh confirmation. This plan states no results.

## Rationale

GBM varies across patients and within tumors. Molecular expression and spatial assays, longitudinal MRI, and treatment simulations describe different aspects of the disease. The repository establishes no patient-level link between the atlas and MRI cohorts, so results will not be fused into a patient-specific model. The project studies (A) donor-held-out cell-label prediction and matched spatial transfer where data permit, (B) MRI growth/volume forecasts and a separate recurrence-location comparison, and (C) adaptive and optimal-control/RL policy simulations plus a separate observational MGMT association. Simulated effects will not be described as clinical treatment effects.

## Research Question/Hypothesis(es)/Engineering Goal(s)

**Track A hypothesis:** Core/Periphery models tested on held-out donors will have lower performance and greater donor-to-donor variation than random-cell splits suggest.

**Track B hypothesis:** On held-out patients, reaction-diffusion forecasts will predict later tumor masks more accurately than persistence and geometric-growth baselines.

**Track C hypothesis:** At equal modeled treatment budget, adaptive schedules will delay simulated progression relative to fixed schedules across prespecified growth and drug-effect settings.

The engineering goal is to build and compare reproducible methods across all three tracks while keeping measured, fitted, assumed, simulated, and observational quantities distinct.

## Materials

- UCSC multiomic-gbm atlas: GBM single-cell and spatial transcriptomic data.
- MU-Glioma-Post v1: longitudinal postoperative glioma MRI, segmentations, and clinical data.
- LUMIERE: longitudinal glioblastoma MRI and expert RANO evaluations; noncommercial-use scope.
- RHUH-GBM: preoperative, early postoperative, and recurrence MRI; use brain-extracted NIfTI derivatives, not restricted raw DICOM.
- Desktop computer: 16 CPU cores, 34 GB RAM, no local GPU.
- Python and the project’s existing analysis environment and scripts.
- No laboratory chemicals, biological agents, patient contact, or new human-subject data collection.

## Procedures

### Track A — Molecular and spatial prediction

1. Download UCSC multiomic-gbm files. Record release and donor IDs. Define Core versus Periphery; exclude NormalBrain and unresolved donors.
2. Split by donor. Fit preprocessing and tuning on training donors only. Compare attention, inductive expression-graph, PCA/logistic, and source-panel models. Build graph edges from expression, not IDs or test labels.
3. Test panels on held-out spatial specimens with measured labels. Report donor/specimen metrics, uncertainty, failures. Mark reused donors exploratory.

### Track B — MRI growth, volume, and recurrence

1. Build a MU patient/visit table. Verify diagnosis, dates, voxel spacing, masks, and volume definitions; convert voxels to physical volume.
2. Keep each patient in one fold. Fit reaction-diffusion models on earlier scans; forecast later masks and compare with persistence and geometric growth.
3. Fit the measurement-aware volume model on training patients. Compare with persistence, median-slope, and time-gap variance baselines. Check calibration, intervals, fit failures, and identifiability after simulation recovery checks.
4. Register LUMIERE/RHUH MRI and recurrence masks. Compare model-risk regions with equal-volume 15 mm expansions. Score held-out enhancing-tumor and all-recurrence coverage; record exclusions and registration quality.

### Track C — Treatment policies and MGMT association

1. Freeze simulator, assumptions, thresholds, parameter grid, seeds, horizon, and budget. Resolve code-setting conflicts first.
2. Pair initial tumor/clone states in adaptive and fixed five-days-on/23-days-off simulations. Separately compare Stupp-like, heuristic, PPO, and oracle-reference policies on held-out RL cases with equal budget.
3. Record burden, modeled progression, resistant fraction, and exposure. Separately estimate MU MGMT/progression association after verifying codebook, time origin, event, censoring, and missingness. Report uncertainty; make no causal claim.

## Risk Assessment

Use de-identified data under each source’s terms. Keep files local; do not upload controlled or unclear-license data. Do not identify participants or join patients across cohorts. Simulations and observational results are not medical advice.

## Data Analysis

Track A uses donor-held-out folds, donor-wise macro-F1/class metrics, and donor-level uncertainty. Track B reports patient-paired Dice for mask forecasting; volume forecasts use patient-level predictive log scores, calibration, intervals, and fit/identifiability diagnostics. The recurrence endpoint reports enhancing-tumor and all-recurrence coverage separately against the equal-volume 15 mm comparator. Track C reports every simulation setting and policy comparison; MGMT associations remain separate observational estimates. Use common evaluation units within comparisons, report failures and exclusions, and do not present reused cohorts as fresh confirmation. Keep all three tracks distinct.

## Bibliography

Abu Khalaf, N., Desjardins, A., Vredenburgh, J. J., & Barboriak, D. P. (2021). Repeatability of automated image segmentation with BraTumIA in patients with recurrent glioblastoma. *AJNR American Journal of Neuroradiology, 42*(6), 1080–1086. https://doi.org/10.3174/ajnr.A7071

Calabrese, E., Villanueva-Meyer, J., Rudie, J. D., Rauschecker, A. M., Baid, U., Bakas, S., Cha, S., Mongan, J., & Hess, C. P. (2022). The University of California San Francisco Preoperative Diffuse Glioma MRI (UCSF-PDGM) dataset. *Radiology: Artificial Intelligence, 4*(6), e220058. https://doi.org/10.1148/ryai.220058

Hormuth, D. A., II, Al Feghali, K. A., Elliott, A. M., Yankeelov, T. E., & Chung, C. (2021). Image-based personalization of computational models for predicting response of high-grade glioma to chemoradiation. *Scientific Reports, 11*, Article 8520. https://doi.org/10.1038/s41598-021-87887-4

Sonpatki, P., Park, H. J., Xing, Y. L., Han, K. Y., Schroeder, B. A., Yu, H. J., Kim, H. J., Chowdhury, T., Hwang, J. H., Nam, S. M., Byun, Y. H., Kang, H., Lee, J. H., Lee, S. T., Won, J. K., Kim, T. M., Choi, S. H., Ku, J.-L., Lee, S., … Shah, N. (2026). A spatially resolved human glioblastoma atlas reveals distinct cellular and molecular patterns of anatomical niches. *Nature Communications, 17*, Article 2951. https://doi.org/10.1038/s41467-026-69716-2

Zimmer, L., Weidner, J., Balcerak, M., Kofler, F., Krupa, M., Ezhov, I., Cepeda, S., Zhang, R. Z., Lowengrub, J., Menze, B., & Wiestler, B. (2026). PREDICT-GBM: A multicenter platform advancing personalized glioblastoma radiotherapy planning. *NPJ Digital Medicine, 9*, Article 686. https://doi.org/10.1038/s41746-026-03194-0

Dataset-specific citations and download links are in the reproduction appendix.
