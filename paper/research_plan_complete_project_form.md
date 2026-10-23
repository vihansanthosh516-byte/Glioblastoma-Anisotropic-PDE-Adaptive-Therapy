# GBM Tracks A, B, and C: research plan

Tracks A, B, and C are equally important parts of this GBM project. Each track uses its own data and endpoint. Earlier exploratory work is disclosed, and previously used test groups are not presented as fresh confirmation. This plan states no results.

## Rationale

GBM varies across patients and within tumors. Molecular expression and spatial assays, longitudinal MRI, and treatment simulations describe different aspects of the disease. The repository establishes no patient-level link between the atlas and MRI cohorts, so results will not be fused into a patient-specific model. The project studies (A) donor-held-out cell-label prediction and matched spatial transfer where data permit, (B) MRI growth/volume forecasts and a separate recurrence-location comparison, and (C) adaptive and optimal-control/RL policy simulations plus a separate observational MGMT association. Simulated effects will not be described as clinical treatment effects.

## Research Question and Hypotheses

**Research question:** In separate GBM datasets and simulations, do donor-held-out molecular models, longitudinal MRI forecasts, and equal-budget treatment policies improve their prespecified evaluation measures over appropriate baselines?

**Track A hypothesis:** I hypothesize that random-cell splits overstate generalization: Core/Periphery models tested on held-out donors will have lower performance and greater donor-to-donor variation.

**Track B hypothesis:** I hypothesize that reaction-diffusion models will achieve higher patient-level Dice scores for later MRI tumor masks than persistence and geometric-growth baselines.

**Track C hypothesis:** I hypothesize that adaptive dosing and reinforcement-learning policies will lengthen simulated time to progression compared with fixed schedules under equal modeled budgets and prespecified settings.

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

1. Obtain UCSC multiomic-gbm single-cell and spatial files. Make a data inventory with release, assay, labels, donor/sample IDs, gene identifiers, and file terms. Resolve donor aliases from source records; do not count unresolved cell IDs as people.
2. Set the primary cell target to GBM Core versus GBM Periphery. Exclude NormalBrain and cells without a source-supported donor or target label. Summarize eligible donors and class counts before modeling.
3. Make donor-grouped training, validation, and test partitions. Keep all cells from each donor in one partition. Fit normalization, feature selection, dimensionality reduction, and model tuning using training donors only. Run a random-cell split only as a leakage diagnostic, not as the main performance estimate.
4. Compare the attention model, inductive expression-graph model, PCA/logistic-regression baseline, and fixed 40-gene anatomical and 348-gene Xenium panels on identical held-out donors. Map each panel to measured assay genes. Build graph edges from expression features only; exclude donor IDs, test labels, anatomical labels, and test-to-test edges.
5. Evaluate spatial transfer as a separate analysis. Apply the fixed panel to spatial specimens with measured Core/Periphery labels, hold out each specimen from fitting, and use only genes measured by that assay. Report donor/specimen metrics, uncertainty, class support, exclusions, and fit failures. Keep previously evaluated donors labeled exploratory.

### Track B — MRI growth, volume, and recurrence

1. Build a MU-Glioma-Post table with patient, diagnosis, visit date, scan interval, image spacing, mask source, and tumor-region labels. Check duplicates, missing visits, scan order, and GBM diagnosis. Convert voxel counts to cubic millimeters using image spacing and apply one written volume definition.
2. Keep each patient in one evaluation fold. For each eligible forecast, use earlier scans for fitting and reserve a later scan for scoring. Fit the reaction-diffusion model and make a later-mask forecast without using the reserved mask. Compare the same patient and time gap with persistence and geometric-growth forecasts.
3. Fit the measurement-aware log-volume model using training patients. Compare next-visit predictions with persistence, median-slope, and time-gap variance baselines. Before patient fitting, test parameter recovery on simulated data using the observed visit spacing. Record parameter bounds, failed fits, predictive intervals, calibration, and identifiability; stop inferential claims if recovery fails.
4. For recurrence location, use LUMIERE and RHUH-GBM MRI and recurrence masks as separate cohorts. Align image and mask coordinates in physical space, define the starting core and brain mask, make a 15 mm expansion, and threshold the model-risk map to the same volume. Score held-out coverage for enhancing recurrence and all recurrence separately. Record registration quality, exclusions, patient overlap, and cohort terms.

### Track C — Treatment policies and MGMT association

1. Keep the two-population simulator and optimal-control/RL simulator as separate experiments. Before either run, record the simulator version, equations, assumptions, parameter grid, thresholds, random seeds, treatment budget, and time horizon. Resolve differences between source comments and executable settings before generating comparisons.
2. In the two-population model, initialize paired runs with the same tumor and sensitive/resistant clone state. Compare adaptive dosing with the fixed five-days-on/23-days-off schedule across the prespecified growth and drug-effect settings. Save modeled burden over time, progression time, resistant fraction, and cumulative modeled dose for each run.
3. In the separate RL model, use the same held-out virtual cases and equal modeled budget to compare the Stupp-like schedule, heuristic policy, PPO policy, and oracle-reference policy. Keep training cases separate from evaluation cases. Treat the oracle as a reference bound, not as a usable treatment policy.
4. Analyze MGMT and progression in MU clinical data as a separate observational analysis. Verify MGMT code meanings, diagnosis, time origin, progression event, censoring, and missingness. Define the association and survival method before estimation; report eligible counts, exclusions, and uncertainty. Do not infer causal treatment benefit or use MGMT to set simulated resistance.

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
