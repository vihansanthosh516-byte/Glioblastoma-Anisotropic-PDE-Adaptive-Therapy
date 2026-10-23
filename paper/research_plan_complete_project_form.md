# GBM Tracks A, B, and C: research plan

This plan treats Tracks A, B, and C as equally important core parts of one GBM research project, while keeping their datasets and endpoints separate. Exploratory analyses and model development have already occurred; disclose that history without presenting earlier test sets as fresh confirmation. This plan states no project result values and makes no directional outcome prediction.

## Rationale

GBM varies across patients and within tumors. Molecular expression and spatial assays, longitudinal MRI, and treatment simulations describe different aspects of the disease. The repository establishes no patient-level link between the atlas and MRI cohorts, so results will not be fused into a patient-specific model. The project studies (A) donor-held-out cell-label prediction and matched spatial transfer where data permit, (B) MRI growth/volume forecasts and a separate recurrence-location comparison, and (C) adaptive and optimal-control/RL policy simulations plus a separate observational MGMT association. Simulated effects will not be described as clinical treatment effects.

## Research Question/Hypothesis(es)/Engineering Goal(s)

**Track A:** How do attention and bounded graph models compare with simple classifiers and an exact source-identified gene panel for donor-held-out GBM Core/Periphery labels? Where matched source-labeled spatial data and reuse terms permit, how do those panels transfer to held-out spatial specimens?

**Track B:** How do reaction-diffusion forecasts compare with persistence and geometric-growth baselines for later MRI masks, and how do measurement-aware models compare with fixed baselines for next-scan core volume? Separately, on rights-cleared source data, how does a model-based recurrence-location map compare with an equal-volume 15 mm geometric expansion?

**Track C:** Under specified model assumptions, how do sensitive/resistant adaptive schedules and optimal-control/RL policies compare with fixed and simple-rule schedules? Separately, what observational association exists between verified MGMT status and progression-time measures in the eligible clinical cohort?

No directional hypothesis is made. The engineering goal is an auditable comparison across three core tracks, with measured, fitted, assumed, simulated, and observational quantities distinguished.

## Materials

Use the UCSC multiomic-gbm atlas only after recording source release, files, donor mapping, and file-level terms. MU-Glioma-Post v1 has listed CC BY 4.0 terms for images, segmentations, and clinical files and will be used locally with citation. LUMIERE item records list CC0 while its paper limits use to noncommercial purposes; resolve exact file terms and use only within clearly permitted noncommercial scope. RHUH terms differ by file type: use only documented brain-extracted NIfTI derivatives if verified; raw DICOM is restricted. Exclude PREDICT-GBM and ALPTDG-derived packs unless exact upstream provenance and permission are documented; ALPTDG also requires its DUA. No patient contact or lab work is planned.

## Procedures

1. **Track A — Molecular and spatial analyses**
   1. Resolve atlas release, files, labels, aliases, and terms. Compare attention-gated, label-safe inductive graph, PCA/logistic, and an exact source-identified panel on donor-held-out Core/Periphery labels; fit preprocessing within training folds.
   2. If exact source panels, matched specimens, spatial labels, and terms permit, evaluate panels on patient-held-out spatial specimens. Report donor/specimen-level metrics and limits; otherwise stop that component.

2. **Track B — MRI growth and recurrence endpoints**
   1. For MU mask/core-volume forecasts, validate visits, dates, labels, and voxel sizes; compare the documented PDE solver with persistence and training-fold geometric growth, and compare a measurement-aware volume model with persistence, median slope, and a time-gap variance baseline.
   2. Separately compare recurrence coverage from source-built MRI risk maps with an equal-volume 15 mm geometric expansion. Use only rights-cleared files; exclude unverified PREDICT/ALPTDG derivatives. The expansion is a geometric comparator, not a clinical treatment plan.

3. **Track C — Policy simulation and clinical association**
   1. Compare the two-population adaptive schedule with the source-defined 5-days-on/23-days-off schedule across the documented effect/growth sweep. In the separate equal-budget RL model compare Stupp-like, heuristic, PPO, and oracle-reference policies under the documented budget and horizon. Keep all effects simulated.
   2. Separately estimate an observational MGMT/progression association only after verifying field coding, time origin, event/censor definitions, and missingness. Do not map MGMT to simulated resistance or infer treatment benefit.

## Risk Assessment

Use de-identified records and comply with exact file terms and agreements. Keep data local unless the source explicitly permits the intended hosting. The repository records prior cloud packaging that includes ALPTDG-derived labels; its DUA and transfer permission have not been established here. Do not re-identify participants or infer patient links across datasets. Simulations and observational associations do not guide clinical care.

## Data Analysis

Track A uses donor-held-out folds, donor-wise macro-F1/class metrics, and donor-level uncertainty. Track B reports patient-paired Dice for mask forecasting; volume forecasts use patient-level predictive log scores, calibration, intervals, and fit/identifiability diagnostics. The recurrence endpoint reports enhancing-tumor and all-recurrence coverage separately against the equal-volume 15 mm comparator. Track C reports every simulation setting and policy comparison; MGMT associations remain separate observational estimates. Use common evaluation units within comparisons, report failures and exclusions, and do not present reused cohorts as fresh confirmation. Keep all three tracks distinct.

## Bibliography

Abu Khalaf, N., Desjardins, A., Vredenburgh, J. J., & Barboriak, D. P. (2021). Repeatability of automated image segmentation with BraTumIA in patients with recurrent glioblastoma. *AJNR American Journal of Neuroradiology, 42*(6), 1080–1086. https://doi.org/10.3174/ajnr.A7071

Calabrese, E., Villanueva-Meyer, J., Rudie, J. D., Rauschecker, A. M., Baid, U., Bakas, S., Cha, S., Mongan, J., & Hess, C. P. (2022). The University of California San Francisco Preoperative Diffuse Glioma MRI (UCSF-PDGM) dataset. *Radiology: Artificial Intelligence, 4*(6), e220058. https://doi.org/10.1148/ryai.220058

Hormuth, D. A., II, Al Feghali, K. A., Elliott, A. M., Yankeelov, T. E., & Chung, C. (2021). Image-based personalization of computational models for predicting response of high-grade glioma to chemoradiation. *Scientific Reports, 11*, Article 8520. https://doi.org/10.1038/s41598-021-87887-4

Sonpatki, P., Park, H. J., Xing, Y. L., Han, K. Y., Schroeder, B. A., Yu, H. J., Kim, H. J., Chowdhury, T., Hwang, J. H., Nam, S. M., Byun, Y. H., Kang, H., Lee, J. H., Lee, S. T., Won, J. K., Kim, T. M., Choi, S. H., Ku, J.-L., Lee, S., … Shah, N. (2026). A spatially resolved human glioblastoma atlas reveals distinct cellular and molecular patterns of anatomical niches. *Nature Communications, 17*, Article 2951. https://doi.org/10.1038/s41467-026-69716-2

Zimmer, L., Weidner, J., Balcerak, M., Kofler, F., Krupa, M., Ezhov, I., Cepeda, S., Zhang, R. Z., Lowengrub, J., Menze, B., & Wiestler, B. (2026). PREDICT-GBM: A multicenter platform advancing personalized glioblastoma radiotherapy planning. *NPJ Digital Medicine, 9*, Article 686. https://doi.org/10.1038/s41746-026-03194-0

Dataset-specific citations and download links are in the reproduction appendix.
