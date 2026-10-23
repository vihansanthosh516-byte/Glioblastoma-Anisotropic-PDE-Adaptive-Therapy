# GBM Tracks A, B, and C: research plan

This document states the selected methods and planned measurements. It contains no project findings and assumes no direction or size of any effect. If any analysis was completed before this plan was submitted, disclose that work as the fair requires.

## Rationale

Glioblastoma (GBM) varies between patients and within a tumor. Molecular data, MRI, and treatment models describe different parts of this disease (Hormuth et al., 2021). The available datasets do not share verified patient identifiers, so we cannot combine them into one patient-level model. We need to study them in three related but separate tracks. The mathematical focus is a state-space model that represents latent postoperative volume change and MRI measurement error separately. Better measurement models may help research teams interpret longitudinal imaging; this retrospective study will not guide patient care.

## Research Question/Hypothesis(es)/Engineering Goal(s)

**Primary question:** How does a hierarchical state-space model of latent postoperative log-volume and MRI measurement error perform for next-scan prediction on held-out patients, compared with population-median log-volume slope and persistence baselines?

**Supporting questions:** How well do molecular features predict tumor-zone labels for donors excluded from training? How do selected reaction-diffusion and recurrence-risk methods compare with prespecified MRI baselines? How do adaptive and standard treatment schedules compare in simulations with assumed drug effects?

**Hypothesis:** We make no directional prediction. We will estimate differences and uncertainty using prespecified two-sided comparisons.

**Engineering goal:** Build an auditable method for longitudinal GBM analysis that reports assumptions, uncertainty, and conditions where the method fails. The tracks remain separate because their patient cohorts cannot be linked.

## Materials

We need a computer with 16 CPU cores and 34 GB RAM, Python, the project’s scientific-computing packages, and the listed public data. Sources are the UCSC multiomic-gbm atlas; MU-Glioma-Post v1; UCSF-PDGM; RHUH-GBM brain-extracted NIfTI files; and LUMIERE. Dataset sizes, citations, URLs, and use terms appear in the appendix. LUMIERE use terms conflict between its article and file record, and file-level reuse terms for atlas exports need checking. We must verify the exact files and permitted use before analysis or transfer. PREDICT-GBM and the Abu Khalaf repeatability study are literature references, not project data. No chemicals, organisms, laboratory equipment, or human contact are involved.

## Procedures

1. **Track A — Molecular tumor-zone prediction**
   1. Manifest the UCSC multiomic-gbm atlas; resolve donors from source metadata and exclude unresolved identities.
   2. Split by donor, keeping each donor's cells in one group. Fit the attention-gated classifier on training and validation donors.
   3. Compare with fixed PCA/logistic-regression and 200-gene logistic-regression baselines. Measure test-donor macro-F1 and donor-level uncertainty.

2. **Track B — MRI growth, recurrence, and measurement-aware modeling**
   1. Manifest eligible MRI visits, scan day, mask source, voxel size, labels, exclusions, and file terms. Do not link patients across cohorts.
   2. Compare patient-fold reaction-diffusion forecasts with persistence and geometric-growth forecasts. Compare selected preoperative recurrence-risk maps with an equal-volume 15 mm expansion; score enhancing recurrence first and all recurrence second.
   3. Fit the hierarchical state-space model to positive postoperative log-volumes, separating latent net trend, process noise, and scan error. Check parameter recovery on simulated visit schedules; compare rolling next-scan predictive scores with population-median-slope and persistence baselines.
   4. Report estimates, two-sided intervals, calibration, exclusions, sensitivity analyses, and prespecified fit failures.

3. **Track C — Treatment-policy simulation**
   1. Run the sensitive/resistant-cell simulator with MRI-derived growth inputs and assumed drug-kill rates. Compare adaptive and maximum-tolerated-dose schedules across the fixed ratio grid.
   2. Record modeled progression time, burden, resistance, dose, and censoring. Analyze MGMT-associated patient outcomes separately as observational data, not treatment effects.

## Risk Assessment

We use de-identified records and follow the terms for each file. We must avoid re-identification, guessed patient links, unauthorized redistribution, and transfer to a cloud GPU unless the exact data terms allow it. Public derivatives and raw DICOM can have different access rules. No participant is contacted and no new human or animal data are collected; we need to ask the fair review committee whether secondary analysis needs prior review.

The main ethical risk is readers treating computational results as medical advice. State that the work is retrospective, uses public data, and does not estimate patient-specific treatment benefit or guide care.

## Data Analysis

Use donor-level units for Track A and patient-level units for Tracks B and C. Report macro-F1 and donor-level uncertainty; paired Dice and equal-volume recurrence coverage with patient-level uncertainty; simulation summaries for every assumed drug-effect setting; and observational MGMT associations separately. For the state-space model, report parameter-recovery coverage, profile-likelihood diagnostics, predictive calibration, and paired patient-level log-score differences against both baselines. Use two-sided intervals and report estimates regardless of direction. Do not infer clinical benefit, causal treatment effects, or a linked molecular-to-MRI model. Report failed fits and eligibility shortfalls without relaxing rules after examining results.

## Bibliography

Abu Khalaf, N., Desjardins, A., Vredenburgh, J. J., & Barboriak, D. P. (2021). Repeatability of automated image segmentation with BraTumIA in patients with recurrent glioblastoma. *AJNR American Journal of Neuroradiology, 42*(6), 1080–1086. https://doi.org/10.3174/ajnr.A7071

Calabrese, E., Villanueva-Meyer, J., Rudie, J. D., Rauschecker, A. M., Baid, U., Bakas, S., Cha, S., Mongan, J., & Hess, C. P. (2022). The University of California San Francisco Preoperative Diffuse Glioma MRI (UCSF-PDGM) dataset. *Radiology: Artificial Intelligence, 4*(6), e220058. https://doi.org/10.1148/ryai.220058

Hormuth, D. A., II, Al Feghali, K. A., Elliott, A. M., Yankeelov, T. E., & Chung, C. (2021). Image-based personalization of computational models for predicting response of high-grade glioma to chemoradiation. *Scientific Reports, 11*, Article 8520. https://doi.org/10.1038/s41598-021-87887-4

Sonpatki, P., Park, H. J., Xing, Y. L., Han, K. Y., Schroeder, B. A., Yu, H. J., Kim, H. J., Chowdhury, T., Hwang, J. H., Nam, S. M., Byun, Y. H., Kang, H., Lee, J. H., Lee, S. T., Won, J. K., Kim, T. M., Choi, S. H., Ku, J.-L., Lee, S., … Shah, N. (2026). A spatially resolved human glioblastoma atlas reveals distinct cellular and molecular patterns of anatomical niches. *Nature Communications, 17*, Article 2951. https://doi.org/10.1038/s41467-026-69716-2

Zimmer, L., Weidner, J., Balcerak, M., Kofler, F., Krupa, M., Ezhov, I., Cepeda, S., Zhang, R. Z., Lowengrub, J., Menze, B., & Wiestler, B. (2026). PREDICT-GBM: A multicenter platform advancing personalized glioblastoma radiotherapy planning. *NPJ Digital Medicine, 9*, Article 686. https://doi.org/10.1038/s41746-026-03194-0

Dataset-specific citations and download links are in the reproduction appendix.
