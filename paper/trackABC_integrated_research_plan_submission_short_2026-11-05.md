# Research Plan

## Project Title

**How Will Multimodal Information Change Simulated Glioblastoma Treatment Scheduling? A Mathematical Study of Molecular Signals, Longitudinal MRI, and Drug Suppression**

## Rationale

Glioblastoma varies across tumor regions and can change during treatment. MRI can show changes in tumor-associated tissue over time, while molecular measurements can describe differences in cell states. Mathematical models can represent tumor growth, treatment suppression, and drug-resistant cells. Prior research has studied MRI-based recurrence risk and mathematical treatment scheduling, so this project will not claim that any one of those ingredients is new by itself.

This project will connect three research tracks without treating their datasets as if they came from the same people. Track A will supply frozen donor-level molecular prediction summaries. Track B will examine longitudinal MRI transitions and ADC measurements in the restricted Duke dataset. Track C will model drug-sensitive and drug-resistant tumor populations and simulated treatment schedules. The central mathematical contribution will be a held-out test of when combining the molecular and MRI/ADC signals changes simulated progression time and treatment-course use, under different assumptions about how those signals relate to drug resistance. The patient MRI analysis and the simulated treatment results will remain separate. The project will not provide a treatment recommendation for any person.

## Research Question, Hypotheses, Engineering Goals, and Expected Outcomes

### Research Questions

**Central mathematical question:** Across four predeclared assumptions about how molecular and MRI/ADC signals relate to latent drug resistance, how will a fused policy’s 730-day restricted mean time to simulated progression and treatment-course use compare with the best development-selected fixed schedule and the best development-selected single-modality policy on held-out virtual patients?

1. **Patient MRI question:** Will newly appearing nonenhancing peritumoral FLAIR regions receive a next-visit necrotic or enhancing core label more often than persistent FLAIR regions at the same distance from the current core?
2. **ADC question:** Will current ADC and its change from the prior visit improve patient-held-out prediction of that transition?
3. **Treatment-simulation question:** Under the moderate-concordance assumption, will the fused policy improve simulated progression time without greater treatment-course use than both selected comparators?
4. **Robustness question:** How will the paired differences change under independent, weak, moderate, and reversed signal-dependence assumptions?

### Hypotheses

1. Newly appearing FLAIR regions will have a higher next-visit core-label conversion rate than distance-matched persistent regions.
2. Adding ADC features will lower held-out patient-macro Brier score compared with the MRI and clinical baseline.
3. Under moderate concordance, the fused policy will have greater 730-day restricted mean time to simulated progression than both the best fixed schedule and best single-modality policy, without using more courses, under both progression definitions and both calibration families.
4. I hypothesize that the simulated advantage of fusion will be smaller under independent or weak dependence than under moderate concordance; I will test whether it will disappear or reverse under the reversed molecular-to-resistance assumption.

The study will report results that support, fail to support, or contradict these hypotheses.

### Engineering Goals and Expected Outcomes

The project will build an auditable patient-level MRI pipeline, quantify ADC’s additional predictive value, and connect the resulting signal distributions to a calibrated two-population drug-suppression simulator. It will compare fused, fixed, single-modality, and no-information policies on the same virtual patients. The analysis will test whether fusion improves, leaves unchanged, or worsens performance under different assumptions. I will label each finding as patient MRI evidence or simulation evidence.

## Materials

- I will use a Dell 16 Plus DB16250 computer running Windows 11 with 32 GB RAM.
- I will use Python 3.14.4 and NumPy, SciPy, pandas, scikit-learn, NiBabel, and CPU-only PyTorch. Before each locked analysis, I will record and freeze the software versions, package files, hardware settings, and code hashes.
- **Track A:** I will use frozen donor-level A5 summaries derived from Yu et al.’s single-cell RNA-sequencing data (GEO GSE117891), Table S1, and GBmap annotation metadata. I will not reopen expression values or retune or rescore A5.
- **Track B:** the restricted Duke longitudinal glioblastoma MRI dataset (Zenodo record 20548528), including permitted visit metadata, contrast-enhanced T1 MRI, FLAIR MRI, ADC maps, segmentation masks, and vessel landmarks where available. I will use these data only under the accepted Duke terms and required school or competition approvals.
- **Track C:** I will use the project’s two-population tumor-growth and drug-suppression simulator, plus published aggregate outcomes from RTOG 0525 and GEINO 14-01. Trial summaries will constrain simulated populations; they will not represent individual patient responses.
- No physical materials will be used. I will not use chemicals, bacteria, cell cultures, animals, or experimental drugs.

## Procedures

1. **Confirm permission and protect data.** Before analysis, I will verify that the Duke agreement and school or competition conditions allow the planned processing and reporting. If terms do not permit a step, I will omit it before data access. I will verify the archive checksum and NIfTI headers before reading voxel arrays, and I will stop on a failed hash or invalid header. I will store restricted data and patient-level outputs locally in access-controlled storage. I will not identify patients, contact them, or upload restricted files or predictions.

2. **Select the Duke cohort from metadata.** Before image or mask reads, I will use permitted metadata to select glioblastoma records whose released IDH field is “Negative.” I will report this only as the dataset-coded field, not as comprehensive IDH1/2 wildtype status. I will form windows from every eligible set of three adjacent visits with T1-c, FLAIR, and segmentation data; I will not skip visits. I will sort by relative scan date, require dates to be finite and increasing, record scan intervals and treatment events, and exclude unresolved windows before voxel reads. I will check prior-use overlap where permitted and will not claim independent confirmation if overlap cannot be established. I will retain surgery-between-visits windows in the primary analysis and report their exclusion as a sensitivity analysis. The study will estimate changes under observed care and will be retrospective.

3. **Register images and define labels and support.** I will rigidly align prior- and next-visit T1-c images to the middle-visit T1-c using image intensities only, then apply each transform to the other modalities and masks. I will accept registrations only if normalized cross-correlation is at least 0.80, the largest absolute Euler angle is at most 10 degrees, and translation length is at most 10 mm; I will not retry or manually adjust failures. I will test registration on synthetic known transforms and report released vessel-landmark errors where available. I will resample to a 2 mm isotropic grid. Labels 1 and 3 will define the operational core, label 2 the FLAIR abnormality, and label 4 the cavity; I will not treat label 0 as normal tissue. I will define support from the prior and middle visits only, using their MRI foreground and released labels. The next-visit image will be used only for evaluation alignment, not for support, features, tuning, or forecast generation. Target-core voxels outside support will count as misses.

4. **Measure the primary MRI transition.** Within the middle-visit FLAIR region, I will exclude a 2 mm margin around the core and cavity. I will classify tissue as newly appearing if it lies outside a 2 mm dilation of prior-visit FLAIR, and persistent if it lies inside. I will use 2–10, >10–20, >20–30, and >30 mm distance bands from the middle-visit core; next-visit labels 1 or 3 will count as recorded core. A window will qualify only if at least one band contains 10 or more voxels in each history-defined region; I will apply this rule without viewing the next-visit labels. I will calculate the new-minus-persistent conversion difference by band, weight each band by the smaller region count, average bands and windows within each patient, and weight patients equally. I will use 10,000 patient-cluster bootstrap resamples. I will stop inference below 20 eligible patients; with 20–34 patients or inadequate precision, I will report feasibility only. I will consider the hypothesis supported only with at least 35 patients, interval half-width no greater than 0.10, and a lower confidence bound above zero.

5. **Test ADC and run secondary forecasts.** I will compare a prespecified MRI/clinical logistic model with the same model plus normalized ADC and ADC change, using five patient-grouped folds and training patients only for feature processing and model selection. I will compare held-out patient-macro Brier scores and will not describe this retrospective analysis as a prospective tool. Separately, I will compare the reaction-diffusion PDE with copy-last, fixed growth, volume-matched persistence, FLAIR-front, two-history shape, PDE-only, and a time-conditioned FiLM U-Net. All arms will use the same patients, folds, support, and targets; fitted settings will use training patients only. A candidate will count as better only if every adjusted lower confidence bound for its paired Dice difference versus a required mask comparator is above zero. A count-only model will be evaluated only with Poisson deviance, count MAE, and Brier score for any new piece, never as a mask.

6. **Build the molecular and drug-suppression inputs.** I will use frozen Track A A5 donor-level balanced-accuracy summaries only as assumed reliability values for a noisy simulated molecular signal. I will not reopen expression data or call A5 a malignant-cell-specific or drug-response marker. I will model sensitive and resistant tumor populations, growth, treatment suppression, and assumed switching or reversion. I will calibrate the Track C simulator against published aggregate trial outcomes, not individual patient responses. I will require at least 30 worlds per calibration family to pass the predeclared calibration checks before policy comparisons.

7. **Construct virtual patients and compare policies.** I will combine the separate Track A, Track B, and Track C distributions without matching people across datasets. I will use independent, weak, moderate, and reversed molecular-to-resistance dependence scenarios while keeping the marginal inputs fixed. The fused policy will use the moderate-concordance prior in every scenario and will not receive the hidden scenario or tumor parameters. I will compare the fused belief-state controller with the best development-selected fixed schedule and single-modality policy, plus Track-C-only, late-fusion, Q-learning, no-information, and signal-shuffle controls. All policies will receive the same virtual patients, scan schedules, noise, and a maximum of 12 courses. I will choose settings and comparators on development simulations, freeze the protocol and seeds, and evaluate 1,000 held-out worlds per dependence scenario once.

8. **Apply stop rules and report results.** I will use the predeclared fallback for expected model failures and record each reason; I will not silently remove patients. Missing required outputs, nonfinite values, or changed protocol hashes will make the comparison incomplete or stop it. Before scoring, I will run a full-size synthetic timing pilot and stop or amend the plan if projected runtime exceeds 1,000 hours. I will save predictions and settings before opening held-out targets, hash the protocol, code, folds, inputs allowed by the agreement, and cohort, and report patient MRI findings separately from simulated treatment outcomes. I will report positive, null, or negative results and will not present simulation as treatment evidence for real people.

## Risk Assessment

The project will use de-identified MRI records under restricted access. Risks will include accidental disclosure and use beyond the data agreement. I will use access-controlled local storage, will not attempt re-identification, and will share only outputs allowed by the agreement and school or competition rules. No person will receive a diagnosis or treatment recommendation, and no drug will be administered. I will describe automated labels as MRI labels, not verified tumor boundaries.

## Data Analysis

For the patient transition analysis, I will report the equal-patient difference in next-visit core-label conversion between newly appearing and persistent FLAIR regions, with a 95% patient-cluster bootstrap interval. For ADC, I will report paired held-out patient-macro Brier scores. For the secondary mask forecasts, I will report patient-macro Dice and simultaneous paired intervals for each required comparator.

For simulation, I will calculate 730-day restricted mean time to progression under two outcomes: first reaching 1.2 times starting burden, and a scan-based RANO-like rule applied every 84 days. The second outcome will be a model approximation, not a formal RANO assessment. For each comparator, I will define ΔRMST as fused-policy RMST minus comparator RMST, and Δcourses as fused-policy course use minus comparator course use. Positive ΔRMST will favor fusion; Δcourses at or below zero will indicate no greater course use.

In the moderate-concordance scenario, I will compare fusion with the development-selected fixed schedule and single-modality policy using simultaneous 95% max-|t| bootstrap intervals across 16 contrasts: 8 ΔRMST contrasts and 8 endpoint-specific Δcourse contrasts across both progression definitions and both calibration families. The primary simulation rule will require every adjusted ΔRMST lower bound to exceed zero, every adjusted Δcourse upper bound to be at or below zero, and successful calibration in both families. I will report the same outcomes separately under all four dependence scenarios. Bootstrap resampling will preserve each simulation world and its related observations as one unit.

## Bibliography

Akbari, H., Macyszyn, L., Da, X., Bilello, M., Wolf, R. L., Martinez-Lage, M., Biros, G., Alonso-Basanta, M., O’Rourke, D. M., & Davatzikos, C. (2016). Imaging surrogates of infiltration obtained via multiparametric imaging pattern analysis predict subsequent location of recurrence of glioblastoma. *Neurosurgery, 78*(4), 572–580. https://doi.org/10.1227/NEU.0000000000001202

Balaña, C., Vaz, M. A., Sepúlveda, J. M., Mesia, C., del Barco, S., Pineda, E., Muñoz-Langa, J., Estival, A., de las Peñas, R., Fuster, J., Gironés, R., Navarro, L. M., Gil-Gil, M., Alonso, M., Herrero, A., Peralta, S., Olier, C., Perez-Segura, P., Covela, M., … Carrato, C. (2020). A phase II randomized, multicenter, open-label trial of continuing adjuvant temozolomide beyond 6 cycles in patients with glioblastoma (GEINO 14-01). *Neuro-Oncology, 22*(12), 1851–1861. https://doi.org/10.1093/neuonc/noaa107

Criscuolo, E., Calabrese, E., Zhang, Z., Wang, Y., Xiong, Y., & Yang, D. (2026). *A longitudinal glioblastoma MRI dataset with anatomical landmark pairs for DIR validation* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.20548528

Ellingson, B. M., Wen, P. Y., & Cloughesy, T. F. (2017). Modified criteria for radiographic response assessment in glioblastoma clinical trials. *Neurotherapeutics, 14*(2), 307–320. https://doi.org/10.1007/s13311-016-0507-6

Gilbert, M. R., Wang, M., Aldape, K. D., Stupp, R., Hegi, M. E., Jaeckle, K. A., Armstrong, T. S., Wefel, J. S., Won, M., Blumenthal, D. T., Mahajan, A., Schultz, L. M., Erridge, S. C., Brown, P. D., Chakravarti, A., Pugh, S., Zadeh, G., Fisher, B., Sulman, E. P., … Mehta, M. P. (2013). Dose-dense temozolomide for newly diagnosed glioblastoma: A randomized phase III clinical trial. *Journal of Clinical Oncology, 31*(32), 4085–4091. https://doi.org/10.1200/JCO.2013.49.6968

Hoff, B. A., Lemasson, B., Chenevert, T. L., Luker, G. D., Tsien, C. I., Amouzandeh, G., Johnson, T. D., & Ross, B. D. (2021). Parametric response mapping of FLAIR MRI provides an early indication of progression risk in glioblastoma. *Academic Radiology, 28*(12), 1711–1720. https://doi.org/10.1016/j.acra.2020.08.015

Liu, Z., Zhang, J., Hong, L., Nie, Q., & Sun, X. (2025). Multiscale mathematical model-informed reinforcement learning optimizes combination treatment scheduling in glioblastoma evolution. *Science Advances, 11*(32), eadv3316. https://doi.org/10.1126/sciadv.adv3316

Moore-Palhares, D., Lawrence, L. S. P., Myrehaug, S., Stewart, J., Detsky, J., Tseng, C.-L., Chen, H., Dinakaran, D., Maralani, P., Ruschin, M., Zhang, B., Perry, J., Lim-Fat, M. J., Sahgal, A., Soliman, H., & Lau, A. Z. (2025). Temporal apparent diffusion coefficient changes during chemoradiation: An imaging biomarker for tumor response monitoring and spatial recurrence prediction in glioblastoma. *International Journal of Radiation Oncology, Biology, Physics, 122*(3), 592–604. https://doi.org/10.1016/j.ijrobp.2025.03.028

Stupp, R., Mason, W. P., van den Bent, M. J., Weller, M., Fisher, B., Taphoorn, M. J. B., Belanger, K., Brandes, A. A., Marosi, C., Bogdahn, U., Curschmann, J., Janzer, R. C., Ludwin, S. K., Gorlia, T., Allgeier, A., Lacombe, D., Cairncross, J. G., Eisenhauer, E., & Mirimanoff, R. O. (2005). Radiotherapy plus concomitant and adjuvant temozolomide for glioblastoma. *New England Journal of Medicine, 352*(10), 987–996. https://doi.org/10.1056/NEJMoa043330

Yu, K., Hu, Y., Wu, F., Guo, Q., Qian, Z., Hu, W., Chen, J., Wang, K., Fan, X., Wu, X., Rasko, J. E., Fan, X., Iavarone, A., Jiang, T., Tang, F., & Su, X.-D. (2020). Surveying brain tumor heterogeneity by single-cell RNA-sequencing of multi-sector biopsies. *National Science Review, 7*(8), 1306–1318. https://doi.org/10.1093/nsr/nwaa099

## Additional Project Information

This project will use de-identified human MRI records, previously generated molecular summaries, and computer simulations. I will follow the accepted Duke data-use terms and school or competition conditions before accessing or analyzing restricted data. I will not recruit participants, contact patients, access direct identifiers, intervene in care, or administer treatment. The project will not be used as a clinical decision tool. I will provide prior-project materials if the competition requires them.
