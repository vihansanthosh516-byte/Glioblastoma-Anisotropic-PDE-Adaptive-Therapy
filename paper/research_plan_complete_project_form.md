# Research Plan: Three Separate Glioblastoma Studies

## Scope and planned-work wording

This project will address three separate computational questions about glioblastoma (GBM). Track A will study donor-held-out molecular classification, Track B will study tumor-mask forecasting from serial MRI, and Track C will compare treatment policies in a mathematical simulation. The tracks will use different data and evaluation units. They will not be joined into a patient-level model or described as a validated digital twin.

The procedures below will be written as future actions. Before submission, I will disclose earlier data access, model development, and exploratory analyses in the fair application’s prior-work section. This plan will not present work already performed as unperformed, and it will not report results or promise that any hypothesis will be supported.

## Rationale

The background will cite single-cell GBM Core/Periphery classification and cancer-cell domain-generalization work (Lopes & Vinga, 2020; Fischer et al., 2024; Zhong et al., 2024), inverse and image-based glioma growth models (Subramanian et al., 2020; Hormuth et al., 2021), and the recent PREDICT-GBM recurrence-prediction and radiotherapy-planning benchmark (Zimmer et al., 2026). It will also cite adaptive-therapy work and GBM treatment-scheduling work that combines multiscale mathematical modeling with reinforcement learning (Zhang et al., 2017; West et al., 2023; Liu et al., 2025). Track A will evaluate generalization to held-out donors, Track B will compare PDE-based shape forecasts with a matched-volume geometric forecast, and Track C will compare policies only under stated simulator assumptions. Since related methods already exist, the project will assess these narrower comparisons and will not claim that a favorable score proves uniqueness or that simulation demonstrates treatment benefit.

## Research questions and hypotheses

**Overall question:** Will the planned method in each track add measurable value over a prespecified, directly comparable baseline on its own evaluation units?

**Track A question:** Will sparse, donor-balanced expression classification improve donor-level Core-versus-Periphery accuracy over ordinary logistic regression when every test donor is excluded from training?

**Track A hypothesis to test:** The study will test whether the sparse, donor-balanced model will have higher mean held-out-donor balanced accuracy than logistic regression. If the data, rights, or donor-count gates fail, the study will report a feasibility stop instead of a model comparison.

**Track B question:** Will an isotropic reaction-diffusion model forecast a later tumor-core mask more accurately than simple forecasts on the same held-out patients and observed scan intervals?

**Track B hypothesis to test:** The study will test whether the PDE forecast will have higher patient-macro Dice than the prespecified matched-volume, no-PDE baseline. The study will treat all MU results as exploratory and will make an external claim only if an untouched cohort and every access, source, and quality gate pass.

**Track C question:** Under stated simulation assumptions, will a resistance-aware policy improve lower-tail simulated time to progression over a simple dose-limited policy?

**Track C hypothesis to test:** The study will test whether the policy will increase the worst-noise-level mean of the lowest 10% of simulated progression times. No simulated result will be described as evidence that a real drug works or that a treatment will help a patient.

## Materials

- The UCSC Cell Browser multiomic-gbm release ([Cell Browser](https://cells.ucsc.edu/?ds=multiomic-gbm)) will be considered for Track A only after I verify the exact downloaded files, donor labels, and reuse terms. Public availability alone will not serve as a reuse license.
- MU-Glioma-Post version 1 ([TCIA collection](https://www.cancerimagingarchive.net/collection/mu-glioma-post/)) will be considered for Track B development. I will cite the TCIA collection page and its listed CC BY 4.0 terms; the prior-exposure manifest will classify all MU evaluation as exploratory.
- The restricted Duke longitudinal GBM collection ([Zenodo record](https://zenodo.org/records/20548528)) will be considered as an external test only after I verify that its visit sequence and tumor masks fit this forecast question and access permission, school-use approval, and protocol gates are complete. I will not download or score it while access remains pending.
- I will cite PREDICT-GBM as related work, but I will not treat its TUM-GBM, LUMIERE, or RHUH source cohorts as new independent test data because the repository will identify those sources as previously examined in this project.
- CFB-GBM version 4 and UCSD-PTGBM will not be counted as Track B PDE validation in this plan. Any annotation-change analysis will remain a separate exploratory side study.
- Track C will use generated virtual scenarios only. Every parameter will be identified as an assumption.
- A Windows computer with 8 CPU cores and 31.6 GiB RAM will be used. No local GPU will be assumed. Patient data will not be uploaded to a cloud service.
- Python and the scientific-computing, machine-learning, and medical-image packages already used by the project will be recorded by version before analysis.
- The project will use no chemicals, biological agents, laboratory procedures, patient contact, or new participant recruitment. I will obtain the school’s required review before any restricted-data use.

## Procedures

### Shared procedure

1. I will create a source and prior-exposure manifest. It will record each dataset version, exact files, checksums, access terms, identifiers, labels, prior use, and planned role. If prior exposure cannot be reconstructed, I will treat the whole source cohort as previously examined.
2. I will keep the three tracks separate and will use donor, patient, or simulated parameter scenario as the independent evaluation unit. Cells, voxels, scans, repeated noise draws, and model particles will not be counted as independent people or scenarios.
3. Before opening any evaluation outcomes, I will freeze each track’s inclusion rules, primary measure, comparator, split, model settings, random seeds, and failure rules. Preprocessing and tuning will use training or development units only.
4. I will report the number of eligible units, exclusions, missing inputs, failed fits, and uncertainty. A result will count as evidence of an advantage only if its prespecified paired 95% interval favors the candidate and all relevant data gates pass. An interval that includes no difference will be reported as inconclusive.
5. I will preserve old outputs and identify prior runs as exploratory, post hoc, or simulated where required. No new split or metric on a previously examined cohort will be called independent confirmation.

### Track A — donor-held-out molecular classification

1. I will inventory the exact multiomic-gbm files and verify the meanings and sources of Core and Periphery labels, sample-to-donor identifiers, assays, measured genes, and permitted reuse. I will not infer donor identity from a filename or merge aliases without source documentation.
2. I will use source-labeled GBM Core and GBM Periphery cells as the only primary target. I will exclude records with unresolved donor identity or target meaning. Spatial labels will not be a scored endpoint unless independent, specimen-linked labels and permitted use are established.
3. I will count eligible donors and class coverage before fitting. Descriptive modeling will require at least six donors with both classes. I will make no inferential superiority claim with fewer than ten eligible donors. Perfect confounding of class with donor, library, or batch will stop the predictive comparison and will be reported as a limitation.
4. I will create outer leave-one-donor-out folds and inner donor-grouped folds. All cells from a donor will stay in one fold. Normalization, gene filtering, feature selection, and regularization selection will be fitted inside training folds only.
5. I will fit one sparse logistic model with equal total training weight per donor. I will compare it with ordinary logistic regression using the same folds and a fixed, small inner-fold tuning budget. Cell-type-only and library/batch-only predictions will be used as leakage diagnostics, not as competing biological models.
6. I will calculate balanced accuracy separately for each held-out donor and use the equal-weight donor mean as the primary outcome. I will calculate a paired donor-level bootstrap interval for the difference from logistic regression. I will show per-donor results and will label all analysis on previously examined atlas records exploratory.
7. If source terms, donor identity, or labels cannot be verified, I will not score the atlas. If too few donors qualify, I will report the data audit and feasibility limit rather than adding more model variants.

### Track B — PDE tumor-growth forecasting

1. I will keep reaction-diffusion tumor-growth forecasting as the primary Track B question. After I verify the release label map, the target will be the source-defined tumor core from labels 1 and 3. I will not call a change in released annotations biological growth, recurrence, or progression.
2. I will use the first two eligible masks as history and the immediate third mask as the target. The patient will remain the unit of analysis, and no patient’s visits will be divided across training and test folds. I will generate each forecast at day 90 after the second history scan and include only patients whose third scan falls 76–104 days after that origin. I will report the actual acquisition interval and describe this as a retrospective nominal 90-day target-window forecast, not as an exact-day or schedule-known-at-origin forecast.
3. MU-Glioma-Post will serve as exploratory development data. The prior-exposure manifest will record its earlier project use. Duke will serve as the only external candidate after access, school approval, overlap, visit chronology, segmentation support, and longitudinal registration checks pass. Before opening target masks, I will freeze the patient list, source hashes, target definition, forecast code, and comparison protocol. Any unresolved access or quality gate will stop external validation.
4. I will use an isotropic Fisher–Kolmogorov PDE with a joint history-only posterior over growth and diffusion parameters. I will combine PDE-predicted shape with a training-only volume-trend model and correlated boundary uncertainty. The forecast will use a 2-mm grid, a day-90 horizon, and 200 joint forecast draws. I will report parameter ridges and bound hits separately; uncertain parameters will not receive a biological interpretation.
5. I will compare the full forecast with a matched-volume, no-PDE signed-distance forecast as the primary comparator. I will average the 200 forecast draws voxel by voxel and threshold inclusion probability at 0.5 to form the hybrid point mask; I will rank the no-PDE inclusion probabilities to match that mask's voxel count. Persistence and geometric-trend forecasts will be secondary controls. I will fit volume and uncertainty terms inside the training patients only and will use the same patients, target masks, and forecast horizon for every method. I will also run the prior-bound and central-gap sensitivity analyses specified in the appendix.
6. I will use patient-macro Dice as the primary geometry measure, with both-empty masks scored as 1 and one-empty masks scored as 0. I will report patient-level paired differences and a patient bootstrap interval. Surface distance and absolute volume error will be secondary measures. I will count every empty mask, invalid forecast, fallback, and technical failure in the cohort flow.
7. Before I score patient masks, source-backed evidence will need to establish longitudinal registration quality and valid target-label support over the entire target grid for every included patient, with label zero documented as a valid negative. I will define the origin-only model domain separately from target scoring support. Matching image headers or a nonzero MRI intensity mask will not count as proof of registration or valid negative labels. Missing evidence will trigger a stop before target-mask scoring and will be reported as a barrier.
8. I will treat all MU patient estimates as exploratory. Duke will be the only possible independent patient test in this plan; UCSD and CFB annotation-change work will remain separate and will not fill that role.

### Track C — simulation-only treatment-policy comparison

1. I will use a two-population mathematical simulator with sensitive burden \(S\), resistant burden \(R\), and total burden \(B=S+R\). Its treatment effects, growth rates, resistance costs, initial resistant fractions, and observation noise will be assumptions, not measurements from a drug or patient.
2. I will use 28-day decision cycles, a five-day treatment course, a maximum of 12 courses, a 52-cycle horizon, and a 0.25-day numerical integration step. The simulated progression event will be the first time \(B\) reaches 1.2 times its initial value; I will label this as a model-defined threshold, not a clinical progression criterion.
3. I will record scenario seeds used earlier in the repository and treat those scenarios as development material. I will generate 300 new development scenarios, 150 filter-calibration scenarios, and 600 held-out parameter scenarios, with 200 held-out scenarios at each observation-noise level (0.05, 0.10, and 0.20). I will generate two independent observation-noise streams per held-out scenario and keep both streams grouped with that scenario when estimating uncertainty. Before evaluating held-out outcomes, I will time a 100-scenario development pilot. If the projected full run exceeds 12 hours, I will reduce the held-out set once to 300 scenarios (100 per noise level), using runtime alone and freezing that size before opening outcomes.
4. I will compare no treatment; periodic schedules with a 1–6-cycle interval and a start at cycle 0 or 2; raw-observation on/off thresholds; raw-observation single thresholds; and a resistance-aware particle-filter policy with 500 particles. I will tune all policy settings on development scenarios only. I will use the separate calibration set to describe 90% interval coverage and width for simulated burden and resistant fraction. I will not interpret simulated calibration as evidence that a patient's hidden resistance can be estimated reliably.
5. I will select the best simple policy from no treatment, periodic schedules, and raw-observation policies using the same development-only lower-tail objective. The primary outcome will be the minimum, across noise levels, of the held-out scenario-level lower-tail mean (CVaR at 10%) of time to progression. For each scenario, I will average its two noise-stream outcomes before calculating CVaR; I will include a fractional contribution from the cutoff scenario when needed to make the lower tail exactly 10%. I will compare the resistance-aware policy with the development-selected simple policy. I will resample scenarios within noise level in a paired bootstrap and recompute the minimum CVaR difference for each of 10,000 replicates. Median progression time, capped fraction, burden integral, resistant fraction, and courses used will be secondary outcomes. I will report solver failures, count them as zero days for the primary outcome, and will not drop failed scenarios.
6. I will run three fresh, prespecified sensitivity sets of 100 scenarios each: 40% missing observations; a 5% decline in treatment strength per completed course; and competition coefficients (a_{sr}=2.0) and (a_{rs}=0.5), compared with the base values of 1.0. I will use new seeds, freeze the settings before testing, and not retune policies on these sets. I will describe robustness only over these named simulated departures.
7. I will report simulation assumptions and failures beside every result. No Track C output will be described as a real-drug comparison, a treatment recommendation, a patient-specific benefit, or clinical evidence.

## Risk assessment

The project will analyze de-identified, previously collected files and generated simulations. I will keep patient-level data local, protect access, follow dataset attribution and use terms, avoid re-identification or cross-dataset patient linkage, and not publish patient-level files. I will obtain school review before using restricted data. I will describe all model outputs as research predictions or simulations, not medical advice or treatment evidence.

## Data analysis and interpretation

For Track A, higher donor-macro balanced accuracy will indicate better classification. For Track B, higher patient-macro Dice and lower surface or volume error will indicate closer agreement with the later released mask. For Track C, a higher lower-tail simulated progression time will indicate better performance under the stated simulator assumptions. Each primary contrast will use the independent unit for that track and a paired 95% bootstrap interval. A favorable estimate with an interval that includes no difference will be called inconclusive. Even an interval that excludes no difference will support only the stated dataset or simulator comparison; it will not establish clinical value or novelty.

## References

- Fischer, F., Fischer, D. S., Mukhin, R., Isaev, A., Biederstedt, E., Villani, A.-C., & Theis, F. J. (2024). scTab: Scaling cross-tissue single-cell annotation models. *Nature Communications, 15*, 6611. https://doi.org/10.1038/s41467-024-51059-5
- Hormuth, D. A., II, Al Feghali, K. A., Elliott, A. M., Yankeelov, T. E., & Chung, C. (2021). Image-based personalization of computational models for predicting response of high-grade glioma to chemoradiation. *Scientific Reports, 11*, 8520. https://doi.org/10.1038/s41598-021-87887-4
- Liu, Z., Zhang, J., Hong, L., Nie, Q., & Sun, X. (2025). Multiscale mathematical model-informed reinforcement learning optimizes combination treatment scheduling in glioblastoma evolution. *Science Advances, 11*(32), eadv3316. https://doi.org/10.1126/sciadv.adv3316
- Lopes, M. B., & Vinga, S. (2020). Tracking intratumoral heterogeneity in glioblastoma via regularized classification of single-cell RNA-Seq data. *BMC Bioinformatics, 21*, 59. https://doi.org/10.1186/s12859-020-3390-4
- Mahmoud, E., et al. (2025). MU-Glioma Post: A comprehensive dataset of automated MR multi-sequence segmentation and clinical features. *Scientific Data, 12*, 1847. https://doi.org/10.1038/s41597-025-06011-7
- Subramanian, S., Scheufele, K., Mehl, M., & Biros, G. (2020). Where did the tumor start? An inverse solver with sparse localization for tumor growth models. *Inverse Problems, 36*(4), 045006. https://doi.org/10.1088/1361-6420/ab649c
- West, J., et al. (2023). A survey of open questions in adaptive therapy: Bridging mathematics and clinical translation. *eLife, 12*, e84263. https://doi.org/10.7554/eLife.84263
- Zhong, Z., Hou, J., Yao, Z., Dong, L., Liu, F., Yue, J., Wu, T., Zheng, J., Ouyang, G., Yang, C., & Song, J. (2024). Domain generalization enables general cancer cell annotation in single-cell and spatial transcriptomics. *Nature Communications, 15*, 1929. https://doi.org/10.1038/s41467-024-46413-6
- Zimmer, L., Weidner, J., Balcerak, M., Kofler, F., Krupa, M., Ezhov, I., Cepeda, S., Zhang, R. Z., Lowengrub, J. S., Menze, B., & Wiestler, B. (2026). PREDICT-GBM: A multicenter platform advancing personalized glioblastoma radiotherapy planning. *npj Digital Medicine, 9*, 686. https://doi.org/10.1038/s41746-026-03194-0
- Zhang, J., Cunningham, J. J., Brown, J. S., & Gatenby, R. A. (2017). Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer. *Nature Communications, 8*, 1816. https://doi.org/10.1038/s41467-017-01968-5

Official dataset pages and detailed parameter definitions will be listed in the companion protocol.
