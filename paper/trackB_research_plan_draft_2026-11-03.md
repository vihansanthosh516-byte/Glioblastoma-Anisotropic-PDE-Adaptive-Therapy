# Track B Research Plan Draft

## Project Title

Testing whether a reaction-diffusion model improves patient-specific forecasting of glioblastoma MRI tumor-core labels

## Rationale

Glioblastoma can spread through brain tissue beyond the region that appears abnormal on an MRI. Reaction-diffusion models represent this spread with parameters for cell dispersal and proliferation. Researchers have fitted such models to patient MRI and have used them to study tumor growth and treatment response (Rockne et al., 2010; Harris et al., 2022). Recent work also combines tumor-growth models with treatment-aware digital twins and compares biophysical forecasts with neural forecasts (Bhattacharya et al., 2025; Zimmer et al., 2026). Therefore, this project will not claim that reaction-diffusion modeling, hybrid modeling, or MRI forecasting is new by itself.

The Duke dataset will provide longitudinal MRI sessions, automated tumor labels, relative treatment dates, and vessel-landmark pairs. The record describes 61 patients and 503 sessions, but those totals do not show how many patients will meet this study’s eligibility rules. The dataset authors state that scans are aligned to a common atlas field of view but that separate sessions still require longitudinal registration. They also report that scan dates were extracted automatically and may be off by several days (Criscuolo et al., 2026; Yang, n.d.).

The scientific question is narrower: when the PDE model and a non-PDE shape forecast receive the same prior scans, forecast interval, support, and history-only volume draws, will the PDE-derived spatial ranking better match the next observed automated tumor-core mask? This controlled comparison will test predictive value beyond volume and interval information. It will not isolate a causal effect of PDE biology. Because the target masks are automated, the result will measure agreement with supplied labels, not accuracy against verified tumor. The work may improve how future tumor-forecasting methods are compared, but it will not guide patient care or make a clinical recommendation.

## Research Question, Hypothesis, Engineering Goals, and Expected Outcomes

### Research Question

Among eligible Duke glioblastoma patients with at least three usable MRI sessions, will a reaction-diffusion PDE-based spatial ranking forecast the next observed automated tumor-core mask more accurately than a two-history non-PDE shape forecast, when both methods use the same history-only volume forecasts, observed forecast interval, and origin-defined support?

### Hypothesis

The PDE-based forecast will have a higher patient-macro Dice score than the two-history non-PDE control. The null expectation will be no improvement or a lower score. A tie, loss, or inconclusive interval will remain plausible.

The primary comparison will be paired by patient. The project will call the hypothesis supported only if the candidate meets the success rule frozen before target-mask scoring. The proposed strict rule will require a positive mean paired Dice difference and a multiplicity-adjusted 95% lower confidence bound above zero against every required mask baseline. The final Duke-specific protocol will name the full comparison family and uncertainty procedure before any target mask is opened.

### Engineering Goals

The project will build a reproducible pipeline that will:

1. Keep every session from one patient in the same validation fold.
2. Fit preprocessing, volume forecasts, and comparator settings using training patients only.
3. Give the focused PDE and two-history forecasts identical volume draws, intervals, support, eligible patients, and target labels.
4. Record exclusions, registration quality, fallbacks, software settings, and file fingerprints.
5. Stop before scoring if access, privacy, cohort, protocol, compute, or integrity requirements fail.

### Expected Outcomes

The candidate may outperform the two-history control, tie it, or perform worse. The project will report each outcome without changing thresholds or model settings after scoring. A positive result will support only a predictive comparison against the named automated labels in this cohort and its observed interval distribution. It will not establish biological causation, true tumor accuracy, clinical usefulness, generalization, or novelty by itself.

## Procedure

1. **Confirm data-use terms and review.** The project will use the restricted Duke archive only under the written access agreement accepted for the dataset. Before analysis, the student will record the applicable terms for analysis and for sharing summary results, figures, predictions, masks, and trained weights. The student will obtain the required school or competition review and written approval. No patient-level output will be shared unless the agreement and review permit it.

2. **Verify the archive.** The project will check that `DATASET_FINAL.zip` has the expected size of 17,043,569,844 bytes and the Zenodo-recorded MD5 `99580fb8ccb3daed186dc6afe5b5bf2f`. It will calculate a local SHA-256 after the expected size and MD5 pass. A mismatch will stop the study. The project will list archive paths and file sizes before opening data members.

3. **Complete a metadata-only cohort audit.** The project will first inspect approved metadata fields only: anonymized patient and session identifiers, relative scan dates, scan availability, diagnosis and IDH fields, treatment dates, and surgery dates. The project will not open MRI images, segmentation masks, brain-support masks, or landmark coordinates during this stage. It will identify complete, chronologically ordered three-session windows with the required images and labels. It will check known patient-level overlap with cohorts already used in the project, where the access terms permit that check. If overlap cannot be resolved, the result will not be called independent external confirmation.

   The planned operational floor will be 20 independent eligible patients. This is a feasibility floor for grouped nested validation, not a power guarantee. Before target masks are opened, a precision analysis using development-only information will either support this floor or raise it. If the eligible count is below the final frozen floor, the patient-scoring branch will stop or be labeled a feasibility study.

4. **Freeze the cohort and analysis rules before target-mask access.** The project will use consecutive three-session windows that pass the metadata rules. All windows from one patient will remain together in one outer fold. The project will calculate a patient-level score by averaging that patient’s eligible-window scores first, then averaging across patients, so patients with more scans will not receive greater weight. Treatment and surgery dates will be recorded. The prediction will not receive future treatment information; results will be interpreted as forecasts of observed label change under the recorded clinical histories, not untreated biological growth.

   The proposed validation design will use five patient-grouped outer folds and five patient-grouped inner folds, with seed 20261103. Every transformation and setting will be fitted inside the corresponding training split. The two-history setting selector will use only the outer-training patients; every inner-training split will contain at least 10 eligible patients. No pooled result across outer folds will change a fold’s prediction. If the final eligible cohort cannot support these rules, the project will stop or revise the design before target scoring.

5. **Register longitudinal scans and define support using MRI only.** After separate approval for image reads, the project will rigidly align each later-session post-contrast T1 image to the origin-session post-contrast T1 image using MRI intensities only. The dataset is resampled to 1 mm isotropic voxels and atlas-aligned for field-of-view standardization, but its author guide says sessions are not registered to each other and recommends rigid or affine alignment. The project will apply each transform to the other modalities and masks, using linear interpolation for intensities and nearest-neighbor interpolation for labels. Before real-patient image QC, synthetic images with known transforms will be used to check the implementation. The final fit region, similarity metric, quality thresholds, retry rule, and support-generation method will be frozen before the real-patient QC pass. Target segmentation masks will remain closed during this QC. If permitted by the dataset agreement and school review, the supplied vessel landmarks will provide an independent registration check for the preoperative-to-recurrence pairs; they will not tune the forecast or select patients.

   The common support will be defined from the origin visit only using a documented brain-support method. It will not be defined from target masks or from a post hoc union of nonzero voxels. The exact support method will be chosen and frozen before patient images are read. Predictions will be restricted to this support. Target-core voxels outside the support will count as misses. The project will freeze the QC-passing patient list and its hash before target labels are opened.

6. **Define labels and eligible forecast windows.** The project will use the dataset’s supplied automated segmentations. The target tumor-core mask will combine label 1 (necrotic tumor core) and label 3 (enhancing tumor). Label 2 (nonenhancing T2-FLAIR signal abnormality) will not be included in this core endpoint. Label 4 (resection cavity) will remain separate and will not count as tumor core (Criscuolo et al., 2026; Zhang et al., 2025). Empty or unusable history masks, missing labels, invalid date order, failed image QC, and missing required modalities will follow rules frozen before target scoring. The project will count each exclusion and failure. It will not treat automated labels as expert ground truth.

7. **Generate forecasts from history only.** The two scans before each target will be the history. The third scan will be the next observed target visit. The actual interval will be calculated from the dataset’s relative dates and supplied to methods that use time. The project will report that these dates may be imprecise by several days. Because the realized interval will be known only after the target visit occurs, this will be a retrospective interval-conditioned evaluation, not a prospective fixed-day forecast.

   The PDE candidate will calibrate a reaction-diffusion model using only the history scans and permitted history MRI. It will produce a spatial ranking for the target interval. A separate history-only volume model will produce target-volume draws using training-patient information only. The candidate will convert its ranking into masks using those draws and the common support.

   The main non-PDE control will extrapolate signed-distance shapes from the two history masks. It will use the same target-volume draws, interval, support, patients, and scoring rule as the PDE candidate. The proposed control grid will be `SHRINK_TAU_DAYS {0, 7, 28} × MAX_FACTOR {1, 4, 8} × MAX_SD_DELTA_MM {10, 30}`, plus the `f = 0` last-mask signed-distance rule. Inner patient-grouped validation will select settings using training patients only. The selector will maximize mean patient-level validation Dice; settings within 0.005 will tie, and the more damped setting will win. The outer-test patient will not affect selection.

   Copy-last will serve as a simple persistence baseline. A fixed-growth baseline will be included only if its time-dependent rate grid and training-only selection rule are frozen before target-mask access. PDE-only or learned-model comparators will be added only through a dated protocol update that defines their inputs, recipes, compute, and fair comparison. The project will not add a comparator after seeing test outcomes.

8. **Apply common failure rules.** The final Duke-specific protocol will name expected recoverable model failures and one common fallback for the focused PDE/control pair. A recoverable failure in either member of the pair will give both methods the same prespecified fallback mask for that patient and will be recorded. Unexpected exceptions, data-integrity errors, or code errors will stop the run as incomplete rather than trigger a fallback. The project will not silently drop a scored patient because a model performed poorly.

9. **Run a synthetic integration and timing pilot.** Before real-patient scoring, the project will test fold isolation, target-read guards, identical patient sets and support, history-only model fitting, common volume draws, fallback accounting, and result hashes on synthetic data. The timing pilot will include every required outer and inner fit, calibration call, comparator, refit, and restart allowance on the actual grid and interval distribution. The student will set a compute ceiling before the pilot. If projected time exceeds that ceiling, the student will stop and revise the plan before target scoring.

10. **Seal and score once.** After written permissions, eligibility, QC, precision, runtime, and protocol checks pass, the project will seal the code, settings, fold map, cohort list, package environment, and permitted input-file hashes. The evaluation will score every required method on the same held-out patients. The project will not tune on outer-test labels or rerun with changed settings after seeing results. Missing arms, non-finite outputs, or a failed seal will produce an incomplete result.

## Materials Required

- **Computer:** Dell 16 Plus DB16250, Windows 11, 32 GB RAM, CPU-only computing environment. The analysis will use the computer’s 8-thread setting for the locked run unless the runtime freeze specifies otherwise.
- **Dataset:** restricted Duke University dataset, Zenodo record 20548528, `DATASET_FINAL.zip`. The record describes 61 patients, 503 MRI sessions, 3–18 sessions per patient, six MRI modalities, relative treatment dates, automated segmentation masks, and preoperative-to-recurrence vessel landmarks (Criscuolo et al., 2026). The project will use only files allowed by the accepted Duke agreement.
- **Software:** Python 3.14.4 (64-bit); NumPy, SciPy, pandas, NIfTI-reading and image-registration software, and project code. PyTorch or scikit-learn will be included only if a learned comparator is added by a dated protocol update. Exact package versions, operating-system build, thread settings, and code commit will be recorded and sealed before scoring.
- **Storage:** local access-controlled storage for the restricted archive and patient-level metadata. Patient-level identifiers, images, masks, landmarks, and outputs will not be committed to the public project repository or uploaded to an external service.
- **Chemicals and biological materials:** none. The project will use existing digital MRI and metadata only. It will not handle bacteria, cell cultures, animals, chemicals, or physical patient samples.

## Risk Assessment

The main risk will be privacy loss from restricted medical images, dates, labels, or derived predictions. The project will use anonymized identifiers, keep the archive and patient-level files in access-controlled local storage, avoid identity linkage, and avoid uploading patient data to external services. It will not publish patient-level images, masks, coordinates, predictions, or trained weights unless the written agreement and school review explicitly permit them. The report will use aggregate results and suppress details that could identify a patient.

The study will not contact patients, intervene in care, diagnose disease, recommend treatment, or return predictions to clinicians or participants. No physical, chemical, or biological laboratory hazards are expected. The computer may become warm during long calculations; it will run on a hard, ventilated surface and will not be left unattended if cooling fails.

## Data Analysis

The primary outcome will be Dice similarity between each forecast mask and the next observed automated core mask:

`Dice = 2 × |forecast ∩ target| / (|forecast| + |target|)`.

If both masks are empty, Dice will equal 1. If only one is empty, Dice will equal 0. Target-core voxels outside origin-defined support will count as false negatives. The project will first average Dice across eligible windows within each patient, then calculate the equal-patient mean. Cells, voxels, and repeated windows will not be treated as independent patients.

The primary contrast will be the paired patient-level difference between the PDE candidate and two-history control. The proposed comparison family will include the PDE candidate versus the two-history control and copy-last. The proposed analysis will resample patients 10,000 times with seed 20261104 and use the same patient resamples for both contrasts. For each bootstrap sample, it will calculate each paired mean difference. It will divide each sample’s absolute deviation from the observed mean by that contrast’s bootstrap standard deviation, take the 95th percentile of the maximum statistic across contrasts, and use that critical value to form simultaneous intervals as observed mean difference plus or minus critical value times bootstrap standard deviation. A zero or non-finite standard deviation will make the result incomplete. The final Duke protocol will freeze the bootstrap method, random seeds, multiplicity adjustment, minimum precision target, and exact pass rule before target masks are opened. Under the proposed strict rule, the candidate will pass only if its mean Dice is higher and every adjusted lower confidence bound is above zero. The project will report per-patient scores, the mean and spread, confidence intervals, interval distribution, exclusions, missing arms, fallbacks, and runtime. Negative or incomplete results will be reported unchanged.

The project will not interpret Dice agreement as agreement with true tumor, and it will not claim that the PDE caused any improvement. Results will be limited to prediction of the dataset’s automated labels for the eligible Duke patients and observed scan intervals.

## Bibliography

Bhattacharya, M., Singh, G., & Prasanna, P. (2025). *SoC-DT: Standard-of-care aligned digital twins for patient-specific tumor dynamics* [Preprint]. arXiv. https://doi.org/10.48550/arXiv.2510.03287

Criscuolo, E., Calabrese, E., Zhang, Z., Wang, Y., Xiong, Y., & Yang, D. (2026). *A longitudinal glioblastoma MRI dataset with anatomical landmark pairs for DIR validation* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.20548528

Harris, D. C., Mignucci-Jiménez, G., Xu, Y., Eikenberry, S. E., Quarles, C. C., Preul, M. C., Kuang, Y., & Kostelich, E. J. (2022). Tracking glioblastoma progression after initial resection with minimal reaction-diffusion models. *Mathematical Biosciences and Engineering, 19*(6), 5446–5481. https://doi.org/10.3934/mbe.2022256

Rockne, R., Rockhill, J. K., Mrugala, M., Spence, A. M., Kalet, I., Hendrickson, K., Lai, A., Cloughesy, T., Alvord, E. C., Jr., & Swanson, K. R. (2010). Predicting the efficacy of radiotherapy in individual glioblastoma patients in vivo: A mathematical modeling approach. *Physics in Medicine & Biology, 55*(12), 3271–3285. https://doi.org/10.1088/0031-9155/55/12/001

Yang, D. (n.d.). *GBM-DIR-QA: Dataset of glioblastoma patient MRI images with blood vessel bifurcation landmarks for deformable image registration quality assurance* [Data documentation]. GitHub. https://github.com/deshanyang/GBM-DIR-QA

Zhang, J., LaBella, D., Zhang, D., Houk, J. L., Rudie, J. D., Zou, H., Warman, P., Mazurowski, M. A., & Calabrese, E. (2025). Development and evaluation of automated artificial intelligence-based brain tumor response assessment in patients with glioblastoma. *American Journal of Neuroradiology, 46*(5), 990–998. https://doi.org/10.3174/ajnr.A8580

Zimmer, L., Weidner, J., Balcerak, M., Kofler, F., Krupa, M., Ezhov, I., Cepeda, S., Zhang, R., Lowengrub, J., Menze, B., & Wiestler, B. (2026). *PREDICT-GBM: A multi-center platform to advance personalized glioblastoma radiotherapy planning* [Preprint]. arXiv. https://doi.org/10.48550/arXiv.2509.13360

## Additional Project Information

This project will analyze pre-existing clinical MRI and related metadata. The student will not recruit or interact with patients. The student will not assume that restricted retrospective data qualify for an exemption. The student will obtain a written determination from the school or competition review body about required human-data review, privacy documentation, and pre-approval before opening MRI, masks, or landmark coordinates.

The dataset authors will have acquired the scans and clinical metadata, identified the vessel landmarks, and generated the automated masks. The student will perform the planned metadata audit, registration quality checks, modeling, validation, and statistical comparisons. The student will not claim to have collected the patient data or created the source labels.

This work will continue the project’s prior glioblastoma modeling work. The student will provide prior abstracts and plans if the competition requests them and will describe the new work as a Duke-cohort, patient-grouped forecast comparison using a matched-volume two-history control. The study will not present earlier exploratory results as independent confirmation.
