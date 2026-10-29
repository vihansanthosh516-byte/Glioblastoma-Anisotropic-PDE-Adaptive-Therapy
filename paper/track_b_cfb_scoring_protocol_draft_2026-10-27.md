# CFB-GBM Track B scoring protocol draft

**Status (updated 2026-10-28):** outcome-blind draft; not authorization to score. All 107 candidates pass file/header checks against the official baseline t0 T1-Gd reference grid and pass t1 origin-MRI normalization QC. These checks establish shared geometry and usable origin intensities, not anatomical registration quality. Mask review lineage and source transform records remain unresolved. This is a separate CFB protocol draft. The locked UCSD protocol remains unchanged.

## 1. Question and scope

Among eligible CFB-GBM patients with three consecutive, file-matched visits (t0, t1, t2), does the spatial pattern of released GTV annotation change from t0 to t1 improve patient-level ranking of newly annotated t2 GTV voxels, compared with a model using only t1 T1-Gd MRI and distance to the t1 GTV?

The target is **future released-annotation addition**, not tumor growth, recurrence, progression, treatment response, or clinical benefit. CFB GTV masks include clinical and nnU-Net-generated annotations; the endpoint can reflect annotation source and review process. No result, positive direction, or uniqueness is guaranteed.

**One primary hypothesis:** the trajectory-augmented model has higher mean patient-level AUC than the frozen t1-only comparator (two-sided alpha .05; direction declared in advance as augmented minus comparator > 0). There is one primary contrast. Dose, alternate modalities, gap strata, and alternate endpoints are secondary or descriptive, not rescue analyses.

## 2. Estimand and endpoint

**Independent unit:** patient. Each patient contributes at most one t1-to-t2 pair, selected by the release-matched metadata rule below, without looking at target voxel values. No voxels are treated as independent observations.

**Patient-level estimand:** for each AUC-evaluable patient i, calculate `d_i = AUC_i(trajectory model) - AUC_i(frozen comparator)`. The primary effect is `Delta = mean_i(d_i)`, giving every AUC-evaluable patient equal weight. Report Delta, its two-sided 95% paired patient-bootstrap confidence interval (10,000 resamples; fixed seed), evaluable n, and the full cohort flow. Bootstrap patients, carrying each patient's paired predictions and labels together.

**Target:** on the shared baseline t0 T1-Gd grid, `Y_i(x)=1` when voxel x is inside the t2 GTV but outside the t1 GTV; otherwise `Y_i(x)=0`. The primary risk set is the intersection of valid t1 brain support and valid t2 annotation support, excluding t1 GTV only. Do not exclude a cavity: the available CFB endpoint is GTV, and a compatible cavity label is not established. Do not infer cavity or edema from the binary GTV mask. Report t2 additions outside the risk set separately; they are not primary positives. The official pipeline fixes baseline t0 T1-Gd and resamples follow-up images to that grid. All MRI and GTV files must share that grid. This header check does not prove anatomical registration quality.

Define valid t2 support from documented annotation coverage and image-acquisition support, using source documentation and image/header information only. Never infer support from t2-mask foreground, connected components, or observed additions. If annotation coverage cannot be established without inspecting target-mask values, stop. Estimate and quality-check spatial transforms from image data and geometry only; never fit, select, or accept a transform using GTV contours or target-derived mask similarity. Freeze the transform method and QC tolerances before applying transforms to masks.

A patient has a defined primary AUC only if the risk set contains at least one positive and one negative voxel. Patients with zero additions, all-positive/all-negative risk sets, or technical failure remain in the flow with reason; never assign a fabricated AUC. Thus the primary estimand is conditional on a patient having an evaluable annotation-addition target. Any all-patient endpoint would require a separate prespecified estimand and threshold and is not added after labels are seen.

## 3. Cohort and visit rule

Use the TCIA current CFB-GBM v4 release and its current linked availability tables. The downloaded package has 264 numeric patient folders and 3,976 NIfTI files, matching every linked availability flag. Current `id_patient` values match package folder IDs and filename prefixes. Every availability flag matches a package file for the listed MRI sequences, GTV, and RTDOSE; all 107 table-defined candidates have all six t0/t1/t2 T1-Gd and GTV files. The linked table and package each contain 112 t2 GTV rows/files. TCIA's collection summary and the paper report 122, and the collection page reports 3,981 images while this package has 3,976 files. These aggregate differences remain unexplained. Define this cohort from current linked table rows with matching files; report aggregate count discrepancies and do not add unlisted cases.

Eligible patient must have: (a) verified identity mapping across t0, t1, t2 in the current v4 release; (b) complete T1-Gd MRI and GTV masks at all three visits; (c) valid, documented spatial transforms and compatible geometry; (d) mask provenance/review metadata at all three visits; (e) t1-to-t2 interval 14–210 days inclusive, matching the currently observed candidate range; and (f) origin-image QC pass at t1. Select the earliest qualifying t1-to-t2 pair for each patient after requiring an immediately preceding eligible t0 visit. If the source metadata does not establish ordered visits or prior-visit relation, stop and resolve it; do not choose based on outcomes.

**Horizon:** next observed eligible scan, variable follow-up. The model must not receive the realized t1-to-t2 interval as a predictor. Report interval distribution and predeclared descriptive strata; do not call this a fixed-horizon forecast. The 45–75-day band (34 table candidates before exclusions) is descriptive sensitivity only and cannot replace the primary cohort after outcome review.

Current metadata joins show 107 file-matched candidate t1-to-t2 pairs with T1-Gd and GTV, median gap 84 days, IQR 56–98 using Tukey hinges, range 14–210. All 107 pass header checks against their t0 T1-Gd grid and t1 T1-Gd robust-normalization QC. Baseline grids use approximately 0.5 × 0.5 × 1.0 mm spacing; 102 are 512 × 512 × 208 and five are 512 × 512 × 224. Shared headers show geometry only. They do not establish anatomical registration quality or eligible scoring n. The source paper lists deceased patients among inclusion criteria; verify whether this applies to these candidates and disclose any resulting survival-related selection.

## 4. Mask provenance gate

Before target-mask voxels are opened, obtain a row-level record of whether every t0, t1, and t2 GTV was clinical/manual, generated, reviewed/corrected, or unknown, plus the evidence source for each classification. In the 107 file-matched candidates, current metadata marks all t1 and t2 GTVs `genSegmentation`; t0 has 93 `humanSegmentation` and 14 `genSegmentation`. The source paper reports that 44 generated masks were sampled and reviewed/corrected by five radiation oncologists, but publishes no patient/timepoint IDs. Current linked TSVs contain no per-mask review field. Review status is therefore unknown for every candidate. For 93/107, the t0-to-t1 trajectory also crosses human-to-generated annotation types and may encode annotation-process change.

Primary confirmatory-style analysis requires one prespecified, consistent provenance and review class across all three visits. Do not combine label-generation classes in the primary set merely to increase n. Current evidence identifies only 14 candidates with generated masks at all three visits; their case-level review status is still unknown, so no provenance-and-review-matched primary set is established. A mixed-source annotation-transition analysis could only be explicitly exploratory under a separate frozen protocol; report provenance strata and the human-to-generated transition. Unknown review status must remain unknown. Do not exclude cases based on future mask appearance or observed AUC.

## 5. Inputs, models, and development

**Primary input modality:** t1 T1-Gd only, with distance to current t1 GTV, both on the shared t0 reference grid. Freeze the exact series name mapping, registration chain, voxel units, interpolation, brain-support rule, robust normalization, and finite-value/zero-scale failure rule before outcome access. Baseline grids use approximately 0.5 × 0.5 × 1.0 mm spacing; compute spatial distances in physical millimetres. Do not reuse the existing 1-mm isotropic four-channel method or resample without a separate frozen plan. T1-Gd is the common candidate modality in the current metadata screen. Do not port the existing four-channel method as-is: its sequence mapping fits only 2 of 107 candidates under one mapping and none under an alternative. Adding FLAIR or any other MRI sequence defines a separate secondary model and requires a new compatible development plan.

**Comparator (frozen baseline):** a patient-level voxel-ranking score trained/developed using only t1 T1-Gd features and distance to the current t1 GTV. Freeze its feature definition, algorithm, coefficients/hyperparameters, and preprocessing before CFB target access.

**Trajectory model:** exactly the same t1-only inputs and score plus an origin-known local t0-to-t1 GTV transition feature, computed only from t0 and t1 annotations after their provenance passes. No t2-derived feature, post-origin scan, realized future interval, or dose enters the primary model. This isolates the incremental value of prior local annotation trajectory.

Before any CFB target labels are accessed, freeze for each arm the complete feature list, model family, preprocessing and scaling, finite hyperparameter/feature-selection search space, tie-breaking rule, random seeds, convergence criteria, and technical-failure rule. Do not expand or revise this search after inspecting any CFB target labels.

**Development path:** preferred path is a compatible independent development cohort with the same image sequence, three-visit GTV endpoint, annotation provenance, grid/registration assumptions, and prediction task. Freeze both models there, then evaluate CFB once. Compatibility must be evidenced, not presumed from a cohort name. The already explored MU cohort and four-modality method are not automatically compatible.

If no independent compatible development cohort exists, use patient-grouped nested cross-validation in CFB: one outer 5-fold grouped split for out-of-fold predictions; all scaling, feature selection, fitting, and tuning inside inner grouped folds of outer-training patients only. Before target labels are opened, freeze the de-identified patient-to-fold mapping, deterministic split seed, inner-fold rule, and complete search grid. Use identical outer folds, patients, voxel risk sets, and model failure policy for both arms. Produce exactly one out-of-fold prediction per patient. This is internal exploratory method development, not external confirmation. Do not use the partially exposed UCSD set as a fresh confirmatory test.

Dose-aware prediction is not primary. Only consider it as a later secondary analysis if delivered-dose timing before t1, dose summation/units, patient mapping, and spatial transform to t1 are all independently verified. Planned dose is not a substitute for delivered dose; any alternative changes the estimand. Do not make causal claims.

## 6. Evaluation and failure handling

The primary metric is patient-specific ROC AUC on the frozen risk set. Both model arms must be scored on exactly the same patients. Report the paired denominator and reasons for every loss from candidate to eligible to AUC-evaluable cohort. No pooled-voxel AUC for inference. Report secondary arms (persistence/no-change, distance-only, T1-Gd-only) descriptively on the same patients if preregistered and technically valid; no cherry-picking the comparator after results.

A technical failure is determined using origin-only information or a predeclared file/geometry rule. Record it without inspecting target values. No method-specific exclusion from the paired primary comparison: if either arm fails for a patient, omit that patient from the paired contrast and report the failure; use one common paired denominator. Do not add missing-modality fallbacks after seeing outcomes. Keep zero-addition and one-class cases visible in cohort flow, even though their AUC is undefined.

Use a patient-cluster bootstrap for the primary interval. Report the paired mean, median, standard deviation, interval, patient-level wins/ties/losses, and per-patient event-evaluable fraction. Any other metric or subgroup is labeled secondary/descriptive and multiplicity is not converted into a new primary claim.

## 7. Power and precision

The final eligible and AUC-evaluable n, event rate, and paired AUC-difference variability are unknown until data gates pass; do not infer power from 107 metadata candidates. Before target access, use independent development/pilot data and the locked cohort/evaluation code to run a patient-cluster simulation for precision/power, preserving event prevalence and patient-level clustering. Freeze a minimum useful effect only if independently justified. The prior UCSD 0.02 AUC margin is not automatically clinically meaningful for CFB.

For scale only, a two-sided normal approximation at alpha .05 and 80% power needs about 49 patients for a 0.02 AUC effect if paired SD is .05, 196 if SD=.10, and 784 if SD=.20. At n=107, approximate detectable effects are .0135, .0271, and .0542 respectively; at n=34, .0240, .0481, and .0961. These are assumption-based calculations, not CFB-derived power. Actual evaluable n will be lower than the table candidate count. If uncertainty is too wide, report an inconclusive exploratory estimate; do not loosen criteria or redefine horizon/target after opening outcomes.

A statistically positive direction requires the two-sided 95% CI for Delta to lie above zero. This does not establish clinical utility or meaningful benefit. A claim that performance exceeds a nonzero practical margin requires an independently justified, preregistered margin and a lower confidence bound above it. No such CFB-specific margin is currently established.

## 8. Stop gates before target-mask access

1. Record current v4 public/CC BY 4.0 terms, required attribution, package version, and prior exposure history. Verify whether the paper’s deceased-patient inclusion rule applies to this cohort and state its generalizability limit.
2. Preserve the verified current-table-to-package mapping: 107 candidate pairs have all six files; the 112 table rows match 112 package t2 GTV files. Keep the 122 and 3,981 aggregate counts as documented, unresolved source discrepancies.
3. All 107 candidates pass six-file header alignment to baseline t0 T1-Gd and t1 origin-only intensity QC. Before mask access, obtain per-case anatomical registration evidence and QC independent of GTV contours; shared shape and affine do not prove anatomical registration quality.
4. Resolve per-mask review status or freeze a separate exploratory plan that explicitly retains unknown review status and mixed t0 provenance. No confirmatory-style provenance-matched subset is currently established.
5. Verify whether an independent development source is compatible. If none exists, freeze CFB-only nested cross-validation as exploratory, not external confirmation.
6. Freeze code, patient and visit rule, physical-distance implementation, feature definitions, fold assignment, failure rules, metrics, bootstrap, power assumptions, and file hashes before any target-mask voxel access.
7. Assess precision using independent development/pilot information. If precision is inadequate, stop confirmatory-style scoring and retain feasibility/exploratory status.

Do not score under this draft. Target-mask review status, per-case registration QC, compatible development, and a complete frozen runner remain open. If a separate exploratory protocol is later approved, open target masks only for that frozen analysis and preserve a one-shot exposure record. No outcome-driven protocol changes.

## 9. Interpretation and Track B boundary

A positive CFB result would support added ranking value for later **released annotation additions** under this cohort and protocol only. It would not establish biological growth, recurrence, treatment response, clinical utility, generalization, or uniqueness. A negative or inconclusive result is reportable and should not be repaired by post hoc endpoint or cohort changes.

CFB is one candidate cohort, not completion of all Track B. Track B also includes the broader reaction-diffusion/PDE parameter estimation and forecasting program, where prior work did not beat simple baselines and fitted diffusion parameters often reached bounds. Completing this CFB protocol or obtaining a positive CFB result does not resolve those limitations. Track B completion requires a separate, coherent evidence package covering the prespecified PDE/forecast estimand, leakage-safe patient-level evaluation, strong simple comparators, parameter identifiability/calibration, external validation, and negative results.

Novelty remains unestablished. Adjacent literature already covers future-mask forecasting, longitudinal tumor-volume changes, and dose-linked recurrence geography. A defensible contribution might be the incremental value of the observed local GTV annotation trajectory for next-visit released-label additions under provenance-aware evaluation; an exact literature and comparator audit is still required before calling it novel or unique.

## 10. Source basis

- Project status and Track B limitations: `START_HERE.md`.
- CFB release, candidate counts, ID and mask-lineage issues: `reports/GBM CFB v4 eligibility audit.md` and `reports/GBM CFB scoring readiness plan.md`.
- Experiment-design conflicts: `research_notes/GBM CFB scoring readiness/experiment.md` and `hypothesis.md`.
- Patient-level estimand and scale calculations: `research_notes/GBM CFB scoring readiness/stats.md`.
- Current CFB collection page: https://www.cancerimagingarchive.net/collection/cfb-gbm/ (DOI 10.7937/v9pn-2f72).
- CFB-GBM v2.0 methods paper (baseline t0 T1-Gd is the fixed registration reference): https://www.melba-journal.org/pdf/2026:038.pdf.
- Adjacent forecasting literature includes TaDiff: https://arxiv.org/abs/2309.05406; longitudinal post-treatment volume change: https://pmc.ncbi.nlm.nih.gov/articles/PMC9530762/; and dose-associated recurrence geography: https://pmc.ncbi.nlm.nih.gov/articles/PMC5628052/.

