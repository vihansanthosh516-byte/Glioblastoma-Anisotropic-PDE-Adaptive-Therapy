# Track B growth-zone plan: mathematical and statistical review

**Scope.** Independent review of `Claude outputs/trackB_growth_zone_plan_claude.md`. I read the plan and searched primary research/methods sources. I did not edit the plan, run project scripts, or analyze patient data. The cohort size, split, segmentation audit, and past-test overlap below are plan-reported, not independently revalidated here.

## 1. What the “steady diffusion reach” field means mathematically

### Takeaway

The proposal is promising as an interpretable *spatial feature*, but “steady PDE reach field” is not yet a defined model. It must specify its source, operator, boundary conditions, tissue coefficients, and decay/scale. It also needs a matched geodesic comparator: diffusion fronts are already known to resemble weighted shortest-path fronts, with smoothing and barrier-related differences.

### Cited findings

- Tumor-front diffusion and weighted shortest-path/eikonal models have closely related front geometry. Bortfeld and Buti compare the equations directly: diffusion fronts are generally smoother, while obstacle effects and complex inhomogeneity can distinguish them. Their paper explicitly treats barriers and anisotropy. **Full text opened** (arXiv manuscript; published DOI: [10.1088/1361-6560/ac8043](https://doi.org/10.1088/1361-6560/ac8043)).
- Engwer and Wenske already derive a stationary approximation to a glioblastoma invasion model, with spatially varying diffusion and a real DTI example. It targets estimation of *current invisible invasion extent*, rather than next-scan new labels, but it means “stationary PDE” by itself is not a novelty claim. **Full text opened**: [10.1007/s00285-021-01563-9](https://doi.org/10.1007/s00285-021-01563-9).
- Akbari et al. used multiparametric MRI and machine learning to map infiltration and predict subsequent recurrence location, with a replication cohort; the model used preoperative imaging and included FLAIR and diffusion/perfusion modalities. **Abstract-only**: [10.1227/NEU.0000000000001202](https://doi.org/10.1227/NEU.0000000000001202).
- A 2023 postoperative radiomics study predicts local recurrence regions from multiparametric MRI and follow-up-defined labels. **Abstract/article-page material only; full text was not opened**: [10.3390/cancers15061894](https://doi.org/10.3390/cancers15061894).
- A 2025 *Radiotherapy and Oncology* article describes DeepGBM-Recure: peri-cavitary FLAIR segmentation, patient-level risk prediction, and voxel-level spatial mapping of recurrence hotspots, tested across four centers. **Abstract-only**: [10.1016/j.radonc.2025.111180](https://doi.org/10.1016/j.radonc.2025.111180).

### Inferences and recommendations

1. **Write down a well-posed equation before implementation.** For a stationary score, one defensible form is
   \[
   -\nabla\!\cdot(D(x)\nabla u(x))+\lambda u(x)=s(x),\qquad \lambda>0,
   \]
   on a specified brain domain, with a specified source `s` on/inside the baseline tumor and specified boundary conditions. The positive decay term makes the screened operator invertible under ordinary boundary choices. A pure steady Neumann Poisson equation with only a nonnegative tumor source is incompatible: integrating `-div(D grad u)=s` over the domain gives zero net boundary flux on the left but positive total source on the right. If the intended PDE is different, state it explicitly.
2. **Use a conservative variable-coefficient operator.** If tissue changes `D(x)`, the intended operator is `div(D grad u)`, not simply `D(x) Δu`; the latter omits flux changes at tissue interfaces. A finite-volume or graph-Laplacian discretization makes the flux and cavity treatment explicit.
3. **Resolve parameter non-identifiability.** If `D(x)=d0 Drel(x)` and `lambda` are both fitted, then
   \[
   (d_0 L+\lambda I)^{-1}s=d_0^{-1}(L+(\lambda/d_0)I)^{-1}s.
   \]
   In a logistic model with a free coefficient on `u`, the global factor `1/d0` is absorbed by that coefficient. The shape is governed by relative tissue conductivities and `lambda/d0`, not separately identified absolute `d0` and `lambda`. Fit dimensionless ratios/length scales or fix one scale; do not report two absolute physical values as independently learned unless the model identifies them.
4. **Do not call a scalar tissue map “fiber anisotropy.”** A higher scalar diffusivity in white-matter atlas voxels creates spatially varying but locally isotropic spread. Directional along-tract spread requires a tensor aligned to fiber directions (typically from DTI), which this plan does not list as an input. Use “white-matter-weighted” or “tissue-weighted” unless a directional tensor is actually available.
5. **Interpret no-flow correctly.** A no-flux cavity boundary `n·D∇u=0` is a reflecting boundary for the field: it blocks flux through the cavity, but permits routes around it. It does not establish that tumor cannot recur at the cavity edge. Keep the hard-barrier ablation, but add a soft-permeability/attenuation version and a cavity-distance baseline. Confirm the segmentation label really denotes surgical cavity/CSF for every included case.
6. **Give the comparator the same information.** As written, H uses white matter and a cavity barrier, while B1 has only Euclidean tumor distance and FLAIR inputs. H beating B1 could therefore be due to extra tissue/cavity information rather than the PDE. Make the principal non-PDE comparator feature-matched: distance + label 2 + FLAIR intensity + white-matter prior + cavity information. Then add a weighted fast-marching/geodesic feature using the *same* tissue costs and barriers. The clean claim is whether the PDE’s path-averaging/smoothing adds value over matched covariates and a weighted shortest path.
7. **A stationary score has no forecast horizon.** Follow-up intervals can differ substantially. Include elapsed days as a prespecified covariate in every learned arm, stratify/restrict to a defined interval window, or use a time-dependent heat field `u(x, Δt)` and make the horizon explicit. If H remains steady, describe it as a static spatial prior for the next observed scan, not a time-calibrated tumor-growth forecast.
8. **Treat synthetic tests as implementation checks, not proof.** If synthetic data are generated by the exact H equation, H is favored by construction (“inverse crime”). Use an independent generator, such as a discrete random walk/agent process on the same tissue/cavity graph, to test whether a separately implemented reach model recovers known spatial ordering. Include isotropic/no-barrier and tissue/barrier scenarios. Keep synthetic results separate from real-patient evidence.

### Gaps

The plan does not define the PDE, source mask, domain, boundary conditions at cavity/CSF/skull, field normalization, physical units, or the two shared values. It is also unclear whether tissue inputs are an atlas, patient-specific segmentation, or DTI tensors. These are specification blockers, not coding details.

## 2. Endpoint, registration, and uncertainty

### Takeaway

The current label subtraction can support a narrowly worded prediction task, but it cannot by itself identify biological growth or clinical recurrence. Registration error, postoperative segmentation variability, scan interval, and a restricted candidate region can dominate a voxelwise result. The analysis unit must remain the patient.

### Cited findings

- RANO 2.0 notes that post-radiotherapy changes can mimic progression; in the first 12 weeks after chemoradiotherapy, a repeat MRI or pathology is used to confirm suspected progression. Thus one scan-to-scan change is not synonymous with true progression. **Abstract/full-text excerpt accessed; full article not independently read here**: [10.1200/JCO.23.01059](https://doi.org/10.1200/JCO.23.01059).
- Visser et al. studied manual longitudinal glioma segmentations at preoperative, postoperative, and progression stages. For enhancing GBM tumor, spatial agreement was poor postoperatively (GCI 0.32) despite high volume ICC; they caution that postoperative spatial overlap is a concern. **Abstract-only**: [10.1016/j.nicl.2019.101727](https://doi.org/10.1016/j.nicl.2019.101727).
- Canalini et al. evaluated longitudinal brain-tumor registration. They report that resection, mass effect, and regrowth deform near-tumor anatomy; in their setting linear registration was not accurate enough, pathology masking improved surrounding-brain registration, and T1-CE was more informative than FLAIR for their registration task. **Full text opened**: [10.3389/fnimg.2022.977491](https://doi.org/10.3389/fnimg.2022.977491).
- Han et al. proposed pathology-aware registration specifically for preoperative and post-recurrence GBM images, noting the missing-correspondence problem caused by pathology/deformation. **Abstract-only**: [10.1007/978-3-030-11723-8_10](https://doi.org/10.1007/978-3-030-11723-8_10).
- Saito and Rehmsmeier show why precision-recall plots can be more informative than ROC plots under severe class imbalance. **Abstract-only**: [10.1371/journal.pone.0118432](https://doi.org/10.1371/journal.pone.0118432).
- ROC methods for clustered/repeated measurements require the sampling/inference to reflect patient clustering; Li and Zhou develop paired clustered ROC comparisons. **Abstract-only**: [10.1198/016214508000000364](https://doi.org/10.1198/016214508000000364). Obuchowski’s earlier paper likewise explains that multiple results from one patient violate ordinary independence assumptions. **Abstract-only**: [PubMed 9192452](https://pubmed.ncbi.nlm.nih.gov/9192452/).

### Inferences and recommendations

1. **Name the outcome literally.** Define `T_k(x)=1` when the dataset’s validated label map is in `{1,3}` at visit `k`, and
   \[
   Y_{ij}(x)=1\{T_{k+1}(x)=1,\ T_k(x)=0\}.
   \]
   Call this “newly observed label-1/3 tissue between scans,” not tumor growth, recurrence, progression, or a confirmed biological event. State that treatment effect/pseudoprogression and segmentation/registration changes are not separated.
2. **Make the ROI conditionality visible.** Define the candidate set using visit-k information only, e.g. baseline brain voxels within 20 mm of `T_k ∪ cavity_k` and outside `T_k`. Report, for every patient and overall, the fraction of all new label-1/3 voxels captured by that ROI. A distant future positive outside the ROI is a missed event, not a no-growth patient.
3. **Do not use “same 1-mm grid” as proof of alignment.** The plan’s cohort-level rule (“if median shift exceeds 2 mm, register all”) can miss individual bad pairs and does not measure local residual error. Set a patient-level QC rule before outcome scoring; inspect landmarks/anatomical correspondence near and away from the tumor. Compare unregistered, rigid, and (only if validated) pathology-aware deformable alignment as a sensitivity analysis. Do not deform the target merely to maximize overlap. Use future images only to define/evaluate the target, not to create prediction features.
4. **Test label-boundary robustness.** Registration and segmentation uncertainty create uncertain boundary voxels. Report primary metrics plus a prespecified tolerance analysis (e.g., derive a boundary band from measured landmark-registration error; separately report the error band rather than silently dilating/eroding the labels). If no landmarks or repeated expert masks are available, say that the boundary uncertainty is unmeasured; use small erosion/dilation perturbations only as a stress test, not as a noise estimate.
5. **Define what counts as a no-growth case.** AUC is undefined if an ROI has no positive or no negative voxels. Report the number excluded from AUC and why. Distinguish no positive inside the 20-mm ROI from no positive anywhere in the follow-up brain. Keep no-positive cases for false-alert burden using a threshold chosen only on development patients.
6. **Keep the independent unit as patient.** The plan says to use all eligible pairs per patient. Pairs from one patient are correlated and horizons/treatment stages differ. Prefer one prespecified eligible pair per patient for the held-out primary estimand. If all pairs remain, first calculate pair-level metrics, aggregate them to one patient-level value by a predefined rule, and bootstrap whole patients (retaining all their pairs and paired H/B1 results). Never bootstrap voxels or treat each scan pair as an independent patient.
7. **Use more than ROC-AUC.** Equal-weight mean patient AUC is a sensible primary ranking metric among patients with both classes, but it can hide poor precision when positives are rare and says nothing about calibration or spatial contiguity. Add per-patient average precision/PR-AUC and a prespecified fixed alert-volume metric. Keep the current `top-N`, `N = true new volume`, only as an explicitly *oracle-size localization diagnostic*: it uses the answer’s volume and is not an operational prediction. With `N=P` positive voxels, precision@N is also recall@N, so it is not an independent second line of evidence.
8. **Make the confidence rule match the claim.** The current rule “point estimate ≥0.02 and 95% CI lower bound >0” supports a positive effect whose point estimate reaches 0.02; it does not support that the true improvement is at least 0.02. If 0.02 is a justified minimum worthwhile gain, require the lower confidence bound to exceed 0.02. Otherwise call 0.02 a prespecified point-estimate benchmark and report the paired effect and CI without claiming a guaranteed useful effect. A paired patient bootstrap should resample patient IDs and preserve paired H/B1 scores; if all visits are retained, resample patient clusters. With 32 planned held-out patients, the achievable precision is unknown until patient-level variance and evaluable counts are known; do not promise power from the nominal cohort count.
9. **Account for model-selection uncertainty.** All transforms, feature standardization, tissue-speed/decay values, logistic coefficients, B2 settings, and any alert threshold must be fit/selected within development folds. Voxel count does not create millions of independent training samples: use patient-balanced fitting weights and report patients as the effective independent sample. The untouched holdout is scored once after the pipeline is frozen. If prior work already inspected these same patients/results, call it a held-out split, not a fresh external test.

### Gaps

No pair interval distribution, per-patient positive/no-growth count, ROI coverage, target-label inter-rater estimate, registration residual, treatment timeline, or held-out patient AUC variance has been checked in this review. Therefore I cannot confirm that 50/32 is adequately powered or that the 20-mm restriction captures the relevant outcomes.

## 3. A feasible and more defensible comparison

### Takeaway

The strongest version is a patient-held-out test of whether a *matched, tissue-weighted path-ensemble feature* adds value beyond both feature-matched logistic prediction and weighted shortest paths. Describe the target as next-scan observed labels. Treat the PDE as a candidate model, not a novel biological mechanism. Existing work already maps recurrence spatially, including postoperative and peri-cavity approaches.

### Cited findings

- Prior spatial recurrence prediction is real, not hypothetical: Akbari 2016 predicts subsequent recurrence location from preoperative multiparametric MRI; the 2023 postoperative radiomics article predicts future local recurrence regions; DeepGBM-Recure (2025) maps peri-cavitary spatial recurrence risk from postoperative MRI. These differ in cohort, timing, features, and models from this plan, but rule out novelty claims based only on “predict where rather than how big,” FLAIR, a cavity region, or a spatial map. **Access:** Akbari and DeepGBM-Recure abstract-only; 2023 article full text not opened (links above).
- Tumor-front literature already compares diffusion and shortest-path models and formalizes the anisotropic eikonal relation. **Full text opened:** Bortfeld and Buti, DOI above.

### Inferences and recommendations

**Define the primary scientific question as:** “Among eligible patients with longitudinal post-treatment MRI, does adding a shared tissue-weighted diffusion/reach score to a feature-matched spatial baseline improve patient-level ranking of next-scan newly observed core/enhancing voxels inside a prespecified ROI, and is that improvement robust to path-model and alignment/label uncertainty?” This says what data can test without asserting tumor biology.

**Recommended frozen comparator ladder:**

1. Euclidean distance-only baseline.
2. Feature-matched logistic baseline: Euclidean distance, visit-k label 2, visit-k FLAIR intensity, white-matter prior, cavity distance/indicator, and interval days.
3. Weighted shortest-path/fast-marching model using the same tissue costs and cavity assumptions.
4. Screened-diffusion/graph-resolvent H using the same inputs; add it to the feature-matched logistic model.
5. H ablations: remove label-2 speed-up, remove white matter, remove cavity barrier, and compare hard vs soft cavity permeability. Include a simple matched cavity-distance baseline so “cavity helps” is not mistaken for “PDE helps.”
6. B2 radiomics as a secondary comparator only if patient-level nested tuning and patient-balanced sampling are feasible with 50 development patients. Do not let millions of voxels be presented as a large independent sample.

**Primary analysis:** one prespecified scan pair per patient; patient-level AUC over the fixed candidate ROI; mean paired `AUC(H) − AUC(feature-matched baseline)` over evaluable patients; a paired patient-resampling CI; report evaluable n, pair interval, class prevalence, ROI coverage, per-patient differences, and leave-one-patient-out influence. Add PR-AUC and fixed-budget precision/recall; evaluate no-positive cases separately at a threshold selected only in development. Keep `top-N true volume` clearly labeled oracle-only. Report 0.02 as a point-estimate benchmark unless its lower confidence bound is also above 0.02.

**Must-pass checks before any held-out scoring:** verify eligible-patient/pair counts and label ontology; check prior use of the proposed holdout; freeze a patient hash split; quantify interval and ROI coverage; define registration QC and label-tolerance sensitivity; freeze the PDE and matched baselines; test synthetic behavior using an independent generator; then freeze code, features, tuning bounds, metrics, and success language.

**What could be genuinely new here (not yet established):** a careful empirical answer, on this already-used public longitudinal cohort, to whether a shared tissue/cavity-aware path-ensemble score adds reproducible localization value beyond equivalent input features and weighted geodesic distance, after accounting for patient clustering and boundary uncertainty. That is a bounded methodological contribution. The research checked here does not establish that this combination is unprecedented, that it will win, or that a 0.02 margin is meaningful.

### Gaps and search record

Focused searches covered: “glioblastoma recurrence spatial prediction FLAIR diffusion MRI,” “voxel-based postoperative recurrence MRI radiomics,” “shortest path and diffusion tumor fronts,” “stationary GBM invasion PDE,” “longitudinal brain tumor registration pathology mask,” “glioma longitudinal segmentation inter-rater,” and “clustered/longitudinal ROC AUC paired bootstrap.” This is a focused prior-art screen, not a systematic review. The search found close prior spatial-localization work but did not complete forward/backward citation tracking, dataset-specific overlap checks, or a full-text comparison of all candidate studies. Novelty remains unproven.

Source access log (status refers to what I personally opened for this note):

- Bortfeld & Buti (2022), *Modeling the propagation of tumor fronts with shortest path and diffusion models*: [DOI](https://doi.org/10.1088/1361-6560/ac8043) — **full text opened** (arXiv PDF).
- Engwer & Wenske (2021), *Estimating the extent of glioblastoma invasion*: [DOI](https://doi.org/10.1007/s00285-021-01563-9) — **full text opened** (Springer article text).
- Canalini et al. (2022), *Quantitative evaluation … registration of longitudinal data acquired during brain tumor treatment*: [DOI](https://doi.org/10.3389/fnimg.2022.977491) — **full text opened** (Frontiers PDF).
- Akbari et al. (2016), *Imaging Surrogates of Infiltration … Subsequent Location of Recurrence*: [DOI](https://doi.org/10.1227/NEU.0000000000001202) — **abstract-only** (PubMed).
- Azzouz et al. (2023), *Predicting Regions of Local Recurrence … Postoperative MRI*: [DOI](https://doi.org/10.3390/cancers15061894) — **abstract/article-page only; full text not opened** (MDPI PDF request was rate-limited).
- Li et al. (2025), *An intelligent postoperative management system … recurrence spatial mapping* (DeepGBM-Recure): [DOI](https://doi.org/10.1016/j.radonc.2025.111180) — **abstract-only** (ScienceDirect).
- Visser et al. (2019), *Inter-rater agreement in glioma segmentations on longitudinal MRI*: [DOI](https://doi.org/10.1016/j.nicl.2019.101727) — **abstract-only** (PubMed).
- Wen et al. (2023), *RANO 2.0*: [DOI](https://doi.org/10.1200/JCO.23.01059) — **abstract/source excerpt only; full article not read here**.
- Saito & Rehmsmeier (2015), *The Precision-Recall Plot Is More Informative than the ROC Plot…*: [DOI](https://doi.org/10.1371/journal.pone.0118432) — **abstract-only**.
- Li & Zhou (2008), *A Unified Approach to Nonparametric Comparison of ROC Curves for Longitudinal and Clustered Data*: [DOI](https://doi.org/10.1198/016214508000000364) — **abstract-only**.
- Obuchowski (1997), *Nonparametric analysis of clustered ROC curve data*: [PubMed](https://pubmed.ncbi.nlm.nih.gov/9192452/) — **abstract-only**.
- Han et al. (2019), *Patient-Specific Registration of Pre-operative and Post-recurrence Brain Tumor MRI Scans*: [DOI](https://doi.org/10.1007/978-3-030-11723-8_10) — **abstract-only** (author lab page).
