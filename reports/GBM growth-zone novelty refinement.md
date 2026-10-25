# Track B: A narrower method worth testing

## Decision

The current “PDE reach map plus MRI features” is not unique as a general method. Published work already predicts voxelwise recurrence, uses tumor-cavity proximity, models longitudinal future tumor maps, and benchmarks PDE and learned recurrence maps. No novelty claim is justified from the current search.

A more specific candidate is a **transition-budgeted, tissue-aware spatial forecast**. It separates two questions that ordinary mask forecasting entangles:

1. How many new annotated core/enhancing voxels may appear by the next visit?
2. Given that predicted amount, where are those voxels most likely to appear?

The first quantity is estimated only from earlier visits. A spatial field allocates that predicted amount across baseline candidate voxels. This removes the future-mask top-N oracle and directly tests whether tissue-aware diffusion improves localization after the predicted amount is held fixed.

This is a candidate gap, not a proven unique contribution. It needs further citation chasing and a search of methods outside glioma imaging before any “new method” claim.

## Candidate method

For patient $i$, after baseline/follow-up alignment, let $T_{i1}$ and $T_{i2}$ be the earlier and current label-1/3 masks, and let $T_{i3}$ be the later outcome mask. These are annotation labels, not verified viable tumor. Define prior and target additions as

\[
G_{i,12}=T_{i2}\setminus T_{i1},\qquad
G_{i,23}=T_{i3}\setminus T_{i2}.
\]

Let $X_{i1}$ and $X_{i2}$ be candidate regions defined at visits 1 and 2, respectively, each using only information available at that visit and excluding its current tumor, cavity, CSF, and non-brain voxels. Use the previous transition, its elapsed time, current tumor volume, and next-visit horizon to estimate the expected number of additions in $X_{i2}$:

\[
\widehat K_i = g\!\left(|G_{i,12}\cap X_{i1}|,\;\Delta t_{i,12},\;\Delta t_{i,23},\;|T_{i2}|\right),
\]

where $g$ is a low-dimensional count model fit only on development patients, with shrinkage toward the cohort rate. Freeze its form before scoring. The exact count model must be selected after an eligibility audit; a simple rate baseline and a shrinkage model should be compared in development. If prior additions are too unstable under registration or annotation perturbations, stop this factorization rather than treat them as signal.

Construct a tissue-aware field $u_i(x)$ from a screened diffusion equation on the visit-2 tissue domain:

\[
-\nabla\!\cdot(A_i(x)\nabla u_i(x))+\kappa^2u_i(x)=q_i(x),
\]

with $q_i(x)$ fixed from the visit-2 label-1/3 source, $A_i(x)$ positive and defined from visit-2 tissue/FLAIR information, and explicit no-flux handling of the cavity and outer tissue boundary. Fix global scale and choose any small number of coefficients inside development folds. Interpret $u_i$ as a spatial ranking feature; it is not a cell-density estimate or biological growth law.

Fit a patient-balanced spatial score $r_i(x)$ using visit-2 features, then choose a patient-specific intercept $a_i$ so the output probabilities satisfy

\[
\sum_{x\in X_{i2}}\operatorname{logit}^{-1}(a_i+r_i(x))=\widehat K_i.
\]

This is a non-oracle expected-count constraint: the amount comes from prior visits and training patients, never from $T_{i3}$. The count model must be fit inside each training fold; training predictions used to fit the spatial model should be cross-fitted. If $\widehat K_i$ is outside $[0, |X_{i2}|]$, use a prespecified boundary rule and report it. Apply the same $\widehat K_i$ constraint to every comparator so only the spatial ranking differs. Compare:

- a feature-matched logistic baseline using distance, FLAIR, white-matter prior, label-2 indicator, cavity distance, and interval;
- a weighted-geodesic score using the same tissue costs and cavity rule;
- the screened-diffusion feature added to the same logistic baseline;
- simple copy/Euclidean-distance references.

If treatment variables are available and reliable, include them symmetrically in the learned arms. Otherwise state the confounding limitation.

## What could be new—and what is already known

| Related work | What it already does | Remaining distinction to test |
|---|---|---|
| Huang et al., voxelwise recurrence prediction with proximity-estimation-coupled SVM | Uses cavity and stem-cell-niche proximity with multiparametric MRI to produce voxelwise recurrence predictions. | It does not use the proposed prior-interval count forecast to constrain a next-visit addition map. The present project’s target is still an imaging-label transition, not clinical recurrence. [PubMed record](https://pubmed.ncbi.nlm.nih.gov/34963559/) (abstract-only here). |
| Chougule et al. (2023), postoperative voxel-based radiomics | Predicts local recurrence regions from multiparametric postoperative MRI. | It does not test the proposed count/location decomposition with an equal-count spatial comparison. [Article](https://doi.org/10.3390/cancers15061894) (full-text excerpts checked; complete text not re-read for this addendum). |
| GliomaDiff (2025 preprint) | Predicts time-conditioned future tumor probability maps from longitudinal scans, using signed-distance maps and learned deformation. | It is a learned future-mask model; the candidate here isolates a low-dimensional predicted addition budget from a tissue-aware spatial allocation feature. This is a methodological distinction, not proof of priority. [arXiv paper](https://arxiv.org/abs/2509.10824) (full text opened). |
| PREDICT-GBM (2026) | Benchmarks PDE and data-driven recurrence maps; selects top-scoring voxels to match a fixed standard-of-care target volume. | Its volume is fixed to a radiotherapy comparator. The candidate’s $\widehat K_i$ is estimated from earlier patient-specific scan transitions and is used to forecast the next scan’s observed additions. This difference needs careful verification against the complete methods and supplement. [Article](https://doi.org/10.1038/s41746-026-03194-0) (search-indexed article passages checked; direct open was blocked). |
| MambaDiff-Net (2026) | Uses preoperative/postoperative MRI for patient-level early recurrence-risk prediction with internal and external cohorts. | The reported outcome is patient-level recurrence risk, not the proposed per-patient next-scan addition allocation. [Article](https://www.mdpi.com/2306-5354/13/6/668) (abstract and indexed method passages checked; full text not independently opened in this search). |

This focused search found no exact match for “prior-transition-estimated addition count + screened tissue-diffusion allocation + matched geodesic control” in the sources above. That is only a provisional gap statement. Search terms included longitudinal glioma mask prediction, voxelwise recurrence mapping, postoperative cavity proximity, PDE recurrence maps, fixed-volume prediction, and spatial recurrence forecasting. The search was not systematic across all oncology, spatial statistics, or constrained-probability literature.

## Evaluation that could falsify the method

- Use patient-level repeated/nested cross-validation on the existing MU cohort for development. Because the same scan pairs were examined in earlier Track B work, call all MU estimates exploratory; do not call a newly drawn split independent confirmation.
- Predict $G_{i,23}\cap X_{i2}$ only from visits 1 and 2. No visit-3 image, mask, registration adjustment, or positive count may enter a predictor, ROI, or $\widehat K_i$.
- Primary comparison: paired patient-level change in PR-AUC or average precision between diffusion and the feature-matched baseline, with the same predicted count budget. Report patient-level AUC as a secondary ranking metric because it can look good under sparse positives.
- Also report count forecast error, precision/recall at the frozen predicted count, fraction of later additions outside the frozen ROI, per-patient results, and no-positive/no-negative cases. Bootstrap patients, not voxels.
- Compare against matched geodesic and negative controls: isotropic distance; flattened/permuted tissue costs that preserve the feature distribution; and cavity-open versus cavity-blocked field. The spatial method earns a physics interpretation only if the anatomically meaningful field beats these controls, not merely if smoothing beats raw distance.
- Perturb baseline masks and registration within prespecified plausible bounds. Report whether the predicted map and paired performance change materially. If alignment uncertainty is too large to distinguish method differences, stop.
- Do not set $\widehat K_i$ from the true future positive volume. A target count from the future scan is an oracle even when only used to choose top voxels.
- Before external/confirmatory claims, freeze the equations, preprocessing, comparators, thresholds, metrics, split, and code, then evaluate on a compatible cohort never used to choose the method. If none is available, retain the exploratory label.

## Success and stop rules

Do not define “success” as merely a positive point estimate. Choose a minimum useful difference before analysis and require its patient-level uncertainty interval to support that threshold. With roughly 50 development patients and fewer eligible patients after restrictions, the proposed $0.02$ AUC threshold has no established power basis. Reconsider the endpoint and practical margin before using it. A negative result is informative: the field should be dropped if it fails against the matched logistic or geodesic comparator, if gains disappear under mask/registration perturbations, or if the amount model is unstable.

The strongest defensible future claim, if all gates pass, would be: “In an exploratory analysis of previously examined public follow-up data, a prior-transition-budgeted screened-diffusion feature improved [prespecified patient-level metric] for predicting next-visit annotated additions over matched non-PDE and geodesic scores.” It would not show actual tumor-cell migration, clinical recurrence prediction, treatment benefit, or external generalization.

## Sources found in the focused search

- [Huang et al., voxelwise prediction using proximity-estimation-coupled SVM](https://pubmed.ncbi.nlm.nih.gov/34963559/)
- [Chougule et al., postoperative voxel-based radiomics](https://doi.org/10.3390/cancers15061894)
- [Kebaili et al., GliomaDiff longitudinal progression prediction](https://arxiv.org/abs/2509.10824)
- [PREDICT-GBM](https://doi.org/10.1038/s41746-026-03194-0)
- [MambaDiff-Net](https://www.mdpi.com/2306-5354/13/6/668)
- [RANO 2.0, explaining why scan-label change is not synonymous with tumor growth](https://doi.org/10.1200/JCO.23.01059)
