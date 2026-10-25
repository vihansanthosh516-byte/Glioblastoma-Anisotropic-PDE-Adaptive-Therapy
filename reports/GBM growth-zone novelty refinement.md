# Track B: A more specific hybrid candidate

## Decision

The current “PDE reach map plus MRI features” is not unique as a general method. Published work already predicts voxelwise recurrence, uses tumor-cavity proximity, models longitudinal future tumor maps, benchmarks PDE and learned recurrence maps, and combines mechanistic models with machine learning. In particular, BaM³ uses mechanistic tumor-growth outputs to guide machine-learning predictions, GliODIL combines MRI data and a reaction-diffusion PDE in one inverse model, and a 2025 preprint combines a mechanistic tumor-burden model with guided diffusion generation. “Hybrid” alone is not a novelty claim.

A more specific candidate is an **uncertainty-gated, residual-corrected, transition-budgeted diffusion forecast** (working description, not a claimed named invention). It fuses four parts into one output:

1. A screened anisotropic diffusion PDE supplies an interpretable tissue-aware spatial prior.
2. A small supervised model learns a residual correction from MRI and prior-visit features; it learns where the PDE prior is systematically wrong.
3. The correction uses the previous interval's local addition map and is shrunk where either that map or the PDE map changes sharply under plausible mask, registration, or tissue-map perturbations.
4. A time-aware count model predicts how many new annotated voxels may appear, and a mass constraint calibrates the final map to that predicted amount without using the future mask.

The components interact: the PDE shapes the initial map, the learner corrects that map, uncertainty controls how much the correction can move it, and the temporal count model sets the map’s total expected positive mass. This is a genuine mechanistic-statistical hybrid design. It is not yet established as unique or useful.

The possible gap is this exact coupling for predicting next-visit annotated additions on longitudinal postoperative MRI, with matched non-PDE, geodesic, and perturbation controls. The search has not established priority. A wider search across mathematical oncology, spatial statistics, and machine-learning calibration is still required before any novelty claim.

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

Generate a fixed, prespecified set of plausible baseline-mask, registration, and tissue-map perturbations. Let $\ell_i^{(b)}(x)=\log(u_i^{(b)}(x)+\epsilon)$ be the log PDE field for perturbation $b$, and let $\bar\ell_i(x)$ be its mean. Also perturb the registered prior-addition map $G_{i,12}$ and measure its local instability. Combine and scale the log-field variance and prior-map instability using development data to obtain $V_i^*(x)\in[0,1]$. These are stability summaries, not calibrated biological uncertainty.

Let $z_i(x)$ contain visit-2 MRI and geometry features available to the non-PDE baseline. Let $m_{i,12}(x)$ encode the locally registered prior-interval additions, computed only from visits 1 and 2. A low-complexity supervised learner estimates a shared residual correction $h_\theta(z_i(x),m_{i,12}(x),\bar\ell_i(x))$ from development patients. Shrink that correction in unstable voxels:

\[
\eta_i(x)=\beta_{\mathrm{PDE}}\bar\ell_i(x)
+\exp(-\lambda V_i^*(x))\,h_\theta(z_i(x),m_{i,12}(x),\bar\ell_i(x)).
\]

Here $\beta_{\mathrm{PDE}}$, $\lambda$, and the residual-model settings are fixed or selected only inside patient-level development folds. Use a small ridge-regularized logistic/GAM residual, not a high-capacity network, unless a later data audit supports more complexity. The learned term is trained to correct the PDE-based log score using MRI and the prior-visit spatial transition, not merely concatenated as another unstructured feature. If the local transition map is unreliable under perturbation, the stability gate reduces the learned correction instead of trusting it.

Convert $\eta_i(x)$ to a patient-specific probability map and choose an intercept $a_i$ so the output satisfies

\[
\sum_{x\in X_{i2}}\operatorname{logit}^{-1}(a_i+\eta_i(x))=\widehat K_i.
\]

This is a non-oracle expected-count constraint: the amount comes from prior visits and training patients, never from $T_{i3}$. The count model must be fit inside each training fold; training predictions used to fit the spatial model should be cross-fitted. If $\widehat K_i$ is outside $[0, |X_{i2}|]$, use a prespecified boundary rule and report it. Apply the same $\widehat K_i$ constraint to every comparator so only the spatial ranking differs. Compare:

- a feature-matched logistic baseline using distance, FLAIR, white-matter prior, label-2 indicator, cavity distance, and interval;
- a weighted-geodesic score using the same tissue costs and cavity rule;
- the PDE-only model;
- the PDE prior plus learned residual correction;
- the hybrid without the previous-transition feature;
- the hybrid without the stability gate;
- the full uncertainty-gated, residual-corrected hybrid;
- simple copy/Euclidean-distance references.

If treatment variables are available and reliable, include them symmetrically in the learned arms. Otherwise state the confounding limitation.

## What could be new—and what is already known

| Related work | What it already does | Remaining distinction to test |
|---|---|---|
| Huang et al., voxelwise recurrence prediction with proximity-estimation-coupled SVM | Uses cavity and stem-cell-niche proximity with multiparametric MRI to produce voxelwise recurrence predictions. | The abstract does not establish whether prior-interval addition counts constrain the spatial map; the full methods need review. The present target is still an imaging-label transition, not clinical recurrence. [PubMed record](https://pubmed.ncbi.nlm.nih.gov/34963559/) (abstract-only here). |
| Chougule et al. (2023), postoperative voxel-based radiomics | Predicts local recurrence regions from multiparametric postoperative MRI. | The full methods need review to determine whether it separates and couples predicted amount with voxelwise location. [Article](https://doi.org/10.3390/cancers15061894) (full-text excerpts checked; complete text not re-read for this addendum). |
| GliomaDiff (2025 preprint) | Predicts time-conditioned future tumor probability maps from longitudinal scans, using signed-distance maps and learned deformation. | It is a learned future-mask model; the candidate here isolates a low-dimensional predicted addition budget from a tissue-aware spatial allocation feature. This is a methodological distinction, not proof of priority. [arXiv paper](https://arxiv.org/abs/2509.10824) (full text opened). |
| PREDICT-GBM (2026) | Benchmarks PDE and data-driven recurrence maps; selects top-scoring voxels to match a fixed standard-of-care target volume. | Its volume is fixed to a radiotherapy comparator. The candidate’s $\widehat K_i$ is estimated from earlier patient-specific scan transitions and is used to forecast the next scan’s observed additions. This difference needs careful verification against the complete methods and supplement. [Article](https://doi.org/10.1038/s41746-026-03194-0) (search-indexed article passages checked; direct open was blocked). |
| MambaDiff-Net (2026) | Uses preoperative/postoperative MRI for patient-level early recurrence-risk prediction with internal and external cohorts. | The reported outcome is patient-level recurrence risk, not the proposed per-patient next-scan addition allocation. [Article](https://www.mdpi.com/2306-5354/13/6/668) (abstract and indexed method passages checked; full text not independently opened in this search). |
| BaM³ (2021) | Uses a mechanistic model's prediction as a Bayesian prior and combines it with a data-derived probability distribution; its glioma demonstration is in silico. | This rules out “mechanistic model + ML correction” as the novelty. The candidate differs only if the voxel-level prior-transition feature, perturbation-stability gate, and next-visit expected-count constraint form a genuinely unreported coupling. [Article](https://www.nature.com/articles/s43856-021-00020-4) (full article opened). |
| GliODIL (2025) | Jointly fits MRI observations and a Fisher–KPP reaction-diffusion model through a physics/data-informed discrete loss; evaluated for recurrence planning. | It is already a physics/data hybrid. The candidate’s possible distinction is an explicit PDE spatial prior plus learned residual correction, perturbation-stability gate, and prior-visit count constraint, rather than a PDE-regularized inverse fit. [Article](https://doi.org/10.1038/s41467-025-60366-4) (article text checked). |
| Laslo et al. (2025 preprint) | Combines mechanistic ODE estimates of tumor burden with guided diffusion generation of future MRI. | This is another direct hybrid precedent. Its output is generated future imaging; the candidate forecasts an annotated transition map with a PDE spatial prior, learned residual, and non-oracle count constraint. [arXiv paper](https://arxiv.org/abs/2509.09610) (abstract checked). |

This focused search found no exact match for “prior-transition count forecast + screened-diffusion prior + learned residual correction gated by perturbation stability + expected-count-constrained voxel probabilities” in the sources above. That is only a provisional gap statement. Search terms included longitudinal glioma mask prediction, voxelwise recurrence mapping, postoperative cavity proximity, PDE recurrence maps, fixed-volume prediction, spatial recurrence forecasting, physics-informed machine learning, mechanistic residual learning, and hybrid tumor-growth prediction. The search was not systematic across all oncology, spatial statistics, uncertainty calibration, or constrained-probability literature.

## Evaluation that could falsify the method

- Use patient-level repeated/nested cross-validation on the existing MU cohort for development. Because the same scan pairs were examined in earlier Track B work, call all MU estimates exploratory; do not call a newly drawn split independent confirmation.
- Predict $G_{i,23}\cap X_{i2}$ only from visits 1 and 2. No visit-3 image, mask, registration adjustment, or positive count may enter a predictor, ROI, or $\widehat K_i$.
- Primary comparison: paired patient-level change in PR-AUC or average precision between diffusion and the feature-matched baseline, with the same predicted count budget. Report patient-level AUC as a secondary ranking metric because it can look good under sparse positives.
- Also report count forecast error, precision/recall at the frozen predicted count, fraction of later additions outside the frozen ROI, per-patient results, and no-positive/no-negative cases. Ablate the amount model too, so the test shows whether the predicted count constraint adds value. Bootstrap patients, not voxels.
- Compare against matched geodesic and negative controls: isotropic distance; flattened/permuted tissue costs that preserve the feature distribution; cavity-open versus cavity-blocked field; PDE-only; no prior-transition feature; and no stability gate. The full hybrid earns a contribution only if each added layer helps beyond its ablation, not merely if the combined score beats raw distance.
- Perturb baseline masks and registration within prespecified plausible bounds. Report whether the predicted map and paired performance change materially. If alignment uncertainty is too large to distinguish method differences, stop.
- Do not set $\widehat K_i$ from the true future positive volume. A target count from the future scan is an oracle even when only used to choose top voxels.
- Before external/confirmatory claims, freeze the equations, preprocessing, comparators, thresholds, metrics, split, and code, then evaluate on a compatible cohort never used to choose the method. If none is available, retain the exploratory label.

## Success and stop rules

Do not define “success” as merely a positive point estimate. Choose a minimum useful difference before analysis and require its patient-level uncertainty interval to support that threshold. With roughly 50 development patients and fewer eligible patients after restrictions, the proposed $0.02$ AUC threshold has no established power basis. Reconsider the endpoint and practical margin before using it. A negative result is informative: the field should be dropped if it fails against the matched logistic or geodesic comparator, if gains disappear under mask/registration perturbations, or if the amount model is unstable.

The strongest defensible future claim, if all gates pass, would be: “In an exploratory analysis of previously examined public follow-up data, a transition-budgeted, uncertainty-gated diffusion/residual hybrid improved [prespecified patient-level metric] for predicting next-visit annotated additions over matched non-PDE, PDE-only, and geodesic scores.” It would not show actual tumor-cell migration, clinical recurrence prediction, treatment benefit, or external generalization.

## Sources found in the focused search

- [Huang et al., voxelwise prediction using proximity-estimation-coupled SVM](https://pubmed.ncbi.nlm.nih.gov/34963559/)
- [Chougule et al., postoperative voxel-based radiomics](https://doi.org/10.3390/cancers15061894)
- [Kebaili et al., GliomaDiff longitudinal progression prediction](https://arxiv.org/abs/2509.10824)
- [Seyedpour et al., BaM³ mechanistic model and machine learning combination](https://www.nature.com/articles/s43856-021-00020-4)
- [Laslo et al., mechanistic learning with guided diffusion models](https://arxiv.org/abs/2509.09610)
- [PREDICT-GBM](https://doi.org/10.1038/s41746-026-03194-0)
- [MambaDiff-Net](https://www.mdpi.com/2306-5354/13/6/668)
- [RANO 2.0, explaining why scan-label change is not synonymous with tumor growth](https://doi.org/10.1200/JCO.23.01059)
