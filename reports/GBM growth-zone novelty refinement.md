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

## Additional prior art found during implementation (2026-10-24)

Rai, Saraswat, and Rai (2026) is a close new match for the **uncertainty-gated residual** component and the dataset. The full 24-page PDF was opened. Their method uses MU-Glioma-Post, makes next-scan patient-level progression predictions from longitudinal radiomics and treatment-history features, adds an uncertainty-gated residual correction and tumor-kinetics expert, and evaluates with nested patient-grouped cross-validation. Their target is a binary increase of at least 25% in tumor volume; they do not describe voxelwise growth-location maps or a screened spatial diffusion prior. This means neither "uncertainty-gated residual" nor patient-grouped next-scan prediction is a new contribution here. The candidate distinction, if any, is the spatial allocation problem and its coupling to a non-oracle addition-count forecast; that distinction still needs a direct methods comparison and a wider search.

- Rai, A. K., Saraswat, S., & Rai, A. K. (2026). *Uncertainty Gated Residual Learning with Patient Level Bagging for Longitudinal Glioma Progression Prediction*. Journal of Intelligent Decision Making and Information Science, 3(8s), 1713–1736. https://doi.org/10.59543/jidmis.v3.2236. **Full text opened.**

The open Track B prototype should therefore be treated as an exploratory methods candidate, not as a novelty-validated or guaranteed-superior method. Its stability gate currently measures sensitivity to chosen one-voxel perturbations; it is not calibrated uncertainty. The code must retain this limitation in reports.

### Tissue-map feasibility check

The MU resource paper says the scans are skull-stripped, resampled, and coregistered to SRI24 space. SRI24 publicly provides 1-mm NIfTI tissue probability maps for CSF, gray matter, and white matter; NITRC lists the atlas under CC-BY-SA and provides a 15-MB tissue-map package. This could support a white-matter-informed PDE and an explicit CSF barrier without running a new tissue segmenter. It is not yet integrated: the normal-control atlas probabilities may be wrong near a distorted tumor or cavity, so compare atlas-masked, atlas-weighted, and no-atlas variants and report failures. The current prototype excludes label-4 resection-cavity voxels, but it does **not** have an explicit CSF mask and must not claim that CSF is blocked.

- SRI24 atlas information and license: [NITRC SRI24 resource page](https://www.nitrc.org/projects/sri24) and [SRI24 atlas paper](https://doi.org/10.1002/hbm.20906). **Resource page opened; paper metadata and tissue-map description verified from the indexed primary article.**

## Prior-art refresh and current-source recheck (2026-10-24)

### What the newer literature changes

| Work | Verified method and evaluation | Consequence for this candidate |
|---|---|---|
| Liu, Trimboli & Li (2026), arXiv:2607.13877 | **Full arXiv text opened.** Adds a 3D learned residual to a reaction-diffusion forecast, updates that residual during recursive prediction, and couples the result to treatment scheduling. Evaluation uses a patient-data-informed synthetic testbed, not clinical longitudinal outcome scans. [Paper](https://arxiv.org/html/2607.13877) | “PDE plus learned residual” is not a new general methodological idea. The remaining difference is the low-capacity voxelwise MRI/transition correction, perturbation gate, and non-oracle count forecast for an annotated next-scan transition. That difference is not yet shown to be a meaningful novelty. |
| Gaillochet, Desrosiers & Lombaert (2026), RW-CP | **Full arXiv text opened.** Diffuses segmentation probabilities over a DINOv3-feature k-nearest-neighbor graph, then applies conformal risk control to FNR. Experiments use CAMUS, ACDC, and MSD pancreas; ACDC/MSD volumes are reduced to selected 2D slices. [Paper](https://arxiv.org/html/2601.18997) | Random-walk spatial smoothing plus conformal uncertainty is already proposed for static segmentation. It does not forecast a future tumor transition or use a reaction-diffusion growth prior. |
| Bereska, Karimi & Samavi (MIDL 2025), SACP | **Abstract and indexed paper passages inspected; direct OpenReview PDF access was blocked.** Uses class-conditional calibration and distance-weighted scoring to expand uncertainty near critical interfaces in multicenter pancreatic-tumor segmentation. [MIDL record](https://mlanthology.org/midl/2025/bereska2025midl-sacp/) | “Spatially varying conformal uncertainty near important boundaries” is also not a new general idea. It is a static segmentation method, not a longitudinal tumor-growth forecast. |
| Viti, Karabelas & Holler (ICLR 2026), CONSIGN | **Official abstract inspected.** Uses spatial groupings/decomposition to include spatial correlation in conformal segmentation sets across medical-image datasets. [Official ICLR record](https://proceedings.iclr.cc/paper_files/paper/2026/hash/47678cf4a2d9bce711b6d59337427dea-Abstract-Conference.html) | Spatial correlation in conformal segmentation is actively studied; the abstract does not establish a longitudinal glioma forecast or mechanistic prior. |
| van Garderen et al. (2024), IEEE TMI | **Official university record and abstract opened; article full text not opened.** Tests structural-MRI and DTI-informed diffusion-proliferation models for recurrent shape in 14 untreated post-resection LGG patients and evaluates spatial ranking with average precision. [DOI](https://doi.org/10.1109/TMI.2023.3298637) | Tissue-informed mechanistic shape prediction and AP evaluation already have real-patient precedent, though in a small LGG cohort rather than the proposed GBM annotated-transition target. |
| DeepGBM-Recure (2026) | **Abstract only.** Describes postoperative FLAIR segmentation, patient-level recurrence risk, and voxelwise recurrence hotspot mapping with supervoxels. [DOI](https://doi.org/10.1016/j.radonc.2025.111180) | Integrated patient/voxel recurrence mapping is not itself a gap; its full methods and comparability remain unchecked. |

This targeted screen weakens the novelty case. It found no source that exactly matches every component and this operational endpoint, but the method is a composition of known ingredients, and no search can prove that the exact composition is absent. Do not call the current prototype unique, first, revolutionary, or superior. The paper-scale contribution would have to come from a clearly defined operation plus a fair, reproducible result—not from adding more modules.

### Current working-copy code findings

These points are **static source inspection only**; no model, patient scoring, tests, or scripts were run for this recheck.

- Two blockers reported in the earlier read-only review are no longer present in the current working copy: `fit_count_prior` has a valid indented body, and `_load_raw_samples` restricts each future transition target to the forecast-origin candidate mask before constructing `RawSample`. The earlier report's syntax and target-shape findings are stale for these files.
- A more consequential design issue remains: the count model currently predicts only an expected count, and `calibrate_probability_mass` adds one shared logistic intercept so probabilities sum to that expectation. For any expected count strictly between zero and the number of candidate voxels, this is a strictly increasing transform of the original score. Therefore average precision and ROC AUC are mathematically unchanged by the count prediction. At boundary counts, all probabilities collapse to zero or one and ranking metrics become degenerate. The current inner-fold ridge choice uses average precision, so it does not select settings based on count quality; a claimed spatial-ranking win cannot be credited to the count component. Count must be assessed separately with count error/calibration, or the method must define and score a joint count-and-location predictive distribution.
- The current deterministic perturbation gate is still a sensitivity heuristic: one-voxel source erosion/dilation and shifts of the previous additions map do not provide calibrated uncertainty. It gates the learned residual, while the PDE score remains active.
- The “no-history” ablation zeros the transition-density feature but the count predictor and stability gate still use transition history. It is not a full ablation of history. The model-selection path also chooses the full hybrid's ridge before comparing several ablations; each arm needs its own training-only tuning for a fair comparison.
- The operational target remains newly annotated labels 1/3 within the fixed candidate region. It is not confirmed tumor-cell growth or clinical recurrence. Treatment effects, registration error, and annotation variability remain possible explanations for transitions.

### A candidate redesign to screen before more code

If the goal is a genuinely integrated mathematical method, the smallest substantive change to investigate is a **joint count-and-location law**, not another standalone module. Let K_i be the number of newly annotated voxels, and let Y_i(x) be a binary indicator that voxel x is a new annotated addition in the frozen candidate region X_i. Use only forecast-origin history H_i to predict a count distribution, for example

\[
K_i\mid H_i \sim \operatorname{NB}(\mu_i,\phi)\;\text{truncated to }\{0,\ldots,|X_i|\},\qquad \log\mu_i=\gamma_0+\gamma_1\log(1+K_{i,12})+\gamma_2\log\frac{\Delta t_{i,23}}{\Delta t_{i,12}}+\gamma_3\log\frac{|X_i|}{|X_{i,12}|}.
\]

Then condition spatial allocation on that count using the frozen PDE-plus-residual score \(\eta_i(x)=\beta u_i(x)+\theta^Tz_i(x)+\omega m_i(x)-\lambda V_i(x)\):

\[
P(Y_i=y\mid K_i=k,H_i)\propto \mathbf 1\!\left\{\sum_{x\in X_i}y_x=k\right\}\exp\!\left(\sum_{x\in X_i}y_x\eta_i(x)\right).
\]

This conditional-Bernoulli construction makes the predicted count and spatial ranking parts of one probability law instead of using an expected-count intercept that has no effect on ranking metrics. It is a **candidate formulation, not a claimed new theorem or verified novelty**; fixed-cardinality mask models and structured segmentation are established areas. Before coding, search those literatures directly and estimate the cost of its normalizer at the actual candidate-region sizes. The most likely failure is computational cost or overconfident independence among neighboring voxels; the weakest scientific point remains that the labels measure annotation transitions, not biological growth.

### Patient-analysis decision

Do not treat the current method as patient-ready or novelty-validated. The next gate is a small synthetic feasibility comparison of the joint count-location law against the current expected-count adjustment and a matched spatial-score baseline. It must verify calibration of the count distribution, spatial ranking, and joint mask score separately under known count overdispersion, spatial clustering, registration perturbations, and a misspecified PDE prior. Keep MU patient results for a later, explicitly exploratory analysis because the cohort has already been examined. Start that analysis only after the method and estimand are frozen, every ablation is tuned fairly within patient-grouped folds, and annotation/registration provenance has been reviewed. For confirmation, use a compatible untouched cohort; do not relabel a different recurrence endpoint as next-scan transition validation.

## Critical novelty check before implementation (2026-10-24)

The latest targeted search found stronger prior art than the earlier table captured. It changes the recommendation: do not present probabilistic future growth maps, hybrid size-plus-location prediction, or fixed-count allocation as the new method by themselves.

| Closest source | What the full text says it does | Consequence for this project |
|---|---|---|
| Petersen et al. (2019), *Deep Probabilistic Modeling of Glioma Growth* | **Full arXiv text opened.** Learns a distribution of plausible future multi-class tumor masks from the current and one prior scan. Data: 199 longitudinal MRIs from 38 patients (23 GBM); four MRI contrasts; 3-D evaluation. [Full text](https://arxiv.org/html/1907.04064) | Probabilistic future-mask prediction from longitudinal glioma scans predates this project. It is data-driven rather than a PDE/residual hybrid, and does not define the proposed next-visit new-label count forecast. |
| Laslo et al. (2025), *Mechanistic Learning with Guided Diffusion Models to Predict Spatio-Temporal Brain Tumor Growth* | **Full arXiv text opened.** A mechanistic ODE predicts tumor burden and conditions guided DDIM synthesis of future MRI; it produces tumor-growth probability maps. Trained with BraTS adult/pediatric data and evaluated on 60 axial slices from in-house pediatric DMG cases. [Full text](https://arxiv.org/html/2509.09610) | Mechanistic burden prediction coupled to probabilistic spatial tumor-growth maps is already a direct hybrid precedent. |
| Laslo et al. (2026), *Generative AI for spatial tumor growth on MRI* | **Full open-access article opened.** A guided DDIM uses baseline MRI and a specified target size; 100 generations are aggregated into a per-pixel growth-probability map. It trains with 1,251 adult and 144 pediatric BraTS23 patients and evaluates on 178 paired slices from 13 DMG patients, reporting mean slice-level continuous Dice 0.79. Its discussion says future work should remove the need to specify target size and add the baseline segmentation as input. [Article](https://doi.org/10.1186/s12916-026-04911-y) | This is very close to the spatial-growth-map endpoint. The opening for this project is narrower: forecast the amount without an oracle target size, use the prior segmentation explicitly, and evaluate 3-D postoperative GBM label additions. The 2025 ODE-plus-DDIM work already narrows the “no oracle size” distinction, so it is not sufficient alone to establish novelty. |
| Heng, Jacob & Ju (2020), conditional Bernoulli sampling | **Full arXiv text opened.** Defines independent Bernoulli indicators conditioned on their sum; discusses exact and swap-MCMC sampling costs. [Full text](https://arxiv.org/abs/2012.03103) | The proposed fixed-cardinality conditional-Bernoulli layer is established mathematics. It can be a baseline or implementation tool, not the claimed novelty. It also conditions on a fixed count and does not itself model positive spatial clustering. |

**How I searched:** targeted web searches used combinations of “glioma probabilistic future mask,” “longitudinal tumor growth probability map,” “mechanistic tumor burden guided diffusion,” “predicted tumor size spatial growth MRI,” and “conditional Bernoulli fixed cardinality.” I opened the primary full texts for the three closest tumor-forecasting papers and the conditional-Bernoulli paper, then followed their methods/discussion sections and linked prior-art references. This is a focused novelty screen, not a systematic review; the released code, all supplements, and full forward-citation coverage were not checked.

### Revised decision

The defensible candidate question is now: **Can a low-dimensional tissue-aware mechanistic prior plus a learned, patient-history correction produce a better calibrated joint forecast of the number and locations of next-scan newly annotated voxels than strong non-mechanistic and mechanistic baselines, without using the future scan or its size?** This is more specific than generic tumor-growth mapping, but it remains a candidate gap, not a verified novelty claim. Petersen et al. already forecast distributions of future masks; Laslo et al. already combine mechanistic amount prediction with probabilistic spatial growth maps. A full methods-level comparison of their equations and released code is still needed before claiming the exact formulation differs materially.

The conditional-Bernoulli formula above does not solve the spatial-dependence problem: given the count, it treats voxel scores as a weighted fixed-size subset and has no attractive neighbor interaction. Adding an arbitrary pairwise field would make normalization and inference difficult at voxel scale and would add complexity without current data support. Do not code that extension now. Keep the existing spatial ranker, the count predictor, and a count-location probabilistic model as separately testable components; the mathematical claim to test first is the joint predictive score and calibration, not a novel distribution family.

### Next gate and failure rule

1. Finish a paper-and-code comparison for Petersen 2019, Laslo 2025, and Laslo 2026: exact input scans, whether future size is known or forecast, target mask/transition, dimensionality, patient-level split, and metrics. Check supplements and released code where available.
2. Before any patient scoring, use synthetic 3-D masks with known count distributions and spatial clustering to compare: (a) the current expected-count intercept, (b) a predicted count distribution with conditional-Bernoulli allocation, and (c) the same spatial scores without count information. Score count calibration, conditional location ranking, and joint predictive log score separately. Include registration perturbation and model misspecification. This is a feasibility and falsification test, not evidence of clinical prediction.
3. If the fixed-count allocation adds no joint-score/calibration value, drop it. If the remaining architecture is only a known PDE/learned-residual combination, stop claiming a new method and narrow the project to an exploratory comparative result. If it passes, freeze the method and evaluate only with patient-grouped CV on the already examined MU data, labeled exploratory; external confirmation needs an untouched compatible cohort.

**Weakest point of this recommendation:** even a successful synthetic gate may show only that a standard count-location factorization behaves as designed. The available MU cohort cannot establish independent confirmation, and the target remains a scan-label transition rather than biological growth.
