# Technical Protocol for the Final GBM Research Plan

This is the methods companion to the form plan, not a second project. Tracks A, B, and C are equal, use separate data and endpoints, and do not form a patient-level digital twin. The protocol states the intended final analyses; it does not list abandoned trial-and-error branches as separate research aims. Earlier exploratory analyses exist and must be disclosed.

## 1. Shared design, data history, and decision rule

### 1.1 Separate evidence

Track A uses molecular and spatial atlas records. Track B uses longitudinal MRI masks. Track C uses simulated tumors only. No cross-source patient linkage is planned. The optional MU MGMT association is excluded from the three primary tracks because it does not validate a drug effect or connect the simulator to patient outcomes.

### 1.2 Previous exposure

The repository records prior analyses involving UCSC multiomic-gbm cells/donors, MU-Glioma-Post scans and splits, RHUH-GBM recurrence cases, LUMIERE cases, and MU MGMT analyses. The exact previously scored identifiers, split files, seeds, and output versions must be copied into a frozen exposure manifest before new evaluation. If exact identifiers cannot be recovered, the entire corresponding source cohort is considered previously examined. RHUH, LUMIERE, and the prior MGMT analysis are exposure disclosures only; they are not planned analyses in this three-track protocol.

All evaluation on these same sources is exploratory. New folds, model variants, or metrics do not turn reused cases into independent confirmation. No independent untouched real-data confirmation cohort is currently part of this plan. Track C uses synthetic scenarios; previously explored scenarios are development data and cannot serve as the locked evaluation set.

### 1.3 Comparison and “better” rule

For each track, freeze one primary score, all required ablations and baselines, eligibility gates, and the evaluation split before opening outer-fold predictions. Every comparison must use the same held-out units and available inputs.

For lower-is-better scores, define improvement as baseline score minus hybrid score. For Track C, define improvement as hybrid score minus comparator score. Calculate paired 95% cluster-bootstrap confidence intervals by resampling whole donors, patients, or virtual parameter scenarios, preserving all repeated cells, visits, noise levels, and replicates within each unit. Use a max-T bootstrap adjustment across all required comparator contrasts within a track. Do not use cells, voxels, or repeated noisy runs as independent bootstrap units.

A track supports a within-study advantage only if all eligibility and validation gates pass and every adjusted interval for the full hybrid versus a required ablation or baseline lies strictly on the favorable side of zero. Track A is descriptive unless at least ten eligible two-class donors remain; below ten, report donor-level estimates but make no inferential superiority claim. Otherwise report the point estimates, intervals, failures, and outcome as inconclusive or scenario-dependent. This rule does not establish a clinically meaningful minimum effect, statistical power for clinical translation, or external validation. The small donor and patient counts may make intervals wide.

### 1.4 Source and terms manifest

| Source | Planned role | Public access and scale | Terms and use gate |
|---|---|---|---|
| UCSC Cell Browser, multiomic-gbm | Track A primary molecular/spatial data | Public Cell Browser page and bulk download | The page supports downloads but does not display an explicit license for each export. Record exact files and terms; keep local, do not upload to Kaggle or redistribute until terms are established. |
| MU-Glioma-Post v1, TCIA | Track B primary MRI and masks | Official TCIA page lists 203 subjects and about 11 GB of public NIfTI/clinical files; the dataset paper reports 594 postoperative MRIs | TCIA lists CC BY 4.0 for its public NIfTI, segmentation-volume, and clinical packages. Reconcile the local visit manifest with the published count before freezing eligibility. |
| Track C virtual simulator | Track C only | Generated locally or in Kaggle; no patient data | Synthetic inputs are project assumptions. Kaggle is permitted for Track C simulations because no patient files are transferred. |

The UCSC source is publicly downloadable, but public access alone does not establish a reuse license. The current public repository already contains derived per-cell UCSC files according to the project audit. Do not add or publish additional source-derived files until their terms are clear; repository history or visibility is outside this protocol.

### 1.5 Computer and resource plan

The measured desktop for this project has 8 CPU cores, 31.6 GiB RAM, and no local GPU. Do not describe it as a 16-core, 34-GB system. Use sparse matrices and patient-level image crops; do not load all MRI volumes as one array. Before opening evaluation scores, profile a complete representative fold, record wall time and peak memory, and estimate the full run. Use Kaggle only for synthetic Track C runs or for a dataset whose exact terms explicitly permit cloud processing. If the locked analysis cannot finish on permitted compute, report the feasibility failure; any reduced design must be selected on development/training data and frozen before evaluation.

## 2. Track A: donor-robust molecular and spatial prediction

### 2.1 Question, target, and eligibility

The primary task is to classify source-labeled GBM Core versus GBM Periphery dissociated cells from expression while testing on donors absent from model fitting. The target is the source label, not a claim that the labels exhaust GBM biology. Keep NormalBrain and other diagnoses outside the primary task.

Create a source-backed alias table mapping every cell/sample identifier to a donor. Merge aliases only when source metadata supports the mapping. Keep unresolved identifiers separate and exclude them from the primary scored cohort. Record assay, library, cell type, label, donor, and specimen metadata. Do not use donor identity or outcome as a feature.

Run descriptive nested leave-one-donor-out evaluation only if at least six eligible donors have both target classes and every inner training split retains at least four such donors. If fewer than six qualify, report feasibility only. Do not make an inferential superiority claim unless at least ten eligible two-class donors remain; below ten, show the donor-level estimates and uncertainty as descriptive because few independent donors cannot support a stable comparative claim.

### 2.2 Full hybrid model

For cell \(i\), let \(x_i\in\mathbb{R}^G\) be the training-normalized expression vector, \(y_i\in\{0,1\}\) the Core/Periphery label, and \(d_i\) its donor. Learn a gene gate \(w=\operatorname{sigmoid}(v)\in[0,1]^G\), a two-class weight matrix \(W\), and intercept \(b\). The classifier is
\[
q_i=\operatorname{softmax}\{W(w\odot x_i)+b\}.
\]

Build a 15-nearest-neighbor graph among training cells only, using cosine distance on training-normalized genes. Symmetrize the graph and use unit weights on retained edges. Let \(E\) be the unique undirected edge set and \(Q\) the matrix of training probabilities. For each donor \(d\), let \(L_d\) be its mean cross-entropy. Use the smooth maximum
\[
L_{\mathrm{smooth}}=\tau\log\left(G_d^{-1}\sum_d \exp(L_d/\tau)\right),\qquad \tau=0.1,
\]
where \(G_d\) is the number of training donors. Minimize
\[
J=L_{\mathrm{smooth}}+\lambda_1\|w\|_1+\frac{\lambda_W}{2}\|W\|_F^2+
\frac{\lambda_2}{|E|}\sum_{(i,j)\in E}\|q_i-q_j\|_2^2.
\]
If there are no graph edges, set the graph term to zero. The pairwise term is zero for constant probabilities and its division by edge count prevents its scale from growing simply because the training set has more cells. The ridge term limits classifier weights so they cannot grow without bound to cancel a shrinking gene gate.

Use Adam at learning rate 0.001 for at most 500 epochs and early-stop after 25 epochs without inner-validation improvement. Search \(\lambda_1,\lambda_2\in\{10^{-4},10^{-3},10^{-2},10^{-1},1\}\) and \(\lambda_W\in\{10^{-3},10^{-2},10^{-1}\}\) in inner donor folds. For a query cell, compute its classifier probability and its 15 nearest training-cell probabilities. Blend them as
\[
p_i=(1-\gamma)q_i+\gamma\sum_{j\in N_{15}(i)}a_{ij}q_j,\quad
a_{ij}\propto \exp(\cos(x_i,x_j)/T).
\]
Choose \(T\in\{0.1,0.25,0.5,1\}\) and \(\gamma\in\{0,0.25,0.5,0.75,1\}\) using inner folds only. Do not form test-to-test graph edges.

### 2.3 Comparators and evaluation

Required component ablations are: (1) \(\lambda_1=0\) with ridge retained; (2) \(\lambda_2=0\) and \(\gamma=0\); (3) replace smooth worst-donor loss with mean donor cross-entropy. Required expression baselines are multinomial logistic regression, a 15-nearest-neighbor vote over training labels, a linear SVM, and a gradient-boosted tree model, all on the same outer donor folds. Add scANVI only if it can be trained on outer-training donors, use the same target labels, and pass a CPU runtime pilot. Tune all models only in inner donor folds with a fixed budget; compare all prespecified baselines with max-T adjustment. A cell-type-only model and a library/batch-metadata model are leakage diagnostics; they are not expression competitors. Add the existing attention classifier only if it can be refit inside folds with matched inputs and a fixed, comparable tuning budget. A donor-composition-only predictor is a leakage diagnostic, not a primary competitor.

For every outer held-out donor, rank cells by \(\max_k p_{ik}\), breaking ties with a fixed seed. At each coverage from 50% to 100%, retain the top fraction of all eligible cells. Compute error separately for retained Core and Periphery cells; if no cell of one class is retained at a coverage point, assign that class risk 1. Average the two class-specific risks, integrate over the coverage grid by the trapezoid rule, and average the resulting AURC equally over donors. Lower class-balanced AURC is better. Do not tune an abstention threshold for this endpoint. Report full-coverage balanced accuracy, calibration, and per-donor values secondarily. Use whole donors for uncertainty. If fewer than ten eligible two-class donors remain, these results are descriptive only and do not support an inferential superiority claim.

Check whether performance is driven by cell type, library, or batch. If source metadata support it, report strata with enough donors and compare a metadata-only diagnostic. The primary expression classifier must not use held-out donor ID, clinical outcomes, or target-derived pseudo-labels. A spatial-transfer extension is allowed only when measured spatial labels, specimen IDs, assay-compatible genes, and reuse terms are verified; hold out whole specimens and do not treat pseudozones as ground truth.

**Spatial label-independence gate:** Before spatial evaluation, document how each target label was assigned and which evidence was used, including histology, expression markers, deconvolution, manual review, and spatial propagation. Map every spot, cell, and section to a source specimen and donor. Do not use expression-derived or expression-propagated labels as independent targets for predictors that use the same expression measurements or overlapping marker features. Primary spatial prediction requires labels fixed from an independent measurement, such as source-authored histology reviewed without expression guidance, and a documented check for target-feature overlap. If independent labels or source-backed specimen-to-donor mapping are unavailable, do not run donor-held-out spatial scoring; describe source annotations only. Spatial results from this previously examined atlas remain exploratory.

### 2.4 Failure interpretation

The main risks are too few independent donors, unresolved aliases, donor/class confounding, and weak gains over logistic regression. If the eligibility gate fails, no resampling method can manufacture donor-level evidence. If the model only works in one library or cell type, report that dependence rather than a general donor-robust method. The result remains exploratory because the source donors have already been examined in the project.

## 3. Track B: uncertainty-aware reaction-diffusion forecasting

### 3.1 Cohort, prior exposure, and endpoint

Use the public MU-Glioma-Post v1 NIfTI, segmentation, and clinical files. TCIA lists the public packages under CC BY 4.0; the dataset paper describes 203 patients and 594 postoperative MRI time points. Record the exact local file manifest and reconcile its visit count with the source records before freezing eligibility.

Build an ordered patient-by-visit manifest with diagnosis, scan date and day origin, image-mask alignment, physical voxel spacing, source-label meanings, treatment context, and duplicate or missing visits. The primary tumor mask is the union of source labels 1 and 3 only after confirming their codebook meanings. Require at least three eligible, ordered, nonempty masks. Keep all visits from a patient together in five outer patient folds. Score the first two masks against the third and then use rolling-origin forecasts, always fitting from visits available before the target.

The cohort and its fixed splits have already been examined in this project. Every MU result is exploratory, even if a new model or score is used. Do not call it independent confirmation.

A close prior-art study is Jung and Shin (2026), an arXiv preprint that cites the MU dataset and reports the same source scale: 203 patients and 594 postoperative MRIs. It forms no-new-treatment triplets and tests a learned tumor-state proxy forecast against persistence and dilation on 15 patient-held-out triplets. Its results were comparable to persistence at one proxy threshold and numerically higher at another, with uncertainty that did not establish a difference from persistence. The study predicts a continuous label-derived proxy from MRI inputs; those published scores are not directly comparable with this plan’s binary mask target and spatial score. Reproduce or adapt its public OASA/SegMamba method only if the target, folds, and available inputs can be aligned fairly. Otherwise retain it as the closest same-dataset work and make no claim that the proposed method beats it.

### 3.2 Shape and size components

Model the latent concentration c(x,t) by
\[
\frac{\partial c}{\partial t}=\nabla\cdot(D\nabla c)+\rho c(1-c/K).
\]
Use K=1, initial concentration 1 inside the input mask and 0 outside, isotropic D in the primary fit, a 2 mm grid, a 30 mm padded domain clipped to the image field, no-flux crop boundaries, and visible contour c=0.5. Use nearest-neighbor mask resampling and sub-voxel trilinear interpolation for the contour. Use forward Euler with
\[
\Delta t\leq \min(0.5,\;0.9\Delta x^2/(6D))\ \mathrm{days},
\]
with D≤10 mm²/day. Do not expand parameter bounds after fitting. Mark any boundary-limited fit.

Estimate a joint empirical-Bayes prior over (log rho, log D) from outer-training patients only. Use a 21-by-21 logarithmic grid over rho in [10^-5, 10^-1] day^-1 and D in [10^-5, 10] mm²/day, with a local 11-by-11 refinement around high posterior mass and grid-cell quadrature weights. Retain correlated parameter pairs rather than substituting independent point estimates. Fit the shape component between two available masks by matching the simulated iso-surface to the second mask’s volume and evaluating signed-distance residuals at 256 fixed surface anchors.

For each rolling-origin forecast, estimate a separate log-volume trend from the last two available masks. Shrink the patient slope toward the outer-training population slope using a time constant chosen from {30, 90, 180, 365} days on training patients only. Model predictive log-volume variance as a linear function of forecast horizon, fitted from outer-training forecast residuals. This component can represent shrinking as well as growing volume.

For each predictive draw, sample a joint (rho,D) pair, use the PDE concentration to rank voxels, sample a positive future volume from the size model, and select the highest-ranked voxels until that sampled volume is reached. This separates the PDE’s spatial ranking from the independently forecast size.

### 3.3 Boundary uncertainty and synthetic gates

Represent effective boundary discrepancy at anchor points p_k by a zero-mean spatial Gaussian field with
\[
\Sigma_{kl}=\sigma_b^2\exp(-\|p_k-p_l\|/\ell)+\tau_b^2\mathbf{1}_{k=l}.
\]
Estimate sigma_b, tau_b, and correlation length ell in {5, 10, 20} mm using cross-fitted rolling-origin residuals from outer-training patients only. For each calibration patient, leave that patient out while estimating the growth prior, forecast a later mask from earlier masks, and exclude that target from all fitting. This is combined mask, registration, treatment, and model discrepancy; it is not a direct scanner-noise estimate.

For each draw, generate correlated displacement on 256 anchors, interpolate it to nearby surface vertices using the three nearest anchors with normalized weights proportional to (distance + 0.1 mm)^-2, move the contour along local normals, and rasterize the mask. Use 200 draws per forecast after a runtime pilot confirms feasibility. Apply the same training-derived error procedure to every probabilistic comparator.

Before patient scoring, generate at least 100 synthetic longitudinal series at observed visit gaps. Use boundary displacement SDs of 0, 2, and 4 mm and correlation lengths of 5, 10, and 20 mm. Include anisotropic growth and shrinking-volume cases. Test numerical stability, parameter recovery, predictive interval coverage, calibration, failure rate, runtime, and peak memory. For separate rho and D interpretation, require at least two fitting intervals, median absolute relative error at most 25% for each parameter, 80–97% coverage of nominal 90% intervals, and at most 10% of complete fits hitting any computational parameter bound. For forecast uncertainty, require 80–97% coverage of nominal 90% predictive intervals. These are project gates, not clinical standards. Failure blocks the affected parameter interpretation but does not remove a patient from forecast scoring.

Profile five synthetic and five training-patient fits on the 8-core desktop. If projected full-run wall time exceeds 48 hours, use the prespecified 3 mm grid, repeat all synthetic gates, and freeze the setting before opening outer-fold scores. Stop and report infeasibility if the 3 mm run still exceeds the limit.

### 3.4 Scores and matched baselines

The primary outcome is a spatial variogram score. Before opening a forecast target, select 256 voxels by a fixed deterministic rule within a 10 mm band around the last observed mask boundary. Form a fixed set of voxel pairs, stratified by physical separation into 0–2, 2–5, and 5–10 mm groups. For binary predictive draws X and held-out labels Y, score each pair by
\[
\left(\mathbb{E}[(X_i-X_j)^2]-(Y_i-Y_j)^2\right)^2.
\]
Give each nonempty distance stratum equal total weight; report empty strata and renormalize over the others. Average within forecast, then equally across forecast origins and patients. Lower is better. The held-out target supplies labels for scoring only; it cannot define model inputs or the scored band. Report energy score, marginal Brier score, volume CRPS, Dice, surface distance, volume error, predictive interval coverage, target voxels outside the band, and performance by observed shrinkage as secondary outcomes.

Use these matched baselines and ablations on identical patient folds and score locations:
1. A no-PDE ranker based on signed distance from the last mask, with the same sampled volume and boundary-error wrapper.
2. The PDE with its own iso-surface volume instead of the separate size forecast.
3. The PDE with independent boundary errors at the same marginal variance.
4. Persistence and volume-trend-only forecasts, each with the same volume and uncertainty treatment where applicable.
5. A 3D U-Net ranker trained on the same mask history and visit interval as the proposed method, followed by the same volume and boundary-error wrapper, if a one-fold compute pilot confirms feasibility.

Also attempt an OASA/SegMamba comparison because it is the closest published same-dataset forecaster. Its original use of MRI channels and continuous proxy labels differs from the mask-only inputs and binary target here. A direct superiority claim is permitted only if a retrained comparison uses the same target, patient folds, forecast origins, score, and a clearly declared input set. If this cannot be done within available licensed compute, report the comparison as infeasible and do not claim to beat the strongest published same-dataset method. Do not compare our new scores numerically with the preprint’s published scores as if they were a head-to-head test.

Keep all eligible cases, including observed shrinkage. Count numerical failures and use persistence with the same error wrapper as the fixed fallback for failed PDE forecasts. Report fallback frequency. The endpoint measures agreement with supplied masks, not the true tumor boundary or scanner noise.

### 3.5 Weakest point

A single supplied segmentation per scan does not identify segmentation noise. The spatial field may absorb real growth, treatment effects, registration error, and mask error. The cohort is post-treatment, some tumors shrink, and all data were previously examined. The PDE’s added spatial value may be zero, and the proposed uncertainty model may be too expensive or too weakly identified for these data.
## 4. Track C: robust virtual treatment control

### 4.1 Simulator and scenario design

Use the two-population sensitive/resistant model
\[
\dot S=\rho S\left(1-\frac{S+R}{K}\right)-\kappa uS,\qquad
\dot R=\rho(1-c)R\left(1-\frac{S+R}{K}\right)-\epsilon\kappa uR,
\]
with total burden B=S+R, carrying capacity K=1, and initial state S0=(1-fR,0)B0, R0=fR,0B0. Use daily simulation with 0.1-day integration steps. At each 28-day observation cycle, select either no course or a five-day course beginning at cycle start; u=1 on its five treatment days and 0 otherwise.

Sample the assumed parameters as follows: rho in [0.001, 0.05] day^-1, kappa in [0.02, 0.08] day^-1, initial resistant fraction fR,0 in [0.001, 0.1], resistance cost c in [0, 0.25], residual resistant-cell drug effect epsilon in [0, 0.2], and initial burden B0 in [0.05, 0.5]. Draw rho, kappa, and fR,0 log-uniformly; draw c, epsilon, and B0 uniformly. These are simulated assumptions, not patient or measured drug values.

Use a fixed 1,456-day horizon, with 52 decision cycles and a maximum of 12 five-day courses (60 active dose days) for every policy. Define simulated progression as the first day true burden reaches 1.2B0, capped at the horizon. Verify reachability without treatment at parameter-bound extremes before generating the locked scenarios. This endpoint belongs only to the simulator and is not a RANO or clinical response rule.

Create disjoint parameter-scenario sets: 400 for policy development and five-fold grouped cross-validation, 200 for filter-validation gates only, and 2,000 held out for the final in-simulator comparison. Evaluate each scenario in the three known observation-noise conditions sigma_y in {0.05, 0.10, 0.20}, with two noise replicates per condition. At cycle starts, observe y_j=log B(t_j)+eta_j, eta_j distributed N(0,sigma_y^2). Use sigma_y=0 only for exact simulator checks. Reuse the same scenario IDs and observation streams for every policy; keep all noise conditions and replicates of a parameter scenario together in the grouped folds. No held-out scenario may tune a threshold, horizon, policy, seed, or checkpoint.

Create four separate 1,000-scenario misspecification suites at sigma_y=0.10, with no tuning on them:
1. Autocorrelated observations: eta_j=0.5 eta_(j-1)+zeta_j, with innovation variance chosen to keep marginal SD at 0.10.
2. Heavy-tailed observation errors: standardized Student-t errors with 4 degrees of freedom and marginal SD 0.10.
3. Changing tumor growth: multiply both population growth rates by 0.5 or 1.5 once at a uniformly sampled cycle between 13 and 39.
4. Waning treatment effect: after each completed course, multiply the next course's kappa by exp(-0.10 times the number of prior courses).
Use common scenario parameters and policy-independent random streams within each suite. Do not adjust any policy after observing suite results.

### 4.2 Particle filter and validation

Represent the hidden state at each cycle by total burden B and resistant fraction fR=R/B, with static parameter particles for rho, kappa, fR,0, c, epsilon, and B0. Use 500 particles and the Gaussian observation likelihood at the known sigma_y for the in-distribution scenarios. Systematically resample when effective sample size falls below 250. Apply Liu–West regularization with shrinkage a=0.98 to standardized static-parameter particles, using kernel covariance (1-a^2) times the weighted covariance and reflecting proposals at prior bounds. After a parameter move, preserve the current B and fR state and reconstruct S=(1-fR)B and R=fRB; parameter jitter must not add or remove current burden.

On the separate 200-scenario filter-validation set, check nominal 90% intervals for burden and resistant fraction; require empirical coverage from 80% to 97%. Across five filter seeds, require posterior probabilities P(fR>=b) to differ by no more than 0.10 at 90% or more of decision times. These gates validate the filter only under the stated simulator. If they fail, retain the policy results but make no reliable hidden-resistance estimation claim.

### 4.3 Hybrid policy and training

At cycle zero, set the reference burden to the filter's posterior mean of B0. At each later cycle, compute the filtered burden ratio r_j=E[B(t_j)|observations]/E[B0|observations at cycle zero]. Tune a_on and a_off from {1.05, 1.10, 1.20, 1.30}, resistance threshold b from {0.05, 0.10, 0.20, 0.35, 0.50}, and resistance-mode threshold a_R from the same burden-threshold grid. Require a_off < min(a_on,a_R). When treatment is off, begin a five-day course if r_j reaches a_on, or if P(fR(t_j)>=b)>=0.80 use a_R in place of a_on and begin the course if r_j reaches a_R. When treatment is on, give the next five-day course only while r_j remains at or above a_off. Stop dosing after 12 courses or progression.

Use five-fold cross-validation grouped by parameter-scenario ID within the 400 development scenarios. Fit each fold's learned component on four folds and evaluate threshold settings on the held-out fold. Average the two noise replicates within each scenario and noise level. Select policy thresholds by maximizing the minimum lower-decile mean time to progression across the three noise levels in pooled out-of-fold predictions. Train and select PPO seeds/checkpoints using training-only out-of-fold predictions, with no evaluation-scenario feedback. The separate 200-scenario filter-validation set tests calibration gates and does not select policy thresholds.

### 4.4 Comparisons and outcome

Compare the full hybrid with: (1) filtered burden-only hysteresis; (2) raw-observation burden-only hysteresis; (3) resistance-only thresholding; (4) a point-estimate filter ablation; (5) threshold selection for mean time rather than lower-decile time; (6) a fixed front-loaded schedule; (7) an evenly paced schedule; (8) the best periodic schedule selected on development scenarios; and (9) PPO trained for the same lower-tail objective with the same observation, filter-summary, treatment-history, remaining-budget, and elapsed-cycle inputs. Use five PPO seeds and select each seed/checkpoint using training-only folds. A schedule optimized with access to true hidden parameters is an unattainable upper-bound illustration only, never a primary comparator.

For each held-out scenario and noise level, average time to progression over its two replicates, capped at 1,456 days. The primary outcome is the minimum across the three noise levels of the lower-decile mean (CVaR at 10%) of these scenario-level times. Higher is better. Use paired scenario bootstrap intervals with max-T adjustment across the prespecified comparisons. Report each noise-level CVaR, capped fraction, median time, burden integral, final resistant fraction, and course count secondarily. Report each misspecification suite separately. Call the policy robust across the four tested departures only if its adjusted paired interval favors the strongest development-selected baseline in all four suites; otherwise name which departures change the ranking. Define strata fR,0<0.01 versus >=0.01 and c<0.125 versus >=0.125 before opening held-out outcomes. Report a ranking reversal only when paired intervals favor different methods in opposite directions in at least two strata.

### 4.5 Weakest point

The treatment model and parameter ranges are assumed, not measured from patients or experiments. A filter may appear calibrated because its assumed likelihood and simulator agree, and all policies may inherit the same structural error. More simulations reduce Monte Carlo error but cannot show that a real drug works or that a policy would help a patient. Track C is an algorithm comparison inside one stated model family.
## 5. Prior-art boundary

Prior work already includes GBM Core/Periphery cell classification, donor-held-out single-cell annotation, graph and sparse-feature methods, reaction-diffusion tumor models, probabilistic mask forecasting, spatial uncertainty models, and model-informed adaptive therapy. In particular, Lopes and Vinga (2020) classify GBM Core versus Periphery cells with sparse logistic models; scTab and scANVI provide strong annotation references; Jung and Shin’s 2026 preprint forecasts a label-derived tumor-state proxy on the same MU collection; and Zhang et al. (2017) and Liu et al. (2025) establish adaptive evolutionary treatment methods. The three proposed combinations remain hypotheses for testing, not established novelties. A limited literature scan cannot support a “first” or “beats all current methods” claim. Before such a claim, read the closest full papers, inspect available code, and run fair, input-matched comparisons on the same units and endpoints.

## 6. References and source records

### Data and source records

- UCSC Cell Browser, multiomic-gbm: https://cells.ucsc.edu/?ds=multiomic-gbm
- MU-Glioma-Post v1, TCIA: https://www.cancerimagingarchive.net/collection/mu-glioma-post/

### Literature

Sonpatki, P., Park, H. J., Xing, Y. L., Han, K. Y., Schroeder, B. A., Yu, H. J., Kim, H. J., Chowdhury, T., Hwang, J. H., Nam, S. M., Byun, Y. H., Kang, H., Lee, J. H., Lee, S. T., Won, J. K., Kim, T. M., Choi, S. H., Ku, J.-L., Lee, S., … Shah, N. (2026). A spatially resolved human glioblastoma atlas reveals distinct cellular and molecular patterns of anatomical niches. Nature Communications, 17, 2951. https://doi.org/10.1038/s41467-026-69716-2

Jung, Y., & Shin, M. (2026). Observation-Anchored Selective Assimilation for Longitudinal Tumor-State Proxy Forecasting in Post-Treatment Glioma. arXiv:2609.12435 [preprint]. https://doi.org/10.48550/arXiv.2609.12435

Lopes, M. B., & Vinga, S. (2020). Tracking intratumoral heterogeneity in glioblastoma via regularized classification of single-cell RNA-Seq data. BMC Bioinformatics, 21, 59. https://doi.org/10.1186/s12859-020-3390-4

Fischer, F., Fischer, D. S., Mukhin, R., Isaev, A., Biederstedt, E., Villani, A.-C., & Theis, F. J. (2024). scTab: Scaling cross-tissue single-cell annotation models. Nature Communications, 15, 6611. https://doi.org/10.1038/s41467-024-51059-5

Hormuth, D. A., II, Al Feghali, K. A., Elliott, A. M., Yankeelov, T. E., & Chung, C. (2021). Image-based personalization of computational models for predicting response of high-grade glioma to chemoradiation. Scientific Reports, 11, 8520. https://doi.org/10.1038/s41598-021-87887-4

Liu, Z., Zhang, J., Hong, L., Nie, Q., & Sun, X. (2025). Multiscale mathematical model-informed reinforcement learning optimizes combination treatment scheduling in glioblastoma evolution. Science Advances, 11(32), eadv3316. https://doi.org/10.1126/sciadv.adv3316

Mahmoud, E., Gass, J., Dhemesh, Y., Greaser, J., Pogorzelski, K., Isufi, E., Garrett, F., Thacker, J., Tahon, N. H. M., Sinclair, J., Layfield, L. J., Altes, T., & Nada, A. (2025). MU-Glioma Post: A comprehensive dataset of automated MR multi-sequence segmentation and clinical features. Scientific Data, 12, 1847. https://doi.org/10.1038/s41597-025-06011-7

Monteiro, M., Le Folgoc, L., Coelho de Castro, D., Pawlowski, N., Marques, B., Kamnitsas, K., van der Wilk, M., & Glocker, B. (2020). Stochastic segmentation networks: Modelling spatially correlated aleatoric uncertainty. Advances in Neural Information Processing Systems, 33, 12756–12767. https://arxiv.org/abs/2006.06015

Xu, C., Lopez, R., Mehlman, E., Regier, J., Jordan, M. I., & Yosef, N. (2021). Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models. Molecular Systems Biology, 17(1), e9620. https://doi.org/10.15252/msb.20209620

Subramanian, S., Scheufele, K., Mehl, M., & Biros, G. (2020). Where did the tumor start? An inverse solver with sparse localization for tumor growth models. Inverse Problems, 36(4), 045006. https://doi.org/10.1088/1361-6420/ab649c

West, J., Adler, F., Gallaher, J., Strobl, M., Brady-Nicholls, R., Brown, J., Roberson-Tessi, M., Kim, E., Noble, R., Viossat, Y., Basanta, D., & Anderson, A. R. A. (2023). A survey of open questions in adaptive therapy: Bridging mathematics and clinical translation. eLife, 12, e84263. https://doi.org/10.7554/eLife.84263

Zhang, J., Cunningham, J. J., Brown, J. S., & Gatenby, R. A. (2017). Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer. Nature Communications, 8, 1816. https://doi.org/10.1038/s41467-017-01968-5

Zhong, Z., Hou, J., Yao, Z., Dong, L., Liu, F., Yue, J., Wu, T., Zheng, J., Ouyang, G., & Song, J. (2024). Domain generalization enables general cancer cell annotation in single-cell and spatial transcriptomics. Nature Communications, 15, 1929. https://doi.org/10.1038/s41467-024-46413-6
