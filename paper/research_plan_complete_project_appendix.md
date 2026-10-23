# Technical Protocol for the Final GBM Research Plan

This is the methods companion to the form plan, not a second project. Tracks A, B, and C are equal, use separate data and endpoints, and do not form a patient-level digital twin. The protocol states the intended final analyses; it does not list abandoned trial-and-error branches as separate research aims. Earlier exploratory analyses exist and must be disclosed.

## 1. Shared design, data history, and decision rule

### 1.1 Separate evidence

Track A uses molecular and spatial atlas records. Track B uses longitudinal MRI masks. Track C uses simulated tumors only. No cross-source patient linkage is planned. The optional MU MGMT association is excluded from the three primary tracks because it does not validate a drug effect or connect the simulator to patient outcomes.

### 1.2 Previous exposure

The repository records prior analyses involving UCSC multiomic-gbm cells/donors, MU-Glioma-Post scans and splits, RHUH-GBM recurrence cases, LUMIERE cases, and MU MGMT analyses. The exact previously scored identifiers, split files, seeds, and output versions must be copied into a frozen exposure manifest before new evaluation. If exact identifiers cannot be recovered, the entire corresponding source cohort is considered previously examined.

All evaluation on these same sources is exploratory. New folds, model variants, or metrics do not turn reused cases into independent confirmation. No independent untouched real-data confirmation cohort is currently part of this plan. Track C uses synthetic scenarios; previously explored scenarios are development data and cannot serve as the locked evaluation set.

### 1.3 Comparison and “better” rule

For each track, freeze one primary score, all required ablations and baselines, eligibility gates, and the evaluation split before opening outer-fold predictions. Every comparison must use the same held-out units and available inputs.

For lower-is-better scores, define improvement as baseline score minus hybrid score. For Track C, define improvement as hybrid score minus comparator score. Calculate paired 95% cluster-bootstrap confidence intervals by resampling whole donors, patients, or virtual parameter scenarios, preserving all repeated cells, visits, noise levels, and replicates within each unit. Use a max-T bootstrap adjustment across all required comparator contrasts within a track. Do not use cells, voxels, or repeated noisy runs as independent bootstrap units.

A track supports a within-study advantage only if all eligibility and validation gates pass and every adjusted interval for the full hybrid versus a required ablation or baseline lies strictly on the favorable side of zero. Otherwise report the point estimates, intervals, failures, and outcome as inconclusive or scenario-dependent. This rule does not establish a clinically meaningful minimum effect, statistical power for clinical translation, or external validation. The small donor and patient counts may make intervals wide.

### 1.4 Source and terms manifest

| Source | Planned role | Public access and scale | Terms and use gate |
|---|---|---|---|
| UCSC Cell Browser, multiomic-gbm | Track A primary molecular/spatial data | Public Cell Browser page and bulk download | The page supports downloads but does not display an explicit license for each export. Record exact files and terms; keep local, do not upload to Kaggle or redistribute until terms are established. |
| MU-Glioma-Post v1, TCIA | Track B primary MRI and masks | Official TCIA page lists 203 subjects, 596 follow-up time points, about 11 GB | TCIA page lists CC BY 4.0 for its public NIfTI/clinical files. The article reports a different time-point count; use the downloaded manifest as the analysis count and cite the discrepancy. |
| RHUH-GBM, TCIA | Track B secondary retrospective recurrence localization | Official TCIA page lists 40 subjects and a 2.9 GB public brain-extracted NIfTI/segmentation/clinical package | The public NIfTI package is listed CC BY 4.0. Exclude DICOM, which has a separate controlled-access policy. |
| LUMIERE, Figshare | No core analysis; possible sensitivity source only | Public 30.33 GB package; 91 patients, 638 study dates, 2,487 MRI images | The Figshare item says CC0, while the associated README and article usage notes say non-commercial use. The terms conflict. Exclude unless the data provider clarifies which terms govern the files. |
| Track C virtual simulator | Track C only | Generated locally or in Kaggle; no patient data | Synthetic inputs are project assumptions. Kaggle is permitted for Track C simulations because no patient files are transferred. |

The UCSC source is publicly downloadable, but public access alone does not establish a reuse license. The current public repository already contains derived per-cell UCSC files according to the project audit. Do not add or publish additional source-derived files until their terms are clear; repository history or visibility is outside this protocol.

### 1.5 Computer and resource plan

The measured desktop for this project has 8 CPU cores, 31.6 GiB RAM, and no local GPU. Do not describe it as a 16-core, 34-GB system. Use sparse matrices and patient-level image crops; do not load all MRI volumes as one array. Before opening evaluation scores, profile a complete representative fold, record wall time and peak memory, and estimate the full run. Use Kaggle only for synthetic Track C runs or for a dataset whose exact terms explicitly permit cloud processing. If the locked analysis cannot finish on permitted compute, report the feasibility failure; any reduced design must be selected on development/training data and frozen before evaluation.

## 2. Track A: donor-robust molecular and spatial prediction

### 2.1 Question, target, and eligibility

The primary task is to classify source-labeled GBM Core versus GBM Periphery dissociated cells from expression while testing on donors absent from model fitting. The target is the source label, not a claim that the labels exhaust GBM biology. Keep NormalBrain and other diagnoses outside the primary task.

Create a source-backed alias table mapping every cell/sample identifier to a donor. Merge aliases only when source metadata supports the mapping. Keep unresolved identifiers separate and exclude them from the primary scored cohort. Record assay, library, cell type, label, donor, and specimen metadata. Do not use donor identity or outcome as a feature.

Run nested leave-one-donor-out evaluation only if at least eight eligible donors have both target classes and every inner training split retains at least six such donors. If the gate fails, report a feasibility/descriptive result and do not claim comparative superiority. A stronger “supported within-study advantage” statement also requires at least ten eligible donors; below ten, present cluster intervals as descriptive because resampling so few donors cannot give stable precision.

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

Required component ablations are: (1) \(\lambda_1=0\) with ridge retained; (2) \(\lambda_2=0\) and \(\gamma=0\); (3) replace smooth worst-donor loss with cell-weighted empirical cross-entropy. Required baseline: multinomial logistic regression on the same expression features and same donor folds. Add the existing attention classifier only if it uses matched inputs and a fixed, comparable tuning budget. A donor-composition-only predictor and metadata-only models are leakage diagnostics, not primary competitors.

For every outer held-out donor, rank cells by \(\max_k p_{ik}\). Calculate retained-cell error over coverage from 50% to 100%, interpolate linearly, integrate by the trapezoid rule, then average the AURC equally over donors. Lower AURC is better. Do not tune an abstention threshold for this endpoint. Report macro-F1, balanced accuracy, calibration, and per-donor values secondarily. Use whole donors for uncertainty.

Check whether performance is driven by cell type, library, or batch. If source metadata support it, report strata with enough donors and compare a metadata-only diagnostic. The primary expression classifier must not use held-out donor ID, clinical outcomes, or target-derived pseudo-labels. A spatial-transfer extension is allowed only when measured spatial labels, specimen IDs, assay-compatible genes, and reuse terms are verified; hold out whole specimens and do not treat pseudozones as ground truth.

### 2.4 Failure interpretation

The main risks are too few independent donors, unresolved aliases, donor/class confounding, and weak gains over logistic regression. If the eligibility gate fails, no resampling method can manufacture donor-level evidence. If the model only works in one library or cell type, report that dependence rather than a general donor-robust method. The result remains exploratory because the source donors have already been examined in the project.

## 3. Track B: uncertainty-aware reaction-diffusion forecasting

### 3.1 Cohort, timing, and endpoint

Use the public MU-Glioma-Post v1 NIfTI and clinical files. Confirm GBM diagnosis for a GBM-specific result; retain the exact diagnosis and treatment context per patient. Build an ordered visit manifest with acquisition day, image-mask alignment, voxel geometry, source mask labels, duplicate/missing visits, and treatment status. Define tumor volume as voxel count times physical voxel volume under a fixed source-label rule.

The first forecast requires at least three ordered masks. Fit using the first two and predict the third, then add each later observed mask only after its prior forecast in rolling-origin evaluation. Keep all visits from each patient together in an outer fold. With one observed interval, fit and predict using the training-derived joint prior but do not interpret \(\rho\) and \(D\) separately. Separate parameter interpretation is a secondary gate requiring at least two fitted intervals and successful synthetic recovery.

Use at least 20 eligible patients for a comparative superiority claim. Below 20, report feasibility and patient-level estimates descriptively. This is a project precision gate, not a clinical standard.

### 3.2 Latent tumor PDE and numerical solution

Model latent concentration \(c(x,t)\) by
\[
\frac{\partial c}{\partial t}=\nabla\cdot(D\nabla c)+\rho c(1-c/K).
\]
Use \(K=1\), an initial condition of 1 inside the first observed mask and 0 outside, an isotropic scalar \(D\) in the primary analysis, a 2 mm physical grid, a 30 mm padded baseline bounding box clipped to the image field, no-flux crop boundaries, and visible contour \(c=0.5\). Resample masks by nearest neighbor. Solve the 3-D finite-difference equation using forward Euler with \(\Delta t=0.05\) day. Hold \(D\le10\;\mathrm{mm}^2/\mathrm{day}\), for which the stated six-neighbor diffusion condition \(\Delta t\le \Delta x^2/(6D)\) is satisfied at \(\Delta x=2\) mm. Do not expand the bound after fitting. If a fit reaches a bound, mark its parameter estimate boundary-limited.

### 3.3 Joint parameter uncertainty and boundary discrepancy

Estimate a bivariate training-only empirical-Bayes distribution for \((\log\rho,\log D)\) from all available ordered scans of outer-training patients, using rolling-origin fits within those patients only. Never use held-out-patient scans to estimate this prior. Use a 21-by-21 logarithmic grid over \(\rho\in[10^{-5},10^{-1}]\;\mathrm{day}^{-1}\) and \(D\in[10^{-5},10]\;\mathrm{mm}^2/\mathrm{day}\), with local 11-by-11 refinement around high posterior mass. Normalize likelihood-times-prior weights with grid-cell quadrature weights. Retain correlated parameter combinations across the \(\rho,D\) ridge rather than substituting separate point estimates.

On each training follow-up mask, sample 256 approximately even surface anchors. Let \(r_k\) be the signed distance from observed anchor \(p_k\) to the predicted \(c=0.5\) surface. Model the effective mask discrepancy as a zero-mean spatial Gaussian field with covariance
\[
\Sigma_{kl}=\sigma_b^2\exp(-\|p_k-p_l\|/\ell)+\tau_b^2\mathbf{1}_{k=l}.
\]
Estimate \(\sigma_b,\tau_b\) and select \(\ell\in\{2,5,10\}\) mm from cross-fitted, rolling-origin forecast residuals from outer-training patients only. For each training patient used to generate calibration residuals, estimate the growth prior without that patient, forecast a later target mask from its preceding masks, and exclude that target mask from all fitting. Pool those out-of-fold signed-distance residuals to estimate the covariance settings for the outer fold. Never estimate error scale or correlation from in-sample fitting residuals or from outer-test masks. This models combined mask, registration, treatment, and model discrepancy; it is not a direct scanner-noise measurement.

For a held-out forecast, draw parameter pairs from the joint weighted grid. Independently draw a spatial error field using the same covariance process on 256 fixed anchors of that forecast surface, then interpolate to its other surface vertices from the three nearest anchors using weights proportional to \((d+0.1\;\mathrm{mm})^{-2}\). Move the predicted contour along local normals and rasterize the binary mask. Use 200 predictive masks per forecast after a runtime pilot confirms feasibility. This is a grid-based predictive mixture, not a Hessian/Laplace approximation.

### 3.4 Synthetic gates, score, and fair baselines

Before patient scoring, create at least 100 synthetic longitudinal series at the observed visit intervals. Generate boundary shifts from spatial Gaussian fields with SD 0, 2, and 4 mm and correlation length 2, 5, and 10 mm. Add anisotropic and misspecified-growth cases. Where at least two fitting intervals are available, check recovery of \(\rho,D\). Separate parameter interpretation passes only when the median absolute relative error for each parameter is at most 25%, its nominal 90% parameter-interval coverage is 80–97%, and no more than 10% of complete synthetic fits hit any computational parameter bound. For forecasts, check nominal 90% predictive interval coverage and width separately; the predictive-coverage gate is 80–97%. These are project gates, not clinical standards. Count each complete fit once in the bound-hit denominator. A failed gate blocks the affected parameter interpretation but does not justify excluding that patient from forecast scoring.

Use a spatial variogram score as the primary metric. On a fixed deterministic sample of 256 voxels inside a 10 mm band around each supplied held-out mask boundary, make a fixed set of voxel pairs stratified by physical separation (0–2, 2–5, and 5–10 mm). For each pair \((i,j)\), let \(X_i,X_j\) be binary labels in predictive masks and \(Y_i,Y_j\) the supplied held-out labels. Score
\[
VS=\sum_{(i,j)}w_{ij}\left(\mathbb{E}[(X_i-X_j)^2]-(Y_i-Y_j)^2\right)^2,
\]
with equal total weight for each nonempty distance stratum; if a stratum has no pairs, report it and renormalize over the remaining strata. Average pair contributions within each forecast, forecasts equally within each patient, then patients equally. Lower is better. The held-out mask may define scoring locations and observed pair outcomes, but never model inputs. Report marginal Brier score as a secondary metric; it alone cannot test spatial dependence.

Compare the full method against: (1) the PDE with independent boundary error and matched marginal error variance; (2) the PDE mean forecast with the training-calibrated correlated-error wrapper; (3) persistence, using the last mask; and (4) time-gap-aware volume trend. Make the volume-trend mask by linearly extrapolating volume from the prior two masks and shifting the prior signed-distance surface by bisection until predicted volume matches that extrapolation, permitting shrinkage. Give persistence and volume trend the same training-derived correlated-error wrapper as the PDE baselines, so every method receives probabilistic samples and the spatial score is comparable. Use identical patient folds and score locations.

Keep all patients with eligible scans, including shrinking-volume cases. Growth-only Fisher-Kolmogorov dynamics may fail on shrinkage. Count the failures and use persistence plus the same error wrapper as a prespecified fallback for a failed PDE forecast; report fallback frequency and include those cases in the primary score. Report Dice, surface distance, volume error, Brier score, calibration, interval coverage, and performance by shrinkage stratum secondarily. The primary endpoint measures agreement with supplied masks, not true tumor boundaries or scanner noise.

### 3.5 Secondary recurrence-localization analysis

Only after the primary forecast is frozen, use public RHUH-GBM brain-extracted NIfTI/segmentation files for a separate retrospective spatial-localization analysis, conditional on the observed early-postoperative-to-recurrence interval. Do not call it a prediction of recurrence time. Require at least 20 verified early-postoperative/recurrent mask pairs for comparison; otherwise show descriptive results only.

For each patient, let \(C\) be the early-postoperative cavity/tumor mask, \(P(x)\) the PDE concentration propagated over the observed interval, and \(D(x)=\exp[-d(x,C)/(15\;\mathrm{mm})]\), where \(d\) is physical distance. Form \(R_w(x)=wP(x)+(1-w)D(x)\), \(w\in\{0,0.1,\ldots,1\}\). Within training patients, select \(w\) to maximize recurrence sensitivity when the predicted high-risk region has the same volume as a 15 mm physical expansion of \(C\); break ties by Dice. For a held-out patient, do not use their recurrence mask to fit or select \(w\). Compare the resulting region against that patient’s recurrence mask and against the 15 mm expansion baseline using sensitivity and Dice. This is retrospective localization conditional on a known follow-up interval, not a radiation margin recommendation. All RHUH results remain exploratory because these cases were previously examined.

### 3.6 Weakest point

A single supplied segmentation per scan does not identify segmentation noise. The spatial field can absorb actual growth, therapy effects, registration error, and labeler/model error. The MU data are mostly post-treatment and patients can shrink, while the primary PDE is growth-only. Many fits may therefore be prior-dominated or use the persistence fallback. All results are exploratory and mask-agreement results.

## 4. Track C: robust virtual treatment control

### 4.1 Simulator and distributions

Use the well-mixed two-population model
\[
\dot S=\rho S\left(1-\frac{S+R}{K}\right)-\kappa uS,\qquad
\dot R=\rho(1-c)R\left(1-\frac{S+R}{K}\right)-\epsilon\kappa uR,
\]
with \(B=S+R\), \(K=1\), \(u\in\{0,1\}\), daily decision times, and a 0.1-day numerical step. Rates \(\rho,\kappa\) are in day\(^{-1}\). Use the existing virtual ranges \(\rho\in[0.005,0.05]\), \(\kappa\in[0.02,0.08]\), initial resistant fraction \(f_R\in[0.001,0.1]\), resistance cost \(c\in[0,0.25]\), residual drug effect \(\epsilon\in[0,0.2]\), and initial \(B/K\in[0.08,0.65]\). Draw rates and \(f_R\) log-uniformly; draw \(c,\epsilon\), and \(B/K\) uniformly. These are assumed ranges, not patient estimates. The upper initial burden is capped at 0.65 so a progression threshold at 1.4 times the running nadir remains below \(K\).

Generate 256 Latin-hypercube training scenarios and 512 disjoint evaluation scenarios. Use five-fold cross-validation grouped by parameter-scenario ID within the 256 training scenarios. In each fold, fit and tune learned components on the other four folds, then score every candidate horizon, trigger threshold, PPO seed, and checkpoint on the held-out fold. Keep all noise levels and replicates of a scenario together. Combine held-out predictions across folds to select these settings, then refit the selected learned policy on all 256 training scenarios. Never use the 512 evaluation scenarios for tuning. For every parameter scenario, simulate at all three log-observation noise SDs \(\sigma_y\in\{0.05,0.10,0.20\}\) with five observation-noise replicates each. Let \(y_t=\log B_t+\varepsilon_t\), \(\varepsilon_t\sim N(0,\sigma_y^2)\). Use \(\sigma_y=0\) only to verify exact simulator behavior. Reuse the same parameter scenario IDs and noise streams across every policy. The test scenarios stay unopened until model, thresholds, horizon, and policy-training choices are frozen.

### 4.2 Particle filter and hybrid rule

Represent hidden state by burden \(B\) and resistant fraction \(f_R=R/B\), plus static parameter particles. Use 500 particles and the Gaussian likelihood of observed \(y_t\). When effective sample size falls below 250, systematically resample. Apply Liu–West regularization with \(a=0.98\) to standardized static-parameter particles, kernel covariance \((1-a^2)\) times the weighted covariance, and reflect proposed parameters at prior bounds. After jitter, preserve current \(B\) and \(f_R\), then set \(S=(1-f_R)B\) and \(R=f_RB\); parameter jitter must not create or remove tumor burden.

Validate 90% particle intervals against known simulated states on training-only validation scenarios, with an 80–97% empirical coverage gate. Across five filter seeds, require the range of posterior tail probabilities at a decision time to be at most 0.10 at 90% or more of validation decision times. These are calibration gates for the assumed simulator. Report them even if they fail.

Let the controller observe normalized burden relative to its running observed nadir,
\[
r_t=\exp\{y_t-\min_{s\le t}y_s\}.
\]
The two-trigger rule gives a dose if \(r_t\ge a\) for \(a\in\{1.05,1.10,1.20,1.30\}\), or if posterior probability that \(f_R\ge b\) is at least 0.80 for \(b\in\{0.05,0.10,0.20,0.35,0.50\}\). The trigger uses only noisy observations and filter summaries. True simulated progression is recorded when \(B_t\ge1.4\min_{s\le t}B_s\); it is not a RANO rule. The resistance trigger can fire even when a burden trigger has not.

### 4.3 Dose budget, horizon, and training

Give each policy the same maximum active-dose-day budget \(B_{\max}=5\lfloor H/28\rfloor+\min(5,H\bmod 28)\), the exact active-day count when a repeated 5-on/23-off schedule starts on an on-day at day 1 of integer horizon \(H\). The fixed comparator uses that schedule. Other policies may place the same number of active days differently; after the budget is spent, treatment is off.

Select \(H\) from 500, 1,000, or 1,500 days using out-of-fold predictions from the training-scenario cross-validation: choose the shortest horizon for which fewer than 90% remain progression-free under the fixed schedule. If none qualifies, use 1,500 days and report the ceiling risk. No horizon may be selected using evaluation outcomes.

Use the five scenario-grouped folds specified in §4.1 for all horizon, threshold, seed, and checkpoint selection. For each candidate policy and threshold setting, obtain predictions only from the fold in which each training scenario was held out; no scenario may be scored by a policy fitted or tuned on that scenario. Average the five noise replicates within each scenario and noise level, then select trigger thresholds by maximizing the minimum lower-decile mean time to progression across the three noise levels using the pooled out-of-fold scenario-level averages. Select horizon and learned-policy seed/checkpoint by the same out-of-fold procedure. Train PPO to maximize the same minimum-across-noise lower-decile objective under the same dose budget. The PPO comparator must receive the same noisy burden, observation and treatment history, posterior mean/SD and resistance-tail probabilities, remaining dose budget, and elapsed-horizon fraction as the hybrid. If risk-sensitive PPO fails to train, report that failure and retain it as a failed comparator; do not substitute a policy selected on evaluation data.

### 4.4 Comparisons and outcome

Required comparisons are: full hybrid; burden-only trigger; resistance-only trigger; point-estimate filter ablation; mean-time-optimized threshold policy; fixed 5/23 schedule; and risk-sensitive PPO. The point-estimate ablation uses an extended Kalman filter and triggers on its point estimate of resistant fraction rather than a posterior probability. Every method uses the same simulator, budget, horizon, evaluation scenarios, and noise streams; evaluate learned policies with deterministic actions. Any policy with access to true hidden states is an upper-bound illustration only and cannot be the primary baseline.

For each scenario and noise level, average time to progression across five replicates, capped at \(H\). Calculate the lower-decile mean (CVaR at the lower 10%) across the 512 scenario-level averages, then take the minimum over \(\sigma_y=0.05,0.10,0.20\) as the primary robust outcome. Higher is better. Report each noise-level score, proportion capped at \(H\), median time, burden integral, resistant fraction, and total dose secondarily. Use paired scenario-level max-T bootstrap intervals for all required comparisons. The 512 evaluation scenarios yield about 51 scenario averages in the lower tail at each noise level; the interval must reflect the small tail sample.

The full hybrid supports a within-simulator advantage only if its adjusted paired interval favors it over every required ablation and baseline. Define quartile cutpoints for initial resistant fraction, growth rate, and drug effect from the training distribution before opening evaluation results. If rankings reverse across these prespecified strata, report scenario-dependent results and no universal preference. Results apply only to the stated simulator and sampled ranges. They do not show that an actual drug works or that a treatment policy helps patients.

### 4.5 Weakest point

The treatment model and parameter ranges are assumed, not measured from patients or experiments. More simulations cannot fix structural misspecification. Even a robust result across the stated noise and parameter distributions can fail under different tumor dynamics, drug mechanisms, or resistance evolution. Track C is an algorithm comparison inside one model family, not treatment evidence.

## 5. Prior-art boundary

These method families already exist: donor-based holdout and cross-domain annotation for single-cell data; graph regularization and sparse feature selection; reaction-diffusion tumor models and inverse parameter fitting; probabilistic spatial segmentation; and model-informed reinforcement learning for GBM therapy. The candidate contributions are the exact tested combinations and their explicit failure rules: donor-balanced sparse graph learning with selective evaluation; joint parameter uncertainty and spatial boundary discrepancy in longitudinal mask forecasting; and hidden-resistance particle filtering with a worst-noise lower-tail dose objective under a fixed dose budget. This protocol does not establish that those combinations are novel or superior. Search each exact combination and compare source code and methods before making a novelty claim.

## 6. References and source records

### Data and source records

- UCSC Cell Browser, multiomic-gbm: https://cells.ucsc.edu/?ds=multiomic-gbm
- MU-Glioma-Post v1, TCIA: https://www.cancerimagingarchive.net/collection/mu-glioma-post/
- RHUH-GBM, TCIA: https://www.cancerimagingarchive.net/collection/rhuh-gbm/
- LUMIERE Figshare item (terms unresolved): https://springernature.figshare.com/articles/dataset/LUMIERE_dataset_-_MRI_data_and_automated_segmentations/21249516

### Literature

Sonpatki, P., Park, H. J., Xing, Y. L., Han, K. Y., Schroeder, B. A., Yu, H. J., Kim, H. J., Chowdhury, T., Hwang, J. H., Nam, S. M., Byun, Y. H., Kang, H., Lee, J. H., Lee, S. T., Won, J. K., Kim, T. M., Choi, S. H., Ku, J.-L., Lee, S., … Shah, N. (2026). A spatially resolved human glioblastoma atlas reveals distinct cellular and molecular patterns of anatomical niches. Nature Communications, 17, 2951. https://doi.org/10.1038/s41467-026-69716-2

Cepeda, S., García-García, S., Arrese, I., Herrero, F., Escudero, T., Zamora, T., & Sarabia, R. (2023). The Río Hortega University Hospital Glioblastoma dataset: A comprehensive collection of preoperative, early postoperative and recurrence MRI scans (RHUH-GBM). Data in Brief, 50, 109617. https://doi.org/10.1016/j.dib.2023.109617

Fischer, F., Fischer, D. S., Mukhin, R., Isaev, A., Biederstedt, E., Villani, A.-C., & Theis, F. J. (2024). scTab: Scaling cross-tissue single-cell annotation models. Nature Communications, 15, 6611. https://doi.org/10.1038/s41467-024-51059-5

Hormuth, D. A., II, Al Feghali, K. A., Elliott, A. M., Yankeelov, T. E., & Chung, C. (2021). Image-based personalization of computational models for predicting response of high-grade glioma to chemoradiation. Scientific Reports, 11, 8520. https://doi.org/10.1038/s41598-021-87887-4

Liu, Z., Zhang, J., Hong, L., Nie, Q., & Sun, X. (2025). Multiscale mathematical model-informed reinforcement learning optimizes combination treatment scheduling in glioblastoma evolution. Science Advances, 11(32), eadv3316. https://doi.org/10.1126/sciadv.adv3316

Mahmoud, E., Gass, J., Dhemesh, Y., Greaser, J., Pogorzelski, K., Isufi, E., Garrett, F., Thacker, J., Tahon, N. H. M., Sinclair, J., Layfield, L. J., Altes, T., & Nada, A. (2025). MU-Glioma Post: A comprehensive dataset of automated MR multi-sequence segmentation and clinical features. Scientific Data, 12, 1847. https://doi.org/10.1038/s41597-025-06011-7

Monteiro, M., Le Folgoc, L., Coelho de Castro, D., Pawlowski, N., Marques, B., Kamnitsas, K., van der Wilk, M., & Glocker, B. (2020). Stochastic segmentation networks: Modelling spatially correlated aleatoric uncertainty. Advances in Neural Information Processing Systems, 33, 12756–12767. https://arxiv.org/abs/2006.06015

Suter, Y., Knecht, U., Valenzuela, W., Notter, M., Hewer, E., Schucht, P., Wiest, R., & Reyes, M. (2022). The LUMIERE dataset: Longitudinal glioblastoma MRI with expert RANO evaluation. Scientific Data, 9, 768. https://doi.org/10.1038/s41597-022-01881-7

Subramanian, S., Scheufele, K., Mehl, M., & Biros, G. (2020). Where did the tumor start? An inverse solver with sparse localization for tumor growth models. Inverse Problems, 36(4), 045006. https://doi.org/10.1088/1361-6420/ab649c

West, J., Adler, F., Gallaher, J., Strobl, M., Brady-Nicholls, R., Brown, J., Roberson-Tessi, M., Kim, E., Noble, R., Viossat, Y., Basanta, D., & Anderson, A. R. A. (2023). A survey of open questions in adaptive therapy: Bridging mathematics and clinical translation. eLife, 12, e84263. https://doi.org/10.7554/eLife.84263

Zhong, Z., Hou, J., Yao, Z., Dong, L., Liu, F., Yue, J., Wu, T., Zheng, J., Ouyang, G., & Song, J. (2024). Domain generalization enables general cancer cell annotation in single-cell and spatial transcriptomics. Nature Communications, 15, 1929. https://doi.org/10.1038/s41467-024-46413-6
