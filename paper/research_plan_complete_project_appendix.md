# Technical Protocol for the Final GBM Research Plan

This is the methods companion to the form plan, not a second project. Tracks A, B, and C are equal, use separate data and endpoints, and do not form a patient-level digital twin. Track B's primary aim is reaction-diffusion/PDE tumor-growth forecasting. UCSD and CFB annotation-change analyses are separate exploratory side studies and do not complete that primary aim. The protocol states the intended analyses; earlier exploratory analyses exist and must be disclosed.

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
| MU-Glioma-Post v1, TCIA | Track B PDE development and exploratory patient-level evaluation | Official TCIA page lists 203 subjects and public NIfTI/clinical files | TCIA lists CC BY 4.0 for public packages. The source has been examined in earlier project work; all MU results are exploratory. |
| UCSD-PTGBM v3, TCIA | Separate, partially exposed annotation-change side study; not PDE validation | Official TCIA page lists 178 subjects and 243 timepoints; standard NIfTI images and segmentations total about 44.98 GB across two packages | Current TCIA page lists CC BY 4.0, data citation DOI 10.7937/fwv2-dt74. The locked attempt failed after outcome access; any further analysis is exploratory. |
| CFB-GBM v4, TCIA | Separate exploratory annotation-change side study; not PDE validation | Local package is file-matched to the downloaded release inventory | Confirm exact source terms for the intended use. Target-mask scoring remains unauthorized until registration QC, mask-review lineage, hashes, and an explicit side-study protocol are complete. |
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

## 3. Track B primary: uncertainty-aware reaction-diffusion forecasting

### 3.1 Question and current evidence

The primary Track B question is whether a patient-specific reaction-diffusion/Fisher–Kolmogorov model can forecast future tumor-core geometry better than persistence, volume-trend, and geometric-trend baselines, with uncertainty and parameter identifiability evaluated explicitly. A forecast target is a later tumor-core mask defined before scoring; the patient, not the scan pair or voxel, is the independent evaluation unit.

The saved MU-Glioma-Post patient run is relevant evidence but remains exploratory because this cohort has been used before. Its all-pairs result was inconclusive against persistence and showed no advantage over the geometric-trend comparator. A first-pair sensitivity was positive against persistence but not clearly better than the geometric baseline. The reported runner conditions on the realized future scan interval, so it is a retrospective conditional-horizon evaluation, not a schedule-known-at-origin forecast. The synthetic joint-posterior hybrid is not yet a patient-data scoring pipeline.

### 3.2 Candidate estimand and cohort handling

The candidate primary design uses two eligible historical tumor-core masks to estimate patient-specific state/parameters and forecasts the next eligible mask. Later rolling origins are secondary. Use only patients with sufficient prior masks and a target inside a prespecified scan-gap range; define tumor core consistently as the source-approved label set and record surgery, treatment, registration, and annotation provenance. If the target interval is the realized interval, state clearly that the estimand is conditional on that interval. Do not call it a fixed-horizon or scheduled forecast.

**Duke external-test candidate (2026-10-28):** use each patient's first two eligible visits as history and the immediate third eligible visit as the single target. Forecast at origin day +90; include the patient only when the published next-visit gap is 76–104 days. This is a nominal 90-day target-window forecast, not exact-day-90 validation. It differs from the historical MU realized-gap evaluation. The external protocol remains a candidate; the user reconfirmed on 2026-10-29 that an access request was submitted, but school-use approval and Duke access are not confirmed. A header-only MU screen found 18 metadata/window candidates. An MRI-only screen read 54 T1n volumes and found high median image-support overlap, but these descriptive checks do not certify per-case registration or valid support. No Duke-specific patient-data hybrid runner exists yet. Patient overlap, date quality, label support, and longitudinal registration gates remain open. See the protocol draft, `output/fk_mu_target_window_screen_20261028/metadata_screen.json`, and `output/fk_mu_mri_only_qc_20261028_v2/mri_only_qc.json`.

MU is development/exploratory only. Keep every patient wholly within one outer fold. Fit volume shrinkage, predictive variance, boundary calibration, and any tuned comparator settings using outer-training patients only. Do not estimate a population rho-D prior from per-patient point fits: with one transition, the pair can lie on a long ridge. Use the protocol's fixed log-uniform prior over prespecified computational bounds and report its bounds sensitivity. Fit each patient's posterior only from history masks; never use that patient's target mask to fit its parameters. Hold target data out until the pipeline, hashes, and failure rules are frozen.

### 3.3 Candidate model and required comparisons

The leading candidate is a joint posterior over reaction-diffusion parameters (including the rho-D ridge when present), combined with a training-only patient-level volume-trend model and correlated boundary uncertainty. This is a candidate, not a proven winner. Reconcile the existing patient runner and synthetic hybrid implementation; freeze the exact PDE, solver, parameter bounds, image/mask resolution, posterior draws, nuisance estimates, seeds, and fallback rules before a new patient run.

Compare the full method with persistence, volume-trend, geometric-train-rate, and prespecified component ablations on exactly the same patients, origins, target masks, and horizons. The primary outcome is patient-macro Dice for the forecast tumor core. Report patient-bootstrap intervals for paired differences. Report surface distance, volume error, and proper probabilistic scores as secondary outcomes; count every patient-level failure in cohort flow.

### 3.4 Identifiability, calibration, and external validation

Report parameter identifiability separately from forecast performance. Show posterior/profile uncertainty, parameter-bound hits, and the rho-D correlation. Do not give fitted parameters a biological interpretation if the data do not separate them. Evaluate uncertainty calibration on held-out patients and include uninformative or failed forecasts rather than dropping them silently.

No untouched external cohort is currently ready. MU, LUMIERE, RHUH, UCSF, and other previously examined sources cannot be relabeled as fresh confirmation without an exact exposure audit. Duke is a restricted candidate, not an approved or ready test set: confirm use terms, access, patient overlap, registration, date quality, and mask provenance before download or label access. CFB-GBM and UCSD-PTGBM annotation-change analyses are separate exploratory side studies; neither can complete the PDE aim.

### 3.5 Failure interpretation

If the PDE does not beat the strongest prespecified simple baseline, report no demonstrated added value. If parameter posteriors remain ridged or hit bounds, do not claim patient-specific biological parameters. If registration, treatment history, mask provenance, scan-gap eligibility, or source rights cannot be established, stop before new patient scoring. Synthetic tests establish code behavior only; they do not establish patient forecast accuracy, novelty, or clinical value.

## 4. Separate exploratory side study: next-visit tumor-core annotation-change ranking

### 4.1 Question, sources, and eligibility

The question is whether a screened graph-diffusion prior combined with an MRI residual attenuated by fixed boundary-perturbation sensitivity improves patient-level ranking of new standard tumor-core annotations at the next available visit. The endpoint is change in released labels. It is not tumor-cell growth, recurrence, progression, response, or a fixed-horizon clinical forecast.

Use MU-Glioma-Post v1 only for method fitting and development-only comparator selection. The project previously examined this source, so all MU metrics remain exploratory. The external source is UCSD-PTGBM v3 from [TCIA](https://www.cancerimagingarchive.net/collection/ucsd-ptgbm/), which lists the image packages and clinical workbook as CC BY 4.0 and requires citation of dataset DOI [10.7937/fwv2-dt74](https://doi.org/10.7937/fwv2-dt74). TCIA lists 178 subjects and 243 visits in v3, split across an Images and Segmentations package and a BraTS-GLI 2024 Test package, totaling approximately 44.98 GB. The current clinical workbook is pinned locally by SHA-256 `bbf40bbdd289266924869052806de4d7e8fc9ee73286cdc0cc066672ef5bb51d`; older workbook versions are not accepted. NIfTI bulk downloads use IBM Aspera. Keep all raw files under the Git-ignored `data/external/UCSD-PTGBM-v3/`; do not move them to Kaggle or commit patient-level files.

Use exact primary diagnosis `Glioblastoma`. Reject duplicate patient/timepoint rows. Group by TCIA subject and sort by numeric visit suffix. If both metadata acquisition-to-anchor values are present, require that the later visit suffix is chronologically later; stop on any contradiction. For each subject, choose at most the first adjacent pair with exactly one standard BraTS segmentation and exactly one each of T1-native, T1-contrast, T2, and FLAIR at each visit. Require internal modality-to-segmentation alignment, millimetre spatial units, the same shape and affine across visits within 0.001 mm, and orthogonal 1.000 ± 0.001 mm spacing on each axis. Resolve files by complete TCIA/BraTS ID tokens; do not permit visit suffix substring collisions. Pair selection may use only diagnosis, source chronology, modality/file completeness, and image-header geometry. It must not depend on target label counts or visual target quality. Use the earlier scan as origin and the next available scan as target; report actual interval in prespecified groups ≤90, 91–180, and >180 days, plus unknown. The interval is not a model input.

### 4.2 Target and risk set

Use only standard compartment labels. At visit t, let C_t(v)=1 for NETC or ET (source codes 1 or 3); SNFH (2) is a current-scan predictor and the resection cavity (4) is excluded. Operational brain support B_0 is the union of nonzero voxels in the four origin MRI modalities and any nonzero origin label. Let Q_0 be the origin cavity. The at-risk voxels are

\[
R_0 = B_0 \setminus (C_0 \cup Q_0), \qquad
Y(v)=\mathbf 1\{v\in C_1\setminus C_0\},\;v\in R_0.
\]

This is the full operational brain support, without a hand-picked 20 mm crop. Report additions outside R_0, including additions in current cavity and outside current support. Do not use the dataset's cellular-tumor segmentations because their definitions use prior/subsequent longitudinal information. TCIA describes the standard segmentation as software-initialized and manually corrected by trained radiologists; this does not make a single released mask error-free truth.

### 4.3 Hybrid method

Robust-normalize each origin MRI modality over B_0. On graph nodes B_0\Q_0, solve

\[
(L_A + \ell^{-2}I)u=\mathbf 1_{C_0}, \qquad \ell=10\,\mathrm{mm},
\]

with a six-neighbor graph and edge conductance

\[
A_{vw}=h_{vw}^{-2}\exp\!\left[-0.5\,\frac14\sum_m(z_m(v)-z_m(w))^2\right]
\exp\!\left[\frac{\log 2}{2}\{\mathbf1_{SNFH}(v)+\mathbf1_{SNFH}(w)\}\right].
\]

Here h_vw is physical neighbor spacing and z_m is the robust-normalized modality m. The graph field is normalized to its maximum over R_0 and used as a spatial ranking prior. It is not a cell-density field or reaction-diffusion estimate.

The ridge-logistic residual uses seven origin-only voxel features: distance to current core, distance to current cavity, four normalized MRI intensities, and current SNFH indicator. No transition history, personalized PDE parameter, interval, or future-derived quantity enters a prediction. Recompute the graph prior after fixed one-voxel erosion and dilation of the current core. Let I(v) be the median absolute change in log-prior across those two perturbations, let s be the median across MU patients of the within-patient 90th percentile of I over R_0, and set g(v)=exp(-I(v)/s). This deterministic gate shrinks the learned residual in spatially sensitive regions; it is not calibrated uncertainty.

For voxel feature vector x(v), the full ranking score is

\[
S(v)=\beta\log u(v)+g(v)[\theta_0+x(v)^T\theta],\qquad\beta\ge0.
\]

Fit the ridge-logistic head on MU training patients only, with equal patient/class weights, ridge penalty 1, prior penalty 10^-3, at most 4,000 sampled positive and negative voxels per patient, and base seed 20261026. Development-fold fits use the base seed plus the one-based fold index; the final MU fit uses the base seed. The gate scale is estimated only from the current training fold. Use the source's existing manifest dataset code `MU` (the canonical MU-Glioma-Post name is accepted as an alias); select patients marked eligible and GBM and their earliest primary pair with nonempty input core. Fail if fewer than five unique development patients remain or any fitting fold lacks the two classes needed by the head.

### 4.4 Development and one-shot external attempt (historical)

This subsection records the protocol that preceded the UCSD run. The one-shot attempt exposed outcomes for at least five patients before failing; the later complete-case run is exploratory and negative. These safeguards document the historical method and do not authorize a fresh confirmatory run.

On MU, run five-fold patient-grouped cross-validation. Fit the full hybrid and prespecified alternatives in each training fold. The non-full alternatives are: ungated hybrid, MRI-only ridge head, diffusion-only score, negative distance-to-core, and uniform no-change score. Select exactly one primary comparator by the highest mean patient-specific development AUC; ties follow that fixed order. Fit the final hybrid and comparator heads on all eligible MU development patients. Do not tune any spatial parameter or model setting on UCSD.

Before reading any UCSD segmentation voxel array, the protocol and all scoring modules must be committed and clean in Git, the local v3 workbook must match its pinned SHA-256, and development/source manifests are hashed. Header-only preflight may read paths, NIfTI shapes, and affines, but not image or target label values. Create an exclusive `analysis_started.json` marker after metadata/header preflight and development fitting; record the locked code digests and selected UCSD path/size manifest. Then hash the complete selected MRI and segmentation file contents, write that digest into the marker, and only then load any UCSD label array. Re-hash selected files after scoring and stop without an aggregate result if any source file changed. The scoring command has no alternate output-path option and refuses to run again if that marker or aggregate result exists. If a failure occurs after the marker, stop and document it before any amended rerun; do not silently restart until a favorable result appears. Write aggregate metrics only, with no subject IDs, per-patient scores, masks, or voxel predictions.

For each test patient whose target in R_0 contains at least one positive and one negative voxel, calculate patient-specific ROC AUC for each score. The primary estimand is the equal-patient mean of the paired difference between the full hybrid AUC and the single comparator selected using MU only. Obtain its two-sided 95% percentile interval from 10,000 patient-level paired bootstrap resamples (seed 20261126). The added-value rule requires the interval lower bound to exceed 0.02 AUC; this project-set margin is not a clinically validated minimum important difference. The metadata maximum of 48 repeat-eligible GBM patients is not a power calculation, and no adequate-power claim is made before image- and endpoint-eligibility are known. Do not assign AUC to patients with no additions or no negatives; keep them in cohort flow and report separately. Report all arm-level patient-macro AUCs, evaluable n, interval-stratified descriptive differences, observed intervals, zero-addition pairs, additions outside R_0, and technical exclusions. No inferential tests are run on secondary arms.

### 4.5 Source and interpretation limits

A positive result supports only ranking of changes in source-provided NETC+ET labels for this public cohort under the specified next-visit sampling rule. It does not show biological tumor growth, a clinically useful forecast, a fixed-day forecast, or superiority to all published methods. UCSD is public and includes BraTS challenge material; other methods may have used these same cases. Cross-cohort overlap with MU cannot be ruled out from a documented patient-level crosswalk. The external test is protocol-locked but is not an externally preregistered clinical validation study.

### 4.6 Weakest point

Only 48 GBM patients have at least two visits in the v3 metadata before image completeness checks. Their observed intervals vary widely, and treatment, registration, contour changes, and source-label corrections can all contribute to target changes. The final event-positive, AUC-evaluable cohort may be small. The test can therefore be imprecise or inconclusive even if the code runs. A future study with an independent, prospectively specified cohort and repeated annotation would be needed to establish biological or clinical value.
## 5. Track C: robust virtual treatment control

### 5.1 Simulator and scenario design

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

### 5.2 Particle filter and validation

Represent the hidden state at each cycle by total burden B and resistant fraction fR=R/B, with static parameter particles for rho, kappa, fR,0, c, epsilon, and B0. Use 500 particles and the Gaussian observation likelihood at the known sigma_y for the in-distribution scenarios. Systematically resample when effective sample size falls below 250. Apply Liu–West regularization with shrinkage a=0.98 to standardized static-parameter particles, using kernel covariance (1-a^2) times the weighted covariance and reflecting proposals at prior bounds. After a parameter move, preserve the current B and fR state and reconstruct S=(1-fR)B and R=fRB; parameter jitter must not add or remove current burden.

On the separate 200-scenario filter-validation set, check nominal 90% intervals for burden and resistant fraction; require empirical coverage from 80% to 97%. Across five filter seeds, require posterior probabilities P(fR>=b) to differ by no more than 0.10 at 90% or more of decision times. These gates validate the filter only under the stated simulator. If they fail, retain the policy results but make no reliable hidden-resistance estimation claim.

### 5.3 Hybrid policy and training

At cycle zero, set the reference burden to the filter's posterior mean of B0. At each later cycle, compute the filtered burden ratio r_j=E[B(t_j)|observations]/E[B0|observations at cycle zero]. Tune a_on and a_off from {1.05, 1.10, 1.20, 1.30}, resistance threshold b from {0.05, 0.10, 0.20, 0.35, 0.50}, and resistance-mode threshold a_R from the same burden-threshold grid. Require a_off < min(a_on,a_R). When treatment is off, begin a five-day course if r_j reaches a_on, or if P(fR(t_j)>=b)>=0.80 use a_R in place of a_on and begin the course if r_j reaches a_R. When treatment is on, give the next five-day course only while r_j remains at or above a_off. Stop dosing after 12 courses or progression.

Use five-fold cross-validation grouped by parameter-scenario ID within the 400 development scenarios. Fit each fold's learned component on four folds and evaluate threshold settings on the held-out fold. Average the two noise replicates within each scenario and noise level. Select policy thresholds by maximizing the minimum lower-decile mean time to progression across the three noise levels in pooled out-of-fold predictions. Train and select PPO seeds/checkpoints using training-only out-of-fold predictions, with no evaluation-scenario feedback. The separate 200-scenario filter-validation set tests calibration gates and does not select policy thresholds.

### 5.4 Comparisons and outcome

Compare the full hybrid with: (1) filtered burden-only hysteresis; (2) raw-observation burden-only hysteresis; (3) resistance-only thresholding; (4) a point-estimate filter ablation; (5) threshold selection for mean time rather than lower-decile time; (6) a fixed front-loaded schedule; (7) an evenly paced schedule; (8) the best periodic schedule selected on development scenarios; and (9) PPO trained for the same lower-tail objective with the same observation, filter-summary, treatment-history, remaining-budget, and elapsed-cycle inputs. Use five PPO seeds and select each seed/checkpoint using training-only folds. A schedule optimized with access to true hidden parameters is an unattainable upper-bound illustration only, never a primary comparator.

For each held-out scenario and noise level, average time to progression over its two replicates, capped at 1,456 days. The primary outcome is the minimum across the three noise levels of the lower-decile mean (CVaR at 10%) of these scenario-level times. Higher is better. Use paired scenario bootstrap intervals with max-T adjustment across the prespecified comparisons. Report each noise-level CVaR, capped fraction, median time, burden integral, final resistant fraction, and course count secondarily. Report each misspecification suite separately. Call the policy robust across the four tested departures only if its adjusted paired interval favors the strongest development-selected baseline in all four suites; otherwise name which departures change the ranking. Define strata fR,0<0.01 versus >=0.01 and c<0.125 versus >=0.125 before opening held-out outcomes. Report a ranking reversal only when paired intervals favor different methods in opposite directions in at least two strata.

### 5.5 Weakest point

The treatment model and parameter ranges are assumed, not measured from patients or experiments. A filter may appear calibrated because its assumed likelihood and simulator agree, and all policies may inherit the same structural error. More simulations reduce Monte Carlo error but cannot show that a real drug works or that a policy would help a patient. Track C is an algorithm comparison inside one stated model family.
## 6. Prior-art boundary

Prior work already includes GBM Core/Periphery cell classification, donor-held-out single-cell annotation, graph and sparse-feature methods, patient-specific reaction-diffusion growth models, probabilistic mask forecasting, spatial uncertainty models, and model-informed adaptive therapy. Patient-specific parameter estimation from longitudinal tumor distributions has been explored with reduced-order models and neural networks ([Cerrone et al., arXiv:2412.05330](https://arxiv.org/abs/2412.05330)). The 2026 [PREDICT-GBM platform](https://doi.org/10.1038/s41746-026-03194-0) benchmarks physics-based and data-driven models across 243 patients for geometric coverage of later recurrence; its TUM-GBM, LUMIERE, and RHUH cohorts overlap with sources already examined in this project. That recurrence-coverage endpoint is related to, but not identical with, future-core PDE forecasting. Other close work includes Lopes and Vinga (2020) on GBM Core versus Periphery cells, scTab and scANVI for annotation, Jung and Shin’s 2026 preprint on a label-derived tumor-state proxy using MU, TaDiff and mechanistic/diffusion MRI-generation work, and Zhang et al. (2017) and Liu et al. (2025) on adaptive evolutionary treatment. Therefore, “PDE tumor-growth forecasting” is not itself a novelty claim. A narrower contribution must be established by a current, focused search and fair, input-matched comparisons on the same patients and endpoints. The three proposed combinations remain hypotheses, not established novelties; do not claim “first” or “beats all current methods.”

## 7. References and source records

### Data and source records

- UCSC Cell Browser, multiomic-gbm: https://cells.ucsc.edu/?ds=multiomic-gbm
- MU-Glioma-Post v1, TCIA: https://www.cancerimagingarchive.net/collection/mu-glioma-post/
- UCSD-PTGBM v3, TCIA: https://www.cancerimagingarchive.net/collection/ucsd-ptgbm/

### Literature

Sonpatki, P., Park, H. J., Xing, Y. L., Han, K. Y., Schroeder, B. A., Yu, H. J., Kim, H. J., Chowdhury, T., Hwang, J. H., Nam, S. M., Byun, Y. H., Kang, H., Lee, J. H., Lee, S. T., Won, J. K., Kim, T. M., Choi, S. H., Ku, J.-L., Lee, S., … Shah, N. (2026). A spatially resolved human glioblastoma atlas reveals distinct cellular and molecular patterns of anatomical niches. Nature Communications, 17, 2951. https://doi.org/10.1038/s41467-026-69716-2

Jung, Y., & Shin, M. (2026). Observation-Anchored Selective Assimilation for Longitudinal Tumor-State Proxy Forecasting in Post-Treatment Glioma. arXiv:2609.12435 [preprint]. https://doi.org/10.48550/arXiv.2609.12435

Lopes, M. B., & Vinga, S. (2020). Tracking intratumoral heterogeneity in glioblastoma via regularized classification of single-cell RNA-Seq data. BMC Bioinformatics, 21, 59. https://doi.org/10.1186/s12859-020-3390-4

Fischer, F., Fischer, D. S., Mukhin, R., Isaev, A., Biederstedt, E., Villani, A.-C., & Theis, F. J. (2024). scTab: Scaling cross-tissue single-cell annotation models. Nature Communications, 15, 6611. https://doi.org/10.1038/s41467-024-51059-5

Hormuth, D. A., II, Al Feghali, K. A., Elliott, A. M., Yankeelov, T. E., & Chung, C. (2021). Image-based personalization of computational models for predicting response of high-grade glioma to chemoradiation. Scientific Reports, 11, 8520. https://doi.org/10.1038/s41598-021-87887-4

Liu, Z., Zhang, J., Hong, L., Nie, Q., & Sun, X. (2025). Multiscale mathematical model-informed reinforcement learning optimizes combination treatment scheduling in glioblastoma evolution. Science Advances, 11(32), eadv3316. https://doi.org/10.1126/sciadv.adv3316

Mahmoud, E., Gass, J., Dhemesh, Y., Greaser, J., Pogorzelski, K., Isufi, E., Garrett, F., Thacker, J., Tahon, N. H. M., Sinclair, J., Layfield, L. J., Altes, T., & Nada, A. (2025). MU-Glioma Post: A comprehensive dataset of automated MR multi-sequence segmentation and clinical features. Scientific Data, 12, 1847. https://doi.org/10.1038/s41597-025-06011-7

Gagnon, L., Gupta, D., Mastorakos, G., et al. (2025). The University of California San Diego annotated post-treatment high-grade glioma multimodal MRI dataset (Version 3) [dataset]. The Cancer Imaging Archive. https://doi.org/10.7937/fwv2-dt74

Monteiro, M., Le Folgoc, L., Coelho de Castro, D., Pawlowski, N., Marques, B., Kamnitsas, K., van der Wilk, M., & Glocker, B. (2020). Stochastic segmentation networks: Modelling spatially correlated aleatoric uncertainty. Advances in Neural Information Processing Systems, 33, 12756–12767. https://arxiv.org/abs/2006.06015

Xu, C., Lopez, R., Mehlman, E., Regier, J., Jordan, M. I., & Yosef, N. (2021). Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models. Molecular Systems Biology, 17(1), e9620. https://doi.org/10.15252/msb.20209620

Subramanian, S., Scheufele, K., Mehl, M., & Biros, G. (2020). Where did the tumor start? An inverse solver with sparse localization for tumor growth models. Inverse Problems, 36(4), 045006. https://doi.org/10.1088/1361-6420/ab649c

West, J., Adler, F., Gallaher, J., Strobl, M., Brady-Nicholls, R., Brown, J., Roberson-Tessi, M., Kim, E., Noble, R., Viossat, Y., Basanta, D., & Anderson, A. R. A. (2023). A survey of open questions in adaptive therapy: Bridging mathematics and clinical translation. eLife, 12, e84263. https://doi.org/10.7554/eLife.84263

Zhang, J., Cunningham, J. J., Brown, J. S., & Gatenby, R. A. (2017). Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer. Nature Communications, 8, 1816. https://doi.org/10.1038/s41467-017-01968-5

Zhong, Z., Hou, J., Yao, Z., Dong, L., Liu, F., Yue, J., Wu, T., Zheng, J., Ouyang, G., & Song, J. (2024). Domain generalization enables general cancer cell annotation in single-cell and spatial transcriptomics. Nature Communications, 15, 1929. https://doi.org/10.1038/s41467-024-46413-6
