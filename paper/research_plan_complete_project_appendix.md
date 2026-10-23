# Technical Protocol Appendix: Three Uncertainty-Aware Hybrid Methods for GBM

This appendix specifies the final proposed analyses, not a chronological record of exploratory attempts. Earlier exploratory analysis exists in the repository. Describe that history truthfully in any submission, and label reused cases exploratory. Tracks A, B, and C have equal standing, separate data, and separate endpoints. No verified patient-level link connects their cohorts.

The shared methodological idea is to combine complementary methods and measure when uncertainty makes an answer unreliable. Each track has a task-specific predictor or simulator, an added uncertainty/robustness component, and an explicit test of whether the full hybrid improves on its parts. Track A addresses donor shift, Track B addresses uncertain MRI boundaries and weak parameter identification, and Track C addresses hidden resistance and uncertain drug effects. The common claim is methodological; no patient-level model joins the three datasets.

## Shared design

1. Keep molecular, MRI, clinical-association, and simulation evidence separate. Do not describe them as one patient-level digital twin.
2. Record each source release, exact files, checksums, labels, identifiers, terms, preprocessing, and exclusions before modeling.
3. Use donor, patient, or virtual parameter scenario as the independent evaluation unit. Cells, voxels, and repeated simulations are not independent subjects.
4. Fit preprocessing, feature selection, tuning, uncertainty calibration, and policy selection only on training donors, patients, or scenarios.
5. Freeze primary outcomes, comparisons, uncertainty methods, and failure rules before outer-fold evaluation. Report every exclusion, failed fit, unstable estimate, and non-estimable endpoint.
6. Use the same cases, input information, and evaluation metric for paired method comparisons. A hybrid counts as better only if it improves its prespecified primary criterion over its components and strongest feasible baseline on held-out units.
7. Treat previously used evaluation cases as exploratory. Changing the model or metric does not make a used test group new confirmation.

## Track A: donor-robust sparse graph prediction

### Target and data

Use the publicly downloadable UCSC multiomic-gbm atlas. The primary target is source-labeled GBM Core versus GBM Periphery in dissociated cells. This is a label-prediction task; it does not prove that the labels capture all GBM biology. Keep NormalBrain and non-GBM samples outside this target. Include only cells with source-traceable labels, assays, and donor identity. Do not count unresolved sample IDs as separate donors.

The Cell Browser offers downloads and bulk rsync, but its page does not show an explicit license for each export. Record the exact files and terms. Keep data local until reuse terms are clear. If donor identities or labels cannot be verified, stop scoring rather than guess.

### Proposed method

Let xᵢ be the expression vector, yᵢ the Core/Periphery label, and dᵢ the donor. Normalize genes using training-donor means and standard deviations only. Learn a nonnegative sparse gene gate w∈[0,1]^G and a two-class linear softmax model qᵢ=softmax(W(w⊙xᵢ)+b). Build a weighted 15-nearest-neighbor cosine graph G among training cells from the training-normalized measured-gene vectors, with edge weights exp(cosine/T); symmetrize by averaging reciprocal weights. Define its symmetric normalized Laplacian as L_G=I−D_G^(−1/2)AD_G^(−1/2), where D_G is the graph degree matrix. Let Q contain the two predicted class probabilities for all training cells. Fit W, b, and w by minimizing:

J(W,b,w) = max over training donors d of L_d(W,b,w) + λ₁‖w‖₁ + (λ_W/2)‖W‖_F² + λ₂ tr(Qᵀ L_G Q).

L_d is mean cross-entropy within donor d. The maximum donor loss reduces domination by donors with many cells. The sparse term drives irrelevant gene gates toward zero; the ridge term on W prevents classifier coefficients from growing without bound to cancel shrinking gates. The graph term encourages connected training cells to have similar class probabilities. Train with Adam at learning rate 0.001 for at most 500 epochs, stopping after 25 epochs without inner-validation improvement. Search λ₁ and λ₂ over {10⁻⁴,10⁻³,10⁻²,10⁻¹,1} and λ_W over {10⁻³,10⁻²,10⁻¹}; tune these and all model settings only inside donor-grouped training folds.

For a held-out query cell i, compute qᵢ from the fitted sparse softmax model, find its 15 nearest training cells by cosine distance in the same normalized measured-gene space, and calculate their normalized weights aᵢⱼ∝exp(cos(xᵢ,xⱼ)/T). Define its graph-smoothed prediction pᵢ=(1−γ)qᵢ+γΣⱼaᵢⱼqⱼ. Select T from {0.1,0.25,0.5,1} and γ from {0,0.25,0.5,0.75,1} inside inner donor folds. Never create test-to-test edges or use donor identity, test labels, tissue annotations, or clinical outcomes as features. Rank confidence by max class probability max_k pᵢₖ. Do not claim formal coverage guarantees unless their assumptions are established and checked.

### Procedure and evaluation

1. Inventory source files, assays, cell and spatial labels, donor/sample aliases, measured genes, and terms.
2. Set the target to Core versus Periphery. Exclude cells with unresolved donor or target labels. Report eligible donor and class counts before fitting; stop if there are too few independent donors to form class-supported folds.
3. Build patient-independent donor-grouped outer folds and nested donor-grouped inner folds. Keep all cells from a donor in one fold. Select gene weights, graph settings, regularization, and calibration only inside the inner folds.
4. Fit the full model (J=\max_d L_d+\lambda_1\|w\|_1+(\lambda_W/2)\|W\|_F^2+\lambda_2\mathrm{tr}(Q^T L_G Q)). Tune \(\lambda_W\in\{10^{-3},10^{-2},10^{-1}\}\) with the other settings inside inner folds. Run three matched ablations: (i) set \(\lambda_1=0\) while retaining the ridge penalty, (ii) remove the graph component by setting \(\lambda_2=0\) and \(\gamma=0\), and (iii) replace the worst-donor loss \(\max_d L_d\) with cell-weighted empirical loss. Also run multinomial logistic regression and the existing attention classifier when its input features match. Tune every method using the same inner folds and a fixed search budget. Use a donor-composition-only predictor only as a leakage diagnostic.
5. For every held-out donor, rank predictions by confidence \(\max_k p_{ik}\) and compute risk (misclassification error among retained predictions) at each coverage. The primary outcome is donor-averaged area under this curve from 50% to 100% coverage, using linear interpolation and trapezoidal integration; lower is better. Report donor-level uncertainty, macro-F1, balanced accuracy, calibration, and paired differences as secondary outcomes; never bootstrap individual cells as if independent.
6. Before full nested evaluation, profile one training fold for wall time and peak memory. Keep the 15-neighbor graph sparse. If the expression matrix or fit exceeds the desktop's 34 GB memory, use another compute host only after source terms permit it; otherwise stop before opening outer-fold scores.
7. Evaluate spatial transfer only as a separate endpoint when measured spatial labels, specimen identity, assay genes, and terms are verified. Hold out whole specimens. Never use predicted pseudozones as spatial truth.

### Prior-art boundary and weakest point

scTab uses donor-based test holdouts in a large cross-tissue benchmark. Cancer-Finder applies domain-generalization training to malignant-cell annotation and evaluates spatial spots. Donor splits, domain generalization, graph models, and selective prediction are not new by themselves. The proposed contribution is the exact combination of sparse feature selection, inductive cell-graph regularization, worst-donor optimization, and selective scoring for this source-labeled GBM target. This combination is a candidate gap, not a verified novelty claim; search group-robust learning, graph regularization, selective prediction, and GBM spatial annotations before claiming novelty.

Weakest point: the number of uniquely mapped donors may be too small for stable tuning or calibration. Spatial labels or dataset terms may also fail the gate. The model may not beat simple baselines at useful coverage.

## Track B: probabilistic reaction-diffusion forecasting

### Dataset and cohort

Use MU-Glioma-Post version 1 for longitudinal tumor-mask and volume forecasting. The official TCIA page lists 203 subjects, 596 follow-up study time points, about 11 GB of files, and CC BY 4.0 for its listed NIfTI and clinical data. Verify diagnosis, scan days, MRI-to-mask alignment, voxel spacing, mask labels, treatment context, and missing or duplicate visits. Restrict to source-confirmed GBM for a GBM-specific claim. Define tumor volume as mask voxel count times physical voxel volume, with included labels fixed before analysis.

All scans from one patient stay in one outer fold. Use earlier scans to predict later scans. Determine eligible-patient counts and achievable precision before looking at outer-fold scores. MU includes mixed glioma diagnoses and post-treatment scans; fitted growth describes the observed care pathway unless untreated growth is separately established.

### Proposed hybrid method

The latent tumor concentration \(c(x,t)\) follows the Fisher-Kolmogorov equation:

∂c/∂t = ∇·(D∇c) + ρc(1 − c/K).

Set \(K=1\), initialize \(c(x,t_0)=1\) inside the first observed mask and 0 outside, and use an isotropic scalar \(D\) in the primary fit. Work in physical millimeters and days on a 2 mm isotropic grid, resampling binary masks by nearest neighbor. Define the computational domain as the baseline-mask bounding box expanded by 30 mm in each direction and clipped to the image field of view; use no-flux boundaries at the crop edge. Set the visible contour to \(c=0.5\). Solve with a three-dimensional central finite-difference Laplacian and forward Euler time steps of 0.05 day; at \(D\le10\;\mathrm{mm}^2/\mathrm{day}\), this satisfies the explicit diffusion stability limit \(\Delta t\le\Delta x^2/(6D)\). Search \(\rho\in[10^{-5},10^{-1}]\;\mathrm{day}^{-1}\) and \(D\in[10^{-5},10]\;\mathrm{mm}^2/\mathrm{day}\) on a 21-by-21 logarithmic grid, then refine around the best pair with an 11-by-11 grid. For each \((\rho,D)\) pair, optimize \(\sigma_b\) and \(\tau\) by bounded L-BFGS-B with at most 200 iterations and relative tolerance \(10^{-6}\). Calculate the Hessian by central finite differences in transformed parameter space with relative step \(10^{-4}\). If an optimum hits a bound, expand that bound tenfold using training patients only and report the hit.

The observed mask is a spatially correlated boundary observation, not a set of independent voxel labels. On each follow-up mask, select 256 deterministic, approximately evenly spaced surface points by farthest-point sampling. For point \(p_k\), let \(r_k\) be its signed distance in millimeters to the predicted \(c=0.5\) contour. Model the residual vector as \(r\sim N(0,\Sigma)\), with \(\Sigma_{kl}=\sigma_b^2\exp(-\|p_k-p_l\|/\ell)+\tau^2\delta_{kl}\). Estimate boundary scale \(\sigma_b\in[0,10]\) mm and nugget \(\tau\in[0.1,5]\) mm; these are computational bounds, not measured error limits. Select \(\ell\in\{2,5,10\}\) mm only inside training folds. The likelihood models effective mask-to-mask boundary discrepancy, not direct scanner noise. For uncertainty draws, transform fitted parameters to \(\eta=(\log\rho,\log D,\log(\sigma_b+0.01),\log\tau)\), calculate the observed Hessian of the negative log likelihood at the fitted value, and use a multivariate normal with covariance equal to the inverse Hessian, truncated to fitted parameter bounds. Draw 200 parameter vectors from this Laplace approximation; this is a likelihood-based sampling approximation, not a Bayesian posterior. Select 256 fixed anchor points on each predicted surface. For each selected \(\ell\), precompute the eigendecomposition \(R_\ell=U\Lambda U^T\) of the 256-by-256 exponential correlation matrix at those anchors. For each draw, generate anchor displacement \(\delta=U\,\mathrm{diag}(\sqrt{\sigma_b^2\lambda_j+\tau^2})z\), with \(z\sim N(0,I)\), then interpolate displacement to every other surface vertex from its three nearest anchors using weights proportional to \((d+0.1\;\mathrm{mm})^{-2}\), normalized to sum to one, where \(d\) is Euclidean distance in physical millimeters. Move the contour along local normals and rasterize the mask. This low-rank surface-field approximation avoids a dense covariance over all vertices and smooths between-anchor variation; validate its marginal coverage and forecast score on synthetic fields before patient scoring. If the Hessian is not positive definite or more than 10% of parameter draws hit computational bounds, mark the fit non-estimable for probabilistic forecasting and report the failure.

### Procedure, gates, and evaluation

1. Build a patient-by-visit manifest; verify image geometry, mask lineage, visit order, diagnosis, and volume definition. Keep all visits from one patient in one outer fold. For each held-out patient, fit to the earliest two masks and forecast the next visit; if later visits exist, use rolling-origin forecasts without exposing each target mask to fitting. Require at least three ordered masks for the core forecast endpoint.
2. Before patient fitting, generate 100 synthetic series at observed visit intervals. Create boundary displacement with Gaussian fields at standard deviations 0, 2, and 4 mm and correlation lengths 2, 5, and 10 mm; include anisotropic and misspecified-growth cases. Hide generating parameters during fitting. Benchmark five complete synthetic fits and five training-patient fits, including all 200 predictive draws; record wall time and peak memory. Estimate full-cohort runtime from the measured median. If the planned 2 mm analysis exceeds desktop capacity, use Kaggle only if the exact data terms permit cloud processing and preserve identical frozen settings. If neither route fits available compute, stop and report the feasibility failure before scoring held-out patients; do not silently change grid, solver, or draw count after seeing results.
3. Measure parameter bias, nominal 90% interval coverage and width, forecast score, numerical failures, and bound hits. A separate parameter claim passes the synthetic gate only if median absolute relative error is at most 25%, 90% interval coverage is 80–100%, and no more than 10% of fits hit a computational bound. Treat these as prespecified project gates, not clinical standards.
4. On real data, withhold a separate \(\rho\) or \(D\) interpretation if its synthetic gate fails, its 95% profile interval touches a search bound, or the interval spans more than tenfold. The forecast may still be reported as conditional on the model and parameter sensitivity range.
5. Compare the full method with three required ablations/baselines: (i) the same PDE with independent rather than spatially correlated boundary residuals, (ii) a deterministic \(c\ge0.5\) PDE mask, and (iii) persistence and time-gap-aware volume-trend forecasts. Use the same patient folds and eligible visits.
6. For each held-out visit, estimate voxel tumor probabilities from the 200 Laplace-parameter and correlated-boundary masks. The primary outcome is the mean patient-level Brier score for reproducing supplied held-out segmentation labels inside a 10 mm band around that mask's boundary; lower is better. This evaluates agreement with the supplied mask only. It does not establish the true tumor boundary or measure segmentation uncertainty. Report Dice, surface distance, volume error, interval coverage, calibration, exclusions, and paired patient-level uncertainty secondarily.

### Separate recurrence extension

Use RHUH-GBM only for a separate, optional early-postoperative-to-recurrence analysis. [TCIA collection page](https://www.cancerimagingarchive.net/collection/rhuh-gbm/) lists 40 subjects and a 2.9 GB CC BY 4.0 package of brain-extracted NIfTI images and segmentations. Original DICOM has a distinct NIH controlled-access policy; exclude it.

After verifying chronology, registration, and label meaning, define \(P(x)\) as the training-fitted PDE concentration at the source-defined recurrence interval and \(D(x)=\exp[-d(x,C)/15\;\mathrm{mm}]\), where \(d(x,C)\) is physical distance from the postoperative tumor/cavity mask \(C\). Set \(R_w(x)=wP(x)+(1-w)D(x)\), choosing \(w\in\{0,0.1,\ldots,1}\) by training-patient inner validation to maximize recurrence sensitivity at the fixed 15 mm-expanded volume; break ties by Dice. For each held-out patient, use only the early postoperative scan and training-derived settings; never use that patient's recurrence mask to fit or tune the map. Select the highest-risk brain voxels until their volume equals the physical 15 mm expansion volume. Compare recurrence sensitivity and Dice against that fixed expansion. This is not a radiation plan or clinical margin recommendation. If chronology, registration, labels, or a held-out partition fail, do not score this extension. Do not use TCGA bulk molecular subtype as spatial truth.

### Prior-art boundary and weakest point

PDE-constrained inversion, sparse initial-condition recovery, patient-specific reaction-diffusion fitting, PINN estimation from MRI, and spatially correlated probabilistic segmentation already exist. For example, Monteiro et al. (2020) model joint, spatially coherent uncertainty over medical-image label maps. The candidate gap is to connect a spatial boundary-error distribution to longitudinal PDE fitting and make the forecast and parameter interpretation conditional on explicit identifiability gates. This exact combination remains a candidate gap, not an established novelty claim; compare with Bayesian state-space, data-assimilation, and segmentation-uncertainty methods before claiming novelty.

Weakest point: one segmentation per scan does not directly measure segmentation noise. The fitted spatial discrepancy can absorb treatment effects, registration error, and real growth. The data may not identify \(\rho\) and \(D\) separately even when synthetic recovery passes.

## Track C: robust virtual treatment control

### Simulation scope

Use a well-mixed two-population reduction of the project’s existing sensitive/resistant growth model so policy search is tractable and every comparator can use the same simulator. Let \(S(t)\) and \(R(t)\) be sensitive and resistant burden, \(K=1\) normalized carrying capacity, and \(u(t)\in\{0,1\}\) a model-defined treatment action:

\[
\dot S=\rho S\left(1-\frac{S+R}{K}\right)-\kappa uS,\qquad
\dot R=\rho(1-c)R\left(1-\frac{S+R}{K}\right)-\epsilon\kappa uR.
\]

Here \(c\) is the assumed resistance fitness cost and \(\epsilon\) the assumed residual treatment effect on resistant cells. Make treatment decisions daily and integrate the ODE with a fixed 0.1-day step. No parameter is a measured patient treatment response. Use the existing project virtual-model stress ranges \(\rho\in[0.005,0.05]\) day⁻¹ and \(\kappa\in[0.02,0.08]\) day⁻¹; they are model assumptions, not patient measurements or validated clinical drug effects. Also use \(f_{R,0}\in[0.001,0.1]\), \(c\in[0,0.25]\), \(\epsilon\in[0,0.2]\), and initial burden \(B_0/K=(S_0+R_0)/K\in[0.08,0.8]\). No mutation or phenotype switching is included in the core model. Draw \(\rho\), \(\kappa\), and \(f_{R,0}\) log-uniformly; draw \(c\), \(\epsilon\), and \(B_0/K\) uniformly. Initialize \(S_0=(1-f_{R,0})B_0\) and \(R_0=f_{R,0}B_0\) with \(K=1\). The report must distinguish source-supported parameter values from assumptions.

The controller observes noisy total burden only, \(y_t=\log[(S_t+R_t)/K]+e_t\), with \(e_t\sim N(0,\sigma_y^2)\); resistant fraction is hidden from every policy. Include \(\sigma_y\in\{0,0.05,0.10,0.20\}\) as prespecified noise settings. Use \(\sigma_y=0\) only for an exact-observation simulator check; exclude it from particle-filter fitting, threshold tuning, and the primary policy comparison because a zero-variance Gaussian likelihood is a point mass. Generate 256 Latin-hypercube training scenarios and 256 disjoint evaluation scenarios from the frozen ranges, assigning 85 or 86 scenarios to each positive noise level \(\{0.05,0.10,0.20\}\), with five noise replicates per scenario. Freeze random seeds and ranges before model fitting. Average the five noise replicates within each parameter scenario before computing scenario-level CVaR or uncertainty.

### Proposed method

The proposed hybrid joins a particle filter with a lower-tail robust threshold controller. For positive \(\sigma_y\), initialize 500 particles from the stated parameter and initial-state priors. Map log-uniform parameters linearly in log space to \([0,1]\), and map uniform parameters linearly to \([0,1]\). Propagate each particle through the deterministic ODE with known past actions and update weights using the Gaussian observation likelihood. When effective sample size falls below 250, systematically resample and apply one Liu–West regularization step to standardized static-parameter particles: \(z_i'=a z_i+(1-a)\bar z+\epsilon_i\), where \(a=0.98\), \(\epsilon_i\sim N(0,(1-a^2)V)\), and \(V\) is their weighted covariance; reflect proposals at prior bounds. This approximates static-parameter posterior updates and does not model biological parameter drift or add tumor-state process noise. Define nominal 90% intervals as weighted, equal-tailed 5th and 95th particle quantiles. On a training-only synthetic validation subset, require interval coverage between 80% and 97%; also run five filter seeds and require the range of posterior probabilities \(P(f_R\ge b_j)\) across seeds to be at most 0.10 at at least 90% of validation decision times. These are project stability gates, not standard guarantees. If either gate fails, do not claim reliable resistance estimates or a successful resistance-trigger controller. Report coverage on disjoint evaluation scenarios without recalibration. Let observed normalized burden be \(\tilde B_t=\exp(y_t-y_0)\) and \(f_{R,t}=R_t/(S_t+R_t)\). Give the model-defined dose if \(\tilde B_t\ge a\) or if posterior probability \(P(f_{R,t}\ge b\mid y_{0:t},u_{0:t-1})\ge0.8\); otherwise withhold dose. Search the fixed grid \(a\in\{0.8,1.0,1.2,1.4,1.6\}\) and \(b\in\{0.05,0.10,0.20,0.35,0.50\}\).

Choose \((a^*,b^*)\) to maximize the empirical lower-tail \(\mathrm{CVaR}_{0.10}\) of capped time to progression over training scenarios. Define empirical lower-tail CVaR as the mean of the lowest 10% of scenario values, using fractional weight at the cutoff when needed. Set \(B_{\max}\) to the number of active dose days in the 5-days-on/23-days-off schedule over the locked horizon. Route every daily action from every arm (hybrid, ablations, fixed schedule, and PPO) through one shared cap wrapper: once cumulative active dose days reach \(B_{\max}\), force all later actions to off. Policies may use fewer than \(B_{\max}\) days; call this a common maximum dose budget, not an equal-dose design. Record active dose days and assert that no arm exceeds the cap. Define progression as burden reaching 1.4 times its running nadir; this is a simulation endpoint, not RANO.

Use a training-only horizon pilot over 500, 1,000, and 1,500 days: select the shortest horizon for which fewer than 90% of training scenarios are still progression-free under the fixed schedule. If none meets that rule, retain 1,500 days and declare the primary outcome ceiling-limited if at least 90% of evaluation scenarios are capped. Freeze the horizon, thresholds, ranges, and dose budget before evaluating held-out scenarios.

### Procedure and evaluation

1. Validate positivity, carrying-capacity behavior, zero-dose growth, resistant-cell survival, positive-noise observation likelihood, and dose accounting on analytic or hand-computable cases before policy search. On a training-only synthetic validation subset, require nominal 90% equal-tailed particle-interval coverage between 80% and 97%. Run five filter seeds and require the range of each posterior probability \(P(f_R\ge b_j)\) to be at most 0.10 across seeds at 90% or more of validation decision times. These are project gates, not statistical standards. If either fails, do not claim reliable resistance estimates or a successful resistance-trigger controller. Freeze all model, filter, and policy settings before evaluation.
2. Time the full 25-threshold controller grid on 10 training scenarios and a 100,000-step PPO run for one seed; record wall time and peak memory, then estimate the full training and evaluation cost. Use Kaggle for this simulation-only computation if desktop capacity is insufficient. If neither environment can complete the frozen design, report the feasibility failure before opening evaluation results; do not reduce the grid or simulation count after seeing results.
3. Keep all parameter combinations in training separate from all evaluation scenarios; pair policies on the same scenario and observation-noise seeds.
4. Compare the full hybrid with three ablations: replace the particle filter with an extended Kalman filter point estimate, remove the resistance trigger (the burden-only adaptive baseline), and use mean-time rather than CVaR threshold selection. Compare also with the fixed 5/23 schedule and PPO trained for 1,000,000 environment steps with five fixed seeds and a two-layer 64-unit tanh actor-critic. Give PPO the same information as the hybrid: observed normalized burden, dose history, particle posterior mean and standard deviation, posterior probabilities \(P(f_R\ge b_j)\) for every threshold \(b_j\) in the policy grid, remaining dose budget, and elapsed fraction of the horizon. Use the same dose cap and training scenarios. The PPO reward is one unit per day until the progression threshold; choose its checkpoint only on a training validation split. A true-state policy may be shown only as an oracle ceiling, never as a fair baseline.
5. Use the same locked parameter scenarios, noise replicates, dose cap, and horizon for every method. Report paired lower-decile mean TTP as primary; report median TTP, fraction capped, burden integral, resistant fraction, dose days, and policy-rank changes secondarily. Compute uncertainty by resampling scenario groups, not daily time points.
6. Define four prespecified stress strata by splitting initial resistant fraction at 0.01 and resistance cost at 0.125. If the best policy or the sign of the full-hybrid comparison changes across strata, report scenario-dependent rankings and no universal preference. Do not tune ranges after seeing evaluation outcomes.

### Separate MGMT association

If the MU clinical workbook has a source-verified codebook and timeline, estimate MGMT–progression association as a separate observational analysis. Record diagnosis, event, censoring, time origin, missingness, and covariates. Do not use MGMT to set virtual resistance or treatment response. If a field is ambiguous, report the association as not estimable.

### Prior-art boundary and weakest point

Adaptive therapy, sensitive/resistant models, optimal control, and reinforcement learning are active areas. Liu et al. (2025) combines multiscale GBM modeling with RL for treatment scheduling; Aif et al. (2026) studies cross-environment RL failure as a way to discover missing tumor dynamics; and West et al. (2023) reviews open questions in adaptive therapy. The candidate contribution is narrower: a controller that filters hidden resistance from noisy total-burden observations and selects thresholds for lower-tail performance under parameter uncertainty and a shared dose cap. This exact combination is not yet established as novel; search partially observed adaptive control, robust/CVaR control, and GBM resistance models before making a novelty claim.

Weakest point: all drug effects are assumptions, not patient-measured responses. Any win applies only within the selected simulator and parameter scenarios. Model misspecification can reverse the policy ranking.

## Public data and use terms

| Source | Planned role | Public access and terms | Use rule |
|---|---|---|---|
| UCSC multiomic-gbm | Track A single-cell target; spatial transfer only with measured labels and a source crosswalk | [Cell Browser download page](https://cells.ucsc.edu/?ds=multiomic-gbm) provides downloads and bulk rsync. The page does not state an explicit license for every export. | Keep local and record exact files and terms. No cloud processing or redistribution until terms are clear. Stop scoring if donor IDs or labels cannot be verified. |
| MU-Glioma-Post v1 | Track B longitudinal MRI; separate MGMT association if fields are valid | [TCIA collection page](https://www.cancerimagingarchive.net/collection/mu-glioma-post/) lists 203 subjects, 596 follow-up study time points, about 11 GB, and CC BY 4.0 for listed files. | Use the cited version, follow TCIA policy and attribution, verify GBM eligibility and clinical-field meaning. |
| RHUH-GBM | Track B recurrence extension | [TCIA collection page](https://www.cancerimagingarchive.net/collection/rhuh-gbm/) lists 40 subjects and a 2.9 GB brain-extracted NIfTI/segmentation package under CC BY 4.0. Original DICOM is controlled access. | Use only the public derivative package and eligible public clinical files. Exclude original DICOM. |
| LUMIERE | Optional sensitivity analysis only; not required for core methods | [Figshare item](https://springernature.figshare.com/articles/dataset/LUMIERE_dataset_-_MRI_data_and_automated_segmentations/21249516) is public (about 30.33 GB) and lists CC0; the related paper states noncommercial use. | Reconcile exact file and article terms before use. No cloud upload until permitted. |
| Kaggle | Optional computation | GPU access does not change data terms. | Upload only sources whose terms allow cloud processing. Keep a CPU-reproducible path. |

Public download does not by itself grant unrestricted reuse. An unresolved term is a stop gate.

## Feasibility and limits

No full compute pilot has yet established runtime or memory for Track A nested donor validation, Track B PDE fitting, or Track C particle filtering and PPO. Track A donor counts and data terms also remain gates. Start with synthetic checks and resource pilots, then run one locked evaluation per track only when each track passes its gate. Use Kaggle for Track C's synthetic simulations; for patient data, use Kaggle only when exact dataset terms permit cloud processing. Keep spatial transfer, recurrence, and MGMT secondary to each track's core method. Do not force one implementation across incompatible datasets. Novelty and superiority remain questions to test, not promises.

## Full bibliography

Aif, S., Eiche, M., Appold, N., Fischer, E., Citak, T., & Kayser, J. (2026). Reinforcement Failing guides the discovery of emergent physical dynamics in adaptive tumor therapy. *npj Artificial Intelligence, 2*, Article 76. https://doi.org/10.1038/s44387-026-00157-4

Cepeda, S., García-García, S., Arrese, I., Herrero, F., Escudero, T., Zamora, T., & Sarabia, R. (2023). The Río Hortega University Hospital Glioblastoma dataset: A comprehensive collection of preoperative, early postoperative and recurrence MRI scans (RHUH-GBM). *Data in Brief, 50*, 109617. https://doi.org/10.1016/j.dib.2023.109617

Fischer, F., Fischer, D. S., Mukhin, R., Isaev, A., Biederstedt, E., Villani, A.-C., & Theis, F. J. (2024). scTab: Scaling cross-tissue single-cell annotation models. *Nature Communications, 15*(1), Article 6611. https://doi.org/10.1038/s41467-024-51059-5

Hormuth, D. A., II, Al Feghali, K. A., Elliott, A. M., Yankeelov, T. E., & Chung, C. (2021). Image-based personalization of computational models for predicting response of high-grade glioma to chemoradiation. *Scientific Reports, 11*, Article 8520. https://doi.org/10.1038/s41598-021-87887-4

Liu, Z., Zhang, J., Hong, L., Nie, Q., & Sun, X. (2025). Multiscale mathematical model-informed reinforcement learning optimizes combination treatment scheduling in glioblastoma evolution. *Science Advances, 11*(32), eadv3316. https://doi.org/10.1126/sciadv.adv3316

Mahmoud, E., Gass, J., Dhemesh, Y., Greaser, J., Pogorzelski, K., Isufi, E., Garrett, F., Thacker, J., Tahon, N. H. M., Sinclair, J., Layfield, L. J., Altes, T., & Nada, A. (2025). MU-Glioma Post: A comprehensive dataset of automated MR multi-sequence segmentation and clinical features. *Scientific Data, 12*, Article 1847. https://doi.org/10.1038/s41597-025-06011-7

Monteiro, M., Le Folgoc, L., Coelho de Castro, D., Pawlowski, N., Marques, B., Kamnitsas, K., van der Wilk, M., & Glocker, B. (2020). Stochastic segmentation networks: Modelling spatially correlated aleatoric uncertainty. *Advances in Neural Information Processing Systems, 33*, 12756–12767. https://arxiv.org/abs/2006.06015

Sonpatki, P., Park, H. J., Xing, Y. L., Han, K. Y., Schroeder, B. A., Yu, H. J., Kim, H. J., Chowdhury, T., Hwang, J. H., Nam, S. M., Byun, Y. H., Kang, H., Lee, J. H., Lee, S. T., Won, J. K., Kim, T. M., Choi, S. H., Ku, J.-L., Lee, S., … Shah, N. (2026). A spatially resolved human glioblastoma atlas reveals distinct cellular and molecular patterns of anatomical niches. *Nature Communications, 17*, Article 2951. https://doi.org/10.1038/s41467-026-69716-2

Subramanian, S., Scheufele, K., Mehl, M., & Biros, G. (2020). Where did the tumor start? An inverse solver with sparse localization for tumor growth models. *Inverse Problems, 36*(4), 045006. https://doi.org/10.1088/1361-6420/ab649c

Suter, Y., Knecht, U., Valenzuela, W., Notter, M., Hewer, E., Schucht, P., Wiest, R., & Reyes, M. (2022). The LUMIERE dataset: Longitudinal glioblastoma MRI with expert RANO evaluation. *Scientific Data, 9*, Article 768. https://doi.org/10.1038/s41597-022-01881-7

Yaseen, D., Garrett, F., Gass, J., Greaser, J., Pogorzelski, K., Isufi, E., Layfield, L. J., Nada, A., Sinclair, J., Tahon, N. H. M., & Thacker, J. (2025). University of Missouri Post-operative Glioma Dataset (MU-Glioma-Post) (Version 1) [Data set]. The Cancer Imaging Archive. https://doi.org/10.7937/7K9K-3C83

West, J., Adler, F., Gallaher, J., Strobl, M., Brady-Nicholls, R., Brown, J., Roberson-Tessi, M., Kim, E., Noble, R., Viossat, Y., Basanta, D., & Anderson, A. R. A. (2023). A survey of open questions in adaptive therapy: Bridging mathematics and clinical translation. *eLife, 12*, e84263. https://doi.org/10.7554/eLife.84263

Zhang, R. Z., Ezhov, I., Balcerak, M., Zhu, A., Wiestler, B., Menze, B. H., & Lowengrub, J. S. (2023). Personalized predictions of glioblastoma infiltration: Mathematical models, physics-informed neural networks and multimodal scans. *arXiv*. https://doi.org/10.48550/arXiv.2311.16536

Zhong, Z., Hou, J., Yao, Z., Dong, L., Liu, F., Yue, J., Wu, T., Zheng, J., Ouyang, G., Yang, C., & Song, J. (2024). Domain generalization enables general cancer cell annotation in single-cell and spatial transcriptomics. *Nature Communications, 15*(1), Article 1929. https://doi.org/10.1038/s41467-024-46413-6
