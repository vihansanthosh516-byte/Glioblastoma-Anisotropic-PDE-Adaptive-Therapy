# Reproduction appendix: GBM Tracks A, B, and C

This appendix specifies selected procedures, inputs, comparisons, and failure rules. It reports no project findings and predicts no direction or size of an effect. If any analysis was completed before submission, disclose it as the fair requires. Use the current source files only after recording their release, URL, file list, and file-level terms.

## Project question and scope

The primary mathematical question is how a hierarchical state-space model represents latent postoperative GBM volume change and MRI measurement error, and how its next-scan predictions compare with fixed baselines on held-out patients. Tracks A, B, and C address molecular tumor-zone prediction, MRI growth and recurrence, and treatment-policy simulation. They use different cohorts and endpoints. No verified patient-level link joins those cohorts, so analyze them separately and do not describe the project as a fused patient-specific model or validated digital twin.

## Shared data rules

1. Use the patient or donor as the split and resampling unit.
2. Link records only with identifiers documented by the data source.
3. Record missing visits, exclusions, image modalities, mask labels, voxel dimensions, release versions, and file-level licenses.
4. Keep patient-data analyses, simulations, external literature values, and assumptions labeled separately.
5. Freeze inclusion rules, partitions, model settings, and comparisons before calculating final evaluation metrics.
6. Use no data on a cloud host unless the exact file terms permit it.
7. Report every prespecified comparison and every model failure. Do not loosen criteria after viewing evaluation results.

## Procedure

### 1. Track A — Patient-level molecular tumor-zone prediction

**Question.** How well do molecular features predict Core, Periphery, or Healthy labels for cells from donors excluded from model fitting?

1. Obtain the UCSC multiomic-gbm atlas files selected for analysis. Record the release, source URL, file names, cell count, labels, donor identifiers, and reuse terms. Do not redistribute files until their specific terms are clear.
2. Define a fixed donor-identifier rule from source metadata. Exclude a donor if its identity cannot be resolved. Keep all cells from one donor in one partition.
3. Freeze donor-grouped training, validation, and test partitions before fitting. If class representation or donor count prevents the prespecified comparison, report the limitation and do not move cells between groups to improve a score.
4. Fit the attention-gated classifier on training donors and use validation donors for model selection. Fit two fixed comparators: PCA-50 followed by logistic regression, and logistic regression using the prespecified 200-gene panel. Apply identical donor partitions and preprocessing safeguards.
5. Calculate macro-F1 on test donors. Calculate paired donor-level differences and uncertainty by resampling donors, not cells. Report per-class metrics, exclusions, partition composition, and any fitting or class-support failure.
6. Treat the single-atlas test as internal evaluation. Do not interpret cell-level scores as evidence of performance on a new institution or population.

**Main limitation.** The number of donors is small, and donor composition can be confounded with tumor-zone labels. Normalized data may limit methods that require raw counts.

### 2. Track B — MRI growth, recurrence, and measurement-aware modeling

**Question.** How do selected MRI growth and recurrence methods compare with fixed baselines, and how does a state-space model represent uncertainty in longitudinal postoperative volumes?

1. Build a manifest for each MRI cohort used: MU-Glioma-Post, LUMIERE, RHUH-GBM, and UCSF-PDGM where the selected analysis requires it. Record patient identifier, visit date, sequence, segmentation source, label definition, voxel size, missingness, exclusions, and exact file terms. Keep cohorts separate.
2. Define eligible scans and patient-level partitions before evaluation. Use only preoperative inputs to construct recurrence-risk maps. Do not train or tune a model on a target cohort's recurrence labels.
3. For growth-mask forecasting, use eligible postoperative masks from MU-Glioma-Post. Generate patient-grouped out-of-fold forecasts with the selected 2 mm positivity-preserving reaction-diffusion solver. Select solver parameters within training folds only. Compare forecasts on the same patient visits with persistence and a population geometric-growth baseline. Report paired patient-level Dice differences and uncertainty. Treat tensor-direction comparisons as secondary.
4. For the recurrence comparison, use the selected H-2 map for LUMIERE/RHUH and H-3 map for the separate UCSF sensitivity cohort. Use preoperative inputs only. Upsample the saved model logits to 1 mm, apply sigmoid, and select highest-scoring voxels inside the brain until the predicted region has the same volume as the 15 mm standard plan. Do not refit on target cohorts. Score enhancing recurrence as primary and all recurrence as secondary. Report paired patient-level coverage differences and uncertainty.
5. Define the new measurement-aware model using positive postoperative segmented volume \(V_{ij}\) for patient \(i\) at visit \(j\), with \(y_{ij}=\log V_{ij}\). Require at least three positive-volume visits with known scan days. Define volume as positive tumor-mask voxels multiplied by physical voxel volume. Calculate \(\Delta t_{ij}\) from the documented MRI day.
6. Let \(x_{ij}\) be latent log-volume. Use an unrestricted patient-specific initial level and exact diffuse Kalman initialization. Specify the observation and transition models as:

   \[
   y_{ij}=x_{ij}+\epsilon_{ij}, \qquad \epsilon_{ij}\sim N(0,\sigma_e^2),
   \]

   \[
   x_{i,j+1}=x_{ij}+g_i\Delta t_{ij}+\eta_{ij}, \qquad
   \eta_{ij}\sim N(0,\sigma_w^2\Delta t_{ij}), \qquad
   g_i\sim N(\mu_g,\tau_g^2).
   \]

   Here \(g_i\) is net postoperative log-volume trend over the observed care pathway. It is not untreated biological growth. Assume scan errors are independent across visits and patients, process innovations are independent across intervals and patients, and both are independent of patient slopes and initial states. Omit treatment-phase terms because treatment timing is incomplete or inconsistently recorded.
7. Augment the Kalman state with the constant patient-specific slope \(g_i\). Let \(z_{ij}=(x_{ij},g_i)^\mathsf T\), with transition

   \[
   z_{i,j+1}=
   \begin{bmatrix}1&\Delta t_{ij}\\0&1\end{bmatrix}z_{ij}
   +\begin{bmatrix}\eta_{ij}\\0\end{bmatrix}.
   \]

   Marginalize the latent volume and slope in the linear-Gaussian likelihood. Estimate \((\mu_g,\tau_g,\sigma_w,s_\Delta)\) by maximum a posteriori optimization and use a Laplace approximation for intervals on population parameters only. Do not include diffuse initial states in those intervals.
8. Set the primary priors to \(\mu_g\sim N(0.002,0.005^2)\), \(\tau_g\sim\mathrm{HalfNormal}(0.01)\), and \(\sigma_w\sim\mathrm{HalfNormal}(0.05)\). Fit prespecified sensitivity priors: \(\mu_g\sim N(0,0.01^2)\), and one-at-a-time half-normal scales of 0.005 and 0.02 for \(\tau_g\), and 0.025 and 0.10 for \(\sigma_w\). Report all fits. Accept a Laplace interval only if optimization converges, the log posterior is finite, and the negative Hessian at the mode is positive definite.
9. Define \(s_\Delta\) as the SD of error in a two-scan log-volume difference. Under independent equal-variance scan errors, set \(\sigma_e=s_\Delta/\sqrt{2}\). Estimate \(s_\Delta\) under \(s_\Delta\sim\mathrm{LogNormal}(\log(0.115),0.4)\). Also fit fixed-scale sensitivity cases at 0.079 and 0.342. These are assumptions based on an ambiguous conversion of external repeatability coefficients, not measurements from the project cohorts. Report the scale convention and uncertainty about that conversion.
10. Before fitting patient outcomes, run 200 parameter-recovery simulations at each of \(s_\Delta\in\{0.079,0.115,0.342\}\). Preserve eligible patients' visit times; use each patient's first observed log-volume as the fixed latent starting value; draw slopes from \(N(0.002,0.01^2)\); set \(\sigma_w=0.05\); and generate independent process and observation errors, including observation error at the first visit. Use seeds 41000–41199, 41200–41399, and 41400–41599 for the three settings. Estimate \(s_\Delta\) under the same central prior in every fit. Calculate 90% interval coverage and exact binomial 95% intervals for each simulated parameter by noise setting. State that this evaluates recovery under the chosen model and simulation settings; it does not establish identifiability in patient data.
11. Profile three variance-parameter pairs: scan-error SD versus process-noise SD, scan-error SD versus slope SD, and process-noise SD versus slope SD. Use a 21-by-21 log-scale grid over the central 95% prior ranges. At each point, maximize likelihood without priors over the population mean slope and remaining variance parameter. When fixing scan-error SD, use \(s_\Delta=\sqrt{2}\sigma_e\). Treat these surfaces as descriptive diagnostics, not formal confidence regions or an identifiability proof.
12. For held-out prediction, create five patient-grouped folds with seed 4200 and reuse them for every model. Fit state-space parameters on training patients only. Estimate the population-median slope from training-patient log-volume slopes. Let persistence predict the previous log-volume and the slope baseline predict previous log-volume plus median training slope times elapsed days. Estimate baseline residual variance from training patients only. For the state-space model, update only with visits available before each prediction and include latent-state and parameter uncertainty in predictive densities.
13. For each held-out patient, average paired one-step log-score differences against the population-median-slope baseline (primary comparison) and persistence (secondary). Summarize patient-level differences with a median and a 2,000-resample patient bootstrap percentile 95% interval using seed 4201. Report two-sided intervals without assuming which model will score higher. Report predictive-interval coverage, calibration, exclusions, and every prespecified prior sensitivity.
14. Use LUMIERE only as a cross-cohort sensitivity analysis because its features were previously handled in this project and its tumor definitions and measurement pipeline differ. If fewer than 50 eligible MU patients remain, do not relax eligibility or claim an inferential five-fold comparison; report the eligible count and limit conclusions to descriptive analysis.

**Main limitation.** The cohorts do not contain a project-specific scan-rescan estimate of segmentation error. Measurement error may remain confounded with true net volume change. The external repeatability conversion is ambiguous and comes from a different cohort, scanner, and tool.

### 3. Track C — Treatment-policy simulation

**Question.** Under explicitly assumed drug-kill rates, how do adaptive and standard maximum-tolerated-dose schedules compare in a sensitive/resistant-cell model?

1. Use the two-population simulator with growth parameters estimated from the MRI cohort. Treat drug-kill rates as assumptions, not measured pharmacology.
2. Compare adaptive and maximum-tolerated-dose schedules at each prespecified E_MAX-to-growth ratio: 500, 1000, and 2000. Use the same initial conditions and simulation settings for both schedules.
3. Record modeled progression time, final tumor burden, resistant fraction, cumulative dose, and censoring for every simulated case and grid value. Report all settings together; do not select a favorable ratio.
4. Analyze MGMT-associated clinical outcomes as observational data only. Report association estimates with uncertainty and do not infer causal treatment benefit.
5. State that assumed kill rates, model structure, and fitted growth inputs limit interpretation. Do not present simulation output as evidence that a drug works in patients.

## Planned synthesis

1. Report each track in a separate row with dataset, population unit, endpoint, method, comparator, estimate, uncertainty, evidence type, and limitations.
2. State whether each prespecified comparison can be estimated, whether model diagnostics pass, and which eligibility or data constraints limit interpretation.
3. Keep the three tracks separate in the final interpretation. Do not claim that a result in one cohort validates another track or establishes a patient-specific digital twin.
4. Report estimates and uncertainty regardless of sign or statistical significance. Do not present a predicted result as a finding.

## Data inventory and use terms

| Source | Public scale and planned use | Access and use terms |
|---|---|---|
| UCSC multiomic-gbm atlas | Published paper reports 223,113 cells from 25 donors; atlas for Track A. Raw download size not recorded. [UCSC collection](https://cells.ucsc.edu/?ds=multiomic-gbm) has browser access and bulk rsync. | Atlas article is CC BY-NC-ND 4.0 and says generated data are publicly available without restricted access. Article license may not cover every exported file. Record exact file terms before reuse or redistribution. [Article](https://doi.org/10.1038/s41467-026-69716-2) |
| MU-Glioma-Post v1 | 203 patients; TCIA lists about 11 GB for images and segmentations; primary postoperative MRI cohort. | TCIA lists CC BY 4.0 for image, segmentation, and clinical files. Dataset DOI 10.7937/7K9K-3C83. [Collection](https://www.cancerimagingarchive.net/collection/mu-glioma-post/) |
| UCSF-PDGM | 495 subjects; TCIA lists about 142 GB at collection level and about 156 GB in image download table; selected separate MRI analysis only. | Public NIfTI and clinical files are listed CC BY 4.0. Record exact release and files because aggregate sizes differ. Dataset DOI 10.7937/tcia.bdgf-8v37; paper DOI 10.1148/ryai.220058. [Collection](https://www.cancerimagingarchive.net/collection/ucsf-pdgm/) |
| RHUH-GBM | 40 patients; brain-extracted NIfTI derivative about 2.9 GB; recurrence analysis. | TCIA marks collection limited overall because original DICOM is controlled. Collection table lists brain-extracted NIfTI images, segmentations, and clinical data under CC BY 4.0. Use only those derivative files. Dataset DOI 10.7937/4545-c905. [Collection](https://www.cancerimagingarchive.net/collection/rhuh-gbm/) |
| LUMIERE | 91 patients, 638 MRI studies, and 2,487 MRI images; Figshare full MRI package 30.33 GB; recurrence and sensitivity analyses. | Article says non-commercial use; Figshare lists CC0. Terms conflict. Verify exact files and permitted use before analysis or transfer. DOI 10.1038/s41597-022-01881-7. [Article](https://doi.org/10.1038/s41597-022-01881-7); [Figshare record](https://springernature.figshare.com/articles/dataset/LUMIERE_dataset_-_MRI_data_and_automated_segmentations/21249516) |
| PREDICT-GBM | Prior-art context only; no data input. | Dataset card lists MIT while paper notes source-data terms continue to apply. No download planned. [Card](https://huggingface.co/datasets/LZimmer/PREDICT-GBM); [paper](https://doi.org/10.1038/s41746-026-03194-0) |
| Abu Khalaf et al. repeatability study | Published study of 20 recurrent-GBM patients scanned about two days apart; external sensitivity reference only. | Published repeatability coefficients are not project images. Conversion to the model's log-volume-error scale is ambiguous. DOI 10.3174/ajnr.A7071. [Article](https://pubmed.ncbi.nlm.nih.gov/33737270/) |

Kaggle is not assumed to be permitted merely because a dataset is public. No restricted QIN-GBM data are included.

## Feasibility and limits

The state-space likelihood uses a Kalman filter over visit sequences, so model fitting scales with visit count rather than image voxel count. The 600 recovery fits, patient folds, profile grids, and image-model comparisons have not been benchmarked as a combined workload on the stated desktop. The full UCSF-PDGM release is large and is not required for the primary state-space analysis. Require at least 50 eligible MU patients for the planned five-fold comparison; this is an operational rule, not a power guarantee. If fewer qualify, retain eligibility rules and report descriptive analysis only.

The main scientific weakness is that scan-rescan error is not measured in the selected postoperative cohort. The state-space model may not separate measurement error from real net volume change. Track A has one atlas and a limited donor count. Track C uses assumed drug effects and cannot establish treatment efficacy.

## Full bibliography

Abu Khalaf, N., Desjardins, A., Vredenburgh, J. J., & Barboriak, D. P. (2021). Repeatability of automated image segmentation with BraTumIA in patients with recurrent glioblastoma. *AJNR American Journal of Neuroradiology, 42*(6), 1080–1086. https://doi.org/10.3174/ajnr.A7071

Calabrese, E., Villanueva-Meyer, J., Rudie, J. D., Rauschecker, A. M., Baid, U., Bakas, S., Cha, S., Mongan, J., & Hess, C. P. (2022). The University of California San Francisco Preoperative Diffuse Glioma MRI (UCSF-PDGM) dataset. *Radiology: Artificial Intelligence, 4*(6), e220058. https://doi.org/10.1148/ryai.220058

Cepeda, S., García-García, S., Arrese, I., Herrero, F., Escudero, T., Zamora, T., & Sarabia, R. (2023). The Río Hortega University Hospital Glioblastoma dataset: A comprehensive collection of preoperative, early postoperative and recurrence MRI scans (RHUH-GBM). *Data in Brief, 50*, 109617. https://doi.org/10.1016/j.dib.2023.109617

Hormuth, D. A., II, Al Feghali, K. A., Elliott, A. M., Yankeelov, T. E., & Chung, C. (2021). Image-based personalization of computational models for predicting response of high-grade glioma to chemoradiation. *Scientific Reports, 11*, Article 8520. https://doi.org/10.1038/s41598-021-87887-4

Mahmoud, E., Gass, J., Dhemesh, Y., Greaser, J., Pogorzelski, K., Isufi, E., Garrett, F., Thacker, J., Tahon, N. H. M., Sinclair, J., Layfield, L. J., Altes, T., & Nada, A. (2025). MU-Glioma Post: A comprehensive dataset of automated MR multi-sequence segmentation and clinical features. *Scientific Data, 12*, Article 1847. https://doi.org/10.1038/s41597-025-06011-7

Sonpatki, P., Park, H. J., Xing, Y. L., Han, K. Y., Schroeder, B. A., Yu, H. J., Kim, H. J., Chowdhury, T., Hwang, J. H., Nam, S. M., Byun, Y. H., Kang, H., Lee, J. H., Lee, S. T., Won, J. K., Kim, T. M., Choi, S. H., Ku, J.-L., Lee, S., … Shah, N. (2026). A spatially resolved human glioblastoma atlas reveals distinct cellular and molecular patterns of anatomical niches. *Nature Communications, 17*, Article 2951. https://doi.org/10.1038/s41467-026-69716-2

Suter, Y., Knecht, U., Valenzuela, W., Notter, M., Hewer, E., Schucht, P., Wiest, R., & Reyes, M. (2022). The LUMIERE dataset: Longitudinal glioblastoma MRI with expert RANO evaluation. *Scientific Data, 9*, Article 768. https://doi.org/10.1038/s41597-022-01881-7

Yaseen, D., Garrett, F., Gass, J., Greaser, J., Isufi, E., Layfield, L. J., Nada, A., Pogorzelski, K., Sinclair, J., Tahon, N. H. M., & Thacker, J. (2025). University of Missouri Post-operative Glioma Dataset (MU-Glioma-Post) (Version 1) [Data set]. The Cancer Imaging Archive. https://doi.org/10.7937/7K9K-3C83

Speir, M. L., Bhaduri, A., Markov, N. S., Moreno, P., Nowakowski, T. J., Papatheodorou, I., Pollen, A. A., Raney, B. J., Seninge, L., Kent, W. J., & Haeussler, M. (2021). UCSC Cell Browser: visualize your single-cell data. *Bioinformatics, 37*(23), 4578–4580. https://doi.org/10.1093/bioinformatics/btab503

Zimmer, L., Weidner, J., Balcerak, M., Kofler, F., Krupa, M., Ezhov, I., Cepeda, S., Zhang, R. Z., Lowengrub, J., Menze, B., & Wiestler, B. (2026). PREDICT-GBM: A multicenter platform advancing personalized glioblastoma radiotherapy planning. *NPJ Digital Medicine, 9*, Article 686. https://doi.org/10.1038/s41746-026-03194-0
