# Track B protocol: locked external test of next-visit tumor-core annotation changes

**Status:** prospective, protocol-locked external patient analysis. The UCSD image and segmentation arrays have not been scored under this protocol. Metadata-only cohort counts and source documentation were inspected to establish feasibility. This is not an externally preregistered clinical study or clinical validation. Freeze this file and all locked scoring code in Git before reading any selected UCSD label arrays.

## Research question

Among patients with glioblastoma and at least two eligible post-treatment MRI visits, does a screened graph-diffusion spatial prior plus an MRI residual attenuated by fixed boundary-perturbation sensitivity improve patient-macro ROC AUC for ranking additions to the next source-provided standard tumor-core annotation, compared with the strongest prespecified non-full model selected using development-only patient-grouped cross-validation?

The endpoint is **change in supplied segmentation labels**. It is not tumor-cell growth, radiographic progression, recurrence, or response to treatment. The evaluation is a project-level locked external test, not clinical validation and not a claim that the method is first or unique.

## Sources and roles

- **Development only:** MU-Glioma-Post v1, with the project’s existing patient and interval manifests. This cohort has been explored in earlier project work. It may train model coefficients and select the primary comparator through patient-grouped development cross-validation; its metrics are not confirmation.
- **External test:** UCSD-PTGBM v3, from the official [TCIA collection page](https://www.cancerimagingarchive.net/collection/ucsd-ptgbm/). TCIA lists CC BY 4.0 for the images, segmentations, and clinical workbook, with 136 subjects in the Images and Segmentations package and 42 subjects in the BraTS-GLI 2024 Test Data package. The two image packages total about 44.98 GB and use IBM Aspera for download. Use the current v3 clinical workbook because TCIA says earlier versions contained incorrectly assigned patient IDs. The locally reviewed workbook is pinned by SHA-256 `bbf40bbdd289266924869052806de4d7e8fc9ee73286cdc0cc066672ef5bb51d`. Cite dataset DOI [10.7937/fwv2-dt74](https://doi.org/10.7937/fwv2-dt74).
- The dataset description reports 178 subjects and 243 timepoints, standard multimodal NIfTI, registration/resampling to a common 1 mm MNI grid, and a FLAIR-derived brain mask applied to the other modalities. Standard compartment segmentations are initialized by software and manually corrected by trained radiologists, according to TCIA; they remain source-provided annotations, not error-free biological truth. TCIA’s page and dataset article are the sources: [TCIA](https://www.cancerimagingarchive.net/collection/ucsd-ptgbm/) and [Gagnon et al.](https://doi.org/10.1038/s41597-025-06499-z).
- The BraTS post-treatment annotation vocabulary distinguishes non-enhancing tumor core (NETC), surrounding non-enhancing FLAIR hyperintensity (SNFH), enhancing tumor (ET), and resection cavity (RC); see the [BraTS 2024 post-treatment challenge description](https://escholarship.org/uc/item/9961p1k9). The analysis uses only the standard BraTS compartment segmentation. It does not use cellular-tumor masks whose definitions use information across multiple timepoints.

### Metadata-only feasibility count

The current v3 clinical workbook has 243 visit rows for 178 unique subjects. Exact `Primary Diagnosis = Glioblastoma` gives 172 unique subjects; 48 have at least two visits. There are 63 adjacent visit intervals among those 48. The metadata field measuring days from acquisition to the initial surgery/treatment/diagnosis is present for 51 intervals; their median absolute gap is 113 days (range 17–838), and only five lie between 76 and 104 days. The numeric visit suffix agrees with later acquisition order in all 36 repeated-patient groups with complete anchor-day values; 12 groups lack complete anchor-day values. These are source-metadata counts, not image-eligibility or outcome results. The actual number eligible after file, modality, and geometry checks will be reported by the preflight before any test scoring.

## Frozen cohort and pair selection

1. Read only the current v3 clinical workbook, pinned by the checksum above. Include rows with exact primary diagnosis `Glioblastoma`; reject duplicate patient/timepoint rows; group by TCIA subject ID and sort visits by the numeric visit suffix. Where both anchor-day values are present, require the later suffix to correspond to a later acquisition (a smaller acquisition-to-anchor value). Stop on a contradiction; do not repair chronology by looking at tumor labels.
2. For each subject, select the **first adjacent pair** for which both visits contain exactly one standard segmentation and exactly one each of T1-native, T1-contrast, T2, and FLAIR, each visit is internally on one grid, and the two visits share shape and affine within 0.001 mm. Continue to the next adjacent pair only when an earlier pair lacks required files or fails this source-grid check. Match complete patient/timepoint and BraTS IDs so suffixes such as visit 1 and visit 10 cannot collide. Never select a pair based on the amount, direction, or apparent quality of the target change.
3. Use at most one pair per subject. A subject with no complete aligned pair remains in the flow table as excluded with a prespecified technical reason. Report counts for all diagnosis-eligible, repeat-eligible, file-complete, aligned, scored, and metric-evaluable subjects. Do not inspect images visually or tune the method on UCSD before the lock.
4. The forecast origin is the earlier scan. The target is the next available scan, regardless of interval length. The interval is reported descriptively in fixed strata (≤90, 91–180, and >180 days) and is not supplied to any model. Do not describe this as a fixed-day forecast.

## Endpoint and mathematical method

Use the standard BraTS compartments: (C_t(v)=1) when voxel (v) is NETC or ET (source label values 1 or 3), and (C_t(v)=0) otherwise. SNFH (label 2) is a predictor compartment, not the target. RC (label 4) is excluded from the addition risk set and blocked in the graph.

Let (B_0) be the current scan’s source skull-stripped support: voxels nonzero in at least one of the four standard modalities, together with current nonzero segmentation voxels. Let (Q_0) be current label-4 resection cavity. Define the forecast-origin risk set

\[
R_0 = B_0 \setminus (C_0 \cup Q_0), \qquad
Y(v)=\mathbf 1\{v\in C_1\setminus C_0\},\quad v\in R_0.
\]

The risk set covers the full current brain support; it has no 20 mm cutoff. Report newly annotated core voxels outside (R_0), including those outside current support and those within current cavity, so missed locations remain visible.

For each origin, robust-normalize its four standard MRI modalities over (B_0). Require the NIfTI spatial units to be millimetres and the source grid to be orthogonal and isotropic at 1.000 ± 0.001 mm per axis for both MU development and UCSD test images; this makes the fixed one-voxel perturbation comparable across subjects. On graph nodes (B_0\setminus Q_0), solve the screened graph-diffusion equation

\[
 (L_A + \ell^{-2}I)u = \mathbf 1_{C_0},\qquad \ell=10\text{ mm},
\]

where six-neighbor edge conductance is

\[
A_{vw}=h_{vw}^{-2}\exp\!\left[-0.5\,\frac14\sum_{m}(z_m(v)-z_m(w))^2\right]
\exp\!\left[\frac{\log 2}{2}\{\mathbf1_{SNFH}(v)+\mathbf1_{SNFH}(w)\}\right].
\]

This is a spatial ranking prior, not a cell-density or biological growth estimate. The MRI residual uses seven origin-only features: distance to current core, distance to current cavity, four robust-normalized MRI intensities, and current SNFH indicator. There is no prior-transition feature, fitted patient-specific PDE parameter, future-volume/count constraint, or realized-interval predictor.

Create two fixed current-mask sensitivity scenarios by one-voxel erosion and dilation, recompute the log prior λu for each, and define

\[
I(v)=\operatorname{median}_{q\in\{\mathrm{erode},\mathrm{dilate}\}}
|\log u_q(v)-\log u_0(v)|,\quad
s=\operatorname{median}_{i\in\mathrm{MU}} Q_{0.90}\{I_i(v):v\in R_i\},\quad
g(v)=\exp[-I(v)/s].
\]

The gate is a deterministic sensitivity heuristic, not calibrated uncertainty. The full model’s voxel score is

\[
S(v)=\beta\log u_0(v)+g(v)\,[\theta_0+x(v)^T\theta],\qquad \beta\ge0.
\]

Fit a ridge-regularized logistic ranking head to the MU development pairs, with equal patient/class weighting, ridge (1.0), prior penalty (10^{-3}), at most 4,000 sampled positive and 4,000 sampled negative voxels per patient, and base seed 20261026. Development-fold fits use base seed plus the one-based fold index; the final MU fit uses the base seed. All spatial parameters are fixed as written above; UCSD is not used to select them, train coefficients, calibrate the gate, or select a comparator.

## Comparators and analysis

Before the external test, use five-fold `GroupKFold` on unique MU development patients to calculate mean patient-specific AUC for these frozen alternatives: ungated hybrid, MRI-only ridge head, diffusion-only prior, negative distance-to-core, and uniform no-change score. Select the one with the highest development mean as the single primary comparator; ties use that listed order. This selection is locked before UCSD scoring. The development cohort is already explored and is used only to train and select.

For each UCSD patient with both at least one positive and one negative voxel in (R_0), calculate its ROC AUC. The primary estimand is the mean paired patient-level difference

\[
\Delta=\frac1{n_E}\sum_{i=1}^{n_E}
\left[\mathrm{AUC}_i(\mathrm{full})-\mathrm{AUC}_i(\mathrm{selected\ comparator})\right].
\]

Use a 10,000-draw paired bootstrap that resamples patients, seed 20261126, to give a percentile 95% confidence interval. Voxels are not inferential units. Patients with no target additions stay in the cohort flow and are reported, but have undefined AUC; do not assign them a fabricated value. Report every arm’s patient-macro AUC, the number evaluable, the event-positive count, observed scan-gap distribution, total additions outside the risk set, and an interval-stratified descriptive table. No extra significance tests are run for secondary arms.

The prespecified added-value criterion is that the paired 95% confidence interval lower bound for Δ exceed 0.02 AUC. This tests evidence beyond the project-set 0.02 margin; the margin is not a clinically validated minimum important difference. If not met, report no demonstrated added value or an imprecise/inconclusive result. A positive result supports only improved ranking of source-annotation additions in this cohort; it does not establish clinical utility, tumor truth, fixed-horizon performance, or superiority to all published systems. The metadata maximum of 48 repeat-eligible GBM patients is not a power calculation; no claim of adequate power is made before image- and endpoint-eligibility are known.

## Source and analysis gates

- Use both official TCIA image packages as needed. The preflight must report package/file completeness and common-grid eligibility before score mode. The image packages have not been downloaded into this repository at protocol freeze.
- Download only to the ignored `data/external/UCSD-PTGBM-v3/` directory. Do not commit images, clinical workbook rows, subject IDs, masks, or voxel predictions. The locked runner verifies that the protocol and scoring code are committed and clean, checks the pinned clinical-workbook checksum, and records a one-shot marker plus the selected file path/size manifest before reading any UCSD segmentation array. It then hashes the selected MRI/segmentation file contents and records the digest in that marker before scoring; it verifies the files did not change during scoring. A failed run keeps the marker and must be investigated before any amended rerun. Aggregate-only results go to ignored `output/track_b_ucsd_locked_test/aggregate_result.json`.
- Require exact source release v3 and standard compartment segmentation. Do not substitute the cellular-tumor segmentations that use future timepoints, and do not silently resample a misaligned pair.
- Any missing modality, ambiguous segmentation filename, unexpected label code, inconsistent visit order, registration mismatch, failed solver, insufficient MU training cases, or empty external cohort is a stop condition to resolve before scoring. Do not edit this protocol after looking at UCSD test outcomes; an amendment would make subsequent analysis exploratory.
- The current test data include public BraTS-GLI challenge material, and other researchers may have used the same collection. The method is not claimed to be novel or to beat the state of the art. Closely related future-mask work includes [TaDiff](https://arxiv.org/abs/2309.05406) and mechanistic/diffusion MRI-generation work ([arXiv:2509.09610](https://arxiv.org/abs/2509.09610)); this test evaluates a narrow cross-cohort annotation-transition question.
- Weakest point before analysis: only 48 GBM patients have repeated timepoints in the clinical metadata, interval lengths vary greatly, and some future annotations may reflect therapy, registration, or contouring rather than biological change. The final image-complete and AUC-evaluable sample may be smaller, so the test may be imprecise or inconclusive. Public availability and inability to rule out overlap with other public cohorts mean this is not clinical confirmation.

## Required data attribution

Gagnon, L., Gupta, D., Mastorakos, G., et al. *The University of California San Diego Post-Treatment Glioblastoma (UCSD-PTGBM) annotated multimodal MRI dataset*. Dataset version 3, The Cancer Imaging Archive. https://doi.org/10.7937/fwv2-dt74.
