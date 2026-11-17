# Imaging and measurement additions

Research date: 2026-11-16. Desk research only; no patient image, mask, reserved outcome, model, or protocol was changed. Proposals below are not established novelty or scoring authorization. The first Firecrawl search/scrape round returned results; a later feedback batch stalled and was interrupted. Two web batches completed as fallback.

## Which addition most directly strengthens the existing forecast paper?

### Takeaway

Prioritize an empirical measurement-error and decision-stability study. Existing project code already propagates assumed spatial boundary error, so the useful addition is to test those assumptions against repeat acquisitions or independent expert labels and measure downstream instability.

### Cited Findings

- The official RIDER-NEURO-MRI record lists 19 recurrent GBM patients with repeat scans approximately two days apart, 17 with repeat DTI and contrast-enhanced FLAIR, and repeat contrast-enhanced 3D FLASH. The current image table points to NIH controlled access; it must not be advertised as unrestricted download. [TCIA RIDER](https://www.cancerimagingarchive.net/collection/rider-neuro-mri/).
- Shiri et al. already studied test-retest and registration effects on GBM MRI radiomic repeatability. Therefore repeatability alone is occupied prior art. [Primary paper](https://pubmed.ncbi.nlm.nih.gov/32615647/).
- Suter et al. already evaluated longitudinal automated measurements on LUMIERE, including response assessment, false positives, and lesion counts. Its two segmentation tools showed limitations in longitudinal consistency. This is direct prior art against generic automated measurement consistency claims. [Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10513769/).
- Local `src/fk_spatial_uncertainty.py:1` explicitly states that its displacement law is assumed or fitted and does not establish scanner-noise measurements. The earlier `reports/Forecastability Bound for Glioma MRI.md` already proposes a narrow observation-aware certificate and records broad prior art. Those ideas cannot be relabeled as a new Track D discovery.

### Inferences

**Candidate D1: Does acquisition and segmentation error change which growth model wins?** Proposed hypothesis: model rankings and predicted growth direction become unstable when observed changes approach empirically estimated repeatability limits. Estimate log-volume paired differences, surface displacement correlation, and growth-sign or progression-threshold flips; propagate those errors through a completely frozen forecaster. Compare with no uncertainty, simple isotropic perturbation, and the project's existing correlated-boundary method. Patient is the sampling unit. Separate acquisition, registration, and segmentation effects where design permits.

Proposed failure rule: if the empirical error model does not improve held-out proper scores/calibration versus a simple error model, or if model rankings remain stable, report that result. Do not tune a noise threshold to create a preferred winner. A 1–2 day interval is not proof of zero biological change; 19 patients cannot support a universal scanner-noise floor. Independent annotators are preferable to treating network agreement as ground truth.

Novelty would be the downstream forecast-comparison and treatment-decision consequence, not uncertainty propagation itself. Best as one substantive addition to Track B, with a potential Track C sensitivity connection. Estimated work: several weeks after data/access and annotation gates, potentially longer with sponsor-dependent access.

### Gaps

- RIDER access has changed from assumptions in older literature; eligible sponsor/access process and usable sequence pairing still need a metadata-only audit.
- No freely reusable multi-rater longitudinal GBM volume dataset was verified in this short search. Repeat segmentation without repeat acquisition measures a different error component.
- Existing MU/UCSD/UCSF/LUMIERE results remain development evidence; CFB/WAW outcomes stay reserved under their locked protocol. A new measurement study does not reset spent forecast tests.

## Which complementary imaging tracks are feasible and distinct?

### Takeaway

Two defensible candidates are a narrow forecast-abstention benchmark and a postoperative-to-recurrence segmentation transport benchmark. Both face substantial prior art; neither justifies another large neural architecture by itself.

### Cited Findings

- CONSeg already applies conformal uncertainty to glioma segmentation, with internal and external calibration/test datasets, and investigates uncertain-case triage. Generic conformal glioma segmentation is not new. Voxel coverage is not patient-level trajectory coverage. [Primary paper](https://pubmed.ncbi.nlm.nih.gov/40610235/).
- A 2026 primary preprint evaluates uncertainty, corruption robustness, failure detection, and selective prediction in brain tumor segmentation. The narrower temporal forecast task must be distinguished from these existing evaluations. [Preprint, not peer-reviewed evidence](https://arxiv.org/abs/2607.22727).
- The official RHUH-GBM collection describes MRI at preoperative, early postoperative, and recurrence timepoints, with tumor segmentations. The official page was retrieved by Firecrawl, but fallback web retrieval timed out. Exact label definitions, date availability, license, complete cases and patient overlap therefore need rechecking before selection. [Official TCIA record](https://www.cancerimagingarchive.net/collection/rhuh-gbm/).
- A 2026 study already uses RHUH-GBM and LUMIERE for preoperative/postoperative recurrence risk. A generic recurrence-prediction claim on these datasets is occupied. [Primary study](https://www.mdpi.com/2306-5354/13/6/668).

### Inferences

**Candidate D2: Can a frozen forecaster identify when it should abstain?** Proposed endpoint: patient-level risk–coverage curve and retained-patient error versus persistence at prespecified coverage; secondary interval coverage and width by site and scan interval. Freeze a simple baseline score (recent observed change, tumor size, scan interval) and compare a pre-outcome uncertainty score with random rejection. Abstention thresholds must come from development/calibration patients, never evaluation outcomes. Failure: no incremental benefit over simple scores at fixed coverage, or coverage collapses externally. Marginal conformal coverage requires appropriate exchangeability and does not automatically hold after site shift. Novelty narrow because the project already has forecastability/trust work; an independently evaluated temporal endpoint would add evidence, not establish a novel general method. Effort: roughly 2–4 weeks if predictions and labels already exist, longer for fresh validation.

**Candidate D3: Does postoperative domain shift propagate into false growth?** Use one frozen segmentation model across preoperative, early postoperative and recurrence images; assess change error and growth-direction flips versus expert masks. Compare independent per-scan segmentation, a simple longitudinal consistency rule, and an ensemble. Measure errors separately by timepoint and subregion; cavity and blood products should be explicit error categories. RHUH is a candidate because those distinct timepoints are described, but it is not a general next-visit forecasting cohort: surgery fundamentally alters anatomy and the final visit is recurrence-selected. Failure: no downstream benefit beyond the basic consistency rule, or insufficient expert labels/time metadata. Do not claim recurrence biology from segmentation discrepancy. Estimated work: 3–6 weeks after eligibility and rights audit. Most useful if user wants an imaging-specific side paper; lower priority than validating the existing primary forecast.

### Gaps

- RHUH patient count and exact released labels were not independently revalidated against official download metadata; do not rely on third-party mirrors for eligibility.
- A confirmation cohort for D2 remains unresolved. The locked CFB/WAW study must finish as specified or receive a separate dated, prospective amendment without outcome-dependent changes.
- The focused searches do not establish global novelty for D2 or D3.

## Could treatment effects and scan timing form another strong track?

### Takeaway

Treatment-aware interpretation is scientifically important, but missing dose, timing, steroid and confirmed progression labels can make this less feasible than it sounds. Prioritize an explicit confounding audit before building a treatment-response classifier or scan-scheduling policy.

### Cited Findings

- RANO 2.0 uses post-radiotherapy baseline imaging for newly diagnosed cases and generally requires confirmation of apparent progression in the first 12 weeks after radiation. The criteria distinguish radiographic changes from definitive progression. [Primary consensus](https://ascopubs.org/doi/10.1200/JCO.23.01059).
- Treatment-aware diffusion modeling for longitudinal glioma MRI generation and growth prediction already exists; merely adding treatment indicators is not a new general approach. [Primary preprint](https://arxiv.org/abs/2309.05406).
- A primary study examines how scheduled and unscheduled GBM follow-up MRI influence treatment decisions and uncertainty. This establishes prior art for clinical scan utility and highlights why visit timing cannot simply be interpreted as random. [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10605145/).

### Inferences

**Candidate D4: Does a model forecast biological growth or the observed treated imaging trajectory?** Before scoring, classify intervals using reliably available RT dates, surgery, treatment changes and confirmation status. Compare a simple treatment-stratified baseline with the frozen growth model. Primary target should remain observed image-label change unless independent clinical confirmation exists; MRI alone cannot provide a true-resistance or pseudoprogression label. Report performance and calibration by treatment phase and scheduled/unscheduled visit status when those fields are present. Freeze missingness handling and avoid making inferred treatment dates from future images. Failure: eligibility audit cannot link treatment dates/status for a prespecified adequate patient count, or apparent subgroup advantage disappears relative to a simple stratified baseline. Effort: 2–4 weeks for a metadata audit plus an uncertain clinical-label acquisition period.

This can strengthen interpretation of the negative Track B results even if no new predictive model wins. It is a weaker standalone new track without reliable treatment metadata. Optimizing a scan schedule from retrospective visit frequency does not establish causal clinical benefit and would duplicate parts of Track C's value-of-information work.

### Gaps

- No new untouched open cohort with linked drug dose, RT dates, steroids, scheduled-versus-triggered scans and confirmed progression was established here.
- ACRIN advanced-MRI feasibility already failed locally; it should not be repackaged as ready to analyze.
- These proposals require new endpoint-specific protocols. No existing frozen rule is relaxed.
