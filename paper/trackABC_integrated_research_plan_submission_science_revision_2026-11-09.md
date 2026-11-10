# Integrating Tumor Biology and Longitudinal MRI for Glioblastoma Progression Modeling and In-Silico Therapy Optimization

**Integrated research plan and procedure — science-corrected submission draft**<br>
**Project period:** 2026–2027<br>
**Project type:** computational continuation project

> This plan separates Track A molecular analyses, Track B imaging analyses, Track C clinical-data and simulation analyses, and the synthetic integration. Existing procedures are identified as completed; proposed work is identified as planned or conditional. This document contains no outcome values. The tracks do not have a verified patient-level link. The integrated therapy-policy study therefore uses simulated signals and does not establish patient-specific treatment benefit. This plan is not a data-use agreement, a protocol seal, or authorization to open restricted files or run a new patient analysis.

## Rationale

Glioblastoma (GBM) varies across tumor regions and changes over time. Single-cell expression can describe molecular patterns in sampled tissue, while serial magnetic resonance imaging (MRI) can describe radiographic changes in tumor-associated tissue. Mathematical models can represent tumor growth, treatment exposure, and uncertainty. These sources can inform related questions, but they measure different things: an expression score is not a drug-response assay, an MRI transition is not a direct resistance measurement, and a simulated treatment outcome is not a patient outcome.

Prior research already combines parts of this problem. Studies have used imaging to estimate infiltration and recurrence, patient-specific mathematical models to study radiotherapy, longitudinal MRI models, and reinforcement learning or other optimization methods for treatment schedules (Akbari et al., 2016; Chaudhuri et al., 2023; Ebrahimi Zade et al., 2020; Liu et al., 2025; Zimmer et al., 2026). A recent single-cell framework, scFOCAL, combines drug-response signatures with GBM cell-state expression to predict sensitive and resistant cell populations and includes experimental validation in treated models; a recent spatial atlas also maps GBM cell states across regions (de Jong et al., 2026; Suter et al., 2026). Ex-vivo drug-sensitivity assays measure response more directly than expression or MRI alone (Ntafoulis et al., 2023). These precedents mean that a broad claim of a first multimodal GBM therapy framework would not be justified.

The project will instead examine three connected but separate evidence streams. Track A will describe regional expression patterns and evaluate donor-held-out prediction. Track B will describe MRI-defined spatial transitions and assess forecast methods only where the data and permissions allow. Track C will analyze an observational MGMT–survival association and develop a mathematical treatment-policy simulation. The integrated question will be tested only in synthetic worlds: whether a controller that models dependence between two imperfect signals has different simulated utility from a matched controller that assumes conditional independence. The molecular-like and MRI-like signals in that experiment will be generated from explicit assumptions; they will not be claimed to be measurements from the same patients.

The project’s contribution will depend on transparent data separation, patient- or donor-level validation, shared baselines, reproducible implementation, and a valid uncertainty analysis. The study will report negative, null, and inconclusive findings under the same rules as positive findings. It will not claim clinical benefit, causal treatment effects, validated patient-specific resistance, or a clinically validated digital twin.

## Research Questions, Hypotheses, Engineering Goals, and Expected Outcomes

### Track A — tumor-region expression

**Research question A1.** Under a frozen six-donor, nested leave-one-donor-out design, will the sparse donor-balanced classifier have lower equal-donor class-balanced area under the risk-coverage curve (AURC) than ordinary L2 logistic regression, with lower AURC in at least five of the six held-out donors?

**Hypothesis A1.** The candidate classifier is hypothesized to meet both parts of that prespecified rule. The rule will not be changed in response to the observed result. The analysis is descriptive for this cohort and does not establish malignant-cell-specific prediction or drug resistance.

**Research question A2.** In paired-donor expression analyses, how will the prespecified SenMayo, OXPHOS, and tumor-associated macrophage programs vary by sampled region after expression-matched controls and the declared covariate analyses?

**Research question A3.** How will a fixed, direction-preserving expression score based on genes associated with in-vitro temozolomide (TMZ) viability align with the sampled Core and Periphery regions? This is an exploratory expression proxy, not a test of response in the patients whose tumors were sequenced.

### Track B — longitudinal MRI and spatial transitions

**Research question B1.** Within a patient and a prespecified physical-distance band, how will the next-visit rate of MRI-label transition differ between history-defined newly appearing FLAIR tissue and persistent FLAIR tissue?

**Exploratory hypothesis B1.** Newly appearing FLAIR-labeled tissue is hypothesized to have a higher subsequent core-label transition rate than persistent FLAIR-labeled tissue at the same distance from the prior core. Any analysis on a cohort whose target labels have already been examined will remain exploratory.

**Research question B2.** In a future eligible, unexposed longitudinal cohort, will a prespecified PDE-based forecast improve patient-level forecast performance over persistence, fixed-growth, PDE-only, and learned-image baselines when all methods use the same patients, image support, target labels, and forecast intervals?

An MRI-defined transition will be called TMZ-associated only if patient-level TMZ exposure and scan chronology are verified. It will not be called direct drug resistance. Confirmatory forecasting will require a separate cohort that passes the access, label, registration, support, coverage, eligibility, and prior-use checks before target labels are opened.

### Track C — clinical association and mathematical treatment simulation

**Research question C1.** In eligible CGGA glioblastoma records, how will MGMT methylation be associated with the database’s recorded overall-survival endpoint after adjustment for age, IDH status, radiotherapy, TMZ status, recurrence status, and batch? Will the association differ by recorded TMZ status?

**Hypothesis C1.** MGMT methylation is hypothesized to be associated with recorded overall survival, and its association may differ by recorded TMZ status. This retrospective analysis will be interpreted as an association, not a randomized treatment effect. The survival-time origin and clinical chronology will be documented before any clinical timing interpretation.

**Research question C2.** Across independent simulated GBM-like worlds, will a dependence-aware minimax policy have greater mean simulated utility than an otherwise matched minimax policy that fixes residual signal dependence at zero? How will the contrast vary with assumed signal quality and dependence?

**Hypothesis C2.** The dependence-aware policy is hypothesized to have greater mean simulated utility when the simulated signals are informative and dependent. The no-information and zero-dependence controls are expected to show little or no advantage. These are simulation hypotheses; they do not predict clinical outcomes.

### Integrated engineering goal

The project will build and evaluate a reproducible pipeline that keeps the three evidence streams separate, makes the synthetic bridge explicit, compares policies on common simulated inputs, and refuses confirmatory evaluation unless the sampling, code, uncertainty, precision, and runtime requirements are met. A positive result is possible but not guaranteed. If the uncertainty method or precision target cannot be met within the prespecified computational budget, the project will report a feasibility result and will not run or describe the policy comparison as confirmatory.

**Expected outputs.** The project will produce donor-level expression summaries and prediction metrics, patient-level MRI transition and forecast metrics when a permitted cohort is available, an observational CGGA association analysis, and a reproducible synthetic policy benchmark only if its inference gates pass. Results may support, fail to support, or be inconclusive for each hypothesis. The project will preserve the original rules and report the outcome actually produced.

## Materials

### Computing materials

- Windows 11 personal computer with 32 GB RAM and CPU-based computation.
- Python environment recorded per analysis, including Python, NumPy, SciPy, pandas, scikit-learn, PyTorch, and relevant image or single-cell packages. Exact versions, hardware, thread settings, code hashes, and protocol hashes will be written to each run manifest; an environment will not be described as locked unless the manifest verifies it.
- Project source code, tests, configuration files, and protocol documents stored in the project repository. Tests and development simulations will run in a separate scratch copy as required by the project workflow.
- No laboratory chemicals, biological agents, animals, or physical treatment devices will be used.

### Data sources and roles

- **Track A:** Yu et al. (2020) single-cell RNA-sequencing counts and metadata, [GEO GSE117891](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117891), with the paper’s Table S1 point labels; [GBmap](https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c) cell-type annotation metadata; Ivy GAP bulk expression ([GSE107559](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107559)); Angers microarray ([E-MTAB-5092](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5092)); and [E-MTAB-8743](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8743) for the separately specified macrophage-program analysis.
- **Track A response-associated references:** the published Ntafoulis et al. (2023) in-vitro TMZ-viability gene panel and its original gene directions; the GBM-Space data and author-defined state annotations for a separate exploratory state-score analysis. These references do not provide patient-specific clinical drug-response outcomes for the Track A samples.
- **Track B:** the project’s previously examined [MU-Glioma-Post](https://www.cancerimagingarchive.net/collection/mu-glioma-post/) and [Duke longitudinal MRI record](https://zenodo.org/records/20548528). The historical Duke analysis is retained as exploratory; because its written data-use terms are not filed in the project, no additional Duke processing or derived-output sharing is planned until those terms are verified and saved. The Duke cohort has already been exposed to target-label analysis and cannot serve as untouched confirmation. A future independent cohort will be used only if it passes all prespecified eligibility and data-use gates.
- **Track C observational:** the Chinese Glioma Genome Atlas (CGGA) clinical and molecular records, subject to the source’s access and use conditions.
- **Track C simulation:** code-generated virtual worlds, parameter schedules, and random-number seeds. No patient image or expression value will be treated as a simulated policy outcome or silently joined to a virtual patient.
- **Literature:** peer-reviewed and primary sources listed in the bibliography.

## Procedures

The procedure distinguishes methods already completed in the continuation project from methods that remain planned or conditional. No numerical outcomes are included here. Existing results and their limitations will be retained in separate analysis records and will not be used to alter the original frozen rules.

### 1. Shared data governance and reproducibility

1. Before processing a dataset, record its source, version, access route, license or data-use terms, allowed derived outputs, and any school or competition review conditions. Keep restricted data and patient-level outputs out of public Git unless the written terms explicitly allow sharing.
2. Create a provenance record for each analysis containing input-file hashes, source-code and protocol hashes, software versions, random seeds, exclusions, folds, and the output path. Hashing a patient file requires reading that file and will occur only within the approved data-use scope.
3. Keep Track A donors, Track B patients, CGGA records, and simulated worlds in separate identifiers and tables. No cross-track patient join will be made because a verified common-patient key and matched drug-response reference do not exist.
4. Use donor as the unit for Track A inference, patient as the unit for Track B inference, and independent simulated world as the unit for Track C policy inference. Cells, voxels, scan windows, and particles will not be treated as independent biological subjects.
5. Run code tests and synthetic mechanics checks before any new target-label analysis. Synthetic checks will be reported as software checks, not evidence of patient validity or model performance.

### 2. Track A — regional expression and molecular proxies

#### A5 donor-held-out classifier procedure already completed

1. Use the frozen Yu et al. (2020) cohort and exact cell, sampling-point, donor, and Table S1 joins. The six protocol-designated donors are GS1, GS2, GS3, GS5, GS6, and GS13. Table S1 labels define Tumoral points as Core and Peritumoral points as Periphery; Normal points and invalid exact joins are excluded.
2. Retain cells with at least 1,000 UMIs and at least 500 detected genes. In each training split, retain genes detected in at least 5% of training cells, normalize each cell to 10,000 counts, apply log1p, select the 2,000 most variable genes from the transformed training cells, and fit scaling means and standard deviations on training cells only. Apply the fitted transformations to validation and held-out cells.
3. Use nested leave-one-donor-out validation. Each outer fold holds out one donor; inner folds use only the remaining donors for model selection. Compare the frozen sparse donor-balanced candidate with L2 logistic regression, 15-nearest-neighbor voting, linear SVM, and gradient-boosted trees. Retain the three prespecified candidate ablations and cell-type-only and cells-per-point diagnostics. Label constant-score controls as degenerate controls.
4. Use equal-donor class-balanced AURC as the primary endpoint. The frozen success rule requires both lower mean AURC than L2 logistic regression and lower donor-specific AURC in at least five of six donors. Balanced accuracy and Brier score are secondary. The A5 pipeline and success rule will not be retuned or reinterpreted as a drug-response test.

#### Regional gene-program and response-proxy procedures

1. Keep the paired-donor SenMayo and OXPHOS analyses in Ivy GAP and Angers separate from the TAM-activation analysis in E-MTAB-8743. Summarize at donor level, use expression-matched gene-set controls and the prespecified exact donor sign-flip procedure, and report the declared covariate analyses separately. Do not pool these results with A5 or treat cells as independent donors.
2. For the Ntafoulis in-vitro TMZ-viability panel, use the original 121-gene source panel and original directions of association with cell viability measured after 100 μM TMZ. The Yu expression matrix has exact symbols for 119 panel genes; omit the two unmapped genes without aliases or substitutions. Create donor-by-region pseudobulk summaries, calculate the fixed signed score, and compare with 1,000 expression-matched control panels. Interpret it only as relative expression alignment with published in-vitro viability directions, not as measured sensitivity or resistance in the sequenced patients.
3. In a separate exploratory GBM-Space analysis, use only the 115 exact mapped genes from the fixed in-vitro viability panel (47 with positive and 68 with negative source directions). Compare pseudobulk scores between the authors’ malignant developmental-like and gliosis/hypoxia state groups, using the donor as the unit and the prespecified minimum of 200 nuclei per group. Use raw-count pseudobulk followed by counts-per-million and log2(CPM+1), percentile scoring, and expression-matched controls. This tests association of a published in-vitro viability panel with author-defined states; it does not test actual TMZ response or establish new spatial biology.
4. Treat calibration, within-library, scANVI, and attention analyses as deferred. They will not be run unless a dated amendment fixes their question, split, endpoint, comparison, and success rule before analysis.
5. Any new external expression cohort must pass all gates before expression values are read: at least six eligible independent donors; verified IDH-wildtype status at the donor level; exact sample-to-region mapping; documented malignant-cell labels and validation; usable feature overlap; clear dataset-specific access and derived-output terms; and no overlap with the project’s training or prior-use records. If no cohort passes, no external-confirmation claim will be made.

### 3. Track B — longitudinal MRI and spatial transitions

#### Existing exploratory procedures

1. Preserve the completed MU-Glioma-Post PDE analysis and the completed Duke transition analysis as exploratory records. The Duke cohort has already had target masks read and its main endpoint scored; it is not an untouched test set. Do not reuse it to claim independent confirmation or to tune a new threshold.
2. Record the exact executed Duke estimand separately from the 2-mm transition-estimand prototype. The historical analysis used a 10-mm history dilation, a 10-mm exclusion around the middle-visit core, and three physical-distance bands. The 2-mm prototype has different region definitions and is not a reproduction of that analysis.
3. Preserve the historical limitation that out-of-field target voxels were resampled as zero without a coverage mask. This can classify unobserved tissue as non-conversion. Do not present the old endpoint as coverage-corrected.
4. Call outcomes MRI-defined transitions or recurrence. Use “TMZ-associated” only when treatment exposure and scan chronology are verified; never call radiographic transitions direct drug resistance. Keep the ADC analysis separate unless its source definitions, read scope, preprocessing, quality checks, and missingness rules are independently frozen.

#### Conditional future forecast study

1. Search for and screen a separate longitudinal cohort using metadata and published documents first. Before reading target labels, verify access and derived-output rights, visit mapping, label meanings, scan timing, registration evidence, valid brain support, target-field coverage, prior-use overlap, and eligible patient count. If any gate fails, stop that cohort’s confirmatory analysis.
2. Freeze a patient-level split and a primary horizon before target labels are opened. Use outer patient-grouped folds and inner patient-grouped selection, so no patient appears in both sides of an evaluation fold. If a separate, adequately sized external cohort later becomes available, reserve it for one final untouched test. Use history-only inputs for forecasting and fit all normalization, model selection, calibration, and baseline choices inside training patients. Do not use held-out patients to choose thresholds, radii, support rules, or model settings.
3. Compare on the same patients, image support, target labels, and actual forecast intervals: copy-last; time-aware fixed growth, with rate k selected from 0.00–0.19 mm per nominal day in 0.01 steps using training patients only; PDE-only, with its concentration cutoff selected from 0.05–0.50 using training patients only and ties resolved toward the larger cutoff; the hybrid PDE model; and a time-conditioned FiLM U-Net. Keep the full-grid time-blind U-Net and gap-blind fixed-growth arm as secondary reference baselines. The proposed FiLM model will use a four-level full-grid 3D U-Net backbone with widths 16/32/64/128, two 3×3×3 convolutions per level, instance normalization, LeakyReLU, three binary mask inputs, and interval-conditioned scale and shift parameters at its normalization layers. The proposed training recipe is Dice plus binary cross-entropy, Adam at 0.001 with zero weight decay, batch size 1, 300 epochs with the final epoch retained, three seeds per fold (20261102 + fold index, then +1000 and +2000), and history-only flips and axial 90-degree rotations. Its full recipe and fingerprint, including both history intervals, must be adopted and sealed before scoring. The count model will define a new piece as a 26-connected component outside the visit-2 core grown by 2 mm, retaining components of at least two voxels; it will fit a Poisson model with log visit-2 volume and interval offset and will be evaluated only with count deviance, count MAE, and the Brier score for any new piece. It will not be converted into a Dice mask. A two-history geometric control will remain descriptive unless its settings are selected only within development folds under a dated amendment.
4. Report patient-macro Dice at every prespecified forecast horizon, along with patient-level paired differences, Brier score for probabilistic outputs, event-count metrics for the count arm, missingness, failures, and fallback use. Freeze the forecast window from eligible cohort metadata before target labels are read; carry actual scan intervals into models that use time. Define empty-mask scoring before analysis (Dice 1 if both masks are empty and 0 if exactly one is empty), and require every method to use the same valid support and target coverage. Use a patient-level uncertainty method selected and validated before target-label scoring. Do not treat voxels, visits, or windows as independent patients.
5. Keep data preprocessing and registration rules source-backed and cohort-specific. MRI-only registration quality checks, valid-support construction, and target coverage must be frozen before target labels are read. Do not transfer LUMIERE thresholds or MU/Duke definitions to a new cohort without a dated amendment.

### 4. Track C — CGGA association and treatment-policy simulation

#### CGGA observational analysis

1. Use eligible CGGA GBM records under the source’s data-use terms. Apply the recorded duplicate rule by retaining a repeated patient identifier from the first batch only. Describe the primary historical WHO Grade IV population and report the IDH-wildtype subset as a sensitivity analysis; do not silently redefine the historical cohort as IDH-wildtype only.
2. Fit prespecified Cox models for the database’s recorded overall-survival time, with MGMT status, age, IDH status, radiotherapy, TMZ yes/no, recurrence status, and batch as specified covariates. Test the MGMT-by-recorded-TMZ interaction as an observational association. TMZ dose and timing are not available as treatment detail in the planned field, so do not infer regimen-specific benefit.
3. Verify and document the source time origin and field definitions. If the time origin cannot be verified, describe the endpoint only as the database’s recorded OS field. Use complete-case handling and report the analyzed denominator and missingness. No causal treatment claim will be made.

#### Mechanistic simulation and treatment-policy comparison

1. Represent sensitive-like tumor burden S and resistant-like burden R with the two-compartment system:

   dS/dt = rS[1 − (S + a_rsR)] − m r u S − μuS + ν(1 − u)R<br>
   dR/dt = r(1 − c)R[1 − (a_srS + R)] − εm r u R + μuS − ν(1 − u)R

   Initialize S(0)=(1−f0)B0 and R(0)=f0B0, with a_rs=a_sr=1. Integrate with fourth-order Runge–Kutta at a 0.25-day step and verify numerical convergence at a 0.05-day step. Parameters and treatment effects are simulated assumptions, not patient estimates.

2. Draw ZR, EB, and EI independently from standard normal distributions. Simulate two imperfect signals from latent molecular-like state ZR and biological state ZB using the declared dependence construction:

   zB = a zR + sqrt(1−a²) EB<br>
   zI = b zR + sqrt(1−b²)[d EB + sqrt(1−d²) EI]

   Set the latent binary state to 1 when zB is at least zero and generate the molecular-like report with error probability 1−q. Treat q as simulated report fidelity, not laboratory assay accuracy. Treat a and b as simulated signal loadings and d as residual dependence. Negative loadings are sign-reversal sensitivity scenarios, not negative controls. The full no-information control is (q,a,b,d)=(0.5,0,0,0); single-channel controls will be reported separately.

3. Compare finite-set minimax controllers with the same candidate parameter set, observations, treatment choices, and computational budget. The primary contrast will be a dependence-aware controller versus its matched version fixed at d=0. The model set will use the 135-model configuration, formed from q in {0.5, 0.6, 0.7, 0.8, 0.9} and a, b, d each in {−0.5, 0, 0.5}, after duplicate removal. Before outcome generation, code, configuration, and driver must be reconciled and tested against that exact set; the current prototype has a smaller grid. If this parity check fails, the policy run will not start.

4. At each replan, evaluate all eight treatment/no-treatment patterns over three 28-day blocks, forecast six 28-day cycles without a terminal value, apply the first three-block action, and replan at the next scheduled scan. Use five-day simulated courses and a maximum of 12 courses. Simulate trajectories through at most 1,456 days; cap progression-time utility at 730 days. The primary scan interval will be 84 days, with 28- and 56-day sensitivity intervals. Do not run the 112-day sensitivity unless the forecast horizon and controller are expanded and frozen before outcomes. The observable stopping proxy will require two consecutive visible noisy readings at least 1.4 times the lowest previous visible burden, with a floor of 0.10 times starting burden; a missing scan resets confirmation. This proxy is not formal RANO.

5. Draw independent worlds under a registered sampling contract. For the single-response family, draw log growth shift from [−0.7, 0.9], log treatment-scale factor from [ln(0.05), 0], log initial-resistant-fraction shift from [−0.7, 2.3], and starting-burden center from [0.12, 0.55] using a scrambled one-point Latin hypercube. For the mixture-response family, use log nonresponder treatment-scale factor [ln(0.02), ln(0.30)] and responder fraction pi_response [0.1, 0.6] in addition to the shared axes. Draw μ and ν independently per world from {0, 0.01, 0.04}; draw q once per world from {0.6, 0.7, 0.8, 0.9}; keep pi_response distinct from q. Generate selection, calibration, and evaluation cohorts from separate registered streams. Within each cohort, draw patient r log-uniformly from [0.001, 0.05] per day; draw baseline m log-uniformly from [6, 20]; draw baseline f0 log-uniformly from [0.001, 0.1], multiply by exp(world log-f0 shift), and cap at 0.9; draw B0 uniformly within ±0.10 of its world center; draw c uniformly from [0, 0.25] and ε uniformly from [0, 0.2]. In the single-response family, multiply m by exp(world log-treatment-scale factor). In mixture worlds, assign responders using pi_response; responders use baseline m, while nonresponders use baseline m multiplied by the weak world treatment-scale factor. Assign observation-noise σ values {0.05, 0.10, 0.20} in balanced counts across sorted patient identifiers within each world. Use explicit PCG64DXSM streams for world, cohort, sigma, signal, and uncertainty draws, and document any remaining generator used by the simulator. Share the same virtual patients, latent parameters, signals, scan schedule, and exogenous noise across paired policies. Validate the implemented sampling law and common-random-number coupling before any policy outcomes are generated.

6. Set the particle-filter prior for the correctly specified primary simulation to cover the full data-generating parameter support. Under the declared world shifts, the simulator’s r support is approximately 0.0005–0.123 per day and B0 support is 0.02–0.65; the current prototype prior uses r from 0.001–0.05 and B0 from 0.05–0.5. Those primary-prior mismatches must be corrected before evaluation. Verify the support of every other parameter against the implemented generator as well. A narrower prior may appear only as a separately named misspecification sensitivity, not as the primary analysis. Freeze every sensitivity, penalty, seed role, and controller setting before policy outcomes.

7. Define simulated utility as the minimum of the two progression times, each capped at 730 days, minus five simulated days per treatment course through day 730. Analyze two progression definitions separately: latent first passage at 1.2 times starting burden, and scan-sampled progression at 1.4 times the larger of the prior latent nadir or 0.10 times starting burden. Report course-use sensitivity at zero and ten simulated days per course. These are mathematical endpoints and do not represent clinical survival or formal RANO progression.

8. The population-mean estimand will remain the paired mean utility difference between the dependence-aware and matched d=0 policy, separately for single and mixture world families. The analysis unit will be an independent world. The success rule will require all simultaneous familywise 95% lower confidence bounds for the prespecified primary contrasts to exceed zero and each corresponding half-width to be no greater than 14 days. Do not switch to a fixed scenario-set estimand or lower the precision target after seeing outcomes.

9. The project already has simulator, signal-likelihood, particle-filter, controller, and world-orchestration prototypes, along with development tests and outcome-blind interval screens. These are methods-development records; they do not count as evaluation of the final contrast. Do not reuse development seeds as evaluation seeds. The replacement uncertainty method must receive independent statistical review and pass outcome-blind coverage checks that reproduce the actual world sampling, pairing, generator coupling, and familywise contrast structure. The current GE-betting prototype is a candidate only; it is not yet adopted. Use a maximum future evaluation budget of 336 continuous wall-clock hours on the specified PC. Before policy outcomes, materialize and collision-check the world seed schedule, build a complete run driver that enforces the configuration, select exact world counts using a separate outcome-blind precision and timing plan, and demonstrate that the 14-day interval half-width target is feasible within the 336-hour cap. If no world count meets both criteria, do not generate confirmatory policy outcomes; report the feasibility analysis and the failed gate.

### 5. Synthetic integration of Tracks A, B, and C

1. Keep the integration at the method level. Track A expression scores and Track B MRI transitions will not be joined to CGGA or to one another at the patient level. The simulation will use explicitly generated molecular-like and MRI-like inputs, with scenario parameters motivated by literature and reported separately from patient analyses.
2. Do not claim that Tracks A or B measure patient drug resistance. The 119-of-121 panel is an expression proxy derived from an in-vitro viability association; MRI recurrence is a radiographic endpoint. Neither is a matched clinical drug-response reference for the simulated patients.
3. Report the integrated result only as performance under the declared synthetic data-generating process. A positive result would motivate a future study with linked molecular, imaging, and direct response data; it would not demonstrate treatment benefit or patient-specific personalization.

### 6. Risk Assessment

This is a computational project. It will use no live biological agents, chemicals, animals, physical treatment, or newly recruited human participants. The main risks are privacy or unauthorized disclosure of restricted clinical data, misinterpretation of exploratory findings, and presenting simulated outcomes as clinical evidence.

To reduce these risks, the project will use de-identified records, store restricted files locally in access-controlled folders, avoid committing patient-level data or derived files unless written terms allow them, and retain only permitted aggregate outputs. Written data-use terms and school or competition review requirements will be checked before restricted-data processing. No model will recommend care for an individual. Reports will label simulation, observational association, exploratory analysis, and clinical measurement distinctly. No biological or therapeutic intervention will be performed.

### 7. Data Analysis

#### Track A

The primary classifier analysis will calculate class-balanced AURC per held-out donor and the equal-donor mean. Lower AURC indicates better selective prediction. The frozen success rule requires a lower mean than L2 logistic regression and wins in at least five of six donors. Balanced accuracy and Brier score will be secondary. Regional program analyses will use donor-level paired scores, expression-matched controls, and the declared exact sign-flip and covariate analyses. The TMZ viability panel scores will be reported as expression proxies with donor-level differences and matched-panel controls. No cell-level p-value will be treated as donor-level evidence.

#### Track B

For the exploratory transition analysis, calculate within-patient transition fractions for history-defined tissue and each physical-distance band; compare newly appearing with persistent tissue within the same patient and band. Preserve coverage as a separate observed/unobserved mask. For a future forecast study, calculate patient-macro Dice at each frozen horizon and paired patient-level differences against every required mask comparator; calculate count metrics only for the count model. Use patient-level resampling or another prespecified patient-level method, and report failures, missing patients, coverage, and fallbacks. No inference will treat voxels or visits as independent patients.

#### Track C observational

Report Cox hazard ratios with confidence intervals, model covariates, sample counts, missingness, proportional-hazards checks, and the MGMT-by-recorded-TMZ interaction. Interpret all estimates as retrospective associations. Do not call the interaction a treatment effect because treatment was not randomized and detailed dose and timing are not available.

#### Track C simulation and integration

For each policy pair, use paired world-level utility differences. Report each progression definition and world family separately, the mean difference, the simultaneous uncertainty interval, interval half-width, and the proportion of worlds with missing or failed arms. Include the no-information, single-channel, zero-dependence, and declared misspecification controls. A positive claim requires the complete frozen success rule. If the interval method fails validation, the precision target is infeasible, or a required arm is missing, report the analysis as not run, incomplete, or inconclusive as appropriate; do not select a more favorable endpoint or sampling rule after observing outcomes.

## Bibliography

Akbari, H., Macyszyn, L., Da, X., Bilello, M., Wolf, R. L., Martinez-Lage, M., Biros, G., Alonso-Basanta, M., O’Rourke, D. M., & Davatzikos, C. (2016). Imaging surrogates of infiltration obtained via multiparametric imaging pattern analysis predict subsequent location of recurrence of glioblastoma. *Neurosurgery, 78*(4), 572–580. https://doi.org/10.1227/NEU.0000000000001202

Chaudhuri, A., et al. (2023). Predictive digital twin for optimizing patient-specific radiotherapy regimens under uncertainty in high-grade gliomas. *Frontiers in Artificial Intelligence, 6*, 1222612. https://doi.org/10.3389/frai.2023.1222612

Ebrahimi Zade, A., Shahabi Haghighi, S., & Soltani, M. (2020). Reinforcement learning for optimal scheduling of glioblastoma treatment with temozolomide. *Computer Methods and Programs in Biomedicine, 193*, 105443. https://doi.org/10.1016/j.cmpb.2020.105443

Hegi, M. E., et al. (2005). MGMT gene silencing and benefit from temozolomide in glioblastoma. *The New England Journal of Medicine, 352*(10), 997–1003. https://doi.org/10.1056/NEJMoa043331

Liu, Q., et al. (2025). Treatment-aware diffusion probabilistic model for longitudinal MRI generation and diffuse glioma growth prediction. *IEEE Transactions on Medical Imaging*. https://doi.org/10.1109/TMI.2025.3533038

Liu, Z., et al. (2025). Multiscale mathematical model-informed reinforcement learning optimizes combination treatment scheduling in glioblastoma evolution. *Science Advances, 11*(32), eadv3316. https://doi.org/10.1126/sciadv.adv3316

Ntafoulis, I., et al. (2023). Ex vivo drug sensitivity screening predicts response to temozolomide in glioblastoma patients and identifies candidate biomarkers. *British Journal of Cancer, 129*(8), 1327–1338. https://doi.org/10.1038/s41416-023-02402-y

de Jong, G., Memi, F., Gracia, T., et al. (2026). A spatiotemporal cancer cell trajectory underlies glioblastoma heterogeneity. *Nature*. https://doi.org/10.1038/s41586-026-11234-8

Ruiz-Moreno, C., et al. (2025). Charting the single-cell and spatial landscape of IDH-wild-type glioblastoma with GBmap. *Neuro-Oncology, 27*(9), 2281–2295. https://doi.org/10.1093/neuonc/noaf113

Saul, D., et al. (2022). A new gene set library for analysis of senescent cells: SenMayo. *Nature Communications, 13*, 4827. https://doi.org/10.1038/s41467-022-32552-1

Suter, Y., et al. (2026). Drug and single-cell gene expression integration identifies sensitive and resistant glioblastoma cell populations. *Nature Communications, 17*, 99. https://doi.org/10.1038/s41467-025-67783-5

Wen, P. Y., et al. (2023). RANO 2.0: Update to the response assessment in neuro-oncology criteria for high- and low-grade gliomas in adults. *Journal of Clinical Oncology, 41*(33), 5187–5199. https://doi.org/10.1200/JCO.23.01059

Yu, K., Hu, Y., Wu, F., Guo, Q., Qian, Z., Hu, W., Chen, J., Wang, K., Fan, X., Wu, X., Rasko, J. E. J., Fan, X., Iavarone, A., Jiang, T., Tang, F., & Su, X.-D. (2020). Surveying brain tumor heterogeneity by single-cell RNA-sequencing of multi-sector biopsies. *National Science Review, 7*(8), 1306–1318. https://doi.org/10.1093/nsr/nwaa099

Criscuolo, E., et al. (2026). *A longitudinal glioblastoma MRI dataset with anatomical landmark pairs for DIR validation* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.20548528

National Center for Biotechnology Information. (2020). *GSE117891: Multiregional glioma single-cell transcriptomic sequencing* [Data set]. Gene Expression Omnibus. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117891

The Chinese Glioma Genome Atlas. (n.d.). *CGGA data portal* [Database]. https://www.cgga.org.cn/

Zhao, Z., et al. (2021). Chinese Glioma Genome Atlas (CGGA): A comprehensive resource with functional genomic data from Chinese glioma patients. *Genomics, Proteomics & Bioinformatics, 19*(1), 1–12. https://doi.org/10.1016/j.gpb.2020.10.005

Zimmer, L., et al. (2026). PREDICT-GBM: A multicenter platform advancing personalized glioblastoma radiotherapy planning. *npj Digital Medicine, 9*, 686. https://doi.org/10.1038/s41746-026-03194-0

## Additional Project Information

This is a continuation computational project using previously collected, de-identified datasets and simulated data. It will not recruit participants, contact patients, collect new clinical specimens, or deliver treatment. The project will follow each dataset’s access and reuse terms and any applicable school or competition review. Written permissions and data-use conditions will be retained before any further restricted-data processing or derived-output sharing. Because Track A, Track B, and CGGA do not have a verified common-patient link or matched direct drug-response endpoint, integration will occur only through the explicitly simulated signals and will be described as in-silico method evaluation.
