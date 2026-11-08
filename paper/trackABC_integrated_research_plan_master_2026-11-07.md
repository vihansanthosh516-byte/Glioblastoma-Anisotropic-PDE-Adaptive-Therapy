# Research Plan

## Project Title

**Integrating Tumor Biology and Longitudinal MRI for Glioblastoma Progression Modeling and In-Silico Therapy Optimization**

## Rationale

Glioblastoma changes across tumor regions and over time. Single-cell RNA sequencing can describe regional expression patterns, longitudinal MRI can describe evolving radiographic patterns, and mathematical models can represent tumor growth, treatment suppression, and resistant populations. These methods address different parts of the disease process, but the available Track A, Track B, and Track C records do not describe the same patients and do not directly measure the same biological variables.

This project will therefore keep three track-specific analyses distinct and connect them only in a synthetic virtual-patient experiment. Track A will examine donor-held-out prediction from mixed-cell expression labels. Track B will examine longitudinal MRI-label transitions and forecasting methods. Track C will examine observational MGMT associations and mathematical treatment policies. The integrated phase will generate molecular-like and MRI-like signals for the same simulated virtual patient; it will not join Track A and Track B people or treat either data source as a direct measurement of drug resistance.

Prior studies have combined tumor-growth models, imaging, molecular features, digital twins, and treatment scheduling (Akbari et al., 2016; Bhattacharya et al., 2025; Chaudhuri et al., 2023; Ebrahimi Zade et al., 2020; Q. Liu et al., 2025; Zimmer et al., 2025; Z. Liu et al., 2025). The project's narrow integrated question will be whether a treatment controller that represents uncertainty about dependence between imperfect signals behaves differently from a matched controller that assumes conditional independence. The project will not claim that broad multimodal treatment modeling is new, that one method is best, or that the simulator can choose treatment for a person. Any positive, null, or harmful result will apply only to the specified data analyses or simulated assumptions. Possible societal value is limited to helping researchers identify methods and assumptions that warrant later validation.

## Research Questions, Hypotheses, Engineering Goals, and Expected Outcomes

### Track A — expression and regional prediction

**Research question.** Under the frozen six-donor nested leave-one-donor-out design, will the sparse donor-balanced expression classifier have lower equal-donor class-balanced area under the risk-coverage curve than ordinary L2 logistic regression and lower AURC for at least five of six held-out donors?

**Hypothesis for the completed A5 procedure.** Before this analysis, the candidate was expected to meet both parts of the frozen comparison rule: lower equal-donor mean class-balanced AURC than L2 logistic regression and lower donor-specific AURC in at least five of the six primary donors. Lower AURC indicates lower error among cells retained at each confidence level.

**Secondary question.** In the separately specified paired-donor signature analyses, will the selected senescence, oxidative-phosphorylation, and tumor-associated-macrophage programs show region-associated patterns that remain interpretable under expression-matched gene-set controls and prespecified adjustment?

### Track B — longitudinal MRI and tissue-state transitions

**Research question.** In a longitudinal glioblastoma MRI cohort, how will history-defined peritumoral FLAIR regions and their physical distance from the current core relate to the next-visit core-label transition, and will a prespecified MRI forecast add information beyond history-only mask baselines?

**Exploratory working hypothesis.** In a future eligible independent cohort, newly appearing FLAIR-labeled tissue will have a higher next-visit core-label transition rate than persistent FLAIR-labeled tissue in the same patient and distance band. This is a working hypothesis for a future independent evaluation; the already-examined Duke cohort will not be presented as an untouched confirmatory test.

### Track C — MGMT association and treatment-policy simulation

**Observational question.** Among eligible glioblastoma records in the public CGGA data, how will MGMT methylation status relate to overall survival after adjustment for the prespecified clinical covariates, and will the association differ by recorded temozolomide status?

**Hypothesis for the completed observational analysis.** The analysis was specified to test whether MGMT methylation was associated with longer overall survival among records marked as temozolomide-treated after covariate adjustment. The MGMT-by-temozolomide analysis is an observational association, not a causal treatment effect.

**Simulation question.** Across independently generated glioblastoma-like worlds, will a finite-set dependence-aware minimax controller achieve higher prespecified robust utility than a matched controller with residual signal dependence fixed at zero, when both receive the same simulated signals, virtual patients, treatment choices, and observation schedules?

**Simulation hypothesis.** In correctly specified scenarios with informative signals and nonzero residual dependence, the dependence-aware controller will have higher equal-world mean robust utility than its matched conditional-independence controller. In uninformative or reversed-signal controls, it will not consistently outperform the no-signal or single-signal controls.

### Integrated question and hypothesis

**Integrated question.** When molecular-like and MRI-like observations are generated jointly from the same simulated tumor state, does explicitly modeling their residual dependence improve the treatment policy's robust simulated utility over a matched conditional-independence policy?

**Integrated hypothesis.** The dependence-aware controller will improve the primary utility contrast in the predeclared informative, nonzero-dependence scenarios while showing no consistent advantage in the predeclared null and reversed-signal controls. This is a simulation hypothesis, not a claim about patient benefit.

### Engineering goals

1. Keep completed methods, remaining methods, data permissions, and research outcomes distinguishable in the project records.
2. Preserve Track A's donor-level split, training-only preprocessing, frozen comparators, and success rule.
3. Audit Track B's source terms, cohort definitions, registration, label support, and target-field coverage before any additional patient-level analysis.
4. Keep the observational MGMT analysis separate from the treatment-policy simulator and avoid causal interpretation.
5. Build a reproducible Track C benchmark with shared virtual patients, common random inputs for paired arms, separate development and evaluation worlds, recorded seeds, and a pre-run code and settings seal.
6. Connect the tracks only with explicitly simulated signals; never match people across datasets or treat synthetic signal parameters as measured resistance biology.

### Expected outcomes

The Track A and completed Track C procedures use their specified analyses to estimate donor-level prediction and observational associations. The exploratory Track B procedure describes MRI-label transitions, while a future independent cohort would test its working transition hypothesis. The proposed treatment-policy simulation may favor the dependence-aware controller, show similar performance, or favor a comparator in some scenarios. This plan reports no observed outcome values or findings; any eventual results will be reported separately and interpreted within their track and design.

## Materials

- **Computer:** Dell 16 Plus DB16250 with Windows 11 and 32 GB RAM. Processor, thread count, and runtime environment will be recorded for any locked run.
- **Software:** Python 3.14.4 and the package snapshot in [the Windows Python environment lock](../tracks/shared/requirements-win-py314.lock.txt). The snapshot records the current environment; a clean install on a second computer remains to be checked before a sealed evaluation.
- **Reproducibility files:** the current proposed settings are in [the integrated v2 settings](../configs/trackABC_integrated_v2.json); stream-role templates are in [the v2 seed-stream file](../configs/trackABC_seed_streams_v2.csv); the input-rights and checksum ledger is [the v2 input manifest](../configs/trackABC_input_manifest_v2.csv); and the reproduction sequence is in [the handoff guide](../tracks/shared/HANDOFF_AND_REPRODUCIBILITY.md). These files define a planning draft, not a frozen execution seal.
- **Track A sources:** Yu et al. single-cell RNA data and sample metadata (GEO GSE117891), Table S1 point labels, GBmap annotation metadata, Ivy GAP bulk RNA data (GSE107559), Angers microarray data (E-MTAB-5092), and bulk expression data (E-MTAB-8743). Public metadata records screened for a possible external cohort include PRJNA1337938, HRA009021, GSE286419, GSE318562, GSE318564, and the Harwood Zenodo record 22013127. Screened records will not contribute expression data unless they pass the stated eligibility and rights rules.
- **Track B sources:** the restricted Duke longitudinal glioblastoma MRI record (Zenodo 20548528), MU-Glioma-Post and LUMIERE source records for earlier protocol development, and synthetic images and masks for method-mechanics checks. Zenodo lists CC BY-NC-ND 4.0 as an access option and describes a separate Duke custom-license route; the exact accepted terms for the local copy will govern its use. Duke files and identifiers will remain in access-controlled local storage, and derived outputs will not be shared unless the accepted terms permit it. MU and LUMIERE will not supply patient outcomes to the new integrated simulation.
- **Track C sources:** the project's two-compartment growth and treatment-suppression simulator; CGGA clinical and molecular records for the separate observational analysis; a TCGA-derived published reference used only as context; and published aggregate RTOG 0525 and GEINO 14-01 results for compatible model checks. Aggregate trial data will not be treated as patient-level dose-response records.
- **Physical materials:** none. The project uses computer files and simulations only; it uses no chemicals, bacteria, cell cultures, animals, or experimental drugs.

## Procedures

### Track A — expression and regional biology

#### Completed procedure to date

1. **Freeze and verify the analysis inputs.** The A5 protocol, donor list, appendix, count matrix, cell table, Table S1, and GBmap snapshot were checksummed against the recorded seal before analysis. The frozen primary donor set was GS1, GS2, GS3, GS5, GS6, and GS13.
2. **Map cells and define labels.** The analysis joined each cell to its sampling point by exact identifiers, then joined points to Table S1. Tumoral points were labeled Core, Peritumoral points Periphery, and Normal points were excluded. The code checked consistency between Table S1 and the donor list.
3. **Apply cell and feature rules.** Cells were retained only with at least 1,000 UMIs and at least 500 detected genes. In each training split, genes detected in at least 5% of training cells were retained, the 2,000 most variable genes were selected, each cell was normalized to 10,000 counts and log1p-transformed, and feature means and standard deviations were fitted on training cells only.
4. **Fit and compare models using donor-held-out splits.** Six outer folds held out one donor each. Inner folds held out each of the remaining donors for setting selection. The fixed candidate used a gene gate, worst-donor loss, a neighbor graph, and a neighbor blend. The same folds were used for L2 logistic regression, the 15-nearest-neighbor vote, linear SVM, gradient-boosted trees, three candidate ablations, cell-type-only and cells-per-point leakage diagnostics. The analysis used the recorded seed 20261031. Constant-score diagnostic behavior was retained as a control rather than relabeled as substantive model performance.
5. **Run regional-program and robustness procedures.** The project conducted paired-donor SenMayo and OXPHOS analyses in Ivy GAP and Angers data and a TAM-activation analysis in E-MTAB-8743. The procedures used donor-paired scores, expression-matched gene-set controls, exact sign-flip tests, and the prespecified hypoxia, neuron-related, and myeloid adjustments. Audits examined donor and library composition, same-cell-type subsets, and the separation between donor-level validation and random-cell validation. Calibration curves were listed as deferred rather than added to A5. External-cohort screening was restricted to public metadata and documented gates; no screened dataset is treated here as an external test.

#### Remaining procedure

1. Keep the frozen A5 model, endpoint, donor list, folds, and success rule unchanged. Any change to these requires a dated amendment before analysis.
2. Decide by amendment whether to run the deferred calibration-curve, within-library, scANVI, and attention analyses. They will not be added to A5 after inspecting its held-out outputs.
3. Consider an external transfer analysis only if one cohort has at least six independent matched donors, a defensible region-to-sample map, validated malignant-cell labels where a malignant-cell claim is intended, clear reuse terms, and verified independence from the project's training and prior-use records. Freeze the target features, model, controls, and donor-level validation before reading expression values. If all gates do not pass, no external-confirmation analysis will be performed.

### Track B — longitudinal MRI and spatial transition analysis

#### Completed procedure to date

1. **Review sources and build synthetic-only methods.** The project reviewed source, access, label, support, and registration documentation for MU-Glioma-Post, LUMIERE, and Duke. It developed synthetic-only longitudinal forecasting components and comparators, including persistence, fixed growth, PDE-only and hybrid forecasts, a 3D U-Net and time-conditioned FiLM U-Net, a count-only model, fold-local calibration, a run guard, and provenance recording. Registration, geometry, coverage, and baseline mechanics were also checked on synthetic inputs. These steps developed and audited code; this plan does not report their test counts or timing results.
2. **Complete the exploratory Duke procedure.** After access was obtained, the local Duke longitudinal MRI archive was used for an exploratory patient-mask analysis. The saved analysis freeze and cohort record defined a 10 mm dilation of prior-visit FLAIR, excluded tissue within 10 mm of the middle-visit core, and summarized subsequent core-label transitions in three distance bands. The analysis read target-label masks and generated patient-level summaries. This cohort and endpoint are exposed and will not be described as an untouched validation set.
3. **Audit the saved workflow without repeating the outcome analysis.** The project checked saved-file provenance, cohort identifiers, analysis-script hashes, and arithmetic from saved summaries. Separate synthetic work examined physical-coordinate transforms and whether resampling preserved an explicit coverage field. These audits do not replace review of the source images, registration quality, label definitions, or target support.

#### Remaining procedure

1. Review and retain the exact accepted Duke terms, school and competition conditions, and any rules on derived summaries, images, model weights, or publication. If permission for a proposed derived use is unclear, do not use or share that output.
2. Independently audit the completed Duke protocol and code, including session chronology, registration, target-field coverage, support, and the relation between the saved analysis freeze and the executed pipeline. Record corrections as a new dated amendment; do not retune or rescore the exposed Duke cohort to seek a preferred result.
3. Keep the separate ADC workflow incomplete until its read scope, source file definitions, preprocessing, quality checks, and analysis question are frozen and approved. Do not infer that an incomplete preflight established what was or was not read in other workflows.
4. Use a future patient cohort only if its access terms, label provenance, longitudinal visit map, registration QC, support, prior-use overlap, and eligible sample count are documented before target labels are opened. Freeze the transition estimand and the same-patient, same-distance controls first. If no independent cohort passes, report no external replication.

### Track C — MGMT association and treatment-policy modeling

#### Completed procedure to date

1. **Carry out the separate observational MGMT analysis.** The CGGA analysis combined the two released clinical-data batches after removing duplicate patient identifiers, applied the recorded glioblastoma and survival-data rules, and included a batch indicator. Cox proportional-hazards models examined age, IDH status, radiotherapy, temozolomide status, recurrence status, and MGMT methylation. The MGMT analysis included a model among temozolomide-treated records and an MGMT-by-temozolomide interaction model in the eligible pooled records. The TCGA age estimate was retained only as a published or previously derived reference, not fitted as part of the CGGA analysis. This is an observational OS analysis; treatment assignment was not randomized.
2. **Develop the mathematical tumor and therapy simulator.** The project implemented sensitive and resistant tumor compartments, logistic competition, treatment suppression, fitness cost, and transition between compartments. The model includes fixed schedules, scan-triggered policies, Q-learning, belief-space model-predictive control, finite-set signal-ambiguity control, and oracle policies for diagnostic comparison. The oracle uses hidden simulated parameters and is not treated as a deployable controller.
3. **Develop signal and benchmark methods.** The project implemented molecular-like and MRI-like signal generation, joint signal likelihoods, particle filtering, common policy inputs, scan-noise and scan-loss mechanics, numerical sensitivity checks, shared policy benchmarks, and development-world seed and provenance utilities. Literature review and sensitivity work were used to identify assumptions and controls that need to be explicit. Previously exposed seeds and development outputs will not be used as fresh evaluation data.

#### Remaining procedure

1. Decide and freeze the treatment model's parameter ranges, response families, signal anchors, observation process, progression definitions, policy information boundary, treatment-course budget, and primary contrast before a research-scale evaluation.
2. Build and review a research-scale driver that reserves fresh independent world seeds, separates selection, calibration, and evaluation patients, records all input and model fingerprints, enforces common inputs for paired policies, and stops on missing or failed arms.
3. Use outcome-blind timing and precision planning to set the number of worlds and virtual patients. Materialize and collision-check numeric seeds, set a runtime cap, validate numerical convergence and interval behavior, and seal the settings, source, environment, and analysis code before evaluation. If the precision or resource target cannot be met, report feasibility or an inconclusive comparison without changing the rules after evaluation outcomes are seen.
4. Keep the CGGA observational association separate from simulated policy evaluation. Do not interpret it as causal evidence that temozolomide works differently by MGMT, and do not use it as patient-specific simulator calibration.

### Integrated procedure — Track A to Track B to Track C

1. **Translate, do not merge, the evidence.** Track A's expression work motivates a hypothetical binary molecular-like observation. The integrated model will use a predeclared synthetic report-quality grid of 0.5, 0.6, 0.7, 0.8, and 0.9; it will not use individual A5 validation metrics as estimates of resistance-assay accuracy. Track B motivates a longitudinal MRI-like observation and scan schedule. The integrated simulation will use a standard-normal MRI-like score unless an authorized, history-only aggregate is specifically approved and frozen. No Track A or Track B person will be linked to another dataset or copied into the virtual cohort.
2. **Generate joint synthetic signals.** For each virtual patient, map resistant fraction f0 to a standard-normal latent score zR using fixed log-uniform reference bounds L = 0.001 × exp(−0.7) and U = 0.9. Generate independent standard-normal errors EB and EI and calculate zB = a zR + sqrt(1 − a²) EB and zI = b zR + sqrt(1 − b²) [d EB + sqrt(1 − d²) EI]. The hidden binary state is C = 1[zB ≥ 0]; the observed molecular-like report flips C with probability 1 − q. The continuous MRI-like score is zI. Coefficients a, b, d, and q are simulation assumptions, not patient measurements. The controller will use the exact joint likelihood in the handoff guide and config v2.
3. **Generate tumor-response worlds.** Use sensitive burden S, resistant burden R, total burden B = S + R, and capacity 1. Use dS/dt = rS[1 − (S + a_rsR)] − mruS − μuS + ν(1 − u)R and dR/dt = r(1 − c)R[1 − (a_srS + R)] − εmruR + μuS − ν(1 − u)R, with S(0) = (1 − f0)B0, R(0) = f0B0, and a_rs = a_sr = 1. Generate single-response and mixture-response families using the ranges and Latin-hypercube order in config v2. Apply μ and ν levels of 0, 0.01, or 0.04 per day across the world design. Draw patient-level c from 0 to 0.25 and ε from 0 to 0.2. Use classical fourth-order Runge–Kutta at 0.25-day steps and compare with 0.05-day steps for numerical convergence.
4. **Compare matched treatment controllers.** Compare a finite-set dependence-aware minimax controller with a matched finite-set controller whose residual-dependence parameter d is fixed to zero. Use the same patient parameters, signals, scan schedule, observation noise, and treatment options within each paired comparison. The proposed model set uses q in the synthetic grid and a, b, d in {−0.5, 0, 0.5}; this gives at most 135 model combinations before duplicate removal. The proposed primary anchor is (a, b, d) = (0.5, 0.5, 0.5); include full no-information and molecular-channel-null controls and report out-of-grid anchors separately. Each model uses 500 filter particles and 16 common stratified prediction draws. Score all eight treatment/no-treatment patterns over three 28-day cycles, forecast six cycles without terminal value, apply the first three-cycle block, and replan. The main scan interval is 84 days; scan intervals of 28, 56, and 112 days and visible scan loss are sensitivity scenarios. A modeled treatment course treats for five days, with a maximum of 12 courses over a 1,456-day simulation.
5. **Separate development and evaluation.** Use disjoint selection, calibration, and evaluation virtual-patient cohorts and distinct random streams. The final world count, cohort sizes, anchors, seed table, inference family, and resource limit are not yet locked. Select these using a written outcome-blind precision and timing procedure; keep single-response and mixture-response families separate. Do not use previously exposed seeds or outcomes for evaluation.
6. **Define outcomes and inference before evaluation.** Record two latent progression measures through 730 days: interpolated time to first reach 1.2 times starting burden, and first scheduled scan at which burden reaches 1.4 times the larger of the prior scan nadir or 0.10 times starting burden. These are simulator outcomes; the second is not a clinical RANO assessment. Proposed utility is the smaller of the two progression times, each capped at 730 days, minus five simulated days per course. The planned primary contrast is the paired candidate-minus-control utility difference within each independent world, averaged equally across worlds. The draft analysis proposes 10,000 world-level bootstrap resamples with max-T adjustment across the primary response-family contrasts and a 14-day confidence-interval half-width precision target. Freeze and validate the exact interval procedure before evaluation. A positive result requires the adjusted lower bounds for all primary contrasts to exceed zero and the precision target to be met.
7. **Seal, run, and report.** Before evaluation, validate the driver on synthetic records with known effects, verify model and input hashes, reserve numeric seeds, record package and operating-system versions, and freeze code, settings, and analysis rules. Run only after the seal passes. Preserve missing arms, failures, fallbacks, exclusions, and all prespecified sensitivity results. Report results in the final project report, not in this research plan. The work will not select therapy for a patient or establish clinical treatment benefit.

## Risk Assessment

Tracks A and B use existing human-derived research data, including access-controlled Duke MRI and segmentation files. No participants will be recruited or contacted, and no treatment will be administered. The completed Duke procedure included patient-mask reads. Any further access, analysis, retention, or sharing will follow the exact accepted data-use terms and school or competition conditions. Zenodo lists CC BY-NC-ND 4.0 as an access option and describes a separate Duke custom-license route; derived outputs will remain local unless the accepted terms permit sharing. Any unresolved permission or data-quality issue will stop the affected analysis. Track C's treatment and resistance variables are mathematical assumptions, not medical instructions. The integrated phase uses synthetic virtual patients and will not recommend treatment for an individual. Long computations will be run on a ventilated surface with the computer connected to power.


## Data Analysis

### Track A

The primary classifier comparison uses class-balanced AURC, averaged equally across the six donors; lower is better. The frozen descriptive success rule requires both lower mean AURC and lower AURC in at least five of six donors. Balanced accuracy and Brier score are secondary. Per-donor values will be shown in the final results report. Within-cell-type analysis will be limited to strata meeting the frozen cell and donor-count rules. The separate signature analyses use paired donor-level scores, expression-matched gene-set controls, exact sign-flip tests, and the prespecified covariate adjustments. Cells will not be treated as independent donors.

### Track B

For the completed exploratory Duke analysis, the unit of summary is the patient, with the saved transition definition and three distance bands. Any audit will first verify the source-to-analysis mapping, registration, support, and target-field coverage. An independent future transition analysis will compare new and persistent history-defined FLAIR regions within patient and physical-distance band, using a patient-level uncertainty method frozen before target-label access. Forecasting comparisons will use the same patients, support, target, and time interval for persistence, fixed-growth, PDE-only, hybrid, and learned baselines. ADC will remain a separate secondary analysis with its own completeness and missingness rules. No outcome values are included here.

### Track C

The observational MGMT analysis uses Cox proportional-hazards models for overall survival, with complete-case reporting and the declared clinical covariates. MGMT analyses among temozolomide-treated records and the MGMT-by-temozolomide interaction will be reported as noncausal associations. The synthetic policy comparison uses paired policy differences within each independently generated world, averages virtual patients within a world, and weights worlds equally. It reports both latent progression measures, treatment courses, censoring, model failures, and sensitivity scenarios. The independent world, not the nested patient, scan, or noise draw, is the uncertainty unit.

### Integrated simulation

For each policy, progression time one is the interpolated time when true burden first reaches 1.2 times starting burden. Progression time two is the first scheduled scan at which true burden reaches 1.4 times the larger of the prior true-burden scan nadir or 0.10 times starting burden. The latter is a latent-burden scan-sampled proxy, not formal RANO. Proposed utility is the smaller of those times, each capped at 730 days, minus five simulated days per treatment course delivered through day 730; course-cost sensitivities will use zero and ten days. The primary contrast is candidate utility minus matched d=0 utility, paired within each world and averaged equally across worlds. The draft proposes 10,000 world-level bootstrap resamples and max-T intervals across the primary response-family contrasts, with a 95% interval half-width target of at most 14 simulated days. The inference method, seed, precision rule, anchors, sample sizes, and world count will be fixed and validated before evaluation. The positive decision rule requires every adjusted lower confidence bound to exceed zero and the precision target to be met. If a resource, precision, seed, or required-arm gate fails, the comparison will be reported as feasibility-only or inconclusive. No outcomes are reported in this plan.


## Bibliography

Akbari, H., Macyszyn, L., Da, X., Bilello, M., Wolf, R. L., Martinez-Lage, M., Biros, G., Alonso-Basanta, M., O’Rourke, D. M., & Davatzikos, C. (2016). Imaging surrogates of infiltration obtained via multiparametric imaging pattern analysis predict subsequent location of recurrence of glioblastoma. *Neurosurgery, 78*(4), 572–580. https://doi.org/10.1227/NEU.0000000000001202

Balaña, C., Vaz, M. A., Sepúlveda, J. M., Mesia, C., del Barco, S., Pineda, E., Muñoz-Langa, J., Estival, A., de las Peñas, R., Fuster, J., Gironés, R., Navarro, L. M., Gil-Gil, M., Alonso, M., Herrero, A., Peralta, S., Olier, C., Perez-Segura, P., Covela, M., & Carrato, C. (2020). A phase II randomized, multicenter, open-label trial of continuing adjuvant temozolomide beyond 6 cycles in patients with glioblastoma (GEINO 14-01). *Neuro-Oncology, 22*(12), 1851–1861. https://doi.org/10.1093/neuonc/noaa107

Bhattacharya, M., Singh, G., & Prasanna, P. (2025). *SoC-DT: Standard-of-Care aligned digital twins for patient-specific tumor dynamics* [Preprint]. arXiv. https://arxiv.org/abs/2510.03287

Chaudhuri, A., Pash, G., Hormuth, D. A., II, Lorenzo, G., Kapteyn, M., Wu, C., Lima, E. A. B. F., Yankeelov, T. E., & Willcox, K. (2023). Predictive digital twin for optimizing patient-specific radiotherapy regimens under uncertainty in high-grade gliomas. *Frontiers in Artificial Intelligence, 6*, 1222612. https://doi.org/10.3389/frai.2023.1222612

Criscuolo, E., Calabrese, E., Zhang, Z., Wang, Y., Xiong, Y., & Yang, D. (2026). *A longitudinal glioblastoma MRI dataset with anatomical landmark pairs for DIR validation* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.20548528

Ebrahimi Zade, A., Shahabi Haghighi, S., & Soltani, M. (2020). Reinforcement learning for optimal scheduling of glioblastoma treatment with temozolomide. *Computer Methods and Programs in Biomedicine, 193*, 105443. https://doi.org/10.1016/j.cmpb.2020.105443

Ellingson, B. M., Wen, P. Y., & Cloughesy, T. F. (2017). Modified criteria for radiographic response assessment in glioblastoma clinical trials. *Neurotherapeutics, 14*(2), 307–320. https://doi.org/10.1007/s13311-016-0507-6

Gilbert, M. R., Wang, M., Aldape, K. D., Stupp, R., Hegi, M. E., Jaeckle, K. A., Armstrong, T. S., Wefel, J. S., Won, M., Blumenthal, D. T., Mahajan, A., Schultz, L. M., Erridge, S. C., Brown, P. D., Chakravarti, A., Pugh, S., Zadeh, G., Fisher, B., Sulman, E. P., & Mehta, M. P. (2013). Dose-dense temozolomide for newly diagnosed glioblastoma: A randomized phase III clinical trial. *Journal of Clinical Oncology, 31*(32), 4085–4091. https://doi.org/10.1200/JCO.2013.49.6968

Liu, Q., Fuster-Garcia, E., Hovden, I. T., MacIntosh, B. J., Grødem, E., Brandal, P., Lopez-Mateu, C., Sederevicius, D., Skogen, K., Schellhorn, T., Bjørnerud, A., & Emblem, K. E. (2025). Treatment-aware diffusion probabilistic model for longitudinal MRI generation and diffuse glioma growth prediction. *IEEE Transactions on Medical Imaging*. https://doi.org/10.1109/TMI.2025.3533038

Liu, Z., Zhang, J., Hong, L., Nie, Q., & Sun, X. (2025). Multiscale mathematical model-informed reinforcement learning optimizes combination treatment scheduling in glioblastoma evolution. *Science Advances, 11*(32), eadv3316. https://doi.org/10.1126/sciadv.adv3316

Stupp, R., Mason, W. P., van den Bent, M. J., Weller, M., Fisher, B., Taphoorn, M. J. B., Belanger, K., Brandes, A. A., Marosi, C., Bogdahn, U., Curschmann, J., Janzer, R. C., Ludwin, S. K., Gorlia, T., Allgeier, A., Lacombe, D., Cairncross, J. G., Eisenhauer, E., & Mirimanoff, R. O. (2005). Radiotherapy plus concomitant and adjuvant temozolomide for glioblastoma. *New England Journal of Medicine, 352*(10), 987–996. https://doi.org/10.1056/NEJMoa043330

Yu, K., Hu, Y., Wu, F., Guo, Q., Qian, Z., Hu, W., Chen, J., Wang, K., Fan, X., Wu, X., Rasko, J. E., Fan, X., Iavarone, A., Jiang, T., Tang, F., & Su, X.-D. (2020). Surveying brain tumor heterogeneity by single-cell RNA-sequencing of multi-sector biopsies. *National Science Review, 7*(8), 1306–1318. https://doi.org/10.1093/nsr/nwaa099

Zimmer, L., Weidner, J., Balcerak, M., Kofler, F., Krupa, M., Ezhov, I., Cepeda, S., Zhang, R., Lowengrub, J., Menze, B., & Wiestler, B. (2025). *PREDICT-GBM: A multi-center platform to advance personalized glioblastoma radiotherapy planning* [Preprint]. arXiv. https://arxiv.org/abs/2509.13360

Hegi, M. E., Diserens, A.-C., Gorlia, T., Hamou, M.-F., de Tribolet, N., Weller, M., Kros, J. M., Hainfellner, J. A., Mason, W., Mariani, L., Bromberg, J. E. C., Hau, P., Mirimanoff, R. O., Cairncross, J. G., Janzer, R. C., & Stupp, R. (2005). MGMT gene silencing and benefit from temozolomide in glioblastoma. *New England Journal of Medicine, 352*(10), 997–1003. https://doi.org/10.1056/NEJMoa043331


## Additional Project Information

This is a continuation project. The Procedures section identifies completed method steps and remaining planned steps for Tracks A, B, and C. This research-plan document contains no observed study-result values or outcome conclusions; any required continuation abstract or results report will be supplied separately under the competition's instructions. Track A and Track B data are not patient-linked. Track C treatment policies are simulated and are not clinical recommendations.
