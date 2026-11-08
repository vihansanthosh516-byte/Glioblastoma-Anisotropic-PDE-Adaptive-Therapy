# Research Plan

## Project Title

**Integrating Tumor Biology and Longitudinal MRI for Glioblastoma Progression Modeling and In-Silico Therapy Optimization**

## Rationale

Glioblastoma changes across tumor regions and over time. Single-cell RNA sequencing can describe regional expression patterns, longitudinal MRI can describe evolving radiographic patterns, and mathematical models can represent tumor growth, treatment suppression, and resistant populations. These methods address different parts of the disease process, but the available Track A, Track B, and Track C records do not describe the same patients and do not directly measure the same biological variables.

This project will keep three track-specific analyses distinct and connect them only in a synthetic virtual-patient experiment. Track A will examine donor-held-out prediction from mixed-cell expression labels and, separately, a published temozolomide (TMZ)-response-associated expression signature. Track B will examine longitudinal MRI-label transitions and forecasting methods; its imaging outcomes will be called treatment-associated only if patient-level exposure and outcome timing can be verified. Track C will examine observational MGMT associations and mathematical treatment policies. The integrated phase will generate molecular-like and MRI-like proxy signals for the same simulated virtual patient. The available Track A and Track B cohorts are not matched and do not contain a verified shared drug-response outcome, so neither will be treated as a direct measurement of resistance or used to estimate the joint signal relationship.

Prior studies have combined tumor-growth models, imaging, molecular features, digital twins, and treatment scheduling (Akbari et al., 2016; Bhattacharya et al., 2025; Chaudhuri et al., 2023; Ebrahimi Zade et al., 2020; Q. Liu et al., 2025; Zimmer et al., 2025; Z. Liu et al., 2025). MGMT promoter methylation is associated with TMZ benefit, while patient-derived ex-vivo drug-response assays measure cellular sensitivity more directly (Hegi et al., 2005; Ntafoulis et al., 2023). The project's narrow integrated question will be whether a treatment controller that represents uncertainty about dependence between imperfect resistance-associated proxy signals behaves differently from a matched controller that assumes conditional independence. The project will not claim that the current patient cohorts directly measure TMZ resistance, that broad multimodal treatment modeling is new, that one method is best, or that the simulator can choose treatment for a person. Any positive, null, or harmful result will apply only to the specified data analyses or simulated assumptions. Possible societal value is limited to helping researchers identify methods and assumptions that warrant later validation.

## Research Questions, Hypotheses, Engineering Goals, and Expected Outcomes

### Track A — expression and regional prediction

**Research question.** Under the frozen six-donor nested leave-one-donor-out design, will the sparse donor-balanced expression classifier have lower equal-donor class-balanced area under the risk-coverage curve than ordinary L2 logistic regression and lower AURC for at least five of six held-out donors?

**Hypothesis for the completed A5 procedure.** Before this analysis, the candidate was expected to meet both parts of the frozen comparison rule: lower equal-donor mean class-balanced AURC than L2 logistic regression and lower donor-specific AURC in at least five of the six primary donors. Lower AURC indicates lower error among cells retained at each confidence level.

**Secondary question.** In the separately specified paired-donor signature analyses, will the selected senescence, oxidative-phosphorylation, and tumor-associated-macrophage programs show region-associated patterns that remain interpretable under expression-matched gene-set controls and prespecified adjustment?

**Resistance-associated proxy question.** Will the 21-gene subset of TMZ-sensitivity-associated genes that Ntafoulis et al. (2023) reported as externally validated show a consistent Core-versus-Periphery difference across the six paired donors? This will be an exploratory molecular proxy analysis; it will not test whether the donors' tumors are clinically sensitive or resistant to TMZ.

### Track B — longitudinal MRI and tissue-state transitions

**Research question.** In a longitudinal glioblastoma MRI cohort, how will history-defined peritumoral FLAIR regions and their physical distance from the current core relate to the next-visit core-label transition, and will a prespecified MRI forecast add information beyond history-only mask baselines? Where patient-level TMZ exposure and progression chronology can be verified, will the imaging pattern also support a treatment-associated radiographic proxy analysis?

**Exploratory working hypothesis.** Newly appearing FLAIR-labeled tissue will have a higher next-visit core-label transition rate than persistent FLAIR-labeled tissue in the same patient and distance band. The Duke analysis is exploratory because its cohort and endpoint have been examined. MRI-defined recurrence or core transition will not be called TMZ resistance unless the relevant exposure and clinical outcome chronology are verified; even then, the result will be a radiographic treatment-associated proxy, not a direct drug-sensitivity measurement.

### Track C — MGMT association and treatment-policy simulation

**Observational question.** Among eligible glioblastoma records in the public CGGA data, how will MGMT methylation status, a treatment-associated biomarker, relate to overall survival after adjustment for the prespecified clinical covariates, and will the association differ by recorded temozolomide status?

**Hypothesis for the completed observational analysis.** The analysis was specified to test whether MGMT methylation was associated with longer overall survival among records marked as temozolomide-treated after covariate adjustment. The MGMT-by-temozolomide analysis is an observational association, not a causal treatment effect.

**Simulation question.** Across independently generated glioblastoma-like worlds, will a finite-set dependence-aware minimax controller achieve higher prespecified robust utility than a matched controller with residual signal dependence fixed at zero, when both receive the same simulated signals, virtual patients, treatment choices, and observation schedules?

**Simulation hypothesis.** In correctly specified scenarios with informative signals and nonzero residual dependence, the dependence-aware controller will have higher equal-world mean robust utility than its matched conditional-independence controller. In uninformative or reversed-signal controls, it will not consistently outperform the no-signal or single-signal controls.

### Integrated question and hypothesis

**Integrated question.** When a molecular resistance-associated proxy signal and an MRI progression proxy signal are generated jointly from the same simulated resistant-tumor state, does explicitly modeling their residual dependence improve robust simulated treatment utility over a matched conditional-independence policy?

**Integrated hypothesis.** The dependence-aware controller will improve the primary utility contrast in the predeclared informative, nonzero-dependence scenarios while showing no consistent advantage in the predeclared null and reversed-signal controls. This is a simulation hypothesis, not a claim about patient benefit.

### Engineering goals

1. Keep completed methods, remaining methods, data permissions, and research outcomes distinguishable in the project records.
2. Preserve Track A's donor-level split, training-only preprocessing, frozen comparators, and success rule.
3. Audit Track B's source terms, cohort definitions, registration, label support, and target-field coverage before any additional patient-level analysis.
4. Keep the observational MGMT analysis separate from the treatment-policy simulator and avoid causal interpretation.
5. Build a reproducible Track C benchmark with shared virtual patients, common random inputs for paired arms, separate development and evaluation worlds, recorded seeds, and a pre-run code and settings seal.
6. Connect the tracks only with explicitly simulated proxy signals; never match people across datasets or treat A5 AURC, MRI recurrence, MGMT association, or synthetic signal parameters as measured patient-specific resistance.

### Expected outcomes

The Track A procedures will estimate donor-level regional prediction and, separately, an exploratory molecular resistance-associated proxy score. The Track B procedure will describe MRI-label transitions; verified treatment chronology may permit a limited radiographic treatment-associated proxy analysis. Track C will separate observational MGMT associations from simulated resistant-cell dynamics and treatment suppression. The proposed policy simulation may favor the dependence-aware controller, show similar performance, or favor a comparator in some scenarios. None of these analyses will directly measure TMZ sensitivity in the current cohorts. This plan reports no observed outcome values or findings; any eventual results will be reported separately and interpreted within their track and design.

## Materials

- **Computer:** Dell 16 Plus DB16250 with Windows 11 and 32 GB RAM. Processor, thread count, and runtime environment will be recorded for any locked run.
- **Software:** Python 3.14.4 and the package snapshot in [the Windows Python environment lock](../tracks/shared/requirements-win-py314.lock.txt). The snapshot records the current environment; a clean install on a second computer remains to be checked before a sealed evaluation.
- **Reproducibility files:** the current proposed settings are in [the integrated v2 settings](../configs/trackABC_integrated_v2.json); stream-role templates are in [the v2 seed-stream file](../configs/trackABC_seed_streams_v2.csv); the input-rights and checksum ledger is [the v2 input manifest](../configs/trackABC_input_manifest_v2.csv); and the reproduction sequence is in [the handoff guide](../tracks/shared/HANDOFF_AND_REPRODUCIBILITY.md). These files define a planning draft, not a frozen execution seal.
- **Track A sources:** Yu et al. single-cell RNA data and sample metadata (GEO GSE117891), Table S1 point labels, GBmap annotation metadata, Ivy GAP bulk RNA data (GSE107559), Angers microarray data (E-MTAB-5092), bulk expression data (E-MTAB-8743), and the published TMZ-sensitivity-associated gene set from Ntafoulis et al. (2023). Public metadata records screened for a possible external cohort include PRJNA1337938, HRA009021, GSE286419, GSE318562, GSE318564, and the Harwood Zenodo record 22013127. Screened records will not contribute expression data unless they pass the stated eligibility and rights rules.
- **Track B sources:** the restricted Duke longitudinal glioblastoma MRI record (Zenodo 20548528), MU-Glioma-Post and LUMIERE source records for earlier protocol development, synthetic images and masks for method-mechanics checks, and RANO 2.0 guidance for distinguishing radiographic progression from treatment-related imaging change. Zenodo lists CC BY-NC-ND 4.0 as an access option and describes a separate Duke custom-license route; the exact accepted terms for the local copy will govern its use. Duke files and identifiers will remain in access-controlled local storage, and derived outputs will not be shared unless the accepted terms permit it. MU and LUMIERE will not supply patient outcomes to the new integrated simulation.
- **Track C sources:** the project's two-compartment growth and treatment-suppression simulator; CGGA clinical and molecular records for the separate observational analysis; a TCGA-derived published reference used only as context; and published aggregate RTOG 0525 and GEINO 14-01 results for compatible model checks. Aggregate trial data will not be treated as patient-level dose-response records.
- **Physical materials:** none. The project uses computer files and simulations only; it uses no chemicals, bacteria, cell cultures, animals, or experimental drugs.

## Procedures

### Track A — expression and regional biology

The frozen A5 analysis used Yu et al. single-cell counts, exact donor/point joins, the 1,000-UMI and 500-gene cell filters, training-only preprocessing, and nested leave-one-donor-out validation across the six prespecified donors. It compared the sparse candidate with logistic regression, secondary baselines, ablations, and leakage diagnostics. A5's endpoint and success rule will remain unchanged. Separately, the project will score the 21-gene subset that Ntafoulis et al. reported as externally validated among its TMZ-sensitivity-associated genes only if the exact genes, directions, and mapping can be frozen before scoring; otherwise, it will omit the proxy analysis rather than substitute another set. Donor-by-region pseudobulk scores will be compared within donor with expression-matched controls. Verified MGMT metadata may be summarized descriptively. Neither analysis will label a patient or donor as TMZ-sensitive or resistant. Deferred calibration, within-library, scANVI, and attention analyses will require a dated amendment. Any external transfer analysis will require at least six independent matched donors, valid region mapping, fit-for-purpose malignant-cell labels, clear reuse terms, and verified lack of prior-use overlap.

### Track B — longitudinal MRI and spatial transitions

Synthetic development built registration, support, forecasting, and comparator components. The exploratory Duke analysis used prior-visit FLAIR, excluded tissue within 10 mm of the middle-visit core, and summarized subsequent core-label transitions in three distance bands; target masks were read, so this cohort and endpoint are exposed. The workflow will be independently audited for provenance, registration, support, and target-field coverage. A TMZ-associated interpretation will require verified patient-level treatment and progression chronology and applicable progression confirmation; otherwise the endpoint will remain MRI-defined recurrence or core-label transition. No retuning of the exposed cohort will be used to seek a preferred result. ADC remains a separate, incomplete workflow requiring its own read-scope, preprocessing, QC, and analysis freeze. Any future cohort must pass access, label, registration, support, independence, and sample-count gates before target labels are opened.

### Track C — MGMT association and treatment-policy modeling

The project conducted an observational CGGA Cox-analysis workflow for overall survival, including covariate-adjusted MGMT models and an MGMT-by-temozolomide interaction; results will be interpreted as associations, not causal effects. The mathematical simulator represents sensitive and resistant tumor compartments, treatment suppression, fitness costs, and switching. It includes fixed and adaptive policies and simulated molecular-like and MRI-like signals. Before evaluation, the project will freeze model ranges, response families, observation and progression rules, policy information, treatment-course limits, cohort sizes, seeds, uncertainty method, runtime limit, and a complete research-scale driver. Previously exposed seeds will not be used for evaluation.

### Integrated Track A-to-B-to-C procedure

Track A will motivate a synthetic molecular resistance-associated proxy signal, and Track B will motivate a synthetic longitudinal MRI progression signal. Their patient records will not be linked, and neither track will calibrate the other's accuracy or the simulator's cross-signal dependence. The simulator will generate both signals from a shared latent resistant state under an explicit quality, bias, missingness, and dependence grid. Dependence-aware and matched conditional-independence controllers, plus no-signal and single-signal controls, will receive identical virtual patients, treatments, scan schedules, and random inputs. Disjoint development, calibration, and evaluation worlds will be used. The outcome definitions, utility, precision target, seeds, and code/configuration seal will be fixed before evaluation. Any positive result will support only a simulated policy comparison under the stated assumptions, not measured patient resistance or clinical benefit.

## Risk Assessment

Tracks A and B use existing human-derived research data, including access-controlled Duke MRI and segmentation files. No participants will be recruited or contacted, and no treatment or laboratory drug assay will be administered. The completed Duke procedure included patient-mask reads. Any further access, analysis, retention, or sharing will follow the exact accepted data-use terms and school or competition conditions. Zenodo lists CC BY-NC-ND 4.0 as an access option and describes a separate Duke custom-license route; derived outputs will remain local unless the accepted terms permit sharing. Any unresolved permission or data-quality issue will stop the affected analysis. The project does not directly measure patient-specific TMZ sensitivity in the current cohorts. Track C's treatment and resistance variables are mathematical assumptions, not medical instructions. The integrated phase uses synthetic virtual patients and will not recommend treatment for an individual. Long computations will be run on a ventilated surface with the computer connected to power.


## Data Analysis

### Track A

The primary classifier comparison uses class-balanced AURC, averaged equally across the six donors; lower is better. The frozen descriptive success rule requires both lower mean AURC and lower AURC in at least five of six donors. Balanced accuracy and Brier score are secondary. Per-donor values will be shown in the final results report. The separate TMZ-response-associated expression score will be summarized as six donor-paired Core-minus-Periphery differences with expression-matched controls; it will not be assigned sensitivity/resistance accuracy or interpreted as a drug-response test. Other signature analyses use paired donor-level scores, expression-matched gene-set controls, exact sign-flip tests, and the prespecified covariate adjustments. Cells will not be treated as independent donors.

### Track B

For the completed exploratory Duke analysis, the unit of summary is the patient, with the saved transition definition and three distance bands. Any audit will first verify the source-to-analysis mapping, registration, support, and target-field coverage. An independent future transition analysis will compare new and persistent history-defined FLAIR regions within patient and physical-distance band, using a patient-level uncertainty method frozen before target-label access. Forecasting comparisons will use the same patients, support, target, and time interval for persistence, fixed-growth, PDE-only, hybrid, and learned baselines. A TMZ-associated interpretation will require verified exposure and clinical chronology; otherwise the endpoint will remain MRI-defined recurrence or core-label transition. ADC will remain a separate secondary analysis with its own completeness and missingness rules. No outcome values are included here.

### Track C

The observational MGMT analysis uses Cox proportional-hazards models for overall survival, with complete-case reporting and the declared clinical covariates. MGMT analyses among temozolomide-treated records and the MGMT-by-temozolomide interaction will be reported as noncausal associations. In the simulator, drug suppression acts on modeled sensitive and resistant compartments; resistant fraction, drug sensitivity, and proxy-signal quality remain simulation inputs rather than patient estimates. The synthetic policy comparison uses paired policy differences within each independently generated world, averages virtual patients within a world, and weights worlds equally. It reports both latent progression measures, treatment courses, censoring, model failures, and sensitivity scenarios. The independent world, not the nested patient, scan, or noise draw, is the uncertainty unit.

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

Hegi, M. E., Diserens, A.-C., Gorlia, T., Hamou, M.-F., de Tribolet, N., Weller, M., Kros, J. M., Hainfellner, J. A., Mason, W., Mariani, L., Bromberg, J. E. C., Hau, P., Mirimanoff, R. O., Cairncross, J. G., Janzer, R. C., & Stupp, R. (2005). MGMT gene silencing and benefit from temozolomide in glioblastoma. *New England Journal of Medicine, 352*(10), 997–1003. https://doi.org/10.1056/NEJMoa043331

Liu, Q., Fuster-Garcia, E., Hovden, I. T., MacIntosh, B. J., Grødem, E., Brandal, P., Lopez-Mateu, C., Sederevicius, D., Skogen, K., Schellhorn, T., Bjørnerud, A., & Emblem, K. E. (2025). Treatment-aware diffusion probabilistic model for longitudinal MRI generation and diffuse glioma growth prediction. *IEEE Transactions on Medical Imaging*. https://doi.org/10.1109/TMI.2025.3533038

Liu, Z., Zhang, J., Hong, L., Nie, Q., & Sun, X. (2025). Multiscale mathematical model-informed reinforcement learning optimizes combination treatment scheduling in glioblastoma evolution. *Science Advances, 11*(32), eadv3316. https://doi.org/10.1126/sciadv.adv3316

Ntafoulis, I., Kleijn, A., Ju, J., Jimenez-Cowell, K., Fabro, F., Klein, M., Yen, R. T. C., Balvers, R. K., Li, Y., Stubbs, A. P., Kers, T. V., Kros, J. M., Lawler, S. E., Beerepoot, L. V., Kremer, A., Idbaih, A., Verreault, M., Byrne, A. T., O’Farrell, A. C., ... Lamfers, M. L. (2023). Ex vivo drug sensitivity screening predicts response to temozolomide in glioblastoma patients and identifies candidate biomarkers. *British Journal of Cancer, 129*(8), 1327–1338. https://doi.org/10.1038/s41416-023-02402-y

Stupp, R., Mason, W. P., van den Bent, M. J., Weller, M., Fisher, B., Taphoorn, M. J. B., Belanger, K., Brandes, A. A., Marosi, C., Bogdahn, U., Curschmann, J., Janzer, R. C., Ludwin, S. K., Gorlia, T., Allgeier, A., Lacombe, D., Cairncross, J. G., Eisenhauer, E., & Mirimanoff, R. O. (2005). Radiotherapy plus concomitant and adjuvant temozolomide for glioblastoma. *New England Journal of Medicine, 352*(10), 987–996. https://doi.org/10.1056/NEJMoa043330

Wen, P. Y., van den Bent, M., Youssef, G., Cloughesy, T. F., Ellingson, B. M., Weller, M., Galanis, E., Barboriak, D. P., de Groot, J., Gilbert, M. R., Huang, R., Lassman, A. B., Mehta, M., Molinaro, A. M., Preusser, M., Rahman, R., Shankar, L. K., Stupp, R., Villanueva-Meyer, J. E., ... Chang, S. M. (2023). RANO 2.0: Update to the response assessment in neuro-oncology criteria for high- and low-grade gliomas in adults. *Journal of Clinical Oncology, 41*(33), 5187–5199. https://doi.org/10.1200/JCO.23.01059

Yu, K., Hu, Y., Wu, F., Guo, Q., Qian, Z., Hu, W., Chen, J., Wang, K., Fan, X., Wu, X., Rasko, J. E., Fan, X., Iavarone, A., Jiang, T., Tang, F., & Su, X.-D. (2020). Surveying brain tumor heterogeneity by single-cell RNA-sequencing of multi-sector biopsies. *National Science Review, 7*(8), 1306–1318. https://doi.org/10.1093/nsr/nwaa099

Zimmer, L., Weidner, J., Balcerak, M., Kofler, F., Krupa, M., Ezhov, I., Cepeda, S., Zhang, R., Lowengrub, J., Menze, B., & Wiestler, B. (2025). *PREDICT-GBM: A multi-center platform to advance personalized glioblastoma radiotherapy planning* [Preprint]. arXiv. https://arxiv.org/abs/2509.13360

## Additional Project Information

This is a continuation project. The Procedures section identifies completed method steps and remaining planned steps for Tracks A, B, and C. This research-plan document contains no observed study-result values or outcome conclusions; any required continuation abstract or results report will be supplied separately under the competition's instructions. Track A and Track B data are not patient-linked, and neither provides direct TMZ-response measurements for the analyzed patients. The molecular and MRI signals in the integrated model are resistance-associated proxies or explicitly simulated observations. Track C treatment policies are simulated and are not clinical recommendations.
