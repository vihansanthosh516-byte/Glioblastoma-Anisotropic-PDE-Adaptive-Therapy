# Research Plan

## Project Title

**Integrating Tumor Biology and Longitudinal MRI for Glioblastoma Progression Modeling and In-Silico Therapy Optimization**

## Rationale

Glioblastoma (GBM) varies across tumor regions and changes over time. Single-cell gene-expression data describe molecular patterns in sampled tissue (Yu et al., 2020); serial magnetic resonance imaging (MRI) describes radiographic changes; and mathematical models can represent tumor growth, treatment suppression, and resistant-cell dynamics. Track A uses mixed-cell expression and author-supplied MRI-marked sampling-point labels, so its prediction is not a malignant-cell-only or drug-response test. Each source measures a different part of the problem. The available molecular and MRI cohorts are not linked patient by patient, and neither provides a shared, direct measurement of temozolomide (TMZ) sensitivity. The project will therefore keep the patient-data analyses separate and connect them only through explicitly simulated signals.

Prior studies have mapped MRI patterns associated with later recurrence, built patient-specific tumor and radiotherapy models, predicted longitudinal MRI, and optimized GBM treatment schedules with reinforcement learning (Akbari et al., 2016; Chaudhuri et al., 2023; Ebrahimi Zade et al., 2020; Q. Liu et al., 2025; Z. Liu et al., 2025; Zimmer et al., 2025). MGMT promoter methylation has clinical evidence as a TMZ-associated biomarker, while ex-vivo assays measure drug response more directly than expression or MRI proxies (Hegi et al., 2005; Ntafoulis et al., 2023). The project will not claim that multimodal GBM modeling or treatment optimization is new by itself. Its focused integrated question will test whether explicitly modeling residual dependence between two imperfect, simulated resistance-associated signals changes robust treatment-policy utility compared with a matched model that assumes conditional independence.

This work is important as a methods study: it will show how regional expression analysis, longitudinal MRI analysis, and mathematical treatment simulation can be compared without treating unrelated patient cohorts as one dataset. The integrated comparison may favor the dependence-aware policy, show little difference, or favor the comparator. Any conclusion will apply only to the stated data analyses or simulated assumptions. The work will not select treatment for a patient or establish clinical benefit.

## Research Questions, Hypotheses, Engineering Goals, and Expected Outcomes

### Track A — regional expression and prediction

**Research question A1.** Under the frozen six-donor nested leave-one-donor-out design, does the sparse donor-balanced expression classifier have lower equal-donor class-balanced area under the risk-coverage curve (AURC) than ordinary L2 logistic regression, with lower AURC in at least five of the six held-out donors?

**Hypothesis A1.** The candidate classifier is expected to have lower equal-donor mean AURC than L2 logistic regression and lower donor-specific AURC in at least five of six donors. This is the predeclared hypothesis and success rule for the completed A5 procedure; it will not be changed using the observed results.

**Research question A2.** In paired-donor analyses, do the prespecified senescence, oxidative-phosphorylation, and tumor-associated-macrophage expression programs show regional patterns that remain interpretable with expression-matched controls and the declared covariate adjustments?

**Research question A3.** Do the paired donors show a consistent Core-versus-Periphery difference in the fixed 119-of-121-gene expression score based on published in-vitro TMZ-viability associations? This exploratory score has no clinical sensitivity or resistance threshold.

### Track B — longitudinal MRI and tissue-state transitions

**Research question B1.** Within each patient and physical-distance band, how is next-visit core-label transition related to history-defined newly appearing versus persistent FLAIR-labeled tissue?

**Exploratory hypothesis B1.** Newly appearing FLAIR-labeled tissue is expected to have a higher next-visit core-label transition rate than persistent FLAIR-labeled tissue in the same patient and distance band. The previously analyzed Duke cohort and endpoint are exposed and will support exploratory description only, not independent confirmation.

**Research question B2.** In a separate eligible cohort, does a prespecified longitudinal MRI forecast improve on history-only mask baselines when all methods use the same patients, spatial support, target labels, and forecast intervals?

### Track C — MGMT association and simulated treatment policies

**Research question C1.** Among eligible CGGA GBM records, how is MGMT methylation associated with overall survival (OS) after adjustment for age, IDH status, radiotherapy, recorded TMZ status, recurrence status, and data batch? Does the association differ by recorded TMZ status?

**Hypothesis C1.** MGMT methylation may be associated with OS among records marked as TMZ-treated. The analysis is retrospective and observational; it will not estimate a causal treatment effect.

**Research question C2.** Across independently generated synthetic GBM worlds, does a finite-set dependence-aware minimax controller achieve higher robust utility than a matched controller whose residual signal-dependence parameter is fixed at zero, when both receive the same virtual patients, signals, treatment choices, and observation schedules?

**Hypothesis C2.** In correctly specified scenarios with informative signals and nonzero residual dependence, the dependence-aware controller is expected to have higher equal-world mean robust utility. In predeclared uninformative and reversed-signal controls, it is not expected to show a consistent advantage over no-signal or single-signal controls.

### Integrated question

When molecular-like and MRI-like proxy signals are generated jointly from the same simulated resistant-tumor state, does explicitly modeling their residual dependence improve robust simulated treatment utility over a matched conditional-independence policy?

### Engineering goals and expected outcomes

1. Preserve the original Track A folds, preprocessing rules, comparisons, and success rule.
2. Keep Track A, Track B, the CGGA analysis, and simulation results separate by cohort, endpoint, and interpretation.
3. Make any permitted Track B analysis reproducible from source files, registration, labels, spatial support, and target-field coverage.
4. Complete and independently review a reproducible Track C benchmark with shared virtual patients and inputs, separate development and evaluation worlds, recorded seeds, and a pre-run code and settings seal.
5. Validate a simultaneous uncertainty method and a feasible precision and runtime plan before generating research-scale policy outcomes.
6. Report whether the simulated policy comparison supports the hypothesis, shows little difference, or favors the comparator. A positive result is not assumed.

Expected outputs are a reproducible account of the Track A procedures, an appropriately bounded Track B transition analysis, a separate observational CGGA report, and a synthetic policy comparison only if its statistical and execution gates pass. None will be presented as a clinical recommendation or a direct measurement of patient-specific TMZ response.

## Materials

- **Computer:** Dell 16 Plus DB16250, Windows 11, and 32 GB RAM. Processor, thread settings, software environment, and runtime will be recorded for future sealed runs.
- **Software:** Python 3.14.4 and the package versions in the project's Windows environment lock file. The integrated v2 settings, seed-stream template, input-rights manifest, and reproduction guide will be used as planning and audit records; they are not, by themselves, a frozen execution seal.
- **Track A data:** Yu et al. single-cell RNA counts and sample metadata (GEO GSE117891), Table S1 sampling-point labels, GBmap cell-type annotation metadata, Ivy GAP bulk RNA data (GSE107559), Angers microarray data (E-MTAB-5092), E-MTAB-8743 bulk expression data, and the Ntafoulis et al. published in-vitro TMZ-viability association panel. Public records screened for external-cohort eligibility will not contribute expression data unless they pass the documented donor, label, independence, and reuse checks.
- **Track B data:** MU-Glioma-Post and LUMIERE records used in prior method development; the restricted Duke longitudinal MRI dataset (Zenodo record 20548528); synthetic images and masks for mechanics checks; and RANO 2.0 guidance. Duke files and identifiers will remain access-controlled. The exact accepted data-use terms will govern analysis and sharing of every derived output.
- **Track C observational data:** the two CGGA clinical-data batches and a published TCGA-derived age reference used only as context.
- **Track C simulation materials:** the two-compartment simulator, policy implementations, versioned configuration, seed-stream templates, and synthetic inputs. No patient-level data will enter the integrated simulator. Published aggregate RTOG 0525 and GEINO 14-01 results may be used only for endpoint-compatible model checks, not as patient-level dose-response records.
- **Physical materials:** none. No chemicals, bacteria, cell cultures, animals, experimental drugs, or participant recruitment will be used.

## Procedures

The procedure records completed analysis steps and identifies work that remains planned or conditional. A completed procedure describes methods already performed; it does not imply that its hypothesis was supported. No study outcome values are presented in this plan.

Before any new patient analysis, run synthetic checks for file handling, registration and support mechanics, split isolation, target leakage, missing-data handling, failure fallbacks, and reproducibility. Record the tested code and settings. Synthetic checks verify implementation behavior; they do not establish patient-label accuracy, clinical validity, or treatment benefit.

### Track A — expression and regional biology

#### Completed A5 classifier procedure

1. Use the six protocol-designated donors: GS1, GS2, GS3, GS5, GS6, and GS13. Verify the frozen checksums and exact joins among cells, points, donors, and Table S1.
2. Define Tumoral points as Core and Peritumoral points as Periphery. Exclude Normal points and invalid joins. Stop if any primary donor has fewer than 20 eligible cells in either class.
3. Retain cells with at least 1,000 UMIs and 500 detected genes. In each training split only, retain genes detected in at least 5% of training cells, select the 2,000 most variable genes, normalize each cell to 10,000 counts, apply log1p, and fit feature means and standard deviations. Apply the fitted transformations to validation and held-out cells.
4. Use six outer leave-one-donor-out folds and inner leave-one-donor-out folds formed from the outer-training donors; use the frozen seed 20261031. Compare the sparse donor-balanced candidate with L2 logistic regression, 15-nearest-neighbor voting, linear SVM, and gradient-boosted trees. Retain the three candidate ablations and cell-type-only and cells-per-point diagnostics; identify constant-score controls as such.
5. Use equal-donor class-balanced AURC as the primary metric. The frozen descriptive success rule requires both a lower donor-mean AURC than logistic regression and lower AURC in at least five of six donors. Balanced accuracy and Brier score are secondary. Do not retune the model or rule using the held-out results.

#### Completed regional-program and response-associated score procedures

1. Analyze SenMayo and oxidative-phosphorylation programs in Ivy GAP and Angers as separate paired-donor analyses; analyze the tumor-associated-macrophage program in E-MTAB-8743 separately. Use donor-level summaries, expression-matched gene-set controls, exact donor sign-flip procedures, and the declared covariate analyses. Do not treat cells or probes as independent donors.
2. For the exploratory TMZ-response-associated expression proxy, use the 121-gene `Correlation_cell_Overlap` panel from Ntafoulis et al. (2023), with the original Spearman directions against in-vitro viability at 100 μM TMZ. Use only the 119 genes with exact symbols in the A5 matrix; omit `AF131215.2` and `PPIAP22` without aliases or substitutions.
3. Sum raw UMIs into donor-by-region pseudobulks, calculate transcriptome-wide expression percentiles, and subtract the mean percentile of the 72 negative-direction genes from the mean percentile of the 47 positive-direction genes. Summarize all paired Core-minus-Periphery differences against 1,000 expression-matched control panels using an exploratory donor-level exact sign-flip calculation. Interpret this only as alignment with published in-vitro viability directions, not as a drug-response test.
4. Use MGMT metadata only if its sample join, assay, and provenance are verified. Any future external-transfer analysis requires a separate dated amendment and an independent cohort with valid donor-region mapping, fit-for-purpose malignant-cell labels, clear reuse terms, adequate matched donors, and no prior-use overlap.

### Track B — longitudinal MRI and spatial transitions

#### Completed exploratory analyses

1. The project conducted an exploratory PDE forecast analysis using MU-Glioma-Post. It will be reported separately and will not serve as fresh validation.
2. The project conducted an exploratory Duke transition analysis on a restricted longitudinal MRI cohort. The recorded procedure expanded prior-visit FLAIR by 10 mm, excluded tissue within 10 mm of the middle-visit core, and summarized subsequent core-label transitions in three physical-distance bands. Target masks were read; therefore, this cohort and endpoint are exposed and will not be described as an untouched test set.
3. Keep the later 2-mm transition-estimand prototype distinct from the executed 10-mm procedure; it is not a reproduction. Do not retune the exposed Duke analysis to seek a preferred result.

#### Planned or conditional Track B work

1. Audit source-to-analysis mapping, image geometry, registration, spatial support, and target-field coverage before further patient-level interpretation. Retain an explicit coverage mask so out-of-field tissue is not counted as observed non-conversion.
2. For any permitted transition analysis, define history-based tissue states before target-label scoring. Compare next-visit transitions for newly appearing and persistent FLAIR-labeled tissue within patient and physical-distance band. Treat the patient as the analysis unit and freeze the patient-level uncertainty method before opening target labels.
3. A forecast comparison will use a separate eligible cohort. Compare persistence, fixed growth, PDE-only, hybrid, 3D U-Net, time-conditioned FiLM U-Net, and count-only baselines on the same patients, support, target labels, and forecast interval. Fit every learned or calibrated component using training patients only. Synthetic mechanics checks will not be treated as patient-image validation.
4. Call an endpoint TMZ-associated only if permitted patient-level treatment exposure, scan chronology, and relevant clinical context are verified. Otherwise describe MRI-defined recurrence or core-label transition. Apply RANO 2.0 confirmation only when the available clinical and scan records support it; a simulated scan-sampled burden threshold is not RANO.
5. Keep ADC as a separate analysis with its own source definitions, read scope, preprocessing, quality checks, and missingness rules. Before a separate cohort is used, verify access rights, longitudinal visit mapping, labels, registration QC, spatial support, prior-use overlap, and sample eligibility. Freeze the registration method and QC thresholds, support construction, target-coverage handling, failure rules, and interval definition before opening target labels for a new analysis.

### Track C — CGGA association and tumor-treatment model

#### Completed observational procedure

Use the two CGGA clinical-data batches, remove duplicate patient identifiers under the recorded rule, apply the declared GBM and survival eligibility criteria, and include a batch indicator. Fit Cox proportional-hazards models for OS with age, IDH status, radiotherapy, recorded TMZ status, recurrence status, and MGMT methylation. Analyze the MGMT association among TMZ-treated records and the MGMT-by-TMZ interaction in eligible pooled records. Report complete-case counts and uncertainty, check proportional-hazards assumptions, and interpret all estimates as noncausal associations. The available treatment fields do not represent full dose and timing histories.

#### Planned synthetic tumor and treatment-policy procedure

1. Represent sensitive burden as `S`, resistant burden as `R`, and total burden as `B = S + R`. Use the declared two-compartment model:

   `dS/dt = rS[1 − (S + a_rsR)] − mruS − μuS + ν(1 − u)R`

   `dR/dt = r(1 − c)R[1 − (a_srS + R)] − εmruR + μuS − ν(1 − u)R`

   Set `S(0)=(1−f0)B0`, `R(0)=f0B0`, and `a_rs=a_sr=1`. The terms represent growth, competition, treatment suppression, resistance cost, and switching. Model values are synthetic inputs, not patient-specific estimates.
2. Integrate using fourth-order Runge–Kutta with a 0.25-day step and compare with a 0.05-day step for numerical convergence. Keep single-response and mixture-response simulation families separate.
3. Generate independent synthetic worlds and separate selection, calibration, and evaluation virtual-patient cohorts. Assign registered random-number streams to world parameters, cohorts, signal generation, policies, and uncertainty calculations. Use common virtual patients, latent parameters, signals, scan schedules, and observation-noise inputs for paired policy arms.
4. Transform simulated resistant fraction to a latent normal score `zR`. Generate a molecular-like latent score `zB = a·zR + √(1−a²)·EB` and an MRI-like score `zI = b·zR + √(1−b²)·(d·EB + √(1−d²)·EI)`, where `EB` and `EI` are independent standard-normal draws. Threshold `zB` to a binary report and flip it with probability `1−q`, where `q` is the report's simulated accuracy. The current planned grid uses `a`, `b`, and `d` values `{−0.5, 0, 0.5}` with primary anchor `(0.5, 0.5, 0.5)`; `q` is drawn once per primary world from `{0.6, 0.7, 0.8, 0.9}`. Use `q=0.5` for the predeclared null controls, including the full no-information setting `(q,a,b,d)=(0.5,0,0,0)`. Reconcile these proposed settings with the executable configuration and freeze them before evaluation. These are simulation assumptions; no patient expression or MRI values will estimate their joint dependence.
5. Compare the finite-set dependence-aware minimax controller with the matched controller fixed to residual dependence `d = 0`. Include no-signal, single-signal, and reversed-signal controls; treat oracle policies as diagnostics only. At each replan, evaluate eight treatment/no-treatment patterns over three 28-day cycles, forecast six cycles without terminal value, apply the first block, and replan at the next scheduled scan. A modeled course lasts five days, with a 12-course cap. The primary scan interval is 84 days; 28-, 56-, and 112-day intervals and visible scan loss are sensitivity scenarios. The planned observable stop proxy requires two consecutive visible noisy readings at or above 1.4 times the lowest earlier observed burden, with a floor at 0.10 times observed starting burden; a missing scan resets confirmation. This is not formal RANO. Give paired policies the same treatment options, virtual-patient inputs, and observation schedule.
6. Before any research-scale policy outcomes are generated, freeze the complete driver, cohort and world counts, seed schedule, model and endpoint rules, uncertainty method, precision plan, runtime cap, package environment, and code/configuration seal. The research-scale comparison is currently stopped; it will resume only if an uncertainty method passes independent outcome-blind coverage validation and the design can meet its precision target at feasible world counts and runtime. If those gates cannot be met, report feasibility or an inconclusive analysis, not a positive policy comparison.

### Integrated Track A-to-B-to-C procedure

1. Use Track A to define the biological interpretation of a simulated molecular-like proxy and Track B to define a simulated longitudinal MRI-like proxy. Keep their patient records separate; do not link patients or estimate joint signal accuracy from unrelated cohorts.
2. In each virtual world, derive both signals from the same latent resistant-tumor state, with declared error and dependence parameters. Compare dependence-aware and conditional-independence controllers under identical synthetic inputs. Add no-signal and single-signal controls to determine whether any apparent gain requires both channels.
3. Define the first progression measure as interpolated time when latent burden reaches 1.2 times starting burden. Define the second as the first scheduled scan when latent burden reaches 1.4 times the larger of the prior latent-burden nadir or 0.10 times starting burden. The second is a simulated scan-sampled proxy, not formal RANO.
4. Define utility as the smaller progression time, each capped at 730 days, minus five simulated days per treatment course delivered through day 730. Run course-cost sensitivities at zero and ten simulated days per course. The primary contrast is candidate utility minus matched `d = 0` utility, paired within each independent world and averaged equally across worlds.
5. Keep development, calibration, and evaluation worlds disjoint. Freeze all outcome definitions, policy comparisons, world counts, seeds, uncertainty procedures, and code/settings hashes before evaluation. Preserve null and comparator-favoring results as well as results that support the hypothesis.

## Risk Assessment

Tracks A and B use previously collected human-derived research data. No participants will be recruited or contacted, and no treatment, survey, or laboratory drug assay will be administered. Restricted Duke files and identifiers will remain in access-controlled storage. All analysis, retention, and sharing will follow the exact accepted data-use terms and applicable school or competition requirements. Derived Duke outputs will not be shared unless the accepted terms explicitly allow it. Any unresolved permission, provenance, registration, coverage, or data-quality issue will stop the affected analysis. This project does not directly measure patient-specific TMZ sensitivity. Track C resistance and treatment variables are mathematical assumptions, not medical instructions. The integrated phase uses synthetic virtual patients and will not recommend treatment for an individual. Long computations will be run on a ventilated surface with the computer connected to power.

## Data Analysis

### Track A

Use equal-donor mean class-balanced AURC as the primary A5 metric, with lower values indicating better risk-coverage performance. Apply the frozen descriptive success rule: the candidate must have both lower mean AURC than L2 logistic regression and lower donor-specific AURC in at least five of six donors. Report balanced accuracy and Brier score as secondary measures. For the separate expression analyses, report donor-level paired scores, matched-gene controls, exact sign-flip calculations, and prespecified covariate analyses. Do not treat cells as independent donors. The 119-of-121 score has no sensitivity/resistance accuracy threshold and is not a direct drug-response measurement.

### Track B

Use the patient as the unit for transition summaries and physical-distance bands. Check source mapping, registration, support, and target-field coverage before interpreting any new analysis. Freeze any patient-level uncertainty method before target-label access. Forecast methods must share the eligible patients, spatial support, target, and interval. Use a TMZ-associated description only when exposure and chronology are verified; otherwise use MRI-defined recurrence or core-label transition.

### Track C observational

Report Cox hazard ratios, confidence intervals, p-values, complete-case counts, covariate definitions, batch handling, and proportional-hazards diagnostics. Interpret MGMT estimates and interaction tests as noncausal associations with OS.

### Integrated simulation

Report both latent progression measures, treatment courses, censoring, failures, missing arms, fallbacks, and prespecified sensitivity scenarios. The independent simulated world is the uncertainty unit; patients, scans, and noise draws within a world will not be treated as independent worlds. The primary comparison will keep single-response and mixture-response families separate. Its positive rule will require both simultaneous adjusted lower confidence bounds to exceed zero and both 95% interval half-widths to be no greater than 14 simulated days. The interval method has not passed the required validation gate; no policy-outcome evaluation will proceed until a replacement method is independently reviewed and passes outcome-blind coverage and precision checks. If the precision target is infeasible, the work will be reported as feasibility-only or inconclusive without changing the rule after outcome access.

## Bibliography

Akbari, H., Macyszyn, L., Da, X., Bilello, M., Wolf, R. L., Martinez-Lage, M., Biros, G., Alonso-Basanta, M., O’Rourke, D. M., & Davatzikos, C. (2016). Imaging surrogates of infiltration obtained via multiparametric imaging pattern analysis predict subsequent location of recurrence of glioblastoma. *Neurosurgery, 78*(4), 572–580. https://doi.org/10.1227/NEU.0000000000001202

Balaña, C., Vaz, M. A., Sepúlveda, J. M., Mesia, C., del Barco, S., Pineda, E., Muñoz-Langa, J., Estival, A., de las Peñas, R., Fuster, J., Gironès, R., Navarro, L. M., Gil-Gil, M., Alonso, M., Herrero, A., Peralta, S., Olier, C., Perez-Segura, P., Covela, M., & Carrato, C. (2020). A phase II randomized, multicenter, open-label trial of continuing adjuvant temozolomide beyond 6 cycles in patients with glioblastoma (GEINO 14-01). *Neuro-Oncology, 22*(12), 1851–1861. https://doi.org/10.1093/neuonc/noaa107

Bhattacharya, M., Singh, G., & Prasanna, P. (2025). *SoC-DT: Standard-of-care aligned digital twins for patient-specific tumor dynamics* [Preprint]. arXiv. https://arxiv.org/abs/2510.03287

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

This is a continuation project using previously collected research data and computer simulations. The Procedures section distinguishes methods already performed from planned or conditional work. No new participants will be recruited or contacted, and no treatment or laboratory drug assay will be administered. Track A and Track B records are not linked at the patient level. Neither provides a direct TMZ-response measurement for the integrated model. The expression score and MRI transitions are proxies; resistant-cell fractions, signal properties, and treatment effects in Track C are simulated assumptions. The treatment-policy comparison is not a clinical recommendation. Any competition-required continuation abstract or outcome report will be supplied separately and will report the study findings accurately.
