# Integrating Tumor Biology and Longitudinal MRI for Glioblastoma Progression Modeling and In-Silico Therapy Optimization

**Research Plan and Current Project Record**

**Project year:** 2026–2027
**Plan version:** 2026-11-08

## Rationale

### Project summary

This project studies glioblastoma (GBM) through three separate research tracks and a synthetic treatment-policy simulation. Track A analyzes regional gene expression and donor-held-out prediction. Track B analyzes longitudinal MRI labels and spatial transitions. Track C includes an observational analysis of MGMT methylation and overall survival, plus mathematical simulation of tumor growth and treatment policies. The integrated study uses synthetic virtual patients to test whether modeling dependence between two imperfect signals changes simulated policy utility.

Track A, Track B, and Track C do not share a verified patient-level linkage or a common measured drug-response outcome. The integrated simulation therefore uses no linked patient records. Molecular-like and MRI-like signals are generated from declared assumptions. Simulation outcomes are not patient outcomes and do not establish clinical treatment benefit.

The project aims to identify conditions in which a dependence-aware controller improves robust simulated utility over a matched conditional-independence controller. Its positive-result rule is specified below. The rule will be applied only after the inference procedure is independently reviewed and passes outcome-blind synthetic validation. A positive result is an empirical outcome to test, not a guaranteed result.

### Scientific scope and interpretation

1. **Track A** measures expression patterns in sampled tumor regions. The A5 labels are author-supplied MRI-marked sampling-point labels, and the analyzed expression is mixed-cell. A5 performance is not resistance-prediction accuracy.
2. **Track B** measures MRI-defined tissue and label transitions. A transition or recurrence label is not a direct measure of temozolomide (TMZ) resistance. A treatment-associated interpretation requires verified patient-level exposure and clinical chronology; even then, it is a radiographic proxy.
3. **Track C observational analysis** estimates retrospective associations in clinical records. MGMT-by-TMZ associations are not randomized treatment effects.
4. **Track C simulation and the integrated study** use assumed latent tumor states, treatment effects, signal quality, and signal dependence. Their outcomes apply to the specified simulated worlds.
5. The project will not label a donor or patient TMZ-sensitive or TMZ-resistant from expression or MRI alone. It will not recommend treatment for a person or claim clinical benefit from simulated policy differences.

### Research record to date

These results describe analyses already recorded for this project. They do not change the project’s hypotheses or success rules.

### Track A

The A5 analysis used six strict donors, 4,139 matched target cells, six outer leave-one-donor-out folds, and nine models. The candidate’s mean class-balanced area under the risk-coverage curve (AURC) was 0.1478656 versus 0.1763646 for L2 logistic regression. The candidate had lower donor-specific AURC in four of six donors; the prespecified rule required at least five of six. Thus, the descriptive two-part success rule was not met, although the mean comparison favored the candidate. Mean balanced accuracy was 0.6786 versus 0.6621, while mean Brier score was 0.2335 versus 0.1944. No inferential superiority claim was made. The frozen A5 procedure and result remain unchanged. See the [A5 execution audit](trackA_execution_audit_20261031.md).

Separate regional-program analyses found unadjusted SenMayo and OXPHOS shifts. Their prespecified adjusted tests did not meet their criteria: SenMayo adjusted p=0.37 and OXPHOS adjusted p=0.67. The E-MTAB-8743 TAM-activation test also did not meet its adjusted rule. These results do not establish that the adjusted effects are absent or that measured covariates explain them. See the [Track A result record](trackA_null_result.md).

The separate TMZ-response-associated expression proxy used 119 of 121 exact-mapped genes from the Ntafoulis in-vitro viability-correlation panel, with no substitutions. Two source symbols were unmapped. Across six paired donors, Core-minus-Periphery score differences split three positive and three negative; the mean was −0.001086 and the exact two-sided donor sign-flip p-value was 0.96875. The score is a mixed-cell expression proxy, not a direct drug-response measurement or validated predictor. See the [project status record](../START_HERE.md) and [amendment](trackA_tmz_response_119of121_amendment_2026-11-08.md); the detailed score JSON is local and not part of the shared Git history.

No donor-linked MGMT assay and sample join is established in the reviewed records. The screened external datasets have not been designated an independent expression-validation cohort.

### Track B

An exploratory MU-Glioma-Post PDE analysis on 133 patients is recorded in an accessible project worktree. For the directional anisotropic r=10 arm, the all-pairs mean Dice difference was +0.0059 versus persistence (95% CI −0.0085 to +0.0202) and −0.0077 versus geometric train-rate (95% CI −0.0209 to +0.0055). The saved manifest did not identify one arm as the sole primary comparison. These patient-level results are exploratory. The numerical summary is recorded on Git branch `codex/duke-access-request-confirmed`, in its `START_HERE.md` MU result entry.

A private Duke analysis completed on 2026-11-04 used target-mask labels. Its saved procedure dilated prior-visit FLAIR by 10 mm, excluded tissue within 10 mm of the middle-visit core, and summarized transitions in three distance bands. The cohort and endpoint have been examined, so this analysis is exploratory evidence, not an untouched confirmation set. The current transition-estimand prototype follows an older 2-mm draft and is not a reproduction of the 10-mm run. The saved summary records out-of-field target voxels resampled as zero; the effect of target coverage on its results is unknown. No patient-level TMZ chronology has been verified for this endpoint. See the [Track B project record](../START_HERE.md) and [Track B handoff](../tracks/track-b/README.md).

Historical synthetic checks cover selected geometry, registration mechanics, masks, and coverage calculations. They do not validate patient registration, label accuracy, target coverage, or clinical meaning. The ADC preflight record is incomplete; its scope does not establish whether every possible ADC analysis occurred.

### Track C

The CGGA observational analysis pooled the two released clinical batches after the recorded duplicate and eligibility rules. Among TMZ-treated records, the adjusted MGMT estimate was HR 0.7967 (95% CI 0.6043–1.0503; n=250; p=0.107). In the pooled MGMT-by-TMZ interaction model, the interaction estimate was HR 0.6026 (95% CI 0.3298–1.1009; n=306; p=0.100). The endpoint is overall survival (OS), not progression-free survival; treatment is recorded as yes/no without dose or timing; the analysis is retrospective and complete-case. These estimates do not show a causal treatment effect. See the [CGGA analysis script](../src/118_cgga_adjusted_age_mgmt.py) and [saved result](../output/cgga_adjusted_age_mgmt.json).

The project has simulator, particle-filter, controller, signal-likelihood, and world-orchestration prototypes, plus historical synthetic development outcomes. Earlier Stage 2 runs are development evidence, not a fresh evaluation of the current v2 contrast. One prior trial-matching exercise did not reproduce clinical PFS because it omitted death and the time origin; it is not treated as model calibration or validation. See the [Track C development report](trackc_stage2_report_2026-11-04.md).

The current primary simulated contrast is a finite-set dependence-aware minimax controller versus a matched minimax controller with residual dependence parameter d=0. No research-scale policy outcome comparison has run under this final contrast; small synthetic contract tests are code checks, not study evidence. An outcome-blind synthetic screen of the proposed observed-SE max-T intervals completed on 2026-11-08 and failed its predeclared coverage criterion: all 24 cells missed the Wilson-lower-bound threshold; observed joint coverage ranged from 70.05% to 93.2% for nominal 95% intervals. The helper is not eligible for policy evaluation or a positive claim. See the [coverage-screen protocol](trackc_max_t_coverage_audit_protocol_2026-11-08.md) and [result](trackc_max_t_coverage_audit_result_2026-11-08.md).

### Prior art

GBM changes across tumor regions and over time. Single-cell expression can characterize sampled regional biology; longitudinal MRI can characterize radiographic change; mathematical models can represent tumor growth and treatment assumptions. These measurements answer different questions and are not interchangeable.

Existing studies establish substantial prior art in imaging-based recurrence mapping, patient-specific radiotherapy digital twins, longitudinal MRI prediction, and model-based treatment scheduling. Examples include imaging surrogates of infiltration for recurrence localization (Akbari et al., 2016), patient-specific radiotherapy modeling under uncertainty (Chaudhuri et al., 2023), treatment-aware longitudinal MRI prediction (Q. Liu et al., 2025), TMZ schedule optimization by reinforcement learning (Ebrahimi Zade et al., 2020), and model-informed reinforcement learning for GBM combination therapy (Z. Liu et al., 2025). The SoC-DT preprint describes reaction-diffusion growth with standard-of-care interventions and genomic or demographic personalization (Bhattacharya et al., 2025). MGMT methylation has clinical evidence as a TMZ-associated biomarker (Hegi et al., 2005), while ex-vivo drug-sensitivity assays measure cellular response more directly (Ntafoulis et al., 2023).

The 2026 PREDICT-GBM publication presents a multicenter platform and longitudinal dataset of 243 patients, with a standardized comparison of biophysical and data-driven recurrence-planning methods (Zimmer et al., 2026). Its reported endpoint is geometric recurrence coverage, not survival benefit. A 2026 prospective pilot studied machine-learning-guided radiation dose escalation; it is a radiotherapy feasibility and safety precedent, not validation of a TMZ treatment controller ("Personalized machine learning-guided radiation dose escalation in newly diagnosed glioblastoma," 2026). The 2025 PREDICT-GBM arXiv preprint and 2026 published paper have different titles and are cited separately below (Zimmer et al., 2025, 2026). The 2026 platform includes TUM-GBM, LUMIERE, and RHUH cohorts already examined by this project, so it is not an independent validation cohort here.

The project’s narrow question is whether uncertainty about dependence between two imperfect, simulated proxy signals changes treatment-policy utility in the declared synthetic worlds. The targeted prior-art review has not verified that this exact contrast is new, and it is not a systematic review. The project will avoid “first,” “unique,” or broad novelty claims.

### Expected contribution and limits

The project’s contribution is a carefully specified comparison of dependence-aware and conditional-independence treatment policies under declared synthetic assumptions, informed—but not calibrated—by separate molecular and imaging analyses. Its value depends on transparent assumptions, validated implementation, sound uncertainty estimates, and a comparison that meets its predeclared precision rule.

The project will not claim a new clinical treatment, a patient-level digital twin, validated patient-specific resistance, causal MGMT treatment benefit, or clinical benefit from simulated utility. A positive simulation result would support only the specified policy contrast in the specified synthetic worlds and would motivate later validation.

## Research Questions, Hypotheses, or Engineering Goals

### Research questions and hypotheses

### Track A — regional expression and prediction

**A1 research question.** Under the frozen six-donor nested leave-one-donor-out design, does the sparse donor-balanced classifier have lower equal-donor class-balanced AURC than L2 logistic regression, with lower AURC in at least five of six held-out donors?

**A1 success rule.** Both the mean-AURC comparison and the five-of-six donor rule must hold. Balanced accuracy and Brier score are secondary measures. The already completed A5 result is reported above under this rule.

**A2 research question.** In paired-donor analyses, do the prespecified SenMayo, OXPHOS, and TAM programs show regional patterns that remain interpretable with expression-matched controls and declared covariate adjustment?

**A3 research question.** Do the six paired donors show a consistent Core-versus-Periphery difference in the fixed 119-of-121 Ntafoulis directional expression proxy? This analysis has no preferred direction or threshold and is not a clinical sensitivity test.

### Track B — longitudinal MRI and tissue-state transitions

**B1 research question.** Within patient and physical-distance band, how does next-visit core-label transition relate to history-defined new versus persistent FLAIR-labeled tissue? Does a prespecified forecast add information beyond history-only mask baselines?

**B1 exploratory hypothesis.** Newly appearing FLAIR-labeled tissue will have a higher next-visit core-label transition rate than persistent FLAIR-labeled tissue in the same patient and distance band. The exposed Duke analysis is exploratory. An independent confirmatory evaluation requires a separate eligible cohort.

**B2 interpretation rule.** MRI-defined recurrence or core-label transition will be called TMZ-associated only if permitted patient-level exposure, scan chronology, and relevant clinical context are verified. It will not be called direct TMZ resistance.

### Track C — MGMT association and simulated treatment policies

**C1 observational question.** Among eligible CGGA GBM records, how is MGMT methylation associated with OS after adjustment for age, IDH status, radiotherapy, TMZ status, recurrence status, and data batch? Does the association differ by recorded TMZ status?

**C1 interpretation rule.** These are retrospective observational associations. They are not randomized treatment effects and will not calibrate individual simulator parameters.

**C2 simulation question.** Across independently generated synthetic GBM worlds, does a finite-set dependence-aware minimax controller achieve higher robust utility than a matched controller whose residual signal-dependence parameter is fixed to zero, when both receive the same virtual patients, signals, treatment choices, and observation schedules?

**C2 hypothesis.** In correctly specified scenarios with informative signals and nonzero residual dependence, the dependence-aware controller will have higher equal-world mean robust utility. In the predeclared uninformative and reversed-signal controls, it will show no consistent advantage over no-signal or single-signal controls.

### Integrated question

When molecular-like and MRI-like proxy signals are generated jointly from the same simulated resistant-tumor state, does explicitly modeling their residual dependence improve robust simulated treatment utility over the matched conditional-independence policy? This is a simulation hypothesis, not a claim about patient benefit.

### Engineering goals and expected outcomes

1. Preserve the frozen Track A analysis and its original success rule.
2. Keep Track A, Track B, the CGGA analysis, and simulation evidence on their own cohorts, endpoints, and interpretation scales.
3. Make any permitted Track B use reproducible from source, image geometry, label support, registration, and target-coverage records.
4. Build an executable Track C research-scale benchmark with paired common inputs, independent development and evaluation worlds, recorded seeds, and a pre-run seal.
5. Select and validate a simultaneous uncertainty method before generating policy outcomes. The screened observed-SE max-T helper is excluded.
6. Test the dependence-aware controller under the positive rule below. Report the estimate and uncertainty as produced, whether the result supports the hypothesis, shows little difference, or favors the comparator.

Expected project outputs are: a reproducible account of the completed Track A work; a clearly bounded exploratory Track B transition analysis and any independently eligible forecast comparison; an observational CGGA report; and, if all synthetic validation and execution gates pass, a synthetic policy comparison. No output is a clinical treatment recommendation.

## Procedures

### Materials and data sources

- **Computing:** Windows 11 workstation with Python 3.14.4 and 32 GB RAM. Record processor, thread settings, operating system, package versions, and runtime for each future sealed run.
- **Track A:** Yu et al. single-cell RNA data and sample metadata (GEO GSE117891; Yu et al., 2020); Table S1 point labels; GBmap annotation metadata; Ivy GAP (GSE107559); Angers (E-MTAB-5092); E-MTAB-8743; and the Ntafoulis et al. (2023) supplementary viability-correlation panel.
- **Track B:** MU-Glioma-Post and LUMIERE records used for prior work and method development; the restricted Duke longitudinal MRI dataset (Zenodo 20548528; Criscuolo et al., 2026); synthetic images and masks for geometry checks; and RANO 2.0 guidance (Wen et al., 2023). The Duke source lists restricted access options. Follow the exact accepted agreement and school or competition conditions for each use and derived output.
- **Track C observational:** the two CGGA clinical-data batches and the published TCGA-derived age reference used only as context.
- **Track C simulation:** the versioned [integrated v2 configuration](../configs/trackABC_integrated_v2.json), seed-stream templates, input manifest, simulator, policy implementations, and synthetic inputs. No patient-level input is allowed in the integrated simulator.
- **Aggregate model checks:** published aggregate RTOG 0525 and GEINO 14-01 results may be used for endpoint-compatible checks (Gilbert et al., 2013; Balaña et al., 2020). They are not patient-level dose-response data.
- **Physical materials:** none. The project uses no cells, animals, chemicals, experimental drugs, or participant recruitment.

### Track procedures

### Track A — expression and regional biology

#### A5 classifier procedure

1. Use the six frozen donors: GS1, GS2, GS3, GS5, GS6, and GS13. Verify the frozen A5 checksums and exact cell, point, donor, and Table S1 joins.
2. Define Tumoral points as Core and Peritumoral points as Periphery; exclude Normal points and invalid joins.
3. Retain cells with at least 1,000 UMIs and 500 detected genes. Within each training split, keep genes detected in at least 5% of training cells, select the 2,000 most variable genes, normalize each cell to 10,000 counts, apply log1p, and fit feature means and standard deviations on training cells only.
4. Use six outer leave-one-donor-out folds and inner leave-one-donor-out folds for setting selection. Compare the frozen sparse donor-balanced candidate with L2 logistic regression, 15-nearest-neighbor voting, linear SVM, and gradient-boosted trees. Keep the candidate’s three ablations, cell-type-only and cells-per-point diagnostics, and constant-score controls identified as such.
5. Evaluate class-balanced AURC as the primary metric. The success rule requires lower equal-donor mean AURC than logistic regression and lower AURC in at least five of six donors. Report balanced accuracy and Brier score as secondary metrics. The A5 result above is the result under this fixed procedure; it will not be retuned.

#### Regional programs and response-associated expression score

1. Keep the paired-donor SenMayo and OXPHOS analyses in Ivy GAP and Angers separate from TAM activation in E-MTAB-8743. Use donor-level summaries, expression-matched gene-set controls, exact donor sign-flip tests, and the declared covariate analyses. Do not treat cells as independent donors.
2. For the TMZ-response-associated proxy, use the 121-gene `Correlation_cell_Overlap` source panel from Ntafoulis et al. (2023) and original Spearman-correlation directions against in-vitro cell viability at 100 μM TMZ. The A5 matrix contains exact symbols for 119 genes; omit `AF131215.2` and `PPIAP22` without aliasing or substitution. Keep the source panel and the dated 119-of-121 amendment together.
3. Create donor-by-region pseudobulk summaries and calculate the fixed signed expression score. Report all six paired Core-minus-Periphery differences, 1,000 expression-matched control panels, and the exploratory exact donor sign-flip calculation. Interpret a score only as relative expression alignment with the published in-vitro viability directions.
4. Use donor-linked MGMT status only if the sample join, assay, and provenance are verified. Use a transfer cohort only after it meets the independent-donor, region-map, cell-label, reuse-rights, and prior-exposure requirements. Any additional analysis requires a dated protocol amendment.

### Track B — longitudinal MRI and spatial transitions

1. Preserve the completed 10-mm Duke transition analysis as exploratory and exposed. Its procedure used prior-visit FLAIR dilation, a middle-visit-core exclusion, and three physical-distance bands. Do not characterize that cohort as an untouched validation set or use it to choose new thresholds.
2. Keep the older 2-mm transition-estimand prototype distinct from the executed 10-mm procedure. Before any new analysis, align the selected estimand, physical coordinates, scan dates, labels, registration transform, resampling support, and target-field coverage with the actual pipeline.
3. For a permitted exploratory transition analysis, define history-based tissue states before target-label scoring. Summarize transitions at the patient level within physical-distance bands. Retain an explicit coverage mask so out-of-field tissue is not treated as observed non-conversion. Freeze the patient-level uncertainty method before any new target-label analysis.
4. For a future forecast comparison, compare persistence, fixed growth, PDE-only, hybrid, 3D U-Net, time-conditioned FiLM U-Net, and suitable count-only baselines on the same patients, support, target, and interval. Fit calibration only within training folds. Treat synthetic mechanics checks as code checks, not as patient-image validation.
5. Verify TMZ exposure and clinical chronology before using “TMZ-associated.” Otherwise call the endpoint MRI-defined recurrence or core-label transition. Apply RANO 2.0 confirmation rules only when the available clinical and scan records allow them (Wen et al., 2023); the simulation’s scan-sampled burden threshold is not RANO.
6. Keep ADC as a separate analysis with its own frozen source definitions, read scope, preprocessing, quality checks, and missingness rules.
7. Use a separate cohort for any confirmatory transition analysis. Verify access rights, source labels, longitudinal visit map, registration QC, spatial support, prior-use overlap, and sample eligibility before target labels are opened.

### Track C — CGGA association and tumor-treatment model

#### Observational analysis

Use the two CGGA clinical-data batches, remove duplicate patient identifiers under the recorded rule, apply the declared GBM and survival criteria, and include a batch indicator. Fit Cox proportional-hazards models for OS with age, IDH status, radiotherapy, TMZ status, recurrence status, and MGMT methylation. Report the MGMT model among TMZ-treated records and the MGMT-by-TMZ interaction in eligible pooled records. Report complete-case counts and uncertainty. Assess proportional-hazards assumptions before interpreting model coefficients. Keep all estimates noncausal; treatment was not randomized and treatment timing and dose are not represented in the saved data.

#### Synthetic tumor dynamics

Represent sensitive burden as S, resistant burden as R, and total burden as B=S+R, with capacity 1. The model is:

`dS/dt = rS[1 − (S + a_rsR)] − mruS − μuS + ν(1 − u)R`

`dR/dt = r(1 − c)R[1 − (a_srS + R)] − εmruR + μuS − ν(1 − u)R`

with `S(0)=(1−f0)B0`, `R(0)=f0B0`, and `a_rs=a_sr=1`. Here `u` is treatment state; `r`, `m`, `f0`, `B0`, `c`, `ε`, `μ`, and `ν` are declared simulation inputs. These are not fitted patient-specific treatment effects.

Use the v2 configuration as the authoritative parameter source. World-level continuous axes are log-growth shift U(−0.7, 0.9), log-resistant-fraction shift U(−0.7, 2.3), and B0 center U(0.12, 0.55). The single-response log-m scale is U(ln 0.05, 0); the mixture nonresponder log-m scale is U(ln 0.02, ln 0.30), and the responder fraction is U(0.1, 0.6). The nine `(μ,ν)` combinations use each rate in `{0, 0.01, 0.04}` per day. Patient-level ranges include `r` log-uniform on [0.001, 0.05] per day, `m` log-uniform on [6, 20], `f0` log-uniform on [0.001, 0.1] before its world shift and 0.9 cap, `B0` uniform within ±0.10 of its world center, `c` uniform on [0, 0.25], and `ε` uniform on [0, 0.2]. In mixture worlds, responder `m` uses [6, 20] and nonresponder `m` uses the world-scaled weak range. Observation-noise σ is balanced over {0.05, 0.10, 0.20}. Keep the single-response and mixture-response families separate.

The particle-filter prior uses `r` log-uniform on [0.001, 0.05], `m` log-uniform on [0.1, 20], `f0` log-uniform on [0.001 exp(−0.7), 0.9], `c` uniform on [0, 0.25], `ε` uniform on [0, 0.2], and `B0` uniform on [0.05, 0.5].

Integrate with classical fourth-order Runge–Kutta at 0.25-day steps and compare with 0.05-day steps for numerical convergence. The configured simulation horizon is 1,456 days; the primary progression and utility outcomes are summarized through day 730. The sealed driver must use these windows consistently and record treatment courses through the primary utility horizon.

### Integrated signals and policy comparison

#### Synthetic signal generation

Map resistant fraction `f0` to latent normal score:

`zR = Φ⁻¹(clip((ln(f0) − ln(L))/(ln(U) − ln(L)), 1e−8, 1−1e−8))`

where `L=0.001 exp(−0.7)` and `U=0.9`. Draw independent `EB, EI ~ Normal(0,1)` and calculate:

`zB = a zR + sqrt(1−a²) EB`

`zI = b zR + sqrt(1−b²) [d EB + sqrt(1−d²) EI]`

Set hidden binary state `C=1[zB≥0]`; flip the observed molecular-like report with probability `1−q`. The MRI-like baseline signal is `zI` unless an allowed, frozen history-only aggregate is used. All `a`, `b`, `d`, and `q` values are synthetic scenario factors, not patient estimates or measured assay accuracy.

The full synthetic q grid is `{0.5, 0.6, 0.7, 0.8, 0.9}`. The primary q sampling rule in the current configuration draws from `{0.6, 0.7, 0.8, 0.9}` once per world and holds it constant across that world’s virtual patients; `q=0.5` is used for null controls. The model grid for `a`, `b`, and `d` is `{−0.5, 0, 0.5}` with primary anchor `(0.5, 0.5, 0.5)`. The full no-information control uses `(q,a,b,d)=(0.5,0,0,0)`. The molecular-channel-null control uses `q=0.5` at the primary `(a,b,d)` anchor. The predeclared out-of-grid sensitivity anchors are `(0.75,0.75,0.75)`, `(0.75,0.75,−0.75)`, `(−0.75,0.75,0.5)`, `(0.5,−0.75,−0.5)`, and `(−0.75,−0.75,0.75)`.

For particle score `zR` and observed `zI`, compute `v=(zI−bzR)/sqrt(1−b²)`, then:

`pB = Φ((a zR + sqrt(1−a²)d v)/sqrt((1−a²)(1−d²)))`

`pReportPositive = q pB + (1−q)(1−pB)`.

The joint likelihood is the normal density of `zI` given `zR` multiplied by the probability of the observed report. At `q=0.5`, the binary report contributes no information. The implementation must follow the likelihood in the [independent handoff guide](../tracks/shared/HANDOFF_AND_REPRODUCIBILITY.md).

#### Matched treatment controllers

Compare the finite-set dependence-aware minimax controller with the matched finite-set controller setting residual dependence `d=0`. Within each paired comparison, provide identical virtual patients, latent parameters, signals, scan schedule, observation noise, and treatment choices. The candidate model set has at most 135 combinations before duplicate removal. Use 500 filter particles per model and 16 common stratified predictive draws per model.

At each replan, score all eight treatment/no-treatment patterns over three 28-day cycles, forecast six cycles (168 days) without terminal value, apply the first three-cycle block, and replan at the next scheduled scan. Break ties by the fixed action order favoring fewer courses. A modeled course treats for five days; the total cap is 12 courses. The primary scan interval is 84 days, with 28-, 56-, and 112-day scan intervals and visible scan loss as sensitivity scenarios. The primary visible scan-loss probability is zero; 0.10 is the declared sensitivity. The synthetic observable stop proxy stops after two consecutive visible noisy readings reach at least 1.4 times the lowest earlier observed volume, with a floor of 0.10 times observed starting volume; a missing scan resets confirmation. This is not formal RANO. The configuration’s policy library also includes no-signal, development-selected fixed schedules, simple scan-based controls, and burden-only Q-learning. The latter is not a same-information primary comparator. Oracle policies use hidden simulated parameters and are diagnostic only.

### Simulation sampling, randomness, and reproducibility

1. Draw worlds independently. For each world, draw continuous axes with a scrambled one-point Latin hypercube and that world's registered world-parameter seed. Draw `μ` and `ν` independently and uniformly from `{0, 0.01, 0.04}`; do not balance their nine combinations across a batch, because the estimand and proposed inference treat worlds as independent. Draw `q_signal` once per world from `{0.6, 0.7, 0.8, 0.9}`. In the mixture family, call the responder fraction `pi_response` so it cannot be confused with signal quality.
2. Generate selection, calibration, and evaluation patients independently with their separate cohort seeds and the configured dimension order `(r,m,f0,B0,c,ε)`. Draw resistance cost `c` per patient. After all three cohorts and IDs exist, assign `{0.05,0.10,0.20}` in balanced counts across all sorted IDs in the world, using the PCG64DXSM child stream `SeedSequence([world_seed, 0x5349474D])`.
3. Use explicit `PCG64DXSM` generators for world parameters, patient cohorts, sigma assignment, joint calibration/evaluation signal generation, and bootstrap draws. Use one joint signal stream per cohort to preserve the declared dependence structure. The existing policy simulator currently uses `numpy.default_rng` (PCG64 in the recorded environment); this generator boundary must be recorded or migrated before evaluation.
4. Supply distinct, registered integer seed roles for world, each cohort, each joint signal cohort, and each policy stage. The same evaluation seed is currently passed to paired policy arms, so they share scan-noise draws and policy-reset stream position. Review this coupling before evaluation. The seed-template file is not materialized or reserved; collision-check all numeric schedule seeds against the registry before a run.
5. Share common patients, latent parameters, signals, scan schedule, and observation-noise realizations across paired policies. Record cohort and world identifiers, input and model hashes, parameter fingerprints, package versions, environment, failures, missing arms, fallbacks, exclusions, and all prespecified sensitivities.
6. Complete the research-scale driver, exact cohort/world counts, outcome-blind timing and precision plan, runtime cap, clean-install check, and independent review of a valid replacement interval method before generating policy outcomes. Existing development seeds and outcomes are not evaluation data.
5. Share common random inputs across paired policies. Record cohort and world identifiers, input and model hashes, parameter fingerprints, package versions, environment, failures, missing arms, fallbacks, exclusions, and all prespecified sensitivities.
6. Complete the research-scale driver, exact cohort/world counts, outcome-blind timing and precision plan, runtime cap, clean-install check, and review of the alternative interval method before generating policy outcomes. Existing development seeds and outcomes are not evaluation data.

The current seed-capacity note identifies 89 free scenario-class seeds in the relevant development block, corresponding to at most 14 worlds at six registered scenario seeds per world before policy seeds and other design constraints. The sigma-assignment child stream does not consume another registered root seed. This is not a precision calculation or a reservation. The final count will come from the written outcome-blind planning rule.

### Project sequence

1. Preserve the completed Track A analyses and their original frozen rules.
2. Complete any planned Track A extensions only under dated amendments and with no reuse of the held-out A5 results for model tuning.
3. Keep the completed Duke transition analysis exploratory. Resolve source, registration, coverage, label, and chronology questions before any additional patient-level interpretation; use a separate eligible cohort for confirmation.
4. Report the CGGA association and diagnostics as an observational OS analysis.
5. Finish and review the Track C simulation protocol, inference method, independent-world driver, seed table, precision plan, and runtime gate using synthetic data only.
6. Seal the evaluation configuration and source, then run the predeclared policy comparison once with the positive rule above. Preserve all results, including null or comparator-favoring results.

## Risk and Safety

Tracks A and B use existing human-derived research data. The project recruits or contacts no participants and administers no treatment or laboratory drug assay. Keep identifiers, restricted files, and any disallowed derived outputs in access-controlled storage. The exact accepted Duke terms and school or competition conditions govern every further use and sharing decision. If permitted use of a derived summary is unclear, do not use or share it; the integrated simulation can use its standard-normal MRI-like fallback.

Do not reopen patient images, masks, ADC maps, or patient-level records for the integrated simulation. It receives only allowed aggregate context or synthetic signals. Long computations will be run on a ventilated surface with the computer connected to power.

## Data Analysis

### Track A

The A5 primary metric is class-balanced AURC averaged equally across the six donors. Lower is better. The descriptive success rule requires both lower mean AURC than L2 logistic regression and lower donor-specific AURC in at least five donors. Balanced accuracy and Brier score are secondary. Regional gene programs use paired donor-level scores, matched controls, exact sign-flip tests, and prespecified covariates. The TMZ-associated score is summarized as six donor-paired Core-minus-Periphery values with its descriptive matched-panel reference; it has no sensitivity/resistance accuracy threshold.

### Track B

The patient is the unit of summary for transitions. Report patient-level transition summaries and the declared distance bands. Any future uncertainty procedure will account for patient-level clustering and will be fixed before opening target labels. Forecast comparisons must share patient, support, target label, and forecast interval across methods. Report MRI outcomes without a TMZ-associated label when exposure chronology is not verified.

### Track C observational

Report Cox hazard ratios, confidence intervals, p-values, complete-case counts, covariate definitions, batch handling, and proportional-hazards diagnostics. Interpret MGMT estimates and interaction tests as noncausal associations with OS.

### Integrated simulation

Record two latent progression measures through 730 days:

1. Interpolated time when true burden first reaches 1.2 times starting burden.
2. First scheduled scan when true burden reaches 1.4 times the larger of the previous true-burden scan nadir or 0.10 times starting burden.

The second measure is a scan-sampled latent-burden proxy, not clinical RANO. Define utility as the smaller progression time, each capped at 730 days, minus five simulated days per treatment course delivered through day 730. Run course-cost sensitivities at zero and ten simulated days per course.

The primary contrast is candidate utility minus matched `d=0` utility, paired within each independent world and averaged equally across worlds. Keep single-response and mixture-response families separate. The positive rule requires both simultaneous adjusted lower confidence bounds to exceed zero and both interval half-widths to be no more than 14 simulated days.

The observed-SE max-T proposal is not used: its predeclared coverage screen failed. Before policy outcomes are generated, an alternative simultaneous inference method must be chosen, independently reviewed, frozen in a dated protocol, and validated on fresh outcome-blind synthetic data that reflect the intended world-level dependence structure. If no method meets its coverage and precision requirements at feasible world counts, report descriptive estimates with their stated uncertainty or report the comparison as inconclusive. Do not alter the positive rule after outcome access.

## References

Akbari, H., et al. (2016). Imaging surrogates of infiltration obtained via multiparametric imaging pattern analysis predict subsequent location of recurrence of glioblastoma. *Neurosurgery, 78*(4), 572–580. https://doi.org/10.1227/NEU.0000000000001202

Balaña, C., et al. (2020). Continuing adjuvant temozolomide beyond 6 cycles in patients with glioblastoma (GEINO 14-01). *Neuro-Oncology, 22*(12), 1851–1861. https://doi.org/10.1093/neuonc/noaa107

Bhattacharya, M., Singh, G., & Prasanna, P. (2025). *SoC-DT: Standard-of-Care aligned digital twins for patient-specific tumor dynamics* [Preprint]. arXiv:2510.03287. https://arxiv.org/abs/2510.03287

Chaudhuri, A., et al. (2023). Predictive digital twin for optimizing patient-specific radiotherapy regimens under uncertainty in high-grade gliomas. *Frontiers in Artificial Intelligence, 6*, 1222612. https://doi.org/10.3389/frai.2023.1222612

Criscuolo, E., et al. (2026). *A longitudinal glioblastoma MRI dataset with anatomical landmark pairs for DIR validation* [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.20548528

Ebrahimi Zade, A., Shahabi Haghighi, S., & Soltani, M. (2020). Reinforcement learning for optimal scheduling of glioblastoma treatment with temozolomide. *Computer Methods and Programs in Biomedicine, 193*, 105443. https://doi.org/10.1016/j.cmpb.2020.105443

Gilbert, M. R., et al. (2013). Dose-dense temozolomide for newly diagnosed glioblastoma: A randomized phase III clinical trial. *Journal of Clinical Oncology, 31*(32), 4085–4091. https://doi.org/10.1200/JCO.2013.49.6968

Hegi, M. E., et al. (2005). MGMT gene silencing and benefit from temozolomide in glioblastoma. *New England Journal of Medicine, 352*(10), 997–1003. https://doi.org/10.1056/NEJMoa043331

Liu, Q., et al. (2025). Treatment-aware diffusion probabilistic model for longitudinal MRI generation and diffuse glioma growth prediction. *IEEE Transactions on Medical Imaging*. https://doi.org/10.1109/TMI.2025.3533038

Liu, Z., et al. (2025). Multiscale mathematical model-informed reinforcement learning optimizes combination treatment scheduling in glioblastoma evolution. *Science Advances, 11*(32), eadv3316. https://doi.org/10.1126/sciadv.adv3316

Ntafoulis, I., et al. (2023). Ex vivo drug sensitivity screening predicts response to temozolomide in glioblastoma patients and identifies candidate biomarkers. *British Journal of Cancer, 129*(8), 1327–1338. https://doi.org/10.1038/s41416-023-02402-y

Zimmer, L., Weidner, J., Balcerak, M., Kofler, F., Krupa, M., Ezhov, I., Cepeda, S., Zhang, R. Z., Lowengrub, J. S., Menze, B., & Wiestler, B. (2026). PREDICT-GBM: A multicenter platform advancing personalized glioblastoma radiotherapy planning. *npj Digital Medicine, 9*, 686. https://doi.org/10.1038/s41746-026-03194-0

Zimmer, L., Weidner, J., Balcerak, M., Kofler, F., Ezhov, I., Menze, B., & Wiestler, B. (2025). PREDICT-GBM: Platform for Robust Evaluation and Development of Individualized Computational Tumor Models in Glioblastoma [Preprint]. arXiv:2509.13360. https://arxiv.org/abs/2509.13360

Wen, P. Y., et al. (2023). RANO 2.0: Update to the response assessment in neuro-oncology criteria for high- and low-grade gliomas in adults. *Journal of Clinical Oncology, 41*(33), 5187–5199. https://doi.org/10.1200/JCO.23.01059

Yu, K., et al. (2020). Surveying brain tumor heterogeneity by single-cell RNA-sequencing of multi-sector biopsies. *National Science Review, 7*(8), 1306–1318. https://doi.org/10.1093/nsr/nwaa099

“Personalized machine learning-guided radiation dose escalation in newly diagnosed glioblastoma: prospective pilot study.” (2026). *Nature Communications, 17*, 6413. https://doi.org/10.1038/s41467-026-72545-y
