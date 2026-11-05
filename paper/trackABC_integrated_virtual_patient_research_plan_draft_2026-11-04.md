# Integrated Track A/B/C Research Plan Draft

**Proposed title:** *A Data-Informed Virtual Glioblastoma Therapy Model: Testing Multimodal Value Under Uncertain Molecular, MRI, and Drug-Response Links*

**Status:** Proposal for a new integrated direction. It is not an adopted protocol and does not authorize new patient-data reads, scoring, or a sealed simulation. Existing Track A results and Track C development results will remain unchanged.

## Rationale

Glioblastoma contains different cell states, changes over time, and can respond differently to treatment. Single-cell RNA sequencing can describe some molecular patterns. MRI and apparent diffusion coefficient (ADC) can show radiographic changes across visits. Mathematical models can represent drug suppression and resistant tumor populations. Combining these signals could help test when adaptive treatment schedules might work.

The available project data do not measure all three signals in the same patients. Track A uses mixed-cell expression and MRI-marked sampling-point labels from Yu et al. (2020). Track B uses Duke longitudinal MRI, ADC, automated masks, and relative treatment dates (Criscuolo et al., 2026). The Duke records do not link to Track A donors or establish individual drug dose, cycles, or adherence. Track C represents drug effects with assumed parameters. Therefore, the project will not join unrelated people or call any result patient-specific therapy.

Instead, this project will build one **data-informed virtual-patient simulation**. Track A will inform uncertainty in a simulated expression-derived spatial-state signal. Track B will inform radiographic observation patterns only if its accepted data terms and approvals permit that analysis. Track C will represent treatment suppression, resistance, and tumor burden. Because the source cohorts do not reveal cross-modal relationships, the simulation will vary those relationships across predeclared scenarios, including no relationship and reversed relationships.

Related work already combines molecular tumor models, drug scheduling, and reinforcement learning (Z. Liu et al., 2025), treatment-aware MRI prediction (Q. Liu et al., 2025), interval-conditioned longitudinal MRI classification (Guo et al., 2026), and physics-guided MRI synthesis (Biller et al., 2026). This project will not claim to be the first multimodal therapy model. Its candidate contribution will be a controlled test of how unknown links between molecular, MRI, and drug-response signals change the apparent value of multimodal treatment planning.

## Research Question, Hypothesis, and Engineering Goals

### Research question

Across held-out virtual glioblastoma worlds, when does a treatment policy that combines expression-derived state, serial MRI/ADC observations, and simulated drug suppression improve simulated tumor control without using more drug than the best fixed schedule? Does any advantage remain when the assumed links among these signals are absent, weaker, or reversed?

### Hypothesis

Under the predeclared moderate-concordance scenario (rank correlation 0.50), the fused policy will increase restricted mean time to a prespecified simulated tumor-burden threshold without increasing modeled drug-course use, compared with the best fixed schedule selected on development worlds. Its advantage will shrink or disappear when the cross-modal links are absent or reversed. This is a conditional simulation hypothesis, not a prediction about patient benefit. The result may be null or negative.

### Engineering goals

1. Build one virtual-patient record that carries molecular, MRI/ADC, tumor-burden, and treatment state together.
2. Mark every value as measured, simulated, or assumed.
3. Test one fused uncertainty-aware belief-state model-predictive controller against fixed-schedule and modality-specific controls.
4. Make every policy use the same virtual patients, scan schedule, observation noise, missing scans, and maximum treatment-course budget.
5. Evaluate only on held-out simulation worlds that are not used to choose model settings.

## Materials

- Dell 16 Plus DB16250 computer running Windows 11, with 32 GB RAM and CPU computation.
- Python 3.14.4 and the project environment. NumPy, pandas, SciPy, scikit-learn, PyTorch, NiBabel, and project-specific Track A/B/C code will be version-pinned before the locked evaluation. Exact versions will be recorded in the run manifest.
- Track A source: Yu et al. single-cell count matrix and sample metadata, GEO GSE117891, plus the study's Table S1 and the saved GBmap cell-annotation metadata. These data will supply donor-level estimates of uncertainty for a simulated expression-derived spatial-state signal. A5's Core/Periphery classifier will not be treated as a drug-sensitivity biomarker.
- Track B source: Duke's restricted longitudinal glioblastoma MRI dataset, Zenodo record 20548528. Its public record describes 61 patients and 503 sessions, with 3–18 sessions per patient, serial MRI/ADC, automated segmentation labels, relative scan dates, and treatment-event dates. These totals do not establish the eligible quality-controlled cohort. The student will use the data only under the exact accepted Duke agreement and school/competition approval. No raw scans, masks, patient-level predictions, or other restricted outputs will be shared unless the agreement permits it.
- Track C source: the project's tumor-growth and drug-suppression simulator, plus published aggregate treatment outcomes from the RTOG 0525 and GEINO 14-01 trials. Trial-level results will constrain simulated scenarios; they will not be treated as patient-level dose-response data.
- No chemicals, biological specimens, patient contact, or treatment intervention will be used.

## Procedure

1. **Freeze the scope and rights.** The student will save the exact Duke access grant, accepted license, and any school or competition conditions. The student will confirm which analyses and derived outputs may be created and shared. If the terms do not permit the planned MRI feature extraction or public summary, the student will omit the Duke-derived component and label the project accordingly. Patient images, masks, landmarks, and patient-level outputs will remain local and access-controlled.

2. **Preserve prior results.** The completed Track A A5 result will remain an exploratory within-cohort result. Its frozen five-of-six donor success rule will not be changed. A5 will not be used as proof of malignancy, invasion, drug resistance, or treatment response. The student will not retune A5 or select a new gene program after viewing integrated results.

3. **Create the Track A molecular input.** The student will use the frozen A5 donor-level balanced-accuracy estimates to set a broad uncertainty range for the reliability of a simulated expression-derived spatial-state signal. Mapping classifier accuracy to signal reliability will be an explicit modeling assumption. The six donors will be the resampling units; cells will not be treated as independent donors. Core and Periphery will retain their original meaning as MRI-marked sampling-point labels. The virtual signal will be generated by the simulation; it will not be assigned to a Duke patient or called a drug-sensitivity score. Its association with simulated drug sensitivity will be an explicit scenario assumption, not an estimate from Yu et al. Because A5 failed its frozen success rule, the simulation will also include weak-information and no-information molecular-channel controls.

4. **Create the Track B radiographic input, if approved.** Before opening MRI or mask voxels, the student will freeze visit selection, interval handling, registration, support generation, image-quality thresholds, label handling, exclusions, and output rules. Registration and image checks will use MRI only. Target masks will not guide registration, support, or exclusions. Label 0 is not defined in the current source record, so it will not be called normal; it will be treated as unknown unless an authoritative source resolves its meaning. The cavity label will be analyzed separately from tumor labels. After the image-only quality-control list is frozen, the student will calculate patient-level changes in released MRI labels and ADC. These will be called radiographic observations, not histologic tumor truth. Patients and visits will not be treated as independent samples.

5. **Check the Track C simulation before policy comparison.** The student will test whether the drug-suppression model can reproduce its predeclared published aggregate anchors and negative controls, including standard-arm progression outcomes, dose-dense comparison, and the six-versus-twelve-cycle comparison where applicable. Existing development results show that this calibration remains weak, especially for the MGMT-stratified check. If the model cannot meet the frozen calibration rules, the student will stop the therapy-policy comparison and report the calibration failure. The student will not change calibration windows to obtain a policy win.

6. **Build shared virtual patients.** Each virtual patient will contain a latent tumor burden, sensitive and resistant tumor compartments, a noisy molecular-state observation, an MRI/ADC observation history, scan timing, missingness, and a treatment-course history. Track A and Track B will provide separate marginal distributions or observation models; they will not be joined by patient ID. The student will use valid, predeclared target-rank-correlation matrices across molecular state, MRI/ADC state, and drug resistance: independence (all off-diagonal values 0); weak concordance (all 0.25); moderate concordance (all 0.50); and a reversed molecular-to-resistance link (molecular-to-MRI 0.25, MRI-to-resistance 0.25, molecular-to-resistance −0.25). The moderate-concordance scenario will be the primary conditional scenario. Every correlation is an analyst-selected simulation assumption, not an observed biological estimate.

7. **Implement the fused treatment policy.** The simulator will use 28-day treatment cycles and an 84-day primary scan interval, matching the planned standard schedule and trial-based scan timing. The treatment action will be a simulated TMZ-like course or no course, with no more than six planned courses over the evaluation horizon. These are model actions, not dose instructions. At each scheduled scan, an uncertainty-aware belief-state model-predictive controller will update its estimate of the virtual patient's tumor state from the available molecular, MRI/ADC, and treatment observations. It will choose the next modeled treatment action subject to the same maximum course budget as the comparators. The policy will use no hidden simulator state. An oracle with access to hidden parameters may be shown as a diagnostic ceiling only; it will not count as a practical competitor.

8. **Compare against fixed and modality-specific controls.** The same virtual patients and random streams will be used for: (a) the published standard fixed schedule; (b) the best fixed schedule selected on development worlds; (c) an MRI/ADC-only policy; (d) an expression-only policy; (e) a Track-C-only controller that receives no A or B observation; (f) a simple late-fusion policy; and (g) shuffled-signal and no-link negative controls. A Q-learning policy may be included as a secondary baseline only if its code passes the same synthetic checks and its settings are frozen before evaluation. The project will not claim superiority over every published method unless that method is implemented and compared fairly.

9. **Separate development from evaluation.** The student will use one registered set of random seeds and worlds for model development. A new, non-overlapping set of seeds and worlds will be sealed before policy settings are chosen. Every policy will receive the same patient draws, scan schedule, observation noise, and missing-scan pattern. No held-out world will be used for feature selection, tuning, calibration, or stopping.

10. **Freeze the scoring and report all results.** Before held-out evaluation, the student will freeze the burden endpoint, RANO-like endpoint, course-use measure, uncertainty method, multiplicity rule, Monte Carlo precision target, and success rule. All coupling scenarios, failed calibrations, model failures, and negative results will be reported. No policy or threshold will be changed after held-out results are inspected.

## Risk Assessment

This project will not recruit or contact people, administer drugs, or provide clinical recommendations. Main risks are privacy, restricted-data misuse, and mistaken medical interpretation. The student will keep Duke files and patient-level outputs on a local access-controlled drive, avoid cloud upload or public source control for restricted material, and share only outputs allowed by the written agreement. The student will report simulated drug suppression as modeled behavior, not measured patient response. The student will describe MRI masks as automated labels, not confirmed tumor biology. School and competition review will determine the applicable rules for secondary patient data.

## Data Analysis

- **Primary simulated endpoint:** restricted mean time, through 730 days, until simulated total tumor burden reaches 1.2 times its starting value.
- **Secondary endpoint:** time to a scan-based RANO-like progression call. This will remain secondary because scan timing and the nadir rule can delay a call without reducing simulated tumor burden.
- **Treatment use:** number of modeled courses delivered and the time-versus-course-use frontier. A policy will not pass the primary claim if it improves time only by using more courses.
- **Independent evaluation unit:** held-out simulation world. Virtual patients are Monte Carlo samples nested inside worlds; they will not be counted as independent worlds.
- **Uncertainty:** each policy contrast will be paired within world. The primary interval will be a two-sided 95% paired bootstrap confidence interval over worlds, not cells, visits, or individual simulated time points. The number of held-out worlds will be selected by a precision pilot before the held-out results are opened.
- **Primary comparison:** fused policy versus the best fixed schedule selected on development worlds, subject to no higher average course use in held-out evaluation. The result will support a simulated superiority claim only if the predeclared interval for the primary paired difference excludes zero in the favorable direction and the course-use condition is met.
- **Secondary comparisons:** modality-only and simpler-fusion controls will be reported with a prespecified multiplicity adjustment. A broad claim that fusion beats all required controls will be made only if all adjusted intervals support that claim.
- **Robustness:** effects will be shown separately for independent, weakly concordant, moderately concordant, and reversed coupling scenarios. Only the moderate-concordance comparison is primary; other scenario estimates will be descriptive. A gain limited to that assumed coupling will be described as conditional, not robust.
- **Calibration gate:** no treatment-policy result will be interpreted as evidence of a useful model if the drug-suppression model fails its frozen clinical aggregate anchors or negative controls.

All conclusions will be limited to the tested simulation assumptions and permitted radiographic summaries. The study will not establish patient-specific treatment benefit, causation, true tumor response, clinical usefulness, or a validated digital twin.

## Bibliography

Balana, C., Vaz, M. A., Sepúlveda, J. M., Mesia, C., del Barco, S., Pineda, E., Muñoz-Langa, J., Estival, A., de las Peñas, R., Fuster, J., Gironés, R., Navarro, L. M., Gil-Gil, M., Alonso, M., Herrero, A., Peralta, S., Olier, C., Perez-Segura, P., Covela, M., … Carrato, C. (2020). A phase II randomized, multicenter, open-label trial of continuing adjuvant temozolomide beyond 6 cycles in patients with glioblastoma (GEINO 14-01). *Neuro-Oncology, 22*(12), 1851–1861. https://doi.org/10.1093/neuonc/noaa107

Biller, V., Bubeck, N., Zimmer, L., Erdur, A. C., Nagar, S., Meyer-Baese, A., Rückert, D., Wiestler, B., & Weidner, J. (2026). *TumorFlow: Physics-guided longitudinal MRI synthesis of glioblastoma growth* [Preprint]. arXiv. https://arxiv.org/abs/2603.04058

Criscuolo, E., Calabrese, E., Zhang, Z., Wang, Y., Xiong, Y., & Yang, D. (2026). *A longitudinal glioblastoma MRI dataset with anatomical landmark pairs for DIR validation* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.20548528

Gilbert, M. R., Wang, M., Aldape, K. D., Stupp, R., Hegi, M. E., Jaeckle, K. A., Armstrong, T. S., Wefel, J. S., Won, M., Blumenthal, D. T., Mahajan, A., Schultz, L. M., Erridge, S. C., Brown, P. D., Chakravarti, A., Pugh, S., Zadeh, G., Fisher, B., Sulman, E. P., … Mehta, M. P. (2013). Dose-dense temozolomide for newly diagnosed glioblastoma: A randomized phase III clinical trial. *Journal of Clinical Oncology, 31*(32), 4085–4091. https://doi.org/10.1200/JCO.2013.49.6968

Guo, W., Yin, C., Giglio, P., Cao, W., Zhang, P., & Mirzaei, G. (2026). *STRIDE: Spatial-temporal representation for interval-conditioned disease evolution in longitudinal glioblastoma MRI* [Preprint]. arXiv. https://arxiv.org/abs/2610.08848

Liu, Q., Fuster-Garcia, E., Hovden, I. T., MacIntosh, B. J., Grødem, E., Brandal, P., Lopez-Mateu, C., Sederevicius, D., Skogen, K., Schellhorn, T., Bjørnerud, A., & Emblem, K. E. (2025). Treatment-aware diffusion probabilistic model for longitudinal MRI generation and diffuse glioma growth prediction. *IEEE Transactions on Medical Imaging*. https://doi.org/10.1109/TMI.2025.3533038

Liu, Z., Zhang, J., Hong, L., Nie, Q., & Sun, X. (2025). Multiscale mathematical model-informed reinforcement learning optimizes combination treatment scheduling in glioblastoma evolution. *Science Advances, 11*(32), eadv3316. https://doi.org/10.1126/sciadv.adv3316

Yu, K., Hu, Y., Wu, F., Guo, Q., Qian, Z., Hu, W., Chen, J., Wang, K., Fan, X., Wu, X., Rasko, J. E., Fan, X., Iavarone, A., Jiang, T., Tang, F., & Su, X.-D. (2020). Surveying brain tumor heterogeneity by single-cell RNA-sequencing of multi-sector biopsies. *National Science Review, 7*(8), 1306–1318. https://doi.org/10.1093/nsr/nwaa099

## Additional Project Information

This proposal will use secondary de-identified research data and computer simulation. The student will not recruit, contact, or treat people. Duke data access and any processing or sharing of derived outputs will follow the exact accepted Duke agreement. The student reports school and competition approval; the approval scope and conditions will be retained with the project records. No human-data exemption will be assumed without the applicable review determination. The project will not provide medical advice or recommend treatment for any person.
