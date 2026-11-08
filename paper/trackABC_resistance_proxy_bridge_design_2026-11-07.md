# Track A–B resistance-proxy bridge design

**Status:** Proposed design note. The separate Track A 119/121 expression proxy has been calculated as an exploratory aggregate score; it does not change the frozen A5 analysis. This note does not authorize new patient-data reads or a sealed Track C evaluation.

## Decision from the literature review

The current Track A and Track B datasets do not directly measure drug resistance. Track A measures mixed-cell expression by sampled tumor region. Track B measures longitudinal MRI and dataset-defined recurrence or progression. Neither provides a verified, same-patient temozolomide (TMZ) dose-response outcome joined to the other track.

The project can test **resistance-associated proxies**. It must not call those proxies direct drug-resistance measurements. A direct resistance study needs a drug-specific reference such as ex-vivo TMZ dose-response, or a clinically adjudicated outcome linked to documented TMZ exposure.

## Track A: molecular resistance-associated proxy

1. Preserve the frozen A5 Core-versus-Periphery classifier and its success rule. Do not reinterpret A5 AURC as resistance accuracy.
2. Use the fixed 121-gene overlap panel from Ntafoulis et al. (2023), Supplementary file 1, `Correlation_cell_Overlap`, with the approved exact-symbol mapping rule: score the 119 genes present in the A5 matrix, report 119/121 coverage, and omit `AF131215.2` and `PPIAP22` without substitution. This is a fixed response-associated gene set, not a validated predictor. Use the original Spearman-correlation signs against in-vitro viability at 100 μM TMZ; do not use the paper's separate survival-validation hits as response-score weights. The source panel and frozen score rules are in `paper/trackA_tmz_response_121_gene_panel_2026-11-07.tsv`, `paper/trackA_tmz_response_score_protocol_2026-11-07.md`, and `paper/trackA_tmz_response_119of121_amendment_2026-11-08.md`.
3. Sum raw counts into donor-by-region pseudobulks and calculate the frozen transcriptome-percentile direction score. Report all six paired donor differences with 1,000 expression-matched gene-set controls and the exploratory exact donor-level sign-flip calculation. The completed local, uncommitted exploratory output is `output/trackA_tmz_response_proxy_119of121/result.json`: the six donor directions split 3/3, the mean Core-minus-Periphery score was -0.00109, and the two-sided exact sign-flip p-value was 0.96875. This does not support a consistent regional shift. Do not treat cells as independent donors or label donors clinically TMZ-sensitive or resistant.
4. Use MGMT promoter-status metadata only if the donor-to-sample join and assay provenance are verified. Treat it as a resistance-associated biomarker, not as observed drug response.
5. Label the output **molecular resistance-associated proxy**. It cannot establish which donor’s tumor is TMZ-sensitive or resistant.

## Track B: treatment-associated MRI progression proxy

1. Keep the completed Duke transition analysis exploratory; its patient cohort and endpoint are exposed and cannot serve as untouched confirmation.
2. Before any new analysis, verify the patient-level TMZ regimen, exposure dates and duration, scan dates, surgery, radiation, and other treatment fields from permitted records.
3. Define an MRI progression endpoint before examining new target labels. Apply a prespecified confirmation rule and handle the early post-chemoradiotherapy pseudoprogression window. Use a TMZ-specific endpoint only if the treatment and progression chronology can be verified.
4. If those fields or confirmation evidence are missing, report **MRI-defined recurrence/progression**, not TMZ resistance. Keep registration, scan coverage, and segmentation checks separate from the clinical endpoint.
5. Use a fresh independent cohort for confirmation. Do not retune or relabel the exposed Duke result.

## Track A–B connection

The two current cohorts are not patient-matched, and their signals do not share a measured drug-response outcome. Their data therefore cannot estimate the cross-modal relationship or calibrate how accurately either signal measures resistance.

The integrated simulator may represent one latent resistant fraction and generate two noisy observations: (a) a molecular-like signal motivated by Track A’s published resistance-associated markers and (b) an MRI-like signal motivated by Track B’s longitudinal progression phenotype. Signal accuracy, bias, missingness, and residual dependence must remain explicit scenario assumptions unless a matched validation cohort supplies calibration. Do not set signal accuracy from A5 AURC or treat Duke recurrence as a TMZ-specific label without the exposure checks above.

The primary simulation comparison will remain the dependence-aware controller versus a matched controller with residual dependence fixed at zero. No-signal, molecular-only, MRI-only, and simple schedule controls will show whether any simulated gain requires both channels. The independent simulation world will be the uncertainty unit. A positive result will support only a simulated policy-utility claim under the stated assumptions; it will not show patient-specific resistance measurement or clinical treatment benefit.

## Evidence needed for direct resistance measurement

To make a direct resistance claim, a future cohort must link, for the same patients:

- tumor molecular measurements;
- serial MRI with documented acquisition and registration;
- named drug exposure, dates, dose or regimen, and relevant co-treatments; and
- direct ex-vivo drug-response measurements or a validated, adjudicated treatment-specific clinical response.

The reviewed IvyGAP MRI–RNA work supplies a partial imaging–molecular bridge, but not the required treatment-response endpoint. The reviewed sources did not establish a ready cohort meeting all four conditions. This is a bounded search result, not proof that no such resource exists.

## Reproducible simulation work still required

1. Implement a world-scale driver that reads one frozen configuration and seed registry, creates disjoint selection, calibration, and evaluation cohorts, runs every paired policy on identical inputs, aggregates at the world level, and writes a complete provenance record.
2. Run a timing-only synthetic pilot at the locked policy settings. The 2026-11-06 pilot used two patients per selection, calibration, and evaluation cohort and took 1,373 seconds for one response family and 1,483 seconds for the other. It did not measure research-scale runtime or policy performance.
3. Use development-only worlds to estimate the world-level variance needed for the predeclared 14-day simultaneous-interval precision target. Calculate and record the exact cohort and world counts before reserving evaluation seeds. Never use evaluation outcomes to resize the study.
4. Use a proposed total wall-clock cap of **24 hours on the current PC** for the locked evaluation; time the pilot separately. This cap is a project constraint, not an estimate from the two-patient pilot. If the precision-targeted design exceeds the cap, optimize or amend the design before evaluation, or report feasibility only. Do not lower the precision target after seeing evaluation outcomes.
5. Use development-only worlds to choose and record exact selection, calibration, and evaluation cohort sizes. Use the predeclared 14-day simultaneous-interval target; calculate the smallest evaluation-world count that meets it from development-only variance. The numeric counts remain open until the complete driver and precision pilot exist.
6. Materialize and collision-check numeric seeds under the new namespace, separate development and evaluation streams, seal code, config, environment, and input manifest, then run evaluation once.

## Claim wording

**Allowed:** “The study compares molecular and MRI progression proxies in separate cohorts and evaluates their assumed joint value in a synthetic treatment-policy simulation.”

**Not supported by current Track A/B data:** “The study measures patient drug resistance,” “the model identifies which real patient will resist TMZ,” or “the optimized schedule improves patient outcomes.”

## Sources reviewed

- Yu et al. (2020), [GSE117891 source paper](https://academic.oup.com/nsr/article/7/8/1306/5849029) and [GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117891).
- Feng et al. (2023), [Supplementary Table S5](https://journal.hep.com.cn/fmd/fileup/2095-0217/SUPPL/20230110084242_1.pdf).
- Hegi et al. (2005), [MGMT methylation and temozolomide benefit](https://doi.org/10.1056/NEJMoa043331).
- Duke, [Zenodo dataset record 20548528](https://zenodo.org/records/20548528) and [dataset README](https://github.com/deshanyang/GBM-DIR-QA).
- Wen et al. (2023), [RANO 2.0](https://pmc.ncbi.nlm.nih.gov/articles/PMC10860967/).
- Ntafoulis et al. (2023), [patient-derived glioblastoma TMZ-response study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10575865/).
- Wang et al. (2023), [glioblastoma TMZ-sensitivity study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10010007/).
- [MRI–RNA IvyGAP study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13114458/) and Park et al. (2025), [spatial MRI/pathology/RNA study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11726254/).
- Full evidence synthesis: [Track A/B resistance signal evidence](trackABC_resistance_signal_evidence_2026-11-07.md).
