# Three-track research execution decisions

2026-10-30. This is an operational decision record, not the science-fair submission. The earlier experiments and unsuccessful tests remain part of the project history. No outcome or novelty rating is guaranteed.

## A. Molecular and spatial analysis

**Decision:** retain the original core-versus-invasive spatial question unless the project explicitly adopts a different question. The E-MTAB-8743 analysis is already exposed and the prespecified adjusted test did not pass. It is not an untouched validation dataset. Its adjusted estimate retained the same direction, and the test does not establish why the unadjusted result changed.

**Next work:** verify a publicly downloadable, reusable dataset's donor and region labels against its exact expression files. Record independent paired-donor counts, label provenance, prior exposure and preprocessing leakage. Check the smaller metadata and crosswalks first. Use donor-level splits and statistics, keep all preprocessing and candidate selection inside training data, and compare the proposed method with a simple prespecified baseline. Apply the existing six-donor descriptive and ten-donor inferential feasibility thresholds; these do not guarantee power.

**Alternative, not an automatic substitution:** GBM CARE provides diagnosis/recurrence data, not core/edge labels. Establish permitted use, exact paired samples, prior exposure and a separate clinically meaningful endpoint before proposing it. The spatial six- and ten-donor feasibility thresholds do not automatically apply to this alternative; it needs its own sample-size and precision assessment. Do not present a different cohort or endpoint as independent validation of the original spatial claim. Do not download a large bundle before its metadata and conditions pass.

**Stop rule:** if no source satisfies the necessary gates, report feasibility limitations rather than scoring a convenient but mismatched dataset or presenting a failed result as positive.

## B. PDE tumor-core forecasting

**Decision:** retain the PDE growth-forecast question; the annotation-change study is separate. The 38 CFB metadata candidates have three expected brain-mask paths each, but no candidate is cleared for scoring. Source descriptions of brain-mask construction require reconciliation. Exact acquisition-day eligibility, target definition, registration quality and valid scoring support remain unresolved.

**Next work:** freeze a metadata manifest recording file hashes, visit paths, label provenance and weekly timing, and reconcile the source-versus-package differences in reported t2 GTV counts. Document released-image reference spaces and brain-mask construction. On non-target files, check brain-mask geometry and define common support prospectively if masks differ; common brain support does not prove label coverage. Then define and review a target-blind, MRI-based longitudinal registration assessment using non-tumor landmarks, predefined failure criteria and a case-level review record, before examining target annotations for performance. Separately establish what regions the tumor labels cover, how binary GTV relates to the original target, and how invalid or ambiguous negative voxels are excluded. Check whether exact acquisition dates exist; otherwise propose a prospective week-resolution protocol amendment rather than presenting approximate dates as exact. Record each case's eligibility and exclusions without looking at predictive outcomes.

Only after the source, reuse, support, registration and target-definition gates pass may a frozen manifest authorize scoring. Compare the PDE forecast with the matched-volume no-PDE signed-distance baseline using patients as independent units and the same target, support, horizon and exclusions. A changed support, target, cohort or horizon requires an explicit amendment and cannot be described as the original test. Registration quality alone does not establish valid tumor-label negatives.

**Stop rule:** if source-backed label coverage or defensible support cannot be established, report the analysis as blocked. No fresh download is currently needed for the metadata checks.

## C. Treatment-policy simulation

**Decision:** forced holidays are not the leading successful method. Script 176's development results show a worse worst-10% outcome than front-loaded dosing in all three reported suites. The prior sealed H1 failures remain failures; simulated findings do not establish clinical benefit or a measured biological recovery rate.

**Next work:** a switching/reversion study changes the existing model and needs a documented protocol amendment; it is not another look at sealed H1. First decide whether the new question is mechanism attribution or policy selection. Define any new method and its selection procedure using development scenarios only. Compare it with front-loaded dosing, the old containment method and a competitive simple schedule chosen using the same development budget. Use matched tumor parameters and noise streams across policies and the following mechanistic conditions: neither switching nor reversion, switching without reversion, and switching with reversion. Record each scenario's outcomes, actual courses and treatment dates; a shared 12-course maximum is not equal actual exposure. Include technical failure handling.

Before running an independent test, freeze the selected method, outcome definitions, paired uncertainty intervals, simulation scenario as resampling unit, multiplicity and suite-level decision rule, sample size justification and independent random seeds. Disclose the previous minus-seven-day margin and clearly label a new test as motivated by exploratory work; a different experiment does not pass the old margin or rewrite H1. Estimate any reversion interaction with a paired difference-in-differences; do not infer an interaction from one significant comparison and another nonsignificant one. Report mean and worst-10% outcomes together, even if they disagree.

The reported recovery of drug sensitivity after at least 30 days off temozolomide in two cell cultures supplies a biological motivation only. It does not measure a human recovery rate or validate a treatment schedule. Review relevant prior art before claiming that any proposed method is new.

## Research-plan deadline

Write the full science-fair procedures in future tense, with achievable goals, specified comparisons, contingencies and expected decision rules. Separately disclose completed work and results as required by the fair. Identify the fair's required six sections and approval rules before treating a submission as final. A feasible plan can describe intended experiments while accurately stating present access and evidence uncertainties; it must not promise positive results.
