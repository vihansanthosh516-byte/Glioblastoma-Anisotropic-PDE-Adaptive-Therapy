# Recurrence measurement reliability: locked protocol

Date: 2026-11-16. Status: **locked** when committed together with `configs/recurrence_reliability.json` (`"locked": true`) and before `src/recurrence_reliability/run.py` is first executed. A Git commit is a timestamped freeze, not independent preregistration.

## 1. Question

How often does a patient's inferred primary-to-recurrence change in the malignant-cell mesenchymal-like (MES-like) program depend on a prespecified technical scoring choice, when the program, the cells, and the data are held fixed?

This is a measurement-reliability study in the Wang et al. (2022) reanalysis cohort (GEO GSE174554). It is not a recurrence classifier, a treatment-causality model, or a replication of Wang et al. Closest prior art and novelty limits: `paper/recurrence_reliability_prior_art.md`.

## 2. Cohort

Patients with one untreated primary and one author-described first-recurrence snRNA-seq specimen, both IDH-wildtype, passing the metadata, annotation, batch (`paper/recurrence_reliability_eligibility.md`) and technical QC (`paper/recurrence_reliability_qc_plan.md`) gates. The final list is `final_cohort` in `output/recurrence_reliability/qc_report.json`, whose SHA-256 is in the locked config: **21 patients**, Table 1 pairs 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 21, 26, 29, 30, 31, 32, 33, 34 (QC result appended to the QC plan). The independent unit is the patient. Cells, libraries, bootstrap draws, and procedures never add patients.

Each specimen uses its primary library only. Nuclei: at least 200 detected genes and mitochondrial fraction at most 2.5% (Wang Methods). Malignant nuclei: author CONICSmat-based `Tumor` label. Any nucleus passing filters without a label stops the run.

## 3. Biological construct (fixed)

Neftel et al. (2019) IDH-wildtype meta-modules MESlike1 and MESlike2 (50 genes each; `data/external/neftel2019/IDHwt.GBM.MetaModules.tsv`, SHA-256 `046905b6ba1e336f28c9c4dca9e8a2fa8f05cd011c2508cf2d8a78826545ce25`). Gene mapping by exact symbol to the GEO features (the GEO feature files are symbol-only; MESlike1 50/50 and MESlike2 49/50 mapped, ERO1L absent); genes with zero or several matching features are dropped (QC lists them); a module below 80% mapped stops the study. Per nucleus, MES-like = max(score MESlike1, score MESlike2), the Neftel definition. The program was chosen because it is the published malignant-cell state set in the closest studies, not from any recurrence result. Literature knowledge disclosed: Wang et al. report a cohort-level mesenchymal shift at recurrence.

## 4. Two fixed procedures for the same construct

- **R (reference): Tirosh/Neftel control-gene score.** E = log2(CPM/10 + 1). Reference pool = all filtered malignant nuclei of the final cohort (each nucleus weighted once). Aggregate expression log2(mean CPM + 1) over the pool, genes ranked into 30 equal-size bins. For each module gene, 100 control genes drawn without replacement from its bin; control set = unique union minus module genes; seed fixed. Score = mean relative expression of module genes minus that of control genes, with relative expression centred on the pool mean.
- **A (alternative): UCell rank score.** Within each nucleus, genes ranked by raw count (descending, ties averaged); ranks above 1500 set to 1501; score = 1 − U/(n·1500) with U = Σrank − n(n+1)/2. No normalization or control genes.

Both are monotone summaries of relative MES-like expression in a nucleus; they are on different scales, so only directions and standardized changes are compared, never raw magnitudes. The specimen score is the mean over its malignant nuclei.

## 5. Patient-level change and indeterminate zone

delta[i,k] = mean score(recurrent i, k) − mean score(primary i, k).

SE[i,k] = sqrt(SE_rec² + SE_prim²), each SE the SD of 500 nucleus-bootstrap means of that specimen (per-nucleus scores fixed). Direction: increase if delta > 1.96·SE, decrease if delta < −1.96·SE, otherwise indeterminate (equality is indeterminate). The 1.96 multiplier comes from the declared sampling model, not from recurrence results. Its calibration is checked in S6 using within-specimen split halves, which contain no recurrence information.

Missing data: a patient without a computable score for either specimen stops the run (QC guarantees support). No ties beyond the rule above.

## 6. Primary endpoint

**Opposite-direction discordance** = number of patients whose determinate R direction is opposite to their determinate A direction, divided by all patients in the final cohort (indeterminate patients stay in the denominator). Interval: exact two-sided 95% Clopper-Pearson. Reported with it: the full 3×3 R-by-A table; indeterminate coverage of each procedure; the count where exactly one procedure is indeterminate; and every patient's continuous delta and SE.

Tolerance: 5% opposite-direction discordance, chosen because a directional recurrence claim that reverses in more than 1 patient in 20 under an equally defensible scorer of the same genes is not a usable patient-level measurement. Reading rule, fixed: CI lower bound > 5% → "unreliable beyond tolerance"; CI upper bound < 5% → "reliable within tolerance"; otherwise "inconclusive".

Precision: the QC-frozen cohort has n = 21. With zero discordant patients, the exact one-sided 95% upper bound is 1 − 0.05^(1/21) = 0.133 (two-sided Clopper-Pearson upper 0.161). A confirmatory "reliable within 5%" reading needs at least 59 independent patients with zero discordance (one-sided) or about 72 for the two-sided interval. **This study cannot confirm reliability; it is predeclared as descriptive for that direction.** It can detect gross unreliability: at n = 21 the lower bound exceeds 5% from 4 discordant patients onward (scipy `beta.ppf`). One primary comparison; no multiplicity adjustment.

## 7. Declared secondary analyses (descriptive; no p-values; cannot rescue the primary)

- S0. Continuous patient changes for every procedure: mean and median delta with a 10,000-draw patient bootstrap percentile interval.
- S1. Composition sensitivity: R and A on all filtered nuclei (own pool, bins, controls) versus malignant-only; 3×3 agreement tables.
- S2. Matched cell number: 200 draws of 50 malignant nuclei per specimen without replacement; per patient, fraction of draws whose R delta sign matches the full-data R sign, and fraction where R and A signs agree.
- S3. Library-aware: for SF12243, SF10433 and SF4209, recompute the specimen score from the independent-capture library alone and pooled with the primary library (library kept in cell IDs; pools and controls unchanged); report the patient delta and direction under each.
- S4. Construct sensitivity: Verhaak et al. (2010) Mesenchymal genes (`TCGA.Verhaak.GBM.txt`, SHA-256 `7815db5d...55a0`) under procedure R; agreement with Neftel R. Disagreement here is construct difference, not technical error.
- S5. Equal-patient mean delta versus pooled-nucleus difference (all recurrent minus all primary malignant nuclei) under R. Different estimands; descriptive.
- S6. Calibration of the indeterminate rule: 20 random half-splits of every specimen's malignant nuclei; rate at which half-versus-half differences are called determinate (expected about 5% if the sampling model is adequate).

No other gene panel, predictor, clinical endpoint, or procedure will be added after scoring.

## 8. Seeds, software, runtime, failure rules

Seeds (numpy `default_rng`): control genes 20261116, nucleus bootstrap 20261117, patient bootstrap 20261118, matched depth 20261119, split half 20261120. Python, numpy and scipy versions and all code hashes are written to `result.json`. Runtime limit 3 hours; exceeding it stops the run without partial reporting.

Stop before or during the run if: config not locked; QC report, matrix, feature, barcode or module hash drift; cohort list differs from the config; a module maps below 80%; an unlabelled nucleus passes filters; a run directory already exists. A code bug found after the run requires a dated correction file, regenerated outputs, and disclosure; it does not restore untouched status.

## 9. Outputs

`output/recurrence_reliability/<run_id>/result.json` and `patient_changes.tsv`; results note `paper/recurrence_reliability_result.md`. Raw matrices stay local.

## 10. External validation

None in this protocol. CARE (GSE274546/GSE274548) remains a metadata lead (Phase G gate). Without it, the output is a within-cohort reliability benchmark only.
