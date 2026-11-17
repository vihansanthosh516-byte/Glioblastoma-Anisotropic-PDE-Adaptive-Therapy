# Recurrence measurement reliability: locked result

Date: 2026-11-16. Run `run_20261116_locked_a1` (protocol locked in commit 642f9ac; amendment 1 in commit 103bc0a, made before any score was computed). Outputs: `output/recurrence_reliability/run_20261116_locked_a1/result.json` (SHA-256 `f5ac9aee…30dc`) and `patient_changes.tsv` (SHA-256 `d8eee79a…a0fb`). Runtime 121.9 s; Python 3.14.4, numpy 2.4.4, scipy 1.18.0. Every number below is from `result.json` unless another file is named.

## Primary endpoint (prespecified reading: inconclusive)

21 patients, 42 specimens, 66,510 malignant nuclei in the reference pool. Program: Neftel MES-like (MESlike1 50 genes, MESlike2 49 genes mapped).

| R (Tirosh/Neftel) \ A (UCell) | A increase | A indeterminate | A decrease |
|---|---|---|---|
| R increase | 7 | 1 | **2** |
| R indeterminate | 0 | 0 | 1 |
| R decrease | **1** | 2 | 7 |

- **Opposite-direction discordance: 3 of 21 = 14.3% (exact 95% CI 3.0% to 36.3%).** The CI contains the 5% tolerance, so the fixed reading is "inconclusive relative to tolerance". The point estimate is about three times the tolerance; the data cannot confirm reliability and cannot confirm unreliability at this n.
- Any disagreement (opposite or one indeterminate): 7 of 21 = 33.3% (95% CI 14.6% to 57.0%). Same determinate direction: 14 of 21.
- Indeterminate coverage: R 4.8% (1 patient), A 14.3% (3 patients). Neither procedure looks reliable by withholding calls.
- The three opposite-direction patients (`patient_changes.tsv`): pair 8 (SF1199→SF1343; R −0.197 ± 0.015, A +0.031 ± 0.002), pair 13 (SF3076→SF3243; R +0.034 ± 0.010, A −0.015 ± 0.001), pair 33 (SF10857→SF12008; R +0.216 ± 0.008, A −0.004 ± 0.001). In each, both procedures are far outside their own nucleus-bootstrap zone, so these are not small changes around zero.

Weakest point: n = 21 gives a wide interval; one cohort, one program, and two scorers only.

## Calibration check (S6)

Half-versus-half splits of the same specimen were called determinate in 5.5% (R) and 6.3% (A) of 840 splits each, close to the 5% expected. The nucleus-bootstrap zone is therefore calibrated for **cell-sampling** noise. Yet 3 patients reverse between scorers with both calls determinate. **Interpretation:** the variation from the scoring choice is much larger than nucleus-sampling uncertainty, so a nucleus-level bootstrap interval alone overstates how firm a patient-level recurrence direction is. This is the main methods finding, and it is exploratory in strength because of n.

## Declared secondary analyses (descriptive)

- **S0 continuous change, patient bootstrap 95% CI:** R mean −0.097 (−0.262 to 0.055), median −0.030; A mean −0.005 (−0.022 to 0.011). Neither interval excludes zero. This study was not designed or powered to test the cohort-level mesenchymal shift reported by Wang et al., and this result must not be read as a contradiction of it (malignant-only, different program definition, 21 of 36 patients).
- **S1 composition (whole tissue vs malignant):** under R, 4 of 21 patients reverse direction and 15 agree; under A, 2 reverse and 13 agree. Whole-tissue R mean change −0.149 (−0.281 to −0.026).
- **S2 matched cell number (50 nuclei, 200 draws):** median fraction of draws matching the full-data R sign = 1.0.
- **S3 library-aware:** SF12243 (1,301 malignant nuclei in the second capture): direction unchanged under R and A. SF10433 (132): R becomes indeterminate with the second capture alone (+0.025 vs −0.386 with the primary library); pooled stays a decrease. SF4209 (4 nuclei): uninformative.
- **S4 construct (Verhaak vs Neftel, both under R):** 4 of 21 reverse, 15 agree. This is construct difference, not technical error.
- **S5 weighting:** equal-patient mean R change −0.097 vs pooled-nucleus difference +0.004. The cohort-level sign depends on whether patients or nuclei are weighted equally (different estimands).

## Exclusion flow

27 metadata candidate pairs → 23 after annotation/support gates (pairs 15, 23, 28, 35 out) → 21 after technical QC (pairs 19, 20 out: SF6118 0 and SF9715 4 malignant nuclei after fixed filters). One unlabelled nucleus dropped (SF4810) under amendment 1. See `paper/recurrence_reliability_eligibility.md` and the QC result in `paper/recurrence_reliability_qc_plan.md`.

## Limitations

1. 21 patients; the primary interval is wide; no confirmatory reliability claim is possible.
2. Two scorers only. Other published scorers (AUCell, JASMINE, ANS, ssGSEA/GSVA, z-score) were not run, so this does not rank scorers or cover "all existing methods".
3. One cohort (Wang reanalysis); no external replication (CARE not audited or opened).
4. Malignant labels are the authors' CNV-based operational calls; first recurrence is author-described.
5. GEO feature files are symbol-only; one MESlike2 gene (ERO1L) was not mapped.
6. Results concern measurement agreement, not biology, treatment effect, or causality.

## What this supports

A narrowly stated methods result: in paired GBM snRNA-seq, the direction of a patient's recurrence change in a fixed MES-like program disagreed between two standard scorers in 3 of 21 patients (95% CI 3% to 36%), while a nucleus-bootstrap interval judged each call firm. Patient-level recurrence directions should be reported with scorer sensitivity, not nucleus-sampling intervals alone.

## Next (not started)

1. Prespecified multi-scorer extension (AUCell, JASMINE, ANS, z-score) on the same locked cohort, framed as an exploratory benchmark.
2. Phase G: audit CARE (GSE274546) metadata and annotation linkage; reserve it before any tuning; apply this frozen pipeline once. CARE has 59 patients in the source paper, which would narrow the interval (not yet a verified eligible count).
