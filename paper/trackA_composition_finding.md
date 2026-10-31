# Track A finding: bulk core-vs-edge differences after gene-score adjustment (2026-10-30)

Exploratory; three sealed tests; bulk tissue; no clinical claim. Every number below is from the cited result file.

## Result
In three tests with protocols fixed before scoring, a gene program differed between labelled tumour core and tumour
edge within patients. In all three, the adjusted test was not statistically significant after the prespecified
adjustment for hypoxia, neuron-related scores and, in test 3, a myeloid-related score. The adjusted point estimates
remained in the same direction. The tests did not meet their conjunctive success rules.

| Test | Cohort | n | Raw within-patient test | After gene-score adjustment | Verdict | Source |
|---|---|---|---|---|---|---|
| 1. SenMayo senescence | Ivy GAP 8 + Angers E-MTAB-5092 7 | 15 | core > edge, p = 0.0035; both cohorts same sign | hypoxia-adjusted p = 0.37 | FAIL: adjusted criterion not met | `output/trackA_senmayo/result.json` |
| 2. Hallmark OXPHOS (alpha 0.025) | same 15 | 15 | edge > core, p = 0.0024; both cohorts same sign | hypoxia + neuron p = 0.67 | FAIL: adjusted criterion not met | `output/trackA_oxphos/result.json` |
| 3. TAM activation (script-157 list) | E-MTAB-8743, 5-ALA-guided invasive vs core | 10 | core > invasive, p = 0.0039 | hypoxia + neuron + myeloid p = 0.44 (myeloid coefficient 0.72) | FAIL: adjusted criterion not met | `output/trackA_tam_5ala/result.json` |

Tests: exact sign-flip over patients, with expression-matched random gene-set nulls; protocols were sealed by hash
before scoring (`paper/trackA_senmayo/`, `paper/trackA_oxphos/`, `paper/trackA_tam_5ala/`).

## What this means
- Raw core-vs-edge differences had the prespecified direction and passed the unadjusted test in each analysis. For
  tests 1 and 2, the two contributing cohorts had matching point-estimate directions.
- None passed its adjusted success criterion. A nonsignificant adjusted test does not prove the difference vanished,
  establish that a particular covariate explains it, or demonstrate no biological effect. With small samples and
  correlated expression-derived covariates, loss of precision and model specification are alternative explanations.
- Bulk expression alone does not identify a cell-intrinsic change. Adjustment by gene scores is an observational
  sensitivity analysis, not a direct measurement of tissue composition or a causal decomposition.

## Limits (state these next to the result)
- Small n: 15 and 10 patients. "Not separable" is not "no effect". Adjusted tests have low power, and the covariates
  are partial controls only.
- Covariate lists are fixed published sets (Hallmark hypoxia, Lein neuron markers) plus an assumed 8-gene pan-myeloid set.
- Test 3 uses the authors' count matrix; how it was made is not documented (`paper/open_data_readiness_addendum_2026-10-30.md`).
- Three tests were run; tests 1 and 2 shared one cohort (alpha split). No correction across the three is applied,
  because all three are reported as failures.
- Prior art: composition confounding of bulk tumour profiles is well known (e.g. purity effects; Ivy GAP anatomic
  structures). These three tests involve different gene programs, so they do not replicate the same hypothesis.
  Novelty is not established.
- The frozen protocols' historical instructions to describe an adjusted nonsignificant result as “explained by” a
  covariate are overstatements. This report corrects that interpretation without changing the sealed protocols or
  recorded results.

## What would be needed for a cell-level claim
Single-cell or single-nucleus data, or same-section spatial data, with region labels made without expression, in at
least 10 donors. No open dataset that meets this was found (`paper/open_data_readiness_addendum_2026-10-30.md`).
