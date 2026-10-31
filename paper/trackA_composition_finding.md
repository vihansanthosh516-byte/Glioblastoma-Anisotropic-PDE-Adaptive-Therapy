# Track A finding: bulk core-vs-edge differences are explained by tissue composition (2026-10-30)

Exploratory; three sealed tests; bulk tissue; no clinical claim. Every number below is from the cited result file.

## Result
In three pre-registered, sealed tests on two independent cohort sets, a gene program differed between labelled
tumour core and tumour edge within patients. In all three, the difference disappeared after the pre-set adjustment for
tissue composition (hypoxia, normal-brain neuron content, and, in test 3, myeloid abundance).

| Test | Cohort | n | Raw within-patient test | After composition adjustment | Verdict | Source |
|---|---|---|---|---|---|---|
| 1. SenMayo senescence | Ivy GAP 8 + Angers E-MTAB-5092 7 | 15 | core > edge, p = 0.0035; both cohorts same sign | hypoxia-adjusted p = 0.37 | FAIL: explained by hypoxia | `output/trackA_senmayo/result.json` |
| 2. Hallmark OXPHOS (alpha 0.025) | same 15 | 15 | edge > core, p = 0.0024; both cohorts same sign | hypoxia + neuron p = 0.67 | FAIL: explained by neuron / normal-brain content | `output/trackA_oxphos/result.json` |
| 3. TAM activation (script-157 list) | E-MTAB-8743, 5-ALA-guided invasive vs core | 10 | core > invasive, p = 0.0039 | hypoxia + neuron + myeloid p = 0.44 (myeloid coefficient 0.72) | FAIL: explained by myeloid abundance | `output/trackA_tam_5ala/result.json` |

Tests: exact sign-flip over patients, with expression-matched random gene-set nulls; protocols were sealed by hash
before scoring (`paper/trackA_senmayo/`, `paper/trackA_oxphos/`, `paper/trackA_tam_5ala/`).

## What this means
- Raw core-vs-edge gene differences in bulk GBM tissue were strong and replicated in direction across cohorts, but
  each was accounted for by differences in what the tissue is made of: more hypoxic tumour in the core, more normal
  brain at the edge, more myeloid cells in the core.
- A bulk core-vs-edge signature should not be read as a change inside tumour or immune cells unless it survives
  composition adjustment.

## Limits (state these next to the result)
- Small n: 15 and 10 patients. "Not separable" is not "no effect". Adjusted tests have low power, and the covariates
  are partial controls only.
- Covariate lists are fixed published sets (Hallmark hypoxia, Lein neuron markers) plus an assumed 8-gene pan-myeloid set.
- Test 3 uses the authors' count matrix; how it was made is not documented (`paper/open_data_readiness_addendum_2026-10-30.md`).
- Three tests were run; tests 1 and 2 shared one cohort (alpha split). No correction across the three is applied,
  because all three are reported as failures.
- Prior art: composition confounding of bulk tumour profiles is well known (e.g. purity effects; Ivy GAP anatomic
  structures). The contribution is a sealed, replicated demonstration in labelled GBM core/edge cohorts, not a new
  concept. Novelty not claimed.

## What would be needed for a cell-level claim
Single-cell or single-nucleus data, or same-section spatial data, with region labels made without expression, in at
least 10 donors. No open dataset that meets this was found (`paper/open_data_readiness_addendum_2026-10-30.md`).
