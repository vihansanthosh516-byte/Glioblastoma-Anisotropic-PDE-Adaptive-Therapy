# Track A sealed protocol 2: OXPHOS program, edge vs core (frozen 2026-10-29)

Second locked test on the same 15 patients. The first (SenMayo, commit fc647e4) failed its hypoxia rule.
This is test 2 of 2 on this cohort. Alpha is split: 0.025. No third list will be tested on this cohort.
Expression data were already read for test 1; the OXPHOS score was never computed or looked at.

## Question
Within each patient, is the Hallmark oxidative phosphorylation program different between labeled edge and core,
beyond expression-matched random sets, beyond the hypoxia gradient, and beyond neuron content? Two-sided.

## Lists (frozen, from MSigDB pages read 2026-10-29)
- oxphos_genes.txt: HALLMARK_OXIDATIVE_PHOSPHORYLATION, 200 genes. Overlap with hypoxia list: GPI, LDHA.
- hallmark_hypoxia_genes.txt: same file as test 1.
- neuron_genes.txt: LEIN_NEURON_MARKERS, 69 genes (Lein 2007, PMID 17151600). Covariate for normal-brain neuron content.

## Method
Same cohort, preprocessing, null (1,000 sets, seed 20261029, decile-matched), per-patient E, cohort RMS scaling and
exact 2^15 sign-flip as src/173_trackA_senmayo_edge_core.py.
Adjusted test: E_adj = residual of no-intercept least squares of E on [H, N] across 15 patients
(H = hypoxia E, N = neuron E), then cohort RMS scaling and sign-flip.

## Success (all required)
1. Unadjusted p < 0.025.
2. Ivy GAP and Angers mean E same sign.
3. Adjusted (hypoxia + neuron) p < 0.025, same sign as unadjusted.
If 1-2 pass and 3 fails: report "explained by hypoxia / neuron content", not new biology.

## Known risks
- Neurons are OXPHOS-rich; edge is mostly brain. Neuron covariate is a partial control only (other brain cells too).
- 15 patients; regression with 2 covariates uses 2 degrees of freedom.
- Prior art: 2 web searches found no direct OXPHOS edge-vs-core test in labeled human GBM; CKB/phosphocreatine at the leading edge (Sci Adv 2025) is close in theme. Not exhaustive.
