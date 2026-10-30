# Track A sealed protocol: SenMayo senescence program, edge vs core (frozen 2026-10-29, before any expression download)

## Question
Within each patient, is the SenMayo senescence/SASP program different between the labeled tumor edge and the core, beyond random gene sets matched on expression, and beyond the hypoxia gradient?
Two-sided. No direction is claimed in advance.

## Signature (frozen)
- File: senmayo_genes.txt, 125 symbols, copied from MSigDB gene set SAUL_SEN_MAYO (page read 2026-10-29).
- Source: Saul et al. 2022, Nat Commun, PMID 35974106. Built in bone and other tissues, not in glioma.
- Symbols missing on a platform are dropped and the count is reported. No other edits.

## Cohort (15 patients; labels not from expression)
- Ivy GAP, 8 donors: 10865, 10926, 12108, 12109, 12110, 12112, 14734, 292023109. Edge = LE-reference-histology; core = CT-reference-histology. All other structures are excluded.
- Angers E-MTAB-5092, 7 patients: FT02, FT04, FT05, FT06, FT08, FT09, FT13. Edge = interface; core = tumor zone. MRI and neuronavigation labels, pathologist-reviewed (Aubry 2015).
- Replicates are averaged within each zone.

## Endpoint and test
1. Expression: log2 scale; z-score each gene within each cohort.
2. Score = mean z over signature genes present.
3. Null: 1,000 random sets, seed 20261029, same size, matched on mean-expression deciles (within each cohort), drawn from the genes not in the signature.
4. Per patient: E = (edge - core signature score) - median over null sets of (edge - core).
5. Each E is divided by the uncentered RMS of E within its cohort.
6. Test: exact two-sided sign-flip over all 2^15 patterns, statistic = mean scaled E. p floor = 2/32768.
7. Success needs all three: p < 0.05, Ivy GAP mean E and Angers mean E have the same sign, and the hypoxia-adjusted test (below) also gives p < 0.05 with the same sign.

## Hypoxia control (required)
- Hallmark hypoxia (MSigDB HALLMARK_HYPOXIA) edge-minus-core score, computed the same way, is regressed out: E_adj = E - b * H, with b from a no-intercept least-squares fit across the 15 patients. Then repeat step 5-6 on E_adj.
- If unadjusted is significant but adjusted is not, report "explained by the hypoxia gradient", not new biology.
- Hallmark hypoxia gene list must be frozen in this folder before expression download.

## Pre-declared secondaries (cannot rescue a null primary)
- Angers peripheral brain zone minus tumor zone (9 patients).
- Ivy GAP-only and Angers-only exact tests.

## Known risks (stated before the run)
- Composition: the edge is mostly normal brain. SASP genes are also made by immune and vascular cells. A difference may reflect cell mix, not tumor-cell senescence. The test can not separate these. Results will say "tissue-level senescence/SASP program", not "senescent tumor cells".
- Prior art: limited web search (3 queries, 2026-10-29) found no SenMayo or senescence test by labeled edge vs core in GBM. Salam 2023 (Nat Commun, PMC9883514) saw senescent cells near proliferation and necrosis in mouse tumors only; no patient regional test. Not exhaustive: no PubMed or citing-paper sweep.
- Odds (judgement): positive unadjusted likely (core > edge, driven by SASP + hypoxia overlap); positive after hypoxia adjustment about 0.3.

## Nothing changes after the seal commit.
