# Track A TMZ-response proxy mapping gate and amended run

**Date:** 2026-11-08
**Original protocol commit:** `168af12`
**Current status:** the original 121-gene analysis stopped under its exact-mapping rule. The user approved a separate 119-of-121 exploratory amendment. The amended score completed on the six strict A5 donors; the original A5 classifier, result, and figures remain unchanged.

## Original 121-gene exact-mapping check

The frozen 121-gene panel was compared with the A5 matrix's gene identifiers. Exact mapping found 119 genes and missed two:

- `AF131215.2`
- `PPIAP22`

The second symbol is the official NCBI gene symbol for PPIAP22, Gene ID 653214; the NCBI record lists `PPIA3L` and `PPIAL3` as alternate names. None of those three names occurs in the A5 matrix. The first source symbol is associated with Ensembl identifier `ENSG00000255310`; neither that identifier nor `AF131215.2` occurs in the matrix. These cross-references do not provide counts for the absent genes. The original 121-gene attempt therefore stopped before scoring.

## Amended 119-of-121 analysis

The 2026-11-08 amendment uses exact matrix symbols for the 119 mapped genes, retains their original in-vitro viability directions, omits the two absent symbols without substitution, and reports 119/121 coverage. The score is an exploratory mixed-cell expression proxy, not a validated response predictor and not measured patient-specific TMZ sensitivity or resistance.

The first runner invocation passed the source seals and loaded the frozen A5 data, then stopped before scoring because the loader returned the six strict plus two lenient donors. It read count values for those eight loader donors but wrote no result and calculated no proxy score. The corrected runner explicitly selects only GS1, GS2, GS3, GS5, GS6, and GS13, and rejects unexpected donor IDs. Its 11-file seal verified after commit `d91259e`; all 9 focused synthetic tests passed in a scratch copy.

The completed result is saved locally and remains uncommitted at `output/trackA_tmz_response_proxy_119of121/result.json`. Its six donor-level paired differences split evenly by direction, and the mean is near zero and within the expression-matched control distribution. This shows no consistent regional shift in the proxy. Exact metrics, all six donor values, and control summaries are in the local result file. The main limits are the small six-donor cohort, two omitted panel genes, mixed-cell expression, an in-vitro-derived direction set, and no matched patient drug-response outcome.

The amendment remains byte-frozen as the pre-score protocol. Its status paragraph records the state at freeze time; this mapping-gate record reports the later completed run.

No gene was substituted or imputed. The amended result did not tune A5, change its frozen five-of-six rule, alter any A5 output, or change any figure. The result JSON remains local and uncommitted.

## Sources

- Ntafoulis, I., et al. (2023). *Ex vivo drug sensitivity screening predicts response to temozolomide in glioblastoma patients and identifies candidate biomarkers.* [Article and supplementary information](https://pmc.ncbi.nlm.nih.gov/articles/PMC10575865/).
- National Center for Biotechnology Information. [PPIAP22, Gene ID 653214](https://www.ncbi.nlm.nih.gov/gene/653214).
- Ensembl. [AF131215.2, ENSG00000255310](https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000255310).
