# Track A 121-gene response proxy mapping gate

**Date:** 2026-11-07  
**Protocol commit:** `168af12`  
**Status:** original 121-gene run stopped before scoring. On 2026-11-08, the user approved a separate 119-of-121 mapped-gene amendment. The original protocol remains unchanged; the amended analysis has not been scored.

See [`trackA_tmz_response_119of121_amendment_2026-11-08.md`](trackA_tmz_response_119of121_amendment_2026-11-08.md). The new runner uses the original source panel and exact symbols only. It omits `AF131215.2` and `PPIAP22` without substitution, records 119/121 coverage, and writes to a separate result path. The original A5 classifier, result, and figures remain unchanged.

The first real invocation on 2026-11-08 passed the source seals and loaded the frozen A5 data, then stopped before scoring because the loader returned the six strict plus two lenient donors. No result file was written and no proxy score was calculated; counts for those eight loader donors were read. The corrected runner now explicitly filters to the six strict donors and rejects unexpected IDs. Its 11-file amendment manifest verifies in scratch, where all 9 focused tests pass. Commit the correction before the one exploratory score.

## Check

The frozen proxy runner verified its seven-file SHA-256 manifest, loaded the frozen A5 count matrix, and compared the 121 source-panel symbols with the matrix's 24,990 gene identifiers. Exact mapping found 119 genes and missed two:

- `AF131215.2`
- `PPIAP22`

The second symbol is the official NCBI gene symbol for PPIAP22, Gene ID 653214; the NCBI record lists `PPIA3L` and `PPIAL3` as alternate names. None of those three names occurs in the A5 matrix. The first source symbol is associated with Ensembl identifier `ENSG00000255310`; neither that identifier nor `AF131215.2` occurs in the A5 matrix. These cross-references do not supply counts for absent genes.

## Outcome

The runner stopped under the frozen exact-mapping rule. It did not calculate donor scores, matched-control scores, or a sign-flip result. `output/trackA_tmz_response_proxy/result.json` was not created. No A5 model, endpoint, result, or figure changed. The A5 count matrix was loaded locally for the mapping check; no count values were saved or reported.

No gene was substituted, imputed, or dropped. A 119-gene score would be a different proxy and is not authorized by the frozen 121-gene protocol. Any alternative would need a new dated protocol before scoring. The proxy therefore remains unavailable for this matrix under the current rule.

## Sources

- Ntafoulis, I., et al. (2023). *Ex vivo drug sensitivity screening predicts response to temozolomide in glioblastoma patients and identifies candidate biomarkers.* [Article and supplementary information](https://pmc.ncbi.nlm.nih.gov/articles/PMC10575865/).
- National Center for Biotechnology Information. [PPIAP22, Gene ID 653214](https://www.ncbi.nlm.nih.gov/gene/653214).
- Ensembl. [AF131215.2, ENSG00000255310](https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000255310).
