# Track A 121-gene response proxy mapping gate

**Date:** 2026-11-07  
**Protocol commit:** `168af12`  
**Status:** stopped before scoring; frozen rules unchanged

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
