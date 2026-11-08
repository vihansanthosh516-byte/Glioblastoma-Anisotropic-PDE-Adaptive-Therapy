# Track A TMZ-response proxy amendment: 119 of 121 mapped genes

**Protocol status at freeze:** adopted by the user on 2026-11-08. The first invocation stopped before scoring because the frozen A5 loader returned six strict plus two lenient donors; the amended scorer filters to the exact six strict donors and rejects unexpected IDs. After this correction, 9 focused synthetic tests passed in a scratch copy. At the time of the frozen protocol version, no proxy values had yet been calculated.

**Post-run record (2026-11-08):** The 119/121 exploratory score is recorded locally and remains uncommitted at `output/trackA_tmz_response_proxy_119of121/result.json` for the six strict donors. The mean paired Core-minus-Periphery score was -0.00109; donor directions split 3 positive and 3 negative; the exact two-sided donor sign-flip p-value was 0.96875. The observed mean was at the 36.5th percentile of the descriptive matched-control means. This shows no consistent regional shift in this small exploratory analysis. It remains an expression proxy associated with in-vitro TMZ viability, not a measured or patient-specific drug response. The result record's `amendment_sha256` is `31ba941c679054ab1be86394e2fe9eb10451f4db81e93dc95b8424d527ad46b6`, the hash of the frozen amendment bytes used by the run. This post-run note changes status documentation only; the frozen scoring rules and A5 results remain unchanged.

## Reason for amendment

The original 2026-11-07 protocol required exact matrix mapping for all 121 source-panel genes and correctly stopped before scoring when two source symbols were absent. The source panel remains unchanged. This amendment authorizes a separate 119-of-121 mapped-gene proxy; it does not revise or replace the original protocol or the A5 classifier.

The two excluded source symbols are `AF131215.2` (sensitivity-associated) and `PPIAP22` (resistance-associated). Neither has an exact symbol in the A5 matrix. The reported alternate names and Ensembl cross-reference also have no exact matrix entry. No alias, ortholog, imputation, substitution, or post hoc gene replacement is permitted.

## Frozen score

- Retain the same six A5 donors, cell and point joins, labels, QC filters, minimum cell counts, donor-level pseudobulk calculation, transcriptome-wide expression ranking, paired Core-minus-Periphery contrast, exact donor sign-flip procedure, and 1,000 expression-matched control panels from the original protocol. The frozen A5 loader may return its six strict and two lenient donors together; the amended scorer will explicitly subset the six strict primary donors before any score calculation. GS11 and GS12 will not enter this score.
- Use the original Ntafoulis et al. viability-correlation directions and equal weighting for the 119 exact-mapped genes. The mapped panel contains 47 positive-coefficient and 72 negative-coefficient genes.
- Record source-panel coverage as **119/121 exact matrix symbols** and list both omitted symbols in the result. Matched controls will have 119 genes, with the same 47/72 direction counts.
- If any missing symbol differs from the two named symbols, any duplicate identifier appears, or any original cohort/data gate fails, stop without writing a result.
- Report all six donor-level scores and the descriptive matched-control summary. Keep the exact sign-flip result exploratory. There is no success threshold, preferred direction, or pass/fail rule.

## Interpretation and output

This score is an exploratory mixed-cell expression proxy associated with in-vitro TMZ viability. It is not a measured drug response, validated predictor, clinical resistance label, or matched patient-specific response. Do not use it to tune A5, calibrate Track C, link Track A donors to Duke patients, or claim treatment benefit.

The exploratory result is written to `output/trackA_tmz_response_proxy_119of121/result.json`. The original output path and all A5 results and figures remain unchanged. The runner refuses to overwrite an existing result and records both the original A5 freeze and the amendment hash used for scoring.

## Source

Ntafoulis, I., et al. (2023). *Ex vivo drug sensitivity screening predicts response to temozolomide in glioblastoma patients and identifies candidate biomarkers.* British Journal of Cancer, 129(8), 1327–1338. https://doi.org/10.1038/s41416-023-02402-y. The project uses the 121-gene `Correlation_cell_Overlap` panel and original in-vitro viability-correlation directions described in the original protocol.
