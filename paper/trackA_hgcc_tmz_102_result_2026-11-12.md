# Track A HGCC 102-gene TMZ bridge: result (2026-11-12)

Run once under amendment 1 (`paper/trackA_hgcc_amendment1_102genes_2026-11-12.md`, commit 3343c27, SHA-256 800a0484...6eedad). Code `src/trackA_hgcc_tmz_102.py` (commit 17454d3, SHA-256 b60bfd1e...184e); 5 synthetic tests passed in a scratch copy before the run. Output: `output/trackA_hgcc_tmz_102/result.json` (aggregate statistics only). Runtime 15 s.

## 1. What was tested

Does the fixed Ntafoulis-derived expression score (102 of 119 genes; 45 resistance-associated, 57 sensitivity-associated) rank measured in-vitro temozolomide viability AUC across 83 HGCC patient-derived glioblastoma lines (94 expression profiles)?

## 2. What was found

- Spearman rho = **-0.102**; one-sided permutation p (positive direction) = **0.820** (100,000 permutations); bootstrap 95% interval **-0.327 to 0.128**.
- Verdict under the frozen rule: **no association shown.** The direction was opposite to the prediction, but the interval includes 0, so it is not a significant reversal.
- Secondary, descriptive: MGMT expression rho = 0.108. The score's rho sat at the 31st percentile of 1,000 matched random control panels (control median -0.051, 95% range -0.269 to 0.168). The panel did no better than random gene sets.

## 3. Weakest points

- 17 of 119 panel genes (mostly non-coding) were absent from the HGCC array, so this is a reduced signature.
- The panel came from a different lab, culture protocol, assay, and platform (RNA-seq versus HTA 2.0 arrays).
- n = 83 lines gives about 80% power only for rho near 0.3; a small true effect could be missed.
- Ju et al. (2025) already modelled TMZ response on these data; this is a secondary reanalysis, not a new method.
- HGCC response-data licence is not explicit; only aggregate statistics are reported.

## 4. What this means and what comes next

The 119-gene TMZ proxy used in Track A did not transfer to an independent cell-line panel with measured TMZ response. Together with the earlier null regional result (`paper/trackA_tmz_response_119of121_amendment_2026-11-08.md`: donor split 3/3, p = 0.97), this is evidence that the panel should **not** be used as a drug-response signal in Track A or as an input to Track C. This is a useful negative finding: it closes one proxy path with measured data. A5's failure and the Track B and Track C statuses are unchanged.
