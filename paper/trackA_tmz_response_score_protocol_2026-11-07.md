# Track A 121-Gene TMZ-Response Direction Score Protocol

**Status:** Frozen for one exploratory analysis on 2026-11-07, before calculating any new Track A proxy scores. This protocol does not change the frozen A5 classifier, its outputs, or its five-of-six success rule.

## Question and interpretation

Among the six A5 donors with paired Core and Periphery samples, do the regions differ in a fixed, directional expression score based on genes associated with temozolomide (TMZ) viability response in patient-derived glioma stem-like cell cultures?

This is an exploratory molecular response-associated proxy. It is not a published or validated 121-gene predictor, a direct drug-response measurement, or a classification of any donor as TMZ-sensitive or TMZ-resistant. A positive score means only that expression ranks align more with genes associated with higher in-vitro viability at 100 μM TMZ.

## Source and fixed gene panel

The panel will contain all 121 genes in Ntafoulis et al. (2023), Supplementary file 1, sheet `Correlation_cell_Overlap`. The panel is copied to [`trackA_tmz_response_121_gene_panel_2026-11-07.tsv`](trackA_tmz_response_121_gene_panel_2026-11-07.tsv). Source workbook: `41416_2023_2402_MOESM3_ESM.xlsx`, PMC accession PMC10575865, version 1, DOI 10.1038/s41416-023-02402-y. The PMC metadata reports CC BY and an MD5 of `d4c65e85ad84c441856590f355a591e3`; the retrieved 842,648-byte workbook matched that MD5. Its SHA-256 was `a21740f50dfc096550b2576f97bff90b70eb6ade62eeba8d25acd3f2b092c7ce`.

The workbook contains 121 unique overlap genes with 48 positive and 73 negative Spearman coefficients against in-vitro cell viability at 100 μM TMZ. Positive coefficients define the `resistance_associated` direction; negative coefficients define the `sensitivity_associated` direction. The score will use these source directions, not survival hazard ratios. No gene will be selected, removed, or reweighted using the Yu et al. regional expression results.

The earlier description of a “21-gene externally validated signature” is withdrawn. The paper's 121-gene overlap joins in-vitro response-associated genes to tumor-survival-associated genes. Its separate external survival checks report 16 nominally significant genes in GLIOTRAIN and 6 in TCGA, with only ZNF540 shared; they do not validate a single 21-gene response score. Those survival checks will not set score weights or define donor response.

## Cohort and input rules

The analysis will use the same six strict A5 donors: GS1, GS2, GS3, GS5, GS6, and GS13. The A5 loader will verify its frozen `SHA256SUMS`, exact cell and point joins, Table S1 labels, and raw-count QC. It will retain cells with at least 1,000 UMIs and at least 500 detected genes. Core means Tumoral; Periphery means Peritumoral; Normal and invalid joins remain excluded. Each donor must have at least 20 QC-passed cells in both regions. If any input hash, join, label check, donor gate, count check, or panel mapping fails, the analysis stops without a result file.

The 121 source symbols must map one-to-one to exact gene identifiers in the A5 count matrix. No aliases, orthologs, imputation, or post hoc replacements will be used. If any source symbol is missing or duplicated, the proxy score will not be calculated. The frozen A5 result remains intact and separate.

## Score calculation

1. For each donor and region, sum raw UMI counts across all eligible cells to create one donor-region pseudobulk. The six donors, not the cells, are the analysis units.
2. For each pseudobulk, divide every gene count by that sample's total counts across all measured genes, multiply by 1,000,000, then calculate `log2(CPM + 1)`.
3. Rank all measured genes within each pseudobulk by `log2(CPM + 1)`, using average ranks for ties. Convert ranks to percentiles with `(rank - 1) / (G - 1)`, where `G` is the number of measured genes.
4. Calculate the response-direction score as:

   `mean(percentile rank of the 48 positive-coefficient genes) - mean(percentile rank of the 73 negative-coefficient genes)`.

   The score ranges from -1 to 1. A larger value means the resistance-associated group ranks higher relative to the sensitivity-associated group across the transcriptome. It does not estimate viability, drug effect, or clinical response.
5. For each donor, calculate `Core score - Periphery score`. Report all six paired differences, their mean, median, range, and numbers above, below, or equal to zero.

## Matched gene-set controls

Generate 1,000 control panels. The background will be all measured genes except the fixed 121-gene panel. For every measured gene, calculate mean per-cell `log2(CPM + 1)` and the fraction of eligible cells with a nonzero count, pooling cells without using Core/Periphery labels. Standardize these two features across measured genes, with a zero standard deviation replaced by 1.

For each control panel, process the 121 source genes in a seeded random order. For each source gene, select one unused background gene uniformly from its 20 nearest available genes by standardized Euclidean distance in the two expression features; ties are resolved by exact gene identifier. Assign the matched control gene the source gene's positive or negative direction. Each control panel must contain 121 unique genes, with 48 positive and 73 negative directions. Use NumPy `PCG64DXSM` seed `20261107`. If a complete unique panel cannot be formed, stop and report the control procedure as incomplete.

Apply the same pseudobulk ranking score to every control panel. Report the observed donor-mean paired difference beside the median and 2.5th–97.5th percentiles of the 1,000 control-panel donor-mean differences. This is a descriptive matched reference, not a validation set or a significance test.

## Uncertainty and decision rule

Calculate one two-sided exact sign-flip p-value for the mean of the six paired donor differences by enumerating all 64 donor-level sign assignments. Treat this p-value as exploratory because the cohort has six donors and the score is a proxy. Cells will never be treated as independent replicates. No multiplicity-adjusted or cell-level inference will be made.

There is no pass/fail threshold and no favored direction. The analysis is complete only if all six donor pairs, the exact sign-flip calculation, and all 1,000 matched panels are available. Otherwise report it as incomplete. Report positive, null, or negative patterns unchanged. Do not use the result to tune A5, assign patient drug response, calibrate Track C signal quality or dependence, or link Track A donors to Duke patients.

## Provenance and output

The run will record the A5 input hashes, exact donor and cell counts, gene mapping, source panel hash, code and protocol hashes, Python, operating-system, NumPy, SciPy, pandas, scikit-learn, and PyTorch versions, seed, all donor-level scores, matched-control summary, and exact sign-flip result. It will save summary metrics only, refuse to overwrite an existing result, and leave raw count data local. The output path will be `output/trackA_tmz_response_proxy/result.json`. No figure will change in this analysis.

## Source

Ntafoulis, I., et al. (2023). Ex vivo drug sensitivity screening predicts response to temozolomide in glioblastoma patients and identifies candidate biomarkers. *British Journal of Cancer, 129*(8), 1327–1338. https://doi.org/10.1038/s41416-023-02402-y. Article and supplementary-file metadata are available through [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10575865/) and the [PMC AWS article record](https://pmc-oa-opendata.s3.amazonaws.com/PMC10575865.1/PMC10575865.1.json).
