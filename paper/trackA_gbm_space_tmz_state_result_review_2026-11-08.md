# Track A GBM-Space TMZ-state score: result review

**Run status:** Completed once under the sealed exploratory protocol on 2026-11-08. The aggregate machine-readable output remains local at `output/trackA_gbm_space_tmz_state/result.json`. The source H5AD and cell-level data were not added to the repository. No figure or A5 result changed.

## Frozen analysis and validation

The run used 12 primary IDH-wildtype GBM-Space donors and the predeclared malignant developmental-like versus gliosis/hypoxia state groups. All donors met the 200-nucleus minimum in both groups; the smallest group had 1,215 nuclei. The source count matrix had the frozen shape of 1,025,329 nuclei by 36,601 features. The one-pass checks passed for every stored `raw.X` value: finite, nonnegative, integer-valued within tolerance, with row sums matching `obs.total_counts`.

The score used 115 exact shared genes from the fixed 121-gene in-vitro TMZ-viability panel: 47 positive-viability directions and 68 negative-viability directions. No aliases or substitutions were used. Six genes were omitted because they were not represented in the exact cross-dataset intersection.

The run used the sealed commit `b2c4bd08881dccbf2ba1076c86a1e908796417e8`, input SHA-256 `9f20ab302091b1ea41016c3bee85d41b3e81cb7e6fc15f6bf4a10555bde414a5`, and completed in 176.74 seconds. The hash fixes the local input bytes; it is not an authoritative upstream checksum.

## Primary result

The frozen estimand was the equal-donor mean of `gliosis/hypoxia score − developmental-like score`. The mean difference was **−0.01128** (95% donor-level t interval **−0.01706 to −0.00549**). All **12 of 12** donor differences were negative. The predeclared exact two-sided sign-flip check returned **p = 0.000488** across 4,096 assignments; this is exploratory and is not a treatment-response test or pass/fail rule.

| Donor | Developmental nuclei | Gliosis/hypoxia nuclei | Developmental score | Gliosis/hypoxia score | Difference |
|---|---:|---:|---:|---:|---:|
| AT3 | 21,967 | 36,583 | −0.054549 | −0.073951 | −0.019403 |
| AT4 | 32,066 | 11,077 | −0.059591 | −0.072519 | −0.012927 |
| AT5 | 34,872 | 2,196 | −0.065335 | −0.083624 | −0.018288 |
| AT6 | 49,607 | 27,659 | −0.088754 | −0.094670 | −0.005916 |
| AT7 | 10,182 | 5,917 | −0.064823 | −0.072482 | −0.007659 |
| AT9 | 101,133 | 9,834 | −0.053610 | −0.084950 | −0.031340 |
| AT10 | 39,654 | 8,249 | −0.090319 | −0.109375 | −0.019056 |
| AT11 | 37,020 | 1,215 | −0.123491 | −0.132662 | −0.009172 |
| AT12 | 6,873 | 25,958 | −0.126039 | −0.127534 | −0.001495 |
| AT13 | 14,083 | 3,753 | −0.073289 | −0.076961 | −0.003672 |
| AT14 | 10,426 | 8,793 | −0.059654 | −0.063067 | −0.003413 |
| AT15 | 106,390 | 7,240 | −0.036917 | −0.039879 | −0.002962 |

## Matched-control context

Across 1,000 expression-matched control panels, the median donor-mean difference was **−0.00170**, with a 2.5th–97.5th percentile range of **−0.01287 to 0.00802**. The observed panel result was at the **4.3rd percentile** of this descriptive control distribution. The observed result remains inside the central 95% control range. These controls do not form a second hypothesis test.

## Interpretation and limits

The fixed score is the mean percentile rank of genes associated with higher in-vitro viability minus the mean percentile rank of genes associated with lower in-vitro viability. A negative contrast means that this expression-derived score was lower in the author-defined gliosis/hypoxia group than in the developmental-like group in each of these 12 tumors.

This is an association between two expression-derived quantities. The state labels and score both depend on transcriptomic data, and this atlas has no matched TMZ dose-response or clinical treatment outcomes for these nuclei. The result does **not** establish TMZ sensitivity or resistance, patient-specific response, clinical benefit, causation, or a new drug-response method. The magnitude is small on the score scale, and the observed difference falls within the matched-control 95% range. Treat the result as exploratory; do not use it to alter Track B or Track C parameters.

The GBM-Space portal lists its single-nucleus transcriptomics dataset as available through CELLxGENE. CELLxGENE's public-data policy states that public datasets are offered under CC BY 4.0 with attribution. This report cites the source and shares aggregate donor-level summaries only; it does not redistribute the H5AD or cell-level data.

## Sources

- de Jong, G., Memi, F., Gracia, T., et al. (2026). *A spatiotemporal cancer cell trajectory underlies glioblastoma heterogeneity*. Nature. https://doi.org/10.1038/s41586-026-11234-8
- Ntafoulis, I., et al. (2023). *Ex vivo drug sensitivity screening predicts response to temozolomide in glioblastoma patients and identifies candidate biomarkers*. British Journal of Cancer, 129(8), 1327–1338. https://doi.org/10.1038/s41416-023-02402-y
- Suter, R. K., et al. (2026). *Drug and single-cell gene expression integration identifies sensitive and resistant glioblastoma cell populations*. Nature Communications, 17, 99. https://doi.org/10.1038/s41467-025-67783-5
- GBM-Space portal. https://www.gbmspace.org/
- Chan Zuckerberg Initiative. *CELLxGENE data submission policy*. https://cellxgene.cziscience.com/docs/032__Contribute%20and%20Publish%20Data
