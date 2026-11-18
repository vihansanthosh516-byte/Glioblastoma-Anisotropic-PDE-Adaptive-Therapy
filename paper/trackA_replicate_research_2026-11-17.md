# Track A research leads after the support threshold simulation

The [Artzi benchmark](trackA_artzi_result_2026-11-17.md) and [support threshold simulation](trackA_support_simulation_result_2026-11-17.md) motivate changing the measurement model. Keep Artzi, Yu and HGCC spent.

## Lead candidate: replicate-smoothed gene order

For each region, build pseudobulks for independently sampled tissue blocks. For each gene pair, average its signed order over blocks to get a soft vote for the region. Use equal weight for each training patient's core-to-edge contrast to orient and rank pairs. Average selected soft votes to score each region. At prediction, calculate each region's score from its own blocks. Cells and adjacent fields from one block cannot serve as independent patient observations.

Separate two questions with an ablation: does averaging within-region votes help relative to ordering the pooled counts, and does the effect persist after matching block size? Compare on the same held-out patients and paired concordance against pooled mean-contrast voting, majority rank voting, an average of block scores from conditional logistic regression, and pooled conditional logistic regression. Keep the remaining eligible methods in the [Track A comparator inventory](trackA_prior_art_full_2026-11-16.md) for a future biological benchmark.

Primary prior art includes [Geman et al. on TSP](https://pmc.ncbi.nlm.nih.gov/articles/PMC1989150/), [Rashid et al. on gene-pair effects that vary across studies](https://arxiv.org/abs/1708.05508), and [Smyth et al. on using correlated replicates](https://doi.org/10.1093/bioinformatics/bti270). Averaging block-level votes is a combination of familiar techniques. This targeted search does not establish uniqueness.

The next simulation should generate gene expression with patient effects, correlated blocks and deliberately unequal block depths. Include a condition where the true regional difference is driven by a rare, localized patch, because averaging block votes could erase it. There must be fresh synthetic test draws and a frozen pass rule before running. For real data, require source-backed patient, region and independently sampled block identifiers. Screen Harwood, Gill and Nam from metadata only using the [existing cohort leads](trackA_fresh_data_followup_2026-11-13.md).

## Alternative leads

- Fit a patient-level covariance penalty to avoid selecting several pairs that all fail in the same patients. Compare against a gene-disjoint selector and regularized paired logistic regression. Small patient counts make estimated correlations unstable.
- Use covariate-adjusted gene order as a diagnostic if an untouched study supplies independent pathology estimates of tissue composition. Compare adjusted and unadjusted paired scores. This changes the biological question and needs its own plan. [Covariate-adjusted TSP](https://pmc.ncbi.nlm.nih.gov/articles/PMC9942303/)
