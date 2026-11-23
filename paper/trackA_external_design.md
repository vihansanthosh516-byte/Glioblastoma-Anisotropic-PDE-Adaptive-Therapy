# Track A external design and power screen

## Candidate with a real technical challenge

Test same-section segmentation crossover before seeking cross-platform claims. Start from the same detected molecules and morphology, reconstruct masks using the vendor and a fixed alternative, and recalculate the same spatial statistic using a frozen marker assignment. This controls physical section differences and tests sensitivity to segmentation and transcript assignment. It cannot diagnose upstream acquisition errors shared by the two reconstructions.

A candidate uses spatial information only when its incremental paired prediction is stable across qualified reconstructions. Reserve complete patients and a segmentation condition before development. The external test must compare the candidate with composition plus malignant expression and density, nominal spatial prediction, and the fixed alternative segmentation without qualification. Use identical paired concordance and keep all tissue blocks nested in patients. Include shuffled state marks within geometry and depth strata. Stop at metadata if independent region labels, images, molecule coordinates, or eligible patient count are missing. Petti could help develop technical checks, but cannot meet the existing external validation gate. [Petti audit](trackA_petti_acquisition_audit.md)

[Segmentation benchmarking](https://www.nature.com/articles/s41592-026-03272-x), [Proseg](https://www.nature.com/articles/s41592-025-02697-0), and [Xenium contamination analysis](https://doi.org/10.1038/s41592-026-03089-8) already cover related technical ideas. The candidate contribution is testing whether same-specimen segmentation stability qualifies *incremental*, patient-held-out core/edge prediction. This is a research design and not an established novelty or an executable external protocol.

## Exact power, design assumptions only

A minimum patient count that allows a significant result is not a power calculation. For one one-sided paired sign test at alpha 0.05, suppose each of n independent patients has a non-tied difference with probability d. Among these, let q be the candidate win probability. Sum over k drawn from Binomial(n,d). For each k sum Binomial(k,q) probability over win counts w whose upper null tail under Binomial(k,0.5) is at most 0.05. The resulting powers below assume independent patients and fixed d and q. These are hypothetical design values, not patient findings.

| Discordance d | Win rate q | n=10 | n=20 | n=40 | n=60 | n=80 | n=100 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 0.70 | 0.012 | 0.087 | 0.234 | 0.375 | 0.487 | 0.580 |
| 0.40 | 0.75 | 0.078 | 0.263 | 0.574 | 0.758 | 0.872 | 0.934 |
| 0.60 | 0.80 | 0.257 | 0.596 | 0.903 | 0.981 | 0.996 | 0.999 |

These values do not include Holm adjustment, a minimum concordance gain, or uncertainty across cohorts. Joint power for a gate demanding all of these can be lower. Do not use the Artzi floor as a target sample size. Pick a meaningful effect, estimate tie frequency using permitted development data, and size the new study for its actual frozen pass rule. The prior Artzi analysis remains unchanged and spent.
