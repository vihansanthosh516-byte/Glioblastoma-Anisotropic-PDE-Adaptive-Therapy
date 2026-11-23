# Prior-art check: gated sister-culture contrast for TMZ prediction

2026-11-23. Candidate: one-slope within-origin contrast plus replication-gated fallback ([dev4 result](trackA_dev4_tmz_result_2026-11-23.md)).

## How it was checked

Firecrawl research index (PubMed, bioRxiv, medRxiv, arXiv), 8 query framings: multi-region cultures and drug response; within vs between patient effects in drug prediction; same-individual contrastive learning; GBM patient-derived TMZ predictors with several lines per patient; invariant learning with ID grouping; Mundlak / fixed-effects penalized prediction; multisector GSC TMZ; sibling designs for polygenic prediction. Abstracts and matched passages were read. Full texts were not read. Patents, Google Scholar and Web of Science were not searched. This is a structured search, not a complete systematic review.

## Closest prior art

| Work | What it does | Overlap with our candidate |
|---|---|---|
| Guan et al. 2025, Nat Genet, family-based GWAS "unified estimator" (snipar), [doi](https://doi.org/10.1038/s41588-025-02118-0) | Combines within-family sibling variation with unrelated individuals to remove family confounding; predictors built this way predict better out of sample. | **Very high.** Same core idea: within-group differences remove group-level confounding, pooled with between-group data. Ours swaps siblings for sister cultures. |
| Selzam 2019; Lello 2020; within-family polygenic prediction | Sibling differences separate direct effects from family confounding in prediction. | High (same logic, genetics). |
| Heinze-Deml & Meinshausen, CoRe, [arXiv 1710.11469](https://arxiv.org/abs/1710.11469) | Uses an ID variable (same person) to penalize within-ID variation and resist shifts in nuisance features. | High in spirit (grouping by identity to drop non-transferable features); different mechanism (variance penalty, not outcome contrast). |
| Mundlak / generalized Mundlak (Arkhangelsky & Imbens, arXiv 1807.02099); Pavlou 2021 (confounding by cluster in risk prediction); SHEL 2026 (arXiv 2605.24587); ARMED (Nguyen 2023) | Group-mean terms or cluster-invariant features to remove cluster confounding; some for prediction on new clusters. | High (statistical basis). The one-slope contrast reduces to a within/between weighting of this kind. |
| Lloyd et al. 2021 PLoS Comput Biol | Between-tissue vs within-tissue signal in pan-cancer drug-response prediction. | Medium (same decomposition, tissue level, no culture contrasts). |
| Nam et al. 2023 Genome Med, [PMC10010007](https://pmc.ncbi.nlm.nih.gov/articles/PMC10010007/) | 69 GBM GSC lines with TMZ screen and ML predictor; plus 52 multisector GSCs from 18 tumours with TMZ response. | Medium for method (uses sectors for validation and heterogeneity, not as within-origin supervision, per abstract). **High as a data lead.** |
| Salatino et al. 2026 bioRxiv; Brynjulvsen et al. 2023 Cancers; Kim 2022 (CRC); Xu 2026 (breast, leave-one-patient-out) | Multi-region cultures with drug screens. They describe heterogeneity or hold whole patients out; none (per abstract) trains on within-origin contrasts. | Low for method; data leads. |
| Ntafoulis 2023; Skaga 2022; Ju; PERCEPTION; CaDRReS-Sc | Expression-based or functional TMZ response prediction without sister-culture supervision. | Comparators already listed. |

## Verdict

- **The core statistical idea is not new.** Using within-group (sibling) differences to remove group-level confounding and combining them with between-group data for out-of-sample prediction is established (family GWAS unified estimator, Mundlak, CoRe). Our dev4 gain is an instance of that idea.
- **What may still be new (unverified):** (1) applying within-origin supervision from independently derived sister tumour cultures to drug-response prediction for new patients; (2) the replication-gated fallback that turns the contrast off when it does not repeat across disjoint origins. No abstract found does either. Full texts of Nam 2023 and Salatino 2026 must be read to confirm.
- **Consequence for the Unique box:** a uniqueness claim cannot rest on the estimator. At most it is an application-plus-safeguard claim. It must beat a Mundlak-style comparator (ridge with origin-mean expression added) and a snipar-style unified within/between estimator on the same data, not only ridge and shuffled links.

## Required changes before any lock

1. Add comparators: Mundlak ridge (expression + origin-mean expression), and a fixed-weight within/between estimator analogous to the unified estimator. Re-run development on fresh seeds.
2. Read full texts of Nam 2023 and Salatino 2026 for method overlap and data structure.
3. Add Nam 2023 (52 sectors, 18 tumours) and Salatino 2026 as data leads: request counts of sectors with untreated RNA plus TMZ viability per tumour.
