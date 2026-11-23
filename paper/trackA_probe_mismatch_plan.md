# Track A probe alignment stress plan

Fresh simulation using the [tradeoff design](trackA_probe_tradeoff_plan.md) and seed 20261204. Add a mismatched-probe regime with the same patient and offset-control distributions as blind mixed but positive-control sensitivity one. This models a challenge that passes on its own target while the actual assay offset remains undetected. The earlier qualified classifier is unchanged.

Compare against a deployment rule that additionally requires independently verified probe alignment in both conditions. In this simulator only the mismatched-probe regime lacks alignment evidence. This metadata flag is supplied by the generator and represents an idealized independent audit, not a learned detector. When evidence is missing, use composition plus expression; otherwise use the earlier decision. Hold all previous tradeoff checks fixed for each model and additionally require concordance no more than 0.02 below omission under mismatch. Report all arms and regimes. No patient data.

In the paper, lead with the final candidate and its evaluation. Keep the development record in an appendix or linked repository. Disclose prior outcomes that motivated changes to methods or success criteria, and report all frozen checks on the final candidate. Novelty and generalization remain open.
