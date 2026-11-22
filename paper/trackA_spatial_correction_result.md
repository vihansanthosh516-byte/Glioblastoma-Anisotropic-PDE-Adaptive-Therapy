# Track A spatial correction result

**Simulation. The correction failed the frozen robustness rule.** Under independent noise and correct rates, it scored 1.000 on signal and 0.513 on the null. Underestimated rates moved the independent null to 0.575. With correlated errors even correct marginal rates gave null concordance of 0.295. Depth matching gave 0.5145 on that same null. [Saved results](../output/trackA_spatial_correction/result.json) | [Frozen plan](trackA_spatial_correction_plan.md)

The product correction assumes independent errors at two endpoints. Marginal error rates do not determine the joint distribution, so independent-control data covering paired errors would be needed to support a joint correction. Depth matching remains a useful control, with the documented coverage loss. The synthetic results leave patient performance and the disruption-budget method untested.
