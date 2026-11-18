# Track A adaptive weighting result

**Simulation. The adaptive ensemble passed the frozen robustness screen.** Under shared shifts it scored 0.976 against 0.790 for fixed averaging and 0.9775 for block logistic. It stayed within the frozen tolerance of the best comparator in both reversal regimes and passed the null check. Its mean was slightly below block logistic in every signal regime. Under shared shifts it selected zero order weight in 88 of 100 replicates. [Saved results](../output/trackA_adaptive_simulation/result.json) | [Frozen plan](trackA_adaptive_plan.md)

The practical finding is that learning the weight avoids much of the harm from fixed averaging, while logistic alone remains the simpler baseline. This synthetic experiment used known signal regimes and fixed gene pairs. A fresh test of added value should vary outlier blocks or assay shifts, state what mechanism should favor order votes, and require improvement over block logistic. Keep spent patient cohorts closed.
