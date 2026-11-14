# Track A gene-pair scoring simulation protocol

This simulation tests whether the proposed feature-selection rule behaves as designed when a shared regional effect competes with a stronger state-specific effect. The generator deliberately encodes that premise. The result can check the implementation and illustrate its trade-off. It cannot establish clinical performance, uniqueness, or superiority to published regional signatures.

Before execution, create the generator and freeze it with this file. Use NumPy with seed 7461. Simulate 96 development donors and 256 independent simulated test donors, with four states and 24 genes per specimen. Draw standard normal donor offsets shared across regions and independent standard normal observation noise for each specimen, with noise SD 0.55. Use independent draws for development and test. Apply one of three fixed margin effects, in this order:

- Shared plus state-specific: increase the first gene by 0.6 and decrease the second by 0.6 in all states. Increase the third gene by 3 and decrease the fourth by 3 in the first state.
- State-specific only: apply the same changes to the third and fourth genes in the first state. Leave the other genes unchanged.
- Null: assign no margin shifts.

Fit every arm on the same simulated development donors. Compare the proposed minimum-state selector, its mean-state selection ablation, pooled-state TSP, and a training-only paired nearest-centroid direction. Use one pair for each rank arm to isolate selection. Pooled TSP uses the highest absolute mean paired vote change on expression averaged across states. The centroid direction uses the vector of mean paired expression differences after averaging across states. Use stable gene order for ties. Evaluate every arm on the same simulated test donors within a scenario. Score each pair by margin-over-core concordance, with one-half credit for ties. Report the four arm scores in each scenario and the selected pair for each rank arm.

Pass the mechanics screen only if the candidate scores strictly above the mean-state and pooled-TSP arms in the shared-effect scenario, and if every null score lies between 0.40 and 0.60 inclusive. The nearest-centroid score is descriptive. Report the state-specific-only scenario as a trade-off, including any losses. These values are prospective simulation settings, with no empirical calibration. No retest with a different seed after seeing the result.
