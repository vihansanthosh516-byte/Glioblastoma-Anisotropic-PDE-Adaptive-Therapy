# Track A outlier block simulation result

**Simulation. The combined model failed the frozen superiority rule.** With contaminated blocks, combined median expression and order votes scored 0.9765 versus 0.9735 for median expression alone, a gain of 0.003 against the required 0.01. Mean block logistic scored 0.94475. Clean and null checks passed. In this generator, the larger improvement came from summarizing blocks by the median. [Results](../output/trackA_outlier_simulation/result.json) | [Frozen rule](trackA_outlier_plan.md)

This tests diffuse reversal with the same synthetic corruption mechanism in training and test. Whether the median drops biologically informative localized signal remains open. A fresh experiment can compare median, trimmed mean and a prespecified robust pooling rule across localized signal and train-to-test changes in corruption. No patient data were used.
