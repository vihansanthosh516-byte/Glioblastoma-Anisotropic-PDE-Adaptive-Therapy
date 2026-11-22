# Track A bounded block pooling simulation

Use the generator and sizes from the [order reversal plan](trackA_order_simulation_plan.md), with seed 20261123 and 100 replicates per regime. Train on clean localized reversal and test independently under clean localized reversal, localized reversal with test-only corruption, diffuse reversal with test-only corruption, and null with test-only corruption. Use the existing [outlier plan](trackA_outlier_plan.md) corruption mechanism. All methods use the same draws. No cohort data.

Fit a no-intercept L2 logistic model (C=1) on each training patient's difference between mean regional log expression blocks and its negation, with opposite labels. Use the [base plan](trackA_block_simulation_plan_2026-11-17.md) normalization. Apply the coefficients to individual blocks. For each patient, center all block scores on the median across their paired regions. Divide by the RMS of centered training block scores, with a scale of one if zero. Centering is label blind at prediction and requires both regions.

Candidate: clip scaled block scores to [-2, 2], then pool each region separately via log-mean-exp at temperature one. Compare mean, median, maximum and unclipped log-mean-exp of the same scores and separately fitted median-expression and pooled-count paired logistic regression. Higher edge score earns one, a tie half, and lower zero.

Pass if candidate exceeds every comparator by at least 0.01 under contaminated localized reversal, is within 0.02 of the strongest comparator under both other signal regimes, and is within 0.02 of chance under null. The result tests symmetric single-block corruption and does not cover correlated regional artifacts. Novelty unresolved. See the [recovery plan](trackA_recovery_plan.md).
