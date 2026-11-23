# Track A positive-control qualification simulation

Use all regimes, cohort sizes, replicate counts and baseline arms of the [crossover simulation](trackA_crossover_plan.md) on fresh seed 20261201. Add half-sensitive mixed: offset is as in mixed but offset control means are half the true offsets. Add 40 independent unit-variance positive-control observations in each condition with mean 1 for full sensitivity, 0.5 for half sensitivity, and 0 for blind mixed. These represent independently calibrated standardized assay challenges, separate from the offset estimates and patient data.

Use corrected spatial prediction only if the one-sided normal lower bound mean minus 1.645/sqrt(40) exceeds 0.60 in both conditions. Otherwise use the fitted composition-plus-expression model. Both fixed crossover arms and all prior comparators use the same patients. Record the fraction using spatial.

Pass if the candidate meets the prior primary gates, loses at most 0.02 versus omission in both incomplete-coverage regimes, and qualifies no more than 0.10 of replicates in each incomplete-coverage regime. This is safe fallback, not rescue of spatial signal when a control is blind. Novelty unresolved. The design assumes the standardized challenge tests the actual nuisance direction. No real patient data.
