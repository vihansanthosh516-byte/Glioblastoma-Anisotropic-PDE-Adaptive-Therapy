# Track A adaptive block score simulation

Simulation follow-up to the [fixed ensemble](trackA_order_simulation_result.md). Use its four regimes, depth schedule, replicate counts and train/test sizes with fresh seed 20261120. No patient files.

Fit block logistic and soft order voting as before. Scale each paired margin by its training RMS, with zero scale replaced by one. Split the ordered training patients into four contiguous equal folds. For each fold, fit both methods and scales on its complement and compute scaled margins on its held-out patients. Among order weights 0, 0.125, 0.25 and 0.5, maximize pooled out-of-fold concordance of (1-w) times logistic plus w times order, minus 0.04 times w. Break ties toward zero. Refit components and scales on all training patients, then apply the selected weight to the untouched synthetic patients.

Compare with all previous arms including the fixed ensemble on identical draws. Pass if adaptive beats fixed by at least 0.05 under diffuse shared shifts, trails the best prior arm by at most 0.02 in each signal regime, and scores within 0.02 of 0.50 under null. This is a robustness screen, with no test of superiority to logistic alone. Cross-validated combination follows familiar stacking practice. No novelty claim.
