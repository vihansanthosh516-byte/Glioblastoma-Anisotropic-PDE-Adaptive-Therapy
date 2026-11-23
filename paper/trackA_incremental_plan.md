# Track A incremental evidence simulation

Simulation only. Test whether a spatial feature adds information beyond composition and malignant-cell expression on fresh synthetic patient pairs. This is a harness check with simulated scalar summaries, not a model of cell geometry or clinical validation. Keep all real test cohorts closed.

Use NumPy default_rng seed 20261127. For each of 100 independent replicates per regime draw independent standard normal features for each patient's paired edge-minus-core contrast. The features represent composition, malignant-cell expression, and spatial arrangement, in that order. Use 80 training and 200 test pairs. Add shifts (0.35,0.35,0.80) in the spatial-signal regime, (0.35,0.35,0) in the expression-only regime, and (0,0,0) in the null. Simulate an additional artifact regime with (0,0,0.80) in training and (0,0,-0.80) in testing, representing acquisition associated with region in opposite directions.

Compare paired L2 logistic regression with C=1 and no intercept on nested feature sets: first one, first two, or all three columns. For training, concatenate paired differences and their negations with labels one and zero respectively. Use identical patients for each comparison, and score the proportion of positive held-out margins with ties worth half.

Pass if the full model improves by at least 0.05 over both simpler arms in the spatial-signal regime, loses no more than 0.02 to either under expression only, and is within 0.03 of chance in both null and artifact regimes. Report each condition even if overall failure. The artifact stress tests whether feature inclusion alone protects against study shift.
