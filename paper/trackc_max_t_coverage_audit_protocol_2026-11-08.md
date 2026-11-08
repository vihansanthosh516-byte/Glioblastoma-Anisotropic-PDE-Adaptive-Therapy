# Track C observed-SE max-T coverage audit

**Status:** preregistered synthetic method check; not a policy evaluation, not a frozen analysis, and not evidence of therapy benefit.

## Question

Does the proposed world-resampled observed-SE max-T interval procedure attain acceptable simultaneous coverage for two independent response-family means under several known distributions and small world counts?

## Fixed procedure

- The target is joint coverage of two nominal 95% intervals for the two response-family means.
- Each observation represents one independent world-level paired candidate-minus-control patient-macro difference. Families are generated independently.
- World counts per family: 5, 10, 14, and 20.
- Distribution pairs: Normal/Normal; centered Gamma(shape 2)/centered Gamma(shape 2); standardized Student t(df 3)/standardized Student t(df 3); bounded Uniform/Uniform; Normal/centered Gamma(shape 2); centered Gamma(shape 2)/standardized Student t(df 3). Each marginal has mean zero and variance one.
- Monte Carlo repetitions: 2,000 per distribution pair and world count.
- Bootstrap resamples per repetition: 2,000, using the production helper and its fixed family streams.
- Master seed: 20261108. Data-generation streams use PCG64DXSM and SeedSequence keyed by scenario, world count, and repetition. Bootstrap seeds use SHA-256 of the fixed audit namespace and those identifiers.
- The empirical joint coverage is the fraction of repetitions in which both intervals include their true mean, zero. The audit also records the false-positive rate of the rule that both lower bounds exceed zero.
- A 95% Wilson interval is computed for joint coverage and for the global-null false-positive rate.

## Decision rule

The proposed method passes this **screen only** if every cell has a 95% Wilson lower bound for joint coverage of at least 0.93, and no cell has an inference failure. Otherwise, record failure or incomplete status. Do not drop a distribution, world count, family, or failed repetition. Do not alter the interval method using these results; any replacement needs a new dated protocol and fresh synthetic validation.

This screen allows up to a 2-percentage-point shortfall below nominal coverage at the Wilson lower bound. Passing does not prove nominal coverage outside these distributions or validate the simulator. The 2,000 bootstrap resamples also differ from the proposed 10,000-resample final setting, so passing cannot freeze the method.

## Reproducibility and limits

The output will record all settings, per-cell counts, Wilson intervals, code/config/protocol SHA-256 fingerprints, and runtime. Inference failures count as noncoverage and keep the cell in the denominator; any failure makes the overall audit incomplete. It will contain no patient data and use no Track C policy outcomes. This evaluates only the interval helper under controlled distributions. It does not validate world independence in a future study, and hash checks cannot prove independence.

An adjacent diagnostic-accuracy simulation found small-sample maxT error inflation, but used a different endpoint and procedure; it motivates this check without predicting its result ([Westphal and Zapf, 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11025299/)).
