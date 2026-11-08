# Track C integrated signal-likelihood audit

**Status:** synthetic/code audit only. No patient data, policy outcomes, or treatment-effect claims were used.

## Finding

For the frozen interior parameter grid, the integrated signal likelihood follows the equations in [`HANDOFF_AND_REPRODUCIBILITY.md`](../tracks/shared/HANDOFF_AND_REPRODUCIBILITY.md): it adds the MRI Gaussian density and the molecular-report probability conditional on the observed MRI score. Existing tests compare the joint likelihood with direct numerical integration and check normalization, factorization when residual dependence is zero, dependence behavior, and the uninformative `q = 0.5` case.

The audit found one edge mismatch in the helper's declared coefficient domain. `SignalModelParameters` permits `a = ±1`. At zero conditional variance and mean exactly zero, the signal generator defines the molecular state as positive (`Z_bio >= 0`), but `_normal_positive_probability` returned 0.5. The probability should be 1 for this deterministic boundary under the generator's rule. This edge is outside the current config's `a,b,d ∈ {-0.5, 0, 0.5}` grid, so it does not change the frozen protocol's current parameter cases.

## Change and verification

The deterministic branch now applies `mean >= 0` when variance is exactly zero. Positive variance uses the normal CDF, including very small positive variances. A synthetic regression test was added. It failed before the code fix and passed afterward.

Scratch-copy verification used Python 3.14.4 with `PYTHONPATH` set to the scratch `src` directory. The edited source and test hashes matched the scratch files. The integrated signal-conditioning, policy benchmark, world runner, and robust belief-MPC suites passed: **83 passed**. `git diff --check` passed; Git emitted only expected LF-to-CRLF working-copy warnings.

## Limits

This verifies mathematical mechanics and policy plumbing for the synthetic signal model. It does not validate the assumed molecular/MRI dependence parameters, establish resistance measurement, evaluate clinical benefit, or show that Track A and Track B patient data are linked. No policy outcome comparison ran.
