# Track C finite-sample signal-recovery protocol

**Status: protocol frozen; synthetic development diagnostic complete; no policy or patient evaluation.**

## Question

When `q`, `a`, `b`, and `d` are all unknown, how often will a finite sample of paired molecular-like and MRI-like signals recover the generating parameters under the project's assumed `f0` distributions? This tests parameter recovery in the specified simulator only. It does not test real biology, treatment benefit, or clinical use.

## Data-generating conditions

- **Truth anchors:** `(q,a,b,d)` = `(0.7,0.5,0.5,-0.5)`, `(0.7,0.5,0.5,0)`, `(0.7,0.5,0.5,0.5)`, `(0.9,0.5,0.5,0.5)`, `(0.5,0.5,0.5,0.5)`, `(0.7,0,0.5,0.5)`, and `(0.7,0.5,0,0.5)`. All anchors lie on the already evaluated parameter grid. The `q=0.5` case is a negative control: a coin-flip molecular report should carry no information about `a` or `d`. The last two anchors vary the truth values of nuisance loadings `a` and `b`.
- **`f0` distributions:** the three declared calibration ranges: low shift, default, and high shift. The reference-normal case is omitted here because its exact `q`/combined-correlation alias is already demonstrated by the analytic control and Fisher grid.
- **Sample sizes:** 18, 28, 100, and 500 virtual patients per cohort. The first two cover the scale of the available cohorts; the larger sizes show how the assumed likelihood behaves as information increases.
- **Replicates:** 25 per truth, `f0` range, and sample size (2,100 synthetic cohorts total), using reserved Stage 2 development root seeds `20261880–20261904`. Each root seed is split deterministically by range, truth, and sample-size indices. The run uses no separate noise-key stream.
- **Signal generation:** draw unobserved `f0` values log-uniformly from the assigned range. Generate the paired observations with `generate_integrated_signals`. The estimator will not receive the latent `f0` values or any latent truth.

## Estimator and reporting

- For each synthetic cohort, evaluate the exact same 81-point grid already used in the Fisher stress test: `q ∈ {0.5, 0.7, 0.9}` and `a,b,d ∈ {-0.5, 0, 0.5}`.
- Score each point by the sum of marginal paired-signal log likelihoods, integrating the unobserved `f0` over its known generating range with 256 Gauss-Legendre nodes. This is a discrete-grid maximum-likelihood diagnostic; there will be no post-run optimizer or grid adjustment.
- If several grid points tie within `1e-8` log-likelihood units, retain and report the full tie set. The point estimate for `d` will be the mean `d` across that tie set.
- Also report exact truth-in-MLE-tie rate, absolute `d` error, sign recovery for nonzero `d`, and the range of `d` values within 2 log-likelihood units of the maximum. This near-maximum range is descriptive and is not called a confidence interval.
- Report every truth/range/sample-size group, including failures and negative-control results. No best-performing group will be selected as the project's success result.

## Verification and limits

Tests will check deterministic substreams, the 81-point candidate set, finite likelihoods, and the `q=0.5` invariance to `a` and `d`. The complete run will use only synthetic observations and registered development seeds. It will not inspect Track A, Duke, or any patient image, expression, mask, or outcome. A positive simulation result would establish only recovery under these assumed distributions and this estimator grid; it would not show that the signal parameters can be learned from the project's separate cohorts or that the treatment policy is better.

## Post-run reporting correction (2026-11-06)

The first completed output computed a `d`-sign score for the `q=0.5` negative-control anchor. That score is not meaningful: under this model, a coin-flip molecular report makes `a` and `d` unidentifiable from the observed pair. The final summary therefore sets that anchor's `d_sign_correct` field to null and excludes it from the sign-recovery denominator. This correction changes no generated observations, likelihoods, candidate fits, or other metrics; the completed run is not repeated and the original seed draws remain spent.

The corrected record contains 2,100 rows, 84 groups, the full 81-candidate grid, and the `postprocessing_correction` note. The 25 used roots are recorded as `dev_used` in `src/trackc_seed_registry.py`; the registry check passed. The focused recovery and registry suites passed 16 tests, with one existing evidence test skipped because the scratch copy did not contain `output/trackc_sealed_containment/frozen_protocol.json`.
