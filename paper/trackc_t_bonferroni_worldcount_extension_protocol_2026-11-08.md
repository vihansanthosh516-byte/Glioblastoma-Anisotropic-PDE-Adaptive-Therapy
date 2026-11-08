# Track C Student-t Bonferroni world-count extension

**Status:** proposed outcome-blind feasibility screen only. This does not adopt the method, authorize policy evaluation, or imply therapy benefit.

## Purpose

The first Student-t Bonferroni screen failed at 5, 10, 14, and 20 worlds per response family. This extension asks whether the same fixed interval procedure reaches the same coverage screen at larger world counts. It does not change the method, confidence level, screen threshold, or policy estimand. It will not select a world count for the final study.

## Fixed interval procedure

Use two-sided Student-t intervals with Bonferroni adjustment across the two required response-family means. For each family, calculate its world mean and sample standard error, then use `t_quantile(0.9875, df=n-1)` for a familywise 95% interval. Exact marginal coverage requires independent, normally distributed world differences within each family. Bonferroni does not require dependence between families. The actual simulator's world-difference distribution is not known before policy evaluation, so this remains only a generic method screen.

## Fixed synthetic screen

- World counts per family: 30, 50, 100, and 200.
- Monte Carlo repetitions per cell: 20,000.
- Master seed: `2026110811`; use PCG64DXSM with `SeedSequence([master_seed, pair_index, world_count, family_index])`.
- Distribution pairs, with each marginal standardized to mean zero and variance one:
  1. Normal / Normal.
  2. Centered Gamma(shape 2) / centered Gamma(shape 2).
  3. Standardized Student t(df 3) / standardized Student t(df 3).
  4. Uniform / Uniform.
  5. Normal / centered Gamma(shape 2).
  6. Centered Gamma(shape 2) / standardized Student t(df 3).
  7. Standardized Beta(2, 5) / standardized Beta(2, 5).
  8. Standardized Beta(2, 5) / standardized Beta(5, 2).
  9. Standardized Bernoulli(p = 0.1) / standardized Bernoulli(p = 0.1).
  10. Centered Gamma(shape 2) / standardized Bernoulli(p = 0.1).
  11. Standardized Student t(df 3) / standardized Bernoulli(p = 0.1).
- Joint coverage means both intervals contain their true mean, zero. Also record the false-positive rate when both lower bounds exceed zero.
- Use 95% Wilson intervals for joint coverage and false-positive rate.
- Pass this screen only if every cell's joint-coverage Wilson lower bound is at least 0.93 and inference failures equal zero. Keep all cells and failures.

No policy outcome or patient data may enter this screen. Passing only permits independent review and a separate validation using a justified model for the actual world-level data structure. It does not freeze the interval method or establish that the t assumptions hold for the simulator.

## Compute reference

The timing-only pilot at two patients per selection, calibration, and evaluation cohort took 1,373.3 seconds for one single-response world and 1,482.7 seconds for one mixture-response world. One hundred worlds per family would therefore require about 79.3 wall-clock hours at those tiny cohort sizes, before any increase in patient count. This is a rough extrapolation from one timing world per family, not a bound or a full-study estimate. See `output/trackc_robust_policy_timing_pilot_20261106_followup/manifest.json`.

## Stop rule and interpretation

Report coverage by cell. Do not choose the smallest passing count as the final design from this screen. If any cell fails, the candidate remains ineligible. Even if all cells pass, the next gates remain: independent statistical review, validation against the actual world-generation and utility process without using evaluation outcomes, a precision calculation, seed capacity, and runtime feasibility. Keep research-scale policy comparisons blocked until those gates pass.
