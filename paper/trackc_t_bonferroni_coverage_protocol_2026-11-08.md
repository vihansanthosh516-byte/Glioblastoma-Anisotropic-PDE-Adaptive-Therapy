# Track C Student-t Bonferroni interval candidate: outcome-blind screen

**Status:** proposed methods check only. This protocol does not adopt an interval method, authorize policy-outcome evaluation, or imply therapy benefit.

## Rationale and estimand

The integrated v2 primary estimand is the equal-weight mean of world-level paired candidate-minus-control patient-macro utility differences, reported separately for the single-response and mixture-response families. The earlier world-resampled observed-SE max-T proposal failed its predeclared coverage screen and remains rejected.

Track C generates each world by calling `draw_worlds(1)` with a unique world seed. Patient cohorts within each world use Latin-hypercube sampling. This screen therefore treats independently seeded worlds as the sampling units; it does not treat patients within one world as independent inferential units.

## Candidate procedure

This screen evaluates one fixed candidate: a two-sided Student-t interval for each response-family mean, with Bonferroni adjustment across the two required families.

For family `j` with `n_j` world differences `d_ij`, calculate:

1. `mean_j = mean(d_ij)`.
2. `SE_j = sample_sd(d_ij, ddof=1) / sqrt(n_j)`.
3. `critical_j = t_quantile(1 - 0.05 / (2 * 2), df=n_j - 1)`.
4. Interval `mean_j ± critical_j * SE_j`.

The method requires at least two finite world differences and positive sample variance in each family. A failure in either family makes the joint inference incomplete. The two required family names are exactly `single-response` and `mixture-response`. Bonferroni controls simultaneous coverage if each family-specific interval has its stated marginal coverage; the t interval has exact marginal coverage only when world differences are independent and normally distributed. For other distributions, coverage is approximate and must be checked. No normality claim is made for the Track C outcome distribution.

This candidate is chosen from standard small-sample simulation-analysis practice before the new synthetic screen. The screen will not select among methods, tune a critical value, or use Track C policy outcomes.

## Fixed outcome-blind screen

- Distribution pairs: Normal/Normal; centered Gamma(shape 2)/centered Gamma(shape 2); standardized Student t(df 3)/standardized Student t(df 3); bounded Uniform/Uniform; Normal/centered Gamma(shape 2); centered Gamma(shape 2)/standardized Student t(df 3).
- Each marginal has mean zero and variance one.
- World counts per family: 5, 10, 14, and 20.
- Monte Carlo repetitions per cell: 20,000.
- Master seed: `2026110809`. Use PCG64DXSM streams derived from `[master_seed, pair_index, world_count, family_index]`.
- Joint coverage: both family intervals contain their true mean, zero.
- Also record false positives for the rule that both lower bounds exceed zero.
- Report 95% Wilson intervals for coverage and false-positive rate.
- Pass this screen only if every cell's Wilson lower bound for joint coverage is at least 0.93 and there are zero inference failures. Keep every cell and failure in the denominator.

The 0.93 lower-bound threshold is retained from the earlier screen; it is not weakened. Passing would qualify this candidate for a separate review and a new validation matched to the actual simulator-generated world-difference distribution. Passing this small family of known distributions would not freeze the method or authorize policy evaluation. Failing rejects this candidate for the current study.

## Distribution-free precision reference

Primary utility is bounded by `[-60, 730]` simulated days: progression time is capped at 730 days and the primary course penalty is 5 days for at most 12 courses. A paired candidate-control world difference is therefore in `[-790, 790]`, with range width 1,580 days.

For two two-sided family intervals at 95% simultaneous confidence, Hoeffding's inequality gives a half-width of `1,580 * sqrt(log(80) / (2*n))` under independent worlds. At 14 worlds this is about 625 days. A worst-case half-width of 14 days would require about 27,900 worlds per family. This is a conservative reference, not a proposed run size; it shows that a distribution-free guarantee cannot meet the current precision target at the present seed capacity.

## Limits and stop rule

This screen tests only interval mechanics on known distributions. It does not validate model biology, treatment benefit, the actual world-difference distribution, seed independence, or clinical use. Do not read patient data or score policy outcomes for this check. Do not change the method or pass rule after seeing coverage. If the candidate fails, keep research-scale policy comparison blocked and write a separate, dated protocol before testing another method.

## Primary sources

- L'Ecuyer, P., Nakayama, M. K., Owen, A. B., & Tuffin, B. (2023). Confidence intervals for randomized quasi-Monte Carlo estimators. *Proceedings of the 2023 Winter Simulation Conference*, 445–456. https://doi.org/10.1109/WSC60868.2023.10408613. The authors study independent randomizations and report that Student-t intervals were more reliable than bootstrap-t intervals on their RQMC integration problems. This is motivation, not direct validation for this project.
- Owen, A. B. (1992). A central limit theorem for Latin hypercube sampling. *Journal of the Royal Statistical Society: Series B (Methodological), 54*(2), 541–551. https://doi.org/10.1111/j.2517-6161.1992.tb01895.x.
- Hoeffding, W. (1963). Probability inequalities for sums of bounded random variables. *Journal of the American Statistical Association, 58*(301), 13–30. https://doi.org/10.1080/01621459.1963.10500830.
