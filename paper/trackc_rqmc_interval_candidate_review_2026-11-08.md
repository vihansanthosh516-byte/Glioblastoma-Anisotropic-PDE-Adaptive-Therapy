# Track C RQMC interval candidate review

**Status:** read-only methods review; no interval method adopted and no policy outcomes or patient data used.

## Conclusion

Randomized quasi-Monte Carlo (RQMC) is a research lead for reducing integration error, not a ready replacement for the Track C uncertainty method. It does not repair the failed Track C Student-t coverage screen by citation alone. Keep population-level policy evaluation blocked until a new method passes an outcome-blind coverage and runtime screen matched to the actual world-level estimand.

## Evidence and fit

- L’Ecuyer, Nakayama, Owen, and Tuffin describe confidence intervals across independent RQMC replications. They report that Student-t intervals performed more reliably than the bootstrap-t alternatives on the integration problems they studied. They also explain that observations inside one randomized point set are dependent, so an ordinary interval over those points is not justified. This is encouraging for some integration problems, but is not a coverage guarantee for this project’s policy contrasts. [GERAD report](https://www.gerad.ca/en/papers/G-2023-67) · [full report PDF](https://www.gerad.ca/fr/papers/G-2023-67.pdf)
- Nakayama and Tuffin give sufficient conditions for a central limit theorem for RQMC error. Their conditions depend on regularity of the integrand, and some settings require the number of independent randomizations to grow sufficiently quickly. A simulator outcome involving thresholded treatment decisions, bounded progression time, censoring, and a maximum across models may not satisfy the needed smoothness. This must be checked, not assumed. [Winter Simulation Conference paper](https://doi.org/10.1109/WSC52266.2021.9715427)
- Jain et al. study empirical-Bernstein and betting confidence intervals for bounded RQMC integrands. They provide a finite-sample route under explicit boundedness assumptions, but their interval width still depends on the replication design and can remain large. The Track C project has already found that the standard empirical-Bernstein range correction is far above its 14-day target at current world counts. A betting interval has not been validated here. [Primary preprint](https://arxiv.org/abs/2504.18677)

## Project-specific limits

The integrated handoff specifies independently seeded, scrambled Latin-hypercube parameter draws within worlds. Track C’s inference unit is the independent world-level paired policy difference, not individual simulated patients or points within a Latin-hypercube draw. RQMC literature therefore does not justify treating within-world patients as independent replicates, nor does it automatically supply valid simultaneous intervals for both response-family means.

The project’s fixed coverage screens already rejected observed-SE max-T and Student-t with Bonferroni adjustment under tested skewed and Bernoulli world-level distributions. The Student-t/RQMC literature does not override those project-specific failures: it studies different integrands and replication structures, and its favorable empirical findings are not universal finite-sample guarantees.

## Decision

Do not add RQMC or a betting interval to the frozen protocol yet. If population-mean inference remains the goal, a new proposal must specify the randomization unit, independent replicate structure, familywise procedure, assumptions, coverage criterion, and runtime cap before policy outcomes are generated. Validate coverage on outcome-blind distributions that include the observed structural features (discreteness, skew, censoring, and paired world contrasts). If it fails either coverage or precision, do not relax the criterion after seeing policy results.

No code, settings, seeds, endpoints, or results changed. No tests ran. This review is not an authorization to run a policy comparison.
