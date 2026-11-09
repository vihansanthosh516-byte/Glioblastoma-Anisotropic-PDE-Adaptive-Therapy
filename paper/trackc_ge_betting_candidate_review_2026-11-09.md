# Track C GE-betting inference candidate review (2026-11-09)

## Decision

GE-betting is a credible fixed-horizon confidence-interval candidate for later outcome-blind evaluation. A standalone implementation now exists at `src/trackc_ge_betting.py`; it is **not wired into the runner, not adopted, and not independently reviewed**. It does not make the current seven-world-per-family design capable of guaranteeing the 14-day precision target. No policy outcomes were generated for this review.

## What the method offers

Martinez-Taboada and Ramdas describe a deterministic interval for a bounded mean by inverting two fixed-horizon betting processes. The paper states finite-sample validity for bounded observations with a common conditional mean, interval-valued inversion, and Gaussian benchmark width asymptotically under iid sampling. This matches a fixed number of independent bounded world-level differences more directly than a sequential confidence sequence. The paper is an arXiv preprint, version 2 dated 2026-10-30; the author's page labels it under review.

For this project, each contrast could be transformed by `X = (D + 790) / 1580`, where the current utility-derived support is `D ∈ [-790, 790]` simulated days. A two-sided `delta = 0.025` interval for each of the two primary contrast means would give at least 95% familywise coverage by Bonferroni, if each interval's sampling assumptions hold. The world schedule must be materialized and shown to generate independent worlds (or otherwise satisfy the common conditional-mean condition); distinct integer seeds alone are not proof of this.

## What it cannot fix

The previously derived all-zero-sample lower bound still applies to deterministic, uniformly valid intervals over the declared bounded class. With seven worlds per family and `delta = 0.025`, the half-width must reach at least 323.60 days on the specified all-zero data sequence. At least 207 worlds per family are needed merely for this particular lower bound to fall below 14 days; that count does not guarantee 14-day precision for the actual design. GE-betting may be narrower on favorable realized samples, but that is not a uniform precision guarantee and cannot be assumed before seeing results.

The paper also presents optional external randomization. It is not selected here: it changes the interval procedure and its realized-width behavior, and would require its own predeclared interpretation and independent review. It must not be used as a way to bypass the sample-size bound.

## Prototype checks

The independent implementation passed **10 tests** in a Git-archive scratch copy with Python 3.14. The checks include:

- Exact coverage enumeration for every seven-observation Bernoulli sequence over 41 probabilities from 0 to 1.
- Exact coverage enumeration for two rare-event distributions at seven worlds, with extreme values at either support endpoint and event probability 0.409.
- Connected accepted-mean sets on a 201-point grid for every four-observation binary sequence.
- The project-scale all-zero sample check, which verifies that the deterministic interval includes at least `[-323.5, 323.5]` days.
- Input checks, predictable variance calculations, and mapping to the `[-790, 790]` day support.

These tests check implementation behavior on specified cases. They do not replace independent statistical review, validate the full project world law, establish the 14-day precision target, or compare interval widths under policy outcomes.

## If this candidate is pursued

1. Keep the deterministic fixed-horizon version, `delta = 0.025` per contrast, and the existing population-mean estimand as the candidate specification.
2. The independent prototype is in place; the authors' linked GitHub repository states that it has no open-source license, so its code was not copied.
3. Basic endpoint inversion, geometry, boundary, and exact-coverage checks have passed. Next, obtain an independent statistical review and test coverage under the exact world-level sampling law and paired-stream coupling.
4. Compare width with HBCI and Hoeffding under the frozen law. Keep the 14-day precision gate and report failure unchanged.
5. Do not wire this candidate into the evaluation runner or generate policy evaluation outcomes unless the independent review and seed/runtime/precision gates pass.

## Sources

- Martinez-Taboada, D., & Ramdas, A. (2026). *Gaussian-efficient testing by betting on the mean of bounded data* (arXiv:2608.21694v2). [Paper](https://arxiv.org/abs/2608.21694).
- Martinez-Taboada, D. [Research page](https://dmartinezt.github.io/) (lists the paper as under review and links its code).
- Authors' [GE-betting repository](https://github.com/DMartinezT/betting) (README says no open-source license has yet been selected).
- Existing project calculation: [Track C population-mean feasibility audit](trackc_population_mean_feasibility_audit_2026-11-08.md).
