# Track C integrated world-sampling contract (development draft)

**Status:** proposed implementation contract; not an execution seal. No policy outcomes were used to choose these rules. No policy evaluation is authorized by this draft.

## Purpose

The integrated comparison estimates an average over independently sampled synthetic worlds. Its sampler must therefore generate independent worlds from the declared world distribution. A batch-balanced design would need a different estimand and an interval method that accounts for its blocks.

## Proposed sampling law

1. **Worlds:** the target law draws each world independently. The continuous world settings use one scrambled Latin-hypercube point generated from that world's registered world-parameter seed. The discrete switching rates `mu` and `nu` are drawn independently and uniformly from `{0, 0.01, 0.04}` per world. They are not balanced across the requested batch. The current sampler accepts per-world seeds; a unique seed check is not proof that a future seed table provides independent streams.
2. **Signal quality:** `q_signal` is drawn once per world, uniformly from `{0.6, 0.7, 0.8, 0.9}`. It is held fixed for that world's virtual patients and is a synthetic assumption, not measured assay accuracy. In the mixture-response family, the responder fraction is named `pi_response`; it is distinct from `q_signal`.
3. **Patients:** selection, calibration, and evaluation cohorts are drawn independently within each world, each from its own registered cohort seed. The Latin-hypercube dimension order is `(r, m, f0, B0, c, epsilon)`; `epsilon` maps to the simulator field `eps`. Resistance cost `c` is patient-level uniform on `[0, 0.25]`.
4. **Observation-noise levels:** after all three cohorts and their patient IDs are created, assign `{0.05, 0.10, 0.20}` in balanced counts across the sorted IDs for that world. Use one deterministic child stream of the world-parameter seed for this assignment. The assignment is shared across the cohorts and does not make the worlds dependent.
5. **Signals:** generate the molecular-like and MRI-like signals jointly from one signal stream per calibration or evaluation cohort. The shared latent draw preserves the dependence model. Use the world's `q_signal` and the predeclared `(a, b, d)` scenario settings. No signal is a patient drug-response measurement.
6. **Random generators:** use explicit `PCG64DXSM` generators for world, patient-cohort, sigma-assignment, and joint-signal generation. The existing policy simulator retains its separately documented `PCG64` trajectory generator until that component is migrated; the same evaluation seed is passed to each paired policy arm for common-random-number coupling. The exact schedule roles and seed registry checks remain required before any run.

## Consequences and limits

- The independent-world law matches the equal-world population-mean estimand. It does not guarantee balanced representation of the nine `(mu, nu)` combinations in a finite sample.
- The seed table is not materialized or reserved. Before claiming independent-world inference, review the actual seed-generation schedule; registry uniqueness alone does not establish independence.
- Patient-level balancing of sigma is only a variance-control design feature. It does not increase the independent-world sample size.
- `q_signal`, response-family parameters, and the tumor dynamics remain synthetic assumptions. This sampler cannot establish a molecular assay's accuracy or patient-specific resistance.
- The current simultaneous interval candidates failed outcome-blind coverage checks, and the distribution-free precision target remains infeasible at the available seed capacity. No policy outcomes may be generated until an alternative inference method passes independent review and outcome-blind validation, and the precision/runtime gates are resolved.

## Implementation boundary

The implementation may be tested on synthetic inputs only. It must not launch controller comparisons, reserve evaluation seeds, read patient MRI or masks, or write a sealed result. Existing Track C results and legacy calibration functions remain unchanged.
