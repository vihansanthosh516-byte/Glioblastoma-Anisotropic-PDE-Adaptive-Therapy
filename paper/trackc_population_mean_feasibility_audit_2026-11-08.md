# Track C population-mean feasibility audit (2026-11-08)

**Status: outcome-free methods audit.** This review used source code, protocol files, saved synthetic coverage records, and timing metadata. It used no patient data or Track C policy outcomes. It reserved no seeds and changed no analysis settings.

## Decision

Keep the population-mean estimand and the 14-day precision target unchanged. Keep policy evaluation blocked. The current design has no validated simultaneous interval, its generator does not match the integrated v2 sampling specification, and current seed capacity cannot guarantee the 14-day target over the full bounded distribution class without stronger assumptions or tighter proven bounds.

## Findings

### 1. Existing interval candidates remain rejected

- The observed-SE max-T screen passed **0/24** predeclared cells.
- Student-t with Bonferroni passed **14/24** cells in its first screen and **35/44** in its larger-world extension. Coverage was not monotonic in world count; the worst larger-screen cell had 88.48% joint coverage.
- The distribution-free Hoeffding interval is valid for independent bounded world differences, but its half-width is far above 14 days at current seed capacity.

These are prior saved synthetic screens, not policy results: [max-T screen](trackc_max_t_coverage_audit_result_2026-11-08.md), [Student-t screen](trackc_t_bonferroni_coverage_audit_result_2026-11-08.md), [larger-world extension](trackc_t_bonferroni_worldcount_extension_result_2026-11-08.md), and [distribution-free precision audit](trackc_distribution_free_precision_audit_2026-11-08.md).

### 2. A distribution-free lower bound rules out guaranteed 14-day precision at current capacity

The configured utility is the minimum of two progression times capped at 730 days, minus five days for each course through the cutoff. With at most 12 courses, patient utility lies in `[-60, 730]` days, so a paired candidate-minus-control world difference lies in `[-790, 790]` days.

For two simultaneous two-sided intervals with familywise 95% coverage, each marginal interval must have at least 95% coverage. Under equal Bonferroni allocation, each interval needs 97.5% marginal coverage. Consider independent world differences supported on `[-790, 790]`:

- `P0`: every difference is zero.
- `P+`: a difference is `+790` with probability `p` and zero otherwise.
- `P−`: a difference is `−790` with probability `p` and zero otherwise.

All three distributions can produce the same all-zero sample. Under `P+` and `P−`, that sample occurs with probability `(1−p)^n`. If this probability exceeds 0.025, a deterministic distribution-free interval with the required coverage must include both alternative means, `−790p` and `+790p`, on that sample. Its half-width must therefore be at least `790p`.

Under equal Bonferroni allocation, the all-zero-data lower bound is `790 × [1−0.025^(1/n)]` days:

| Independent worlds per family | Necessary half-width on all-zero data |
| ---: | ---: |
| 7 | 323.60 days |
| 14 | 182.99 days |
| 200 | 14.44 days |

At least **207 independent worlds per family** are needed merely to clear this bound under equal Bonferroni allocation. This is not sufficient to guarantee 14-day precision for the actual design. More generally, the familywise 95% requirement implies at least 95% marginal coverage per family; even using that less stringent marginal level, seven worlds require a half-width of at least 275.05 days on all-zero data. These bounds apply to deterministic, uniformly valid intervals built from independent bounded world differences. A narrower interval would need defensible additional assumptions or a tighter proven support; a method name alone cannot supply either.

This is a worst-case guarantee result. It does not say every interval computed from every possible nonzero sample must exceed 14 days. It says the current design cannot guarantee that precision for every distribution in the stated bounded class without extra assumptions or a tighter proven range.

The current development seed block has 89 free scenario seeds. The proposed schedule uses six scenario seeds per world, giving an upper capacity of 14 total worlds, or seven per family with an equal split. The 207-world threshold projects to about 164.22 serial hours for both families even at the existing two-patient-per-cohort timing reference; larger cohorts and method overhead add cost. The timing source is a small synthetic pilot, not a bound.

### 3. Integrated v2 config and prototype generator disagree

The integrated configuration is not an execution seal, and the development runner is not yet a faithful implementation of it:

- The configuration specifies balanced `(mu, nu)` counts across the requested world set. `prepare_synthetic_world` instead calls `draw_worlds(1)` separately for each world. That samples one pair per call; it does not enforce the stated global balance.
- The configuration puts `c` in the patient axes. `trackc_calibration.draw_worlds` includes `c` among world axes, and `draw_patients` assigns the same `c` value to every patient in that world.
- The configuration specifies a per-world `q_signal` draw from `{0.6, 0.7, 0.8, 0.9}`. The prototype receives fixed `signal_parameters` from its caller and uses that same `q` for signal generation.

Further source inspection found these additional execution mismatches:

- **Patient-axis order:** the configuration orders patient dimensions as `[r, m, f0, B0, c, epsilon]`; `trackc_sim.PARAMS` orders them as `[r, m, f0, c, eps, B0]`. Because the code uses Latin-hypercube coordinates, this changes which coordinate receives each stratified draw, even though marginal bounds remain the same.
- **Observation-noise assignment:** the configuration assigns a balanced `{0.05, 0.10, 0.20}` sigma vector once per world over sorted patient IDs. `trackc_sim.draw_scenarios` instead permutes a resized sigma vector on every call. `prepare_synthetic_world` calls it separately for selection, calibration, and evaluation cohorts, so each cohort gets a separate assignment. With two patients per cohort, resizing assigns only `0.05` and `0.10`; `0.20` is absent.
- **Bit generator and seed derivation:** the configuration names `PCG64DXSM` and SHA-256-derived stream seeds. The prototype uses `np.random.default_rng`, which returned `PCG64` in the recorded Python environment, and accepts preassigned integer seeds. `trackc_world_schedule` checks seed availability and uniqueness; it does not derive or reserve the configured streams.
- **Signal-stream names:** the configuration lists separate molecular and MRI signal streams. `generate_integrated_signals` uses one joint seed for the correlated latent draws and report flip. Separate channel streams are not implemented. A joint signal stream should be named explicitly so future edits do not break the intended dependence structure.
- **Policy and observation streams:** `run_shared_policy_world` passes one evaluation-policy seed to every arm. `simulate_two_endpoints` draws scan noise, then passes the same generator to `policy.reset`; it also derives missing-scan draws from that seed. The configuration instead names observation-noise and candidate, `d=0`, and secondary-policy streams. This is a reproducibility-spec mismatch; it does not by itself show biased scores, but the random-number coupling must be specified before evaluation.
- **Treatment-stop rule:** the configuration specifies the observed two-consecutive-scan rule. `PolicyBenchmarkSettings` still defaults to `stop_endpoint="either"`, and `_evaluate_policy` forwards that value. The integrated world-preparation path now rejects settings other than `observed_two_scan` before generating a world. There is still no config loader, and direct calls to the low-level benchmark can bypass this guard.
- **Penalty sensitivities and bootstrap:** the configuration lists course-penalty sensitivities of 0 and 10 days and a `world_bootstrap` stream. One `run_shared_policy_world` call uses one penalty value (default 5); the development schedule has no bootstrap role. The research driver and sensitivity orchestration remain unimplemented.

The target population, categorical weights, randomization unit, and generator must agree before a confidence method can be validated for the intended estimand. If independent world-level inference is retained, a defensible option is to draw `(mu, nu)` and `q_signal` independently from their declared uniform distributions per world and draw `c` per patient. A balanced design can use complete independently randomized blocks, or another interval that accounts for the actual stratified design. Do not treat patients within one world, or points within one randomized QMC scramble, as independent worlds.

### 4. Literature does not remove these project-specific gates

Independent randomized QMC replications can support interval estimation, but ordinary Student-t intervals are approximate and depend on the integrand and replication regime. Published CLT results state regularity and relative-growth conditions; they are not a guarantee for this simulator's thresholded decisions, censoring, and minimum-of-capped-endpoints utility. Betting intervals for bounded means offer finite-sample guarantees under their sampling assumptions, but boundedness does not remove the rare-event precision limit above. No cited result validates an interval for this project's current generator without a matching outcome-blind check.

Primary sources:

- L’Ecuyer, P., Nakayama, M. K., Owen, A. B., & Tuffin, B. (2023). *Confidence intervals for randomized quasi-Monte Carlo estimators*. [GERAD report](https://www.gerad.ca/en/papers/G-2023-67).
- Nakayama, M. K., & Tuffin, B. (2024). *Sufficient conditions for central limit theorems and confidence intervals for randomized quasi-Monte Carlo methods*. [Primary paper](https://doi.org/10.1145/3643847).
- Waudby-Smith, I., & Ramdas, A. (2024). *Estimating means of bounded random variables by betting*. [JRSS B paper](https://doi.org/10.1093/jrsssb/qkad009).
- Jain, A., Hickernell, F. J., Owen, A. B., & Sorokin, A. G. (2026). *Empirical Bernstein and betting confidence intervals for randomized quasi-Monte Carlo*. [Primary paper](https://doi.org/10.1093/imaiai/iaag003).

## Next work before any policy outcomes

1. Amend the synthetic generator and configuration together. Lock the target population law, `(mu, nu)` weights, `q_signal` distribution, patient-level `c`, patient-axis order, sigma assignment, response-family relationship, and independent sampling unit.
2. Align the implementation contract: choose the bit generator and seed derivation; define joint signal streams, shared observation-noise coupling, per-arm policy streams, and the course-penalty sensitivity runs. Keep the integrated-builder stop-rule guard, then make the eventual driver use the same configured rule and prevent direct-run bypass.
3. Obtain independent statistical review of a fixed-time simultaneous interval under that exact law. Keep 95% familywise coverage and the 14-day target unless a scientific rationale changes them before any outcomes.
4. Preflight precision and runtime analytically first. If a valid method cannot meet both targets with a feasible, registered seed block, stop the confirmatory comparison. Do not lower the target after seeing outcomes.
5. Only after those gates pass, implement the complete driver, test it on synthetic data in a scratch copy, register fresh seeds, freeze code/settings/environment hashes, and run the locked evaluation.

No policy outcome run, patient-data read, seed reservation, or result tuning is authorized by this memo.
