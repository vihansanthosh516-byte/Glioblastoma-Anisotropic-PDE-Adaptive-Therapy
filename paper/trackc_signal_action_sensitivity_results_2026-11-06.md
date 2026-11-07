# Track C paired-signal action-sensitivity result (2026-11-06)

**Status: completed synthetic diagnostic.** It tests whether repeated paired signals change a controller's selected plan. It does not test whether treatment improves tumor outcomes.

## Frozen design and records

The run followed the committed [protocol](trackc_signal_action_sensitivity_protocol_2026-11-06.md), [pre-run manifest](../output/trackc_signal_action_sensitivity_20261106_v2/manifest.json), and Stage 2 development-seed entry in `src/trackc_seed_registry.py`. The project title stayed **“Integrating Tumor Biology and Longitudinal MRI for Glioblastoma Progression Modeling and In-Silico Therapy Optimization.”** The result is [diagnostic.json](../output/trackc_signal_action_sensitivity_20261106_v2/diagnostic.json).

The diagnostic made 128 planning calls: 64 paired signal cases across four synthetic states, with the full 81-model set and its matched 27-model `d=0` set. It used one assumed signal generator, exact baseline `log(B0)`, no burden measurement noise, and no treatment-outcome simulation. No patient data entered the run.

The result contains 128 rows. The manifest hash, six source hashes, and seed-registry hash match their recorded files. Environment: Python 3.14.4, NumPy 2.4.4, SciPy 1.18.0, Windows 11. The per-call planning timers sum to 3,768.29 seconds; this excludes process startup and other run overhead.

## Results

The full model set selected more than one action pattern across repeated signals in three of four states. The matched `d=0` set did so in one of four. The two sets selected the same pattern in 43 of 64 paired cases (67.2%).

| Synthetic state | Full 81: distinct patterns | Matched `d=0`: distinct patterns | Paired agreement |
| --- | ---: | ---: | ---: |
| State 0 | 2 | 2 | 6/16 |
| State 1 | 1 | 1 | 16/16 |
| State 2 | 2 | 1 | 11/16 |
| State 3 | 2 | 1 | 10/16 |

The overall agreement and all pattern counts are descriptive. They show decision sensitivity under these chosen synthetic inputs. They do not show that either plan improves outcomes or that the full set is superior. The controller's planning-value margins are not observed treatment effects.

## Output and input-boundary deviations

The result rows contain the latent simulator value `f0`. The frozen protocol did not list it among the saved outputs. The runner used true `f0` to generate signals, then saved it after planning. The controller's reset path receives the scenario object but [RobustBeliefMPC.reset](../src/208_trackc_robust_belief_mpc.py) reads only `sigma` to initialize its filter; the controller conditions prior particles on the generated signals and does not use the scenario's true `f0` value. This output-field deviation does not establish decision leakage, but it permits unplanned post-hoc stratification. No analysis used `f0`; keep the saved result unchanged and do not stratify or select on that field.

A separate interface audit found that the diagnostic called `policy.reset(scenario, ...)` with the full simulator scenario dictionary. That dictionary contains latent values such as `f0`, `r`, `m`, `c`, and `eps`, even though the protocol said the controller would not receive latent truth. Inspection of the current reset path found that it reads `sigma` and initializes its own prior; no use of the other scenario fields was found. This means there is no demonstrated truth-use in the recorded plan choice, but the protocol's input-isolation boundary was not enforced. Treat this as an interface-level protocol deviation. Do not alter or rerun the historical diagnostic; future policy runs must pass a whitelist of observable fields only.

## Checks and limits

The final focused scratch suites passed 26 tests; one seed-registry evidence test skipped because the scratch copy omitted `output/`. A broader controller/simulator suite passed 127 tests with one such skip before the last manifest-metadata guard change. These tests check code behavior, not scientific benefit.

The result uses four states, one assumed signal-generating parameter set, a dependence value outside the controller grid, and an exact baseline burden. Its reset interface was not truth-blinded, and its output includes an unplanned latent field. It cannot establish biological validity, clinical value, method uniqueness, or superiority over existing therapy policies. Track A and Track B data remain unlinked; the molecular-to-resistance relation remains a Track C simulation assumption.

## Next test

Freeze a separate multi-world outcome benchmark before scoring. Compare the signal-ambiguity controller with matched `d=0`, a no-signal robust controller, standard MPC, and predeclared fixed-schedule controls. Match treatment budgets and evaluate both progression definitions on independent Stage 2 development worlds. Report beneficial and harmful worlds, paired uncertainty intervals, and all failures. Use fresh seeds. Do not claim a best or unique method unless this outcome benchmark and a current prior-art review support that claim.
