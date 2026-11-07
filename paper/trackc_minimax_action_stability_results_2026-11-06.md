# Track C minimax action-stability results (2026-11-06)

**Status: synthetic controller diagnostic completed. No treatment-outcome comparison ran. No patient data were read.** This result does not establish therapy benefit, biological validity, clinical value, or novelty.

## Record

The pre-run settings are frozen in [the protocol](trackc_minimax_action_stability_2026-11-06.md). The run produced 72 plans: four analyst-chosen signal anchors, two model sets, three posterior-prediction sample counts, and three prediction seeds. The 81-model set used `q ∈ {0.5, 0.7, 0.9}` and `a,b,d ∈ {-0.5, 0, 0.5}`. The matched control used the same `q/a/b` grid with `d=0`.

The run used one synthetic patient state, 128 posterior particles, exact `log(B0)` as the baseline burden reading, no observation noise, a three-cycle forecast, and a three-cycle scan interval. It recorded plans and planning values only. It did not simulate realized outcomes or calculate time-to-progression, survival, or treatment benefit.

The saved record is [`output/trackc_minimax_action_stability_20261106/diagnostic.json`](../output/trackc_minimax_action_stability_20261106/diagnostic.json), with its pre-run manifest beside it. All six code/protocol hashes match the worktree. The recorded seed-registry hash matches the file used at run time. The JSON has 72 unique finite rows and states `outcome_metrics_included: false`.

The first launch stopped at controller construction because duplicate dynamic imports created incompatible model classes. It ended before the planner forecast or any action was computed and wrote no result. The loader was corrected before the final run. The focused scratch suites then passed 32 tests; one seed-registry test skipped because its required output artifact was absent from the scratch copy.

## Observed plan stability

| Model set | 16 vs. 32 draws | 32 vs. 64 draws | Plans at 64 draws |
| --- | ---: | ---: | --- |
| Full 81 models | 8/12 same (66.7%) | 12/12 same (100%) | `100` in all anchors and roots |
| Matched `d=0` set | 4/12 same (33.3%) | 8/12 same (66.7%) | `100` in all anchors and roots |

Here, `100` means treat in cycle 1, then do not treat in cycles 2 or 3. The four signal anchors produced the same plan pattern at each sample count. At 64 draws, both model sets selected the same plan for all four anchors and all three prediction seeds. At 32 draws, the full and matched-`d=0` sets agreed in 8/12 paired cases; they agreed in 12/12 cases at 16 and 64 draws.

The full set was stable from 32 to 64 draws in this narrow diagnostic. The matched `d=0` set was not: four of twelve paired plans changed. The shared `100` choice at 64 draws does not show that the signal inputs improve decisions; the choice stayed unchanged across all four signal anchors.

The sum of measured planning-call times was 1,759.3 seconds. This is a single-machine development timing, not a general runtime estimate.

## Interpretation and next test

This is **partial numerical evidence for one synthetic state**, not a pass of the broader action-stability gate. Four hand-chosen signal inputs and three Monte Carlo seeds cannot establish stability across the simulator's state space. The exact-baseline reading without measurement noise leaves stability under measurement error untested and may overstate robustness. The test did not compare treatment outcomes.

The next development test should use multiple independently generated synthetic states, include held-out signal-model settings outside the fitted 81-point grid, add observation noise, and compare stagewise minimax against a persistent-model controller and no-signal MPC. It should keep the same virtual patients and treatment budgets across arms, predeclare world-level outcome metrics, and report harmful as well as helpful regimes. Do not tune the signal map or success rule to make the minimax arm win. Do not use sealed seeds.

## Prior-art update

The broad integrated project theme is already crowded. A 2026 preprint combines reaction-diffusion tumor growth, a learned residual, repeated digital-twin updating, and therapy scheduling by MPC on synthetic brain-tumor trajectories; it reports improved simulated forecasts and treatment scheduling. This overlap means “integrated tumor model plus adaptive therapy” is not a sufficient novelty claim. A separate 2025 GBM study combines multiscale mathematical modeling with reinforcement-learning combination-treatment scheduling. A 2026 preprint develops robust treatment-regime learning under model and measurement uncertainty in a metastatic-prostate-cancer simulator, showing that uncertainty-aware adaptive therapy is also active outside GBM. See [Liu et al. (2026)](https://arxiv.org/abs/2607.13877), [M4RL (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC13155510/), and [Xie et al. (2026)](https://arxiv.org/abs/2610.22336).

The narrow finite-set paired-signal question remains a development hypothesis, not a verified novelty claim. Track A expression, Track B MRI, and Track C virtual patients are not linked by patient identity. The signal-to-resistance mapping is a simulator assumption. The chosen umbrella title can remain, but the current evidence does not show one validated, patient-specific integrated therapy framework or a method that is best.
