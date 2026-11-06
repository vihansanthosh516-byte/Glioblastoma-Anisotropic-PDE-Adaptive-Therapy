# Track C finite-set signal-ambiguity MPC prototype (2026-11-06)

**Status: synthetic prototype only. Not adopted as the primary method. No policy-outcome comparison or patient analysis ran.** The umbrella title remains **“Integrating Tumor Biology and Longitudinal MRI for Glioblastoma Progression Modeling and In-Silico Therapy Optimization.”**

## Why this prototype exists

The current signal-conditioned controller receives fixed values for the paired-signal parameters. The leave-one-out model-averaging predictive gate failed at cohort sizes 18 and 28, and the finite-sample recovery study found weak identification at those sizes. A controller that treats one assumed signal model as known would hide this uncertainty. The prior results and their limits are in [`trackc_method_reconciliation_and_next_steps_2026-11-06.md`](trackc_method_reconciliation_and_next_steps_2026-11-06.md).

## Proposed decision rule

Given a finite, predeclared set of signal models \(\mathcal M\), the prototype keeps one particle filter per model. Every filter receives the same patient’s paired baseline signals and subsequent visible burden scans. For each candidate treatment pattern \(a\), it estimates posterior-predictive utility \(Q_m(a)\) under each model \(m\). It chooses:

\[
a^* = \arg\max_a \min_{m \in \mathcal M} Q_m(a).
\]

The policy propagates the chosen action through every filter. Model order is canonicalized. Ties go to the first pattern in the existing deterministic pattern order. This is a finite-set minimax policy. It does not learn a biological signal model or validate any signal-to-resistance relationship.

## Implementation

- New module: [`src/216_trackc_robust_signal_ambiguity.py`](../src/216_trackc_robust_signal_ambiguity.py). The number is 216 because script 215 is already assigned to the LOO predictive check.
- New tests: [`tests/test_trackc_robust_signal_ambiguity.py`](../tests/test_trackc_robust_signal_ambiguity.py).
- A caller must provide the model set. No evaluation outcomes select models or settings.
- The default latent resistance-fraction support matches the simulator range, 0.001–0.1. A different world must pass its frozen range explicitly.
- A one-model set matches the existing fixed-likelihood controller when both use the same prior range.
- The policy is not wired into `src/211_trackc_policy_benchmark.py`.

## Mechanics verification

The new suite covers worst-case action selection, deterministic ties, invalid and duplicate model sets, simulator-matched default support, separate model-conditioned filters, missing-scan handling, one-model equivalence, course-budget behavior, and model-order invariance. The four focused suites passed **68 tests in 23.86 seconds** in a scratch copy. The tested source and test hashes matched the worktree. This validates mechanics only. It does not show that the controller improves utility, is clinically useful, or is novel.

The first one-model equivalence test failed because a cloned filter used a separate random generator. The first model now keeps the established controller’s random-stream behavior. A second failure came from comparing controllers with different \(f_0\) supports; the test now supplies the same support to both. The final focused run passed.

### Follow-up audit (2026-11-06)

The initial filename used `215`, which already belongs to the LOO predictive-check script. The prototype is now `216_trackc_robust_signal_ambiguity.py`. The seed registry now records controller RNG seeds separately from scenario seeds and observation-noise keys, detects numeric collisions across streams, and records the discarded 20261111/20261112 timing attempt. A follow-up scratch run passed **80 tests, with 1 skipped in 25.12 seconds**, across the ambiguity, signal-conditioning, robust-MPC, LOO, and seed-registry suites. The skip is the registry test that requires a separate sealed-run artifact omitted from the scratch copy. These remain code and bookkeeping checks, not policy outcome evidence.

## Synthetic timing

One timing call measured reset, baseline signal conditioning, and the first plan only. It used the existing 81-model grid, 128 filter particles, 8 prediction particles, a 3-cycle scan interval, and a 3-cycle forecast. The fresh-seed timing record is [`output/trackc_robust_signal_ambiguity_20261106/timing.json`](../output/trackc_robust_signal_ambiguity_20261106/timing.json).

For the 81-model grid, elapsed time was 15.185 seconds for one synthetic patient, 29.860 seconds for four, and 50.208 seconds for eight. These single timing observations do not provide a full-study runtime estimate. Repeated scans, all benchmark arms, world replication, and runtime variation remain unmeasured.

**Seed correction:** an initial timing attempt used scenario seed 20261111, which is already registered for Track C script 166 development, and controller RNG seed 20261112, the same integer as script 166’s sealed scenario seed. That attempt produced timing only. No endpoint, utility, or policy comparison value was saved or examined. Its timings were discarded. The timing pilot was repeated with checked seeds 20262604 and 20262605. The registry records both the discarded attempt and the fresh timing run.

## Prior-art and novelty boundary

The broad title is not itself a novelty claim. M4RL already joins multiscale GBM biology and reinforcement-learning drug scheduling ([Liu et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC13155510/)). A 2026 preprint joins reaction–diffusion growth, residual learning, digital-twin updating, and MPC treatment scheduling ([W. Liu et al., 2026](https://arxiv.org/abs/2607.13877)). A 2026 GBM paper uses multimodal MRI, a latent neural ODE, and proof-of-concept deterministic and reinforcement-learning control, but derives its dynamics from cross-sectional rather than longitudinal scans ([Volovăț et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC12785452/)). A separate 2026 preprint studies robust adaptive treatment policies under model and measurement uncertainty in a metastatic-prostate-cancer simulator, not GBM ([Xie et al., 2026](https://arxiv.org/abs/2610.22336)). These sources make broad multimodal GBM modeling, model-predictive treatment scheduling, reinforcement-learning control, or robust treatment planning alone insufficient novelty claims.

The exact use of a finite ambiguity set over paired molecular-like and MRI-like signal models in this GBM simulator remains unverified as prior art. That is a hypothesis for further searching, not a claim of novelty. The finite-set rule may be overly conservative or may lose to the no-signal controller. Its behavior depends on the supplied model set, which is analyst-defined and not fitted from linked patient data.

## Gates before any outcome comparison

1. Freeze the full and independence-only model sets, including parameter bounds and any out-of-grid truth anchors. Do not choose the set after viewing policy outcomes.
2. Wire the controller into the shared policy benchmark with identical virtual patients, observations, scan schedule, treatment budget, and endpoint handling.
3. Include the matched conditional-independence set, no-signal MPC, failed LOO model average as a secondary comparator, fixed schedules, equal-course frontier, Q-learning, and applicable scan/state-estimation controls.
4. Freeze the world-level estimand, uncertainty intervals, precision target, number of worlds, and fresh development seeds before scoring. Keep selection worlds disjoint from evaluation worlds.
5. Report benefit and harm across every prespecified world. Do not tune the set, policy, endpoint, or success rule after scores appear. Do not use sealed seeds.

The integrated title still exceeds the empirical connection between tracks: Track A expression, Duke MRI, and Track C virtual patients are not patient-linked. The signal-to-resistance mapping remains a simulator assumption. Any result from this prototype would support only a conditional simulation claim.
