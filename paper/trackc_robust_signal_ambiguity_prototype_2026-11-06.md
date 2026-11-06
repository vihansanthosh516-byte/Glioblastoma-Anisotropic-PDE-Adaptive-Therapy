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
- The policy can be enabled as an opt-in arm in `src/211_trackc_policy_benchmark.py`; it remains outside the frozen primary contrast.

## Shared benchmark wiring (2026-11-06)

`run_shared_policy_world` now accepts `include_ambiguity_controller=False`. When enabled, the benchmark loads the finite-set minimax policy and supplies the same evaluation-cohort baseline observations, predeclared inference model set, MPC settings, virtual patients, scan schedule, endpoint rules, and evaluation random seed used by the other arms. It records the arm's inclusion status and model-set size. The existing primary contrast remains dependence-aware LOO MPC minus conditional-independence LOO MPC.

The integration test enabled the new arm on a two-patient synthetic evaluation cohort and a two-model inference set. It verified outcome-array fields, shared evaluation mechanics, and the unchanged primary contrast. This was a mechanics test: it did not save or interpret arm utility differences. No planned policy-outcome analysis has run. The controller remains exploratory and has not shown a benefit.

### Matched d=0 minimax control (2026-11-06)

The opt-in benchmark now also evaluates `finite_set_d0_minimax_mpc`. It uses the same observation data, MPC settings, virtual patients, and random streams as `finite_set_minimax_mpc`; its model set preserves the supplied q/a/b support and fixes d to 0. This tests whether allowing nonzero dependence values changes the finite-set minimax controller, holding the q/a/b support and other policy settings fixed. Both arms remain additional development arms outside the frozen primary contrast. The integration test checks the d=0 model-set size, output arrays, and unchanged primary contrast. It does not compare or interpret their outcome values.

## Mechanics verification

The new suite covers worst-case action selection, deterministic ties, invalid and duplicate model sets, simulator-matched default support, separate model-conditioned filters, missing-scan handling, one-model equivalence, course-budget behavior, and model-order invariance. The four focused suites passed **68 tests in 23.86 seconds** in a scratch copy. The tested source and test hashes matched the worktree. This validates mechanics only. It does not show that the controller improves utility, is clinically useful, or is novel.

The first one-model equivalence test failed because a cloned filter used a separate random generator. The first model now keeps the established controller’s random-stream behavior. A second failure came from comparing controllers with different \(f_0\) supports; the test now supplies the same support to both. The final focused run passed.

### Follow-up audit (2026-11-06)

The initial filename used `215`, which already belongs to the LOO predictive-check script. The prototype is now `216_trackc_robust_signal_ambiguity.py`. The seed registry now records controller RNG seeds separately from scenario seeds and observation-noise keys, detects numeric collisions across streams, and records the discarded 20261111/20261112 timing attempt. A follow-up scratch run passed **80 tests, with 1 skipped in 25.12 seconds**, across the ambiguity, signal-conditioning, robust-MPC, LOO, and seed-registry suites. The skip is the registry test that requires a separate sealed-run artifact omitted from the scratch copy. These remain code and bookkeeping checks, not policy outcome evidence.

After benchmark wiring, six focused suites passed **95 tests, with 1 skipped in 55.66 seconds**, in scratch: policy benchmark, ambiguity MPC, signal conditioning, robust MPC, LOO model averaging, and seed registry. The skipped test needs an `output/` sealed-run artifact omitted from the scratch copy. The benchmark mechanics test exercised a two-model ambiguity arm on two virtual patients; it checked output structure and primary-contrast stability, not comparative performance. All 441 source and test files matched the final worktree after the scratch test copy was synchronized.

### Uninformative-signal negative control (2026-11-06)

A new end-to-end regression test sets `q=0.5` and `a=b=d=0`, so both synthetic signals carry no information about resistance. Under identical settings and random seeds, the finite-set controller exactly matched no-signal MPC for baseline endpoint time, RANO-like endpoint time, and courses used in a two-patient synthetic run. The six focused suites now pass **97 tests, with 1 skipped in 84.64 seconds**, in a scratch copy. The skipped registry test needs a sealed-run artifact omitted from that copy. SHA-256 comparison confirmed that the copied `src` and `tests` matched the worktree. This validates a null-information behavior only; it is not evidence of policy benefit, patient relevance, or novelty.

After adding the matched d=0 minimax control, the same six focused suites passed **97 tests, with 1 skipped in 94.53 seconds**, in a fresh scratch copy. SHA-256 comparison confirmed that all 443 copied Python source and test files matched the worktree, excluding generated `__pycache__` files. The skipped registry test needs a sealed-run artifact omitted from the copy. The new check verifies plumbing only; no policy outcomes were generated or compared.

### Forecast-sampling fairness correction (2026-11-06)

Before any policy-outcome comparison, code review found that each model's action values used 16 random posterior-particle draws from a sequential generator. The 81-model candidate and 27-model d=0 control therefore gave matched d=0 models different forecast draws. Taking a minimum over many noisy estimates can also bias the candidate's estimated worst-case value downward. This was a comparison-quality risk, not evidence that either policy wins.

The multi-model controller now uses keyed, stratified uniforms for each virtual patient and scan cycle. Every model uses the same uniforms, independent of model order and set size. A regression test verifies that a shared d=0 model selects identical posterior particles in the full and reduced sets when its filter state matches. The one-model path remains unchanged. Eight focused suites passed **107 tests, with 2 skipped in 88.21 seconds**, in a scratch copy. The skips required an omitted sealed-run artifact and an omitted timing record. SHA-256 comparison found no differences in copied Python source or test files after excluding generated `__pycache__` files.

This correction improves paired Monte Carlo sampling. It does **not** remove the minimum-of-noisy-estimates bias, establish policy benefit, or validate the signal-to-resistance mapping. Before any outcome comparison, test action stability as `prediction_particles` increases and compare the stagewise minimax rule with a persistent-model uncertainty treatment. The current policy recomputes the worst model at each scan; it does not maintain one posterior probability for a fixed unknown global signal model.

## Synthetic timing

One timing call measured reset, baseline signal conditioning, and the first plan only. It used the existing 81-model grid, 128 filter particles, 8 prediction particles, a 3-cycle scan interval, and a 3-cycle forecast. The fresh-seed timing record is [`output/trackc_robust_signal_ambiguity_20261106/timing.json`](../output/trackc_robust_signal_ambiguity_20261106/timing.json).

For the 81-model grid, elapsed time was 15.185 seconds for one synthetic patient, 29.860 seconds for four, and 50.208 seconds for eight. These single timing observations do not provide a full-study runtime estimate. Repeated scans, all benchmark arms, world replication, and runtime variation remain unmeasured.

**Seed correction:** an initial timing attempt used scenario seed 20261111, which is already registered for Track C script 166 development, and controller RNG seed 20261112, the same integer as script 166’s sealed scenario seed. That attempt produced timing only. No endpoint, utility, or policy comparison value was saved or examined. Its timings were discarded. The timing pilot was repeated with checked seeds 20262604 and 20262605. The registry records both the discarded attempt and the fresh timing run.

## Prior-art and novelty boundary

The broad title is not itself a novelty claim. M4RL already joins multiscale GBM biology and reinforcement-learning drug scheduling ([Liu et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC13155510/)). A 2026 preprint joins reaction–diffusion growth, residual learning, digital-twin updating, and MPC treatment scheduling ([W. Liu et al., 2026](https://arxiv.org/abs/2607.13877)). A 2026 GBM paper uses multimodal MRI, a latent neural ODE, and proof-of-concept deterministic and reinforcement-learning control, but derives its dynamics from cross-sectional rather than longitudinal scans ([Volovăț et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC12785452/)). A separate 2026 preprint studies robust adaptive treatment policies under model and measurement uncertainty in a metastatic-prostate-cancer simulator, not GBM ([Xie et al., 2026](https://arxiv.org/abs/2610.22336)). These sources make broad multimodal GBM modeling, model-predictive treatment scheduling, reinforcement-learning control, or robust treatment planning alone insufficient novelty claims.

The exact use of a finite ambiguity set over paired molecular-like and MRI-like signal models in this GBM simulator remains unverified as prior art. That is a hypothesis for further searching, not a claim of novelty. The finite-set rule may be overly conservative or may lose to the no-signal controller. Its behavior depends on the supplied model set, which is analyst-defined and not fitted from linked patient data.

### Prior-art follow-up (2026-11-06)

Targeted primary-source searches found closer overlap and narrowed the novelty boundary:

- A prospective study of 21 patients with high-grade glioma personalized a family of biology-based tumor models from baseline-to-mid-treatment mpMRI, then forecast later tumor response during and after chemoradiation ([Hormuth et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12815495/)). It directly overlaps patient-specific biology-based longitudinal MRI forecasting. It does not establish the same paired molecular-like/MRI signal uncertainty policy comparison.
- Glioblastoma reinforcement learning for temozolomide scheduling predates this project ([Zade et al., 2020](https://doi.org/10.1016/j.cmpb.2020.105443)). M4RL later combined multiscale GBM biology with reinforcement-learning drug-combination scheduling ([Liu et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC13155510/)).
- A 2026 preprint combines reaction-diffusion, residual learning, adaptive digital-twin updating, and MPC treatment scheduling on synthetic tumor trajectories ([W. Liu et al., 2026](https://arxiv.org/abs/2607.13877)). A multicenter 2026 GBM study compares multimodal MRI fusion architectures for survival prediction ([Frontiers Oncology, 2026](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2026.1867582/full)).
- Worst-case control over uncertain POMDP transition or observation models is established outside oncology. Examples include moment-based distributionally robust POMDPs ([Nakao et al., 2020](https://epubs.siam.org/doi/10.1137/19M1268410)), robust uncertain-POMDP synthesis ([Suilen et al., IJCAI 2020](https://www.ijcai.org/proceedings/2020/0569.pdf)), and 2026 analysis of robustness to observation-model perturbations ([Kraske et al., 2026](https://arxiv.org/abs/2604.21256)). A finite-set minimax policy alone is therefore not a method novelty claim.

The narrow remaining candidate is an application and evaluation question: whether a predeclared uncertainty set over paired biology-like and MRI-like signals changes treatment decisions and improves held-out simulated outcomes in a GBM-motivated model, after comparison at equal drug use with strong controls. The model set is not empirically calibrated from patient-linked Track A and MRI records. Any positive result would remain conditional on simulator assumptions. The targeted search is not exhaustive; novelty remains unverified.

## Gates before any outcome comparison

1. Freeze the full and independence-only model sets, including parameter bounds and any out-of-grid truth anchors. Do not choose the set after viewing policy outcomes.
2. **Complete:** wire the controller into the shared policy benchmark with identical virtual patients, observations, scan schedule, treatment budget, and endpoint handling. Mechanics tests pass; performance remains unknown.
3. Include the matched conditional-independence set, no-signal MPC, failed LOO model average as a secondary comparator, fixed schedules, equal-course frontier, Q-learning, and applicable scan/state-estimation controls.
4. Freeze the world-level estimand, uncertainty intervals, precision target, number of worlds, and fresh development seeds before scoring. Keep selection worlds disjoint from evaluation worlds.
5. Report benefit and harm across every prespecified world. Do not tune the set, policy, endpoint, or success rule after scores appear. Do not use sealed seeds.
6. Before scoring, show that robust action choices are stable across larger posterior-prediction sample counts. Compare stagewise minimax with a persistent-model or model-averaged control; label any remaining model-switching assumption.

The integrated title still exceeds the empirical connection between tracks: Track A expression, Duke MRI, and Track C virtual patients are not patient-linked. The signal-to-resistance mapping remains a simulator assumption. Any result from this prototype would support only a conditional simulation claim.
