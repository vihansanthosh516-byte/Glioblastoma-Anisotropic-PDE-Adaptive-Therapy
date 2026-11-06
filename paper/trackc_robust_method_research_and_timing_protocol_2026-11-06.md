# Track C robust-method research update and timing protocol

**Status:** prior-art update complete; timing-only pilot protocol frozen; pilot not yet run. The umbrella title remains **“Integrating Tumor Biology and Longitudinal MRI for Glioblastoma Progression Modeling and In-Silico Therapy Optimization.”** This note does not amend the patient analyses or authorize a sealed evaluation.

## Prior-art update

The broad title-level combination is not novel by itself. The 2025 M4RL study combines a multiscale glioblastoma model with a physics-informed surrogate and reinforcement learning to schedule drug combinations. The authors use simulated tumor dynamics and spatial-transcriptomic data to initialize virtual tumor microenvironments. This is close prior art for biology-informed GBM simulation plus adaptive drug scheduling ([Liu et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC13155510/)).

A patient-specific high-grade-glioma digital twin has also combined MRI-based calibration, uncertainty analysis, and risk-aware adaptive radiotherapy planning ([predictive digital-twin study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10598726/)). A 2026 preprint combines reaction-diffusion growth, learned residual correction, recursive digital-twin updates, and model-predictive chemotherapy/radiotherapy scheduling on 387 synthetic trajectories ([Liu et al., 2026](https://arxiv.org/abs/2607.13877)). These works make “tumor biology + MRI + adaptive therapy” an unsupported novelty claim.

Uncertainty analysis for MRI-informed glioblastoma response metrics also predates this project ([Hawkins-Daarud et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6633916/)). General robust POMDP and minimax control methods are established outside this exact GBM use. The narrower candidate question remains whether a finite-set minimax controller over paired molecular-like/MRI-like signal-dependence models changes simulated policy value relative to a matched `d=0` controller. This exact combination remains **unverified** as novel; the targeted search is not exhaustive.

## Integration boundary

Track A and Duke MRI records have no shared patient identifiers. Track A's A5 balanced-accuracy values are not a measured assay of drug resistance. The simulator's links from molecular-like signal to resistance and from MRI-like signal to tumor state remain assumptions. Track B's patient-derived MRI analysis remains separate. The proposed integrated result is therefore a synthetic decision-method stress test informed by separate summaries, not a validated patient-specific therapy system.

The LOO signal-prediction gate failed at `n=18` and `n=28`; its spent seeds cannot tune a replacement. The finite-set minimax controller and matched `d=0` arm remain development candidates. No policy-outcome comparison has shown that either improves results. Any positive synthetic result would remain conditional on the stated simulator.

## Timing-only pilot

### Purpose and status

This pilot will estimate elapsed time for the current full `run_shared_policy_world` harness with the finite-set minimax arm and matched `d=0` arm enabled. It will test the current default computation path before any multi-world outcome comparison. It will not select methods, settings, anchors, or hypotheses.

The function returns synthetic policy outcomes in memory. The pilot script will discard each returned object without reading or saving outcome fields. It will save only elapsed seconds, environment details, source hashes, protocol hash, and seed assignments. This is a runtime measurement, not performance evidence.

### Frozen pilot design

- Use one synthetic world each from the `single` and `mixture` response families.
- Use three disjoint cohorts per world: selection, calibration, and evaluation; each cohort will contain two virtual patients. Patient IDs will be unique across cohorts.
- Use the same predeclared signal-generation setting in both worlds: `q=0.7`, `a=0.75`, `b=0.75`, `d=0.75`. This is a timing setting, not a selected evaluation anchor.
- Use the full 81-point signal-model grid from `src/212_trackc_signal_fisher_grid.py`. The matched independence arm will use the 27-point grid with the same `q/a/b` support and `d=0`.
- Enable both finite-set arms in `run_shared_policy_world`. Use `PolicyBenchmarkSettings()` defaults: scan interval 3 cycles, 730-day horizon, 5-day course penalty, 500 MPC particles, 16 prediction particles, six forecast cycles, and 50 Q-learning episodes. Other settings will remain at code defaults.
- Run the two worlds sequentially on an otherwise idle computer. Do not launch another timing job at the same time.
- If a run fails, record the exception type and message, discard any partial outcomes, and stop. Do not retry with changed settings.

### Seed allocation

All seeds come from the reserved Stage 2 development block; none come from the Stage 3 sealed block.

| Response family | World seed | Scenario and signal-generation seeds | Policy seeds |
| --- | ---: | --- | --- |
| `single` | 20261930 | 20261931–20261935 | 20261936–20261940 |
| `mixture` | 20261941 | 20261942–20261946 | 20261947–20261951 |

Within each row, scenario seeds will be assigned in order to selection, calibration, and evaluation cohorts; the next two seeds will generate evaluation and calibration signals. Policy seeds will be assigned in order to schedule selection, Q-learning training, calibration simulation, evaluation simulation, and schedule-mixture assignment. All numeric seeds will be distinct across streams.

### Record and stop rule

The pre-run manifest will include the exact seed map, cohort sizes, response family, signal parameters, grid sizes, code and protocol hashes, environment, and default settings. It will be written before either run starts. The post-run record will add elapsed seconds and terminal status only. It will not contain policy scores, utilities, endpoint times, courses, or rank ordering.

No patient data, MRI, masks, or restricted files will be read. The timing record will not authorize or start a performance comparison. After this pilot, update the seed registry and verify the output before launching any further run. A separate dated protocol must freeze evaluation worlds, world-level estimand, precision target, interval method, and outcome anchors before any policy scores are examined.
