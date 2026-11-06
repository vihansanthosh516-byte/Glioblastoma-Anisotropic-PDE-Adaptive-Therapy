# Track C: peer-only model averaging for paired baseline signals

**Status: synthetic method prototype; mechanics tested; no policy outcomes or patient data evaluated.** The project title remains **“Integrating Tumor Biology and Longitudinal MRI for Glioblastoma Progression Modeling and In-Silico Therapy Optimization.”** Track A and Duke are not linked to this Track C simulator.

## Proposed method

The signal model has four shared parameters, `q`, `a`, `b`, and `d`. A fixed 81-point grid supplies the candidate parameter settings. For each virtual patient *i*, the method calculates model weights from the paired baseline signals of all other virtual patients. It integrates each candidate model over the declared reference distribution for `f0`. It then averages patient *i*'s conditional signal likelihood over those peer-only weights and uses the mixture to update that patient's `f0` particles.

The method uses no latent resistance values, future simulator outcomes, or treatment responses to estimate the peer weights. It assumes all patients in a synthetic cohort share one signal-generating model. It is a cohort-calibrated simulation method, not a validated clinical biomarker workflow.

## Implementation and verification

- Code: `src/214_trackc_loo_model_average.py`.
- Tests: `tests/test_trackc_loo_model_average.py`.
- The focused test file passed 7 tests in the required scratch copy. Combined with the existing signal-conditioning, robust-MPC, and paired-benchmark suites, 62 tests passed in scratch.
- Tests cover leave-one-out isolation (changing a patient's own signal does not change that patient's model weights), row-normalized model weights, exact zero prior mass, equivalence to the direct likelihood for a one-model mixture, invalid inputs, and finite MPC initialization.
- A zero prior weight remains exactly zero. Peer likelihood sums explicitly omit the held-out patient to prevent floating-point cancellation from allowing the held-out value to affect its own weight.
- **Execution note:** one initial 7-test command was accidentally launched from the working tree. The same tests were then rerun in scratch, and the combined 62-test suite passed in scratch. This was synthetic-only and read no patient data; the command-location mistake is recorded here.

## Evidence and limits

This implementation is not a result that the method is better, unique, or biologically correct. The simulator chooses the molecular-to-resistance and MRI-to-resistance relationships. Current finite-sample recovery results show weak recovery at sample sizes comparable to the project's small virtual cohorts; see [`trackc_finite_sample_recovery_protocol_2026-11-06.md`](trackc_finite_sample_recovery_protocol_2026-11-06.md) and [`output/trackc_finite_sample_signal_recovery/diagnostic.json`](../output/trackc_finite_sample_signal_recovery/diagnostic.json). The new method must therefore pass a prospective held-out signal-prediction diagnostic before any treatment-policy claim is tested.

The broader frame is not new by itself. Prior work already describes MRI-personalized digital twins and risk-aware radiotherapy optimization in high-grade glioma ([Chaudhuri et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10598726/)). A 2025 primary study uses longitudinal multiparametric MRI to personalize a family of biology-based high-grade glioma models, averages model weights, and selects the highest-weight model for cohort-level error analysis ([Hormuth et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12815495/)). These overlaps mean the project must not claim that multimodal digital twins, MRI-informed therapy planning, or model averaging are new. A future contribution, if supported, would be the specific controlled simulation result about peer-calibrated dependence-aware baseline signals versus matched controls. Exact prior-art coverage for that intersection is still unverified.

## Next gate, before treatment-outcome evaluation

The synthetic leave-one-out predictive-score diagnostic is now frozen in [`trackc_loo_predictive_check_protocol_2026-11-06.md`](trackc_loo_predictive_check_protocol_2026-11-06.md); it has not run. It compares the peer-weighted mixture with a fixed uniform-grid mixture, a prespecified nominal model, a matched conditional-independence model, and the generating-parameter oracle on the same generated cohorts. It uses fresh registered development seeds, includes the `q=0.5` negative control, reports every predeclared regime and sample size, and uses held-out paired-signal log predictive density. If the peer-weighted method does not pass the frozen gate, do not advance it as the primary policy method.

Only if this gate is met should the project freeze a paired policy comparison. That comparison must include both progression definitions, same virtual worlds and scan-noise streams, equal treatment budgets, uncertainty summaries, and clinically named baseline policies. It must report harms and failures as well as favorable regimes. Any positive result remains an in-silico result under the chosen simulator assumptions, not evidence of clinical benefit.
