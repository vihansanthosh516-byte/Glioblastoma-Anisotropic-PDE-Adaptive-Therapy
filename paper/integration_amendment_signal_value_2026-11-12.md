# Integration amendment: test the added value of combining biology and MRI (2026-11-12)

Written before any Track C evaluation world or development-pilot world has run. No policy outcomes have been seen.

## Overall project

*Integrating Tumor Biology and Longitudinal MRI for Glioblastoma Progression Modeling and In-Silico Therapy Optimization.*

- **Track A (tumor biology):** which biology signals relate to tumor region or measured treatment response, with patients kept separate in validation. Current A5 / TMZ results did not show the planned signal.
- **Track B (MRI over time):** do earlier scans predict later tumor burden better than simple forecasts, with uncertainty? ACRIN failed its rule; LUMIERE stopped at registration QC.
- **Track C (simulator):** tumor growth, drug suppression and treatment choice under uncertainty; robust policy vs matched policies in simulated worlds.
- **Integration:** do biology and longitudinal MRI signals **together** improve simulated treatment decisions beyond either signal alone or no signal, on the same simulated patients with fair comparison policies?

## Change

1. **Integration question (new, primary for the integration).** In the primary cell, the combined controller (`finite_set_minimax_mpc`) is compared with three matched controls that use the same controller, model set, prior, patients and seeds:
   - `combined_vs_molecular_only` (control `finite_set_minimax_molecular_only_mpc`)
   - `combined_vs_mri_only` (control `finite_set_minimax_mri_only_mpc`)
   - `combined_vs_no_signal` (control `no_signal_mpc`)
2. **Rule.** Each contrast uses the frozen Track C rule unchanged (GE-betting, delta 0.025 per family, both family lower bounds > 0 and both half-widths <= 14 days). **Integration positive iff all three contrasts are positive** (intersection-union test: no multiplicity correction is needed because every contrast must pass).
3. **Track C primary unchanged.** `minimax` (link-aware vs matched `d = 0`) stays the Track C primary contrast and the method-novelty test. Its rule is not changed.
4. Code: the three contrasts are now named and saved per world (`src/211_trackc_policy_benchmark.py`); the analysis reports them (`src/trackc_v2_analysis.py`). The four arms were already required by `check_arm_plan` (`src/trackc_integrated_driver.py`).

## Weak points

- Simulation only. No patient-level link exists across tracks. This cannot show that treatment helps real patients.
- Signal strengths are still assumed (config anchor). They become data-grounded only if Track A/B measured values are put into the config before the Spark run; this is a separate open step.
- Novelty: combining biology and MRI for treatment is established in general (`paper/trackc_integrated_control_prior_art_update_2026-11-06.md`). The distinctive parts are the fair added-value test and the link-aware controller; a full prior-art search is still open.
- Requiring all three contrasts to pass is strict; power is lower than for one contrast.
