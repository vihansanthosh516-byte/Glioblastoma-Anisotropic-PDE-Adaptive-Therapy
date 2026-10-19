> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Historical Track C positioning: Track C (Digital Twin Reactor)

Per the ISEF outline, **Track C** is the component the outline called most revolutionary. That is the outline's wording; no literature check supports "first" or "revolutionary".

## Outline assertions (novelty unverified)

~~The outline explicitly lists these as "first work" achievements, all within Track C:~~

~~1. **First** to include anisotropy + stromal feedback in 3D digital twin~~
~~2. **First** to train RL on mechanistic models with biomarker-driven stratification  ~~
~~3. **First** to provide clinical decision rules (ρ > 0.024 day⁻¹ threshold → RL preferred)~~
~~4. The paper's one-liner conclusion: *"First end-to-end framework linking single-cell GRNs to RL adaptive therapy, with real-patient physics validation and actionable stratification rules."*~~

## Unsupported outline positioning

~~- **Integrative capstone**: Combines molecular insights (Track A) + spatial biophysics (Track B) into a clinically actionable pipeline~~
~~- **Real-patient validation credibility**: The 62-patient tensor comparison study (Figure 4) validates the core anisotropic modeling assumption directly on real patient data — this is positioned at the paper's start as the credibility-setting piece: *"Start strong: Real-patient validation study (Figure 4) ← this is your credibility"*~~
~~- **Tangible clinical benefit**: 13% cumulative drug reduction with equivalent tumor control (non-inferior time-to-progression)~~
~~- **Actionable stratification**: The proliferation rate biomarker (ρ > 0.024 day⁻¹) enables patient-specific therapy selection (>80% win rate in high-ρ tumors, vs 36.7% overall)~~
~~- **Mechanistic RL**: Policies trained on validated models outperform fixed schedules (RL 73.8 mm³ vs Stupp 137.8 mm³, 46.5% reduction, p=0.00067)~~

## Supporting Tracks

- **Track A** (Molecular → Spatial → Causal): Provides the GRN landscape, master switches (APOD, S100B, MT3), and invasion dynamics that parameterize Track C's models
- **Track B** (Anisotropic PDE + Therapy): Supplies the anisotropic diffusion tensor geometry and adaptive therapy policy space that Track C's RL operates within

## Summary

~~Track C is the most revolutionary because it **translates molecular and biophysical insights into a clinically viable adaptive therapy framework** — with real-patient validation, biomarker-driven stratification, and quantified drug-sparing benefits. It is the track that makes the entire computational platform actionable in a clinical context.~~

Audit correction (2026-10-19): The struck outline assertions are withdrawn.
No literature review establishes a first/revolutionary contribution or clinical
viability. Track C is a historical assumed-model policy experiment; Track A-to-B
measured-expression linkage is unverified (ledger 2m). The 62-case tensor
comparison is not independent therapy validation; its exact manifest/output
bundle is missing from this note. Dose-saving, non-inferiority and the
73.8/137.8-mm3 p-value claim need a matching run and prespecified endpoint/margin;
no clinical benefit follows. The old 0.024/day, >80% and 36.7% rule is superseded
by the distinct modeled script-62 threshold 0.0751/day, CI [0.0746,0.0820]
(`../output/biomarker_stability_metrics.json`, ledger 2). It is not actionable.
Later RL results remain simulations with assumed kill rates and no TTP win
(ledger C1). Current work: `../START_HERE.md`; corrected phase comparisons:
`docs/PROJECT_REVIEW_SUMMARY.md`. No model is certified as biologically validated.
