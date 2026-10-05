# Paper: Figure Inventory and Plan

Written 2026-10-04. I listed files and read the JSON that feeds them. **I have not opened any PNG.** Check each image by eye before use.
`paper/figures/` exists and is empty.

## Rule
Make every paper figure from a JSON or CSV with one script (suggest `scripts/make_paper_figures.py`). Do not reuse old figures that were drawn when numbers were different.

## Main-text figures

| Planned | Finding | Existing file | Status |
|---|---|---|---|
| Fig 1 (A1) classification leaderboard | C-GAT vs baselines | `output/benchmark_bar_chart.png` | **Unknown content.** The benchmark TSV has only LR, RF, Transformer, Hybrid. C-GAT and scVI are not in it. **Generate new** from `gat_metrics.json` + `benchmark_comparison.tsv` + scVI/NMF numbers. |
| Fig 1 (A2) SPIB saddle | index-1 saddle | `output/spib_landscape_validation.png`, `spib_rc_validation.png`, `energy_potential.png` | Exist (Sep 19). Probably usable. Check. |
| Fig 1 (A3) TCGA | age HR | `output/forest_plot.png`, `km_survival_curves.png` | Exist. Check they match HR 1.031 and n=518. |
| Fig 2 (B1) adaptive vs MTD | resistance separation | `output/adaptive_therapy_comparison.png`, `adaptive_therapy_metrics_summary.png`, `adaptive_therapy_dynamics.png` | Exist (Sep 23). Must show **final mass** too, not only resistance. Likely generate new. |
| Fig 2 (B2) forecast | core-target Dice | `output/forecast_validation/forecast_results.png` is the **old `mask>0` target**. No PNG in `forecast_labels_core/`. | **Generate new** from `output/forecast_labels_core/results.json`. Show all-patients, grew, shrank. |
| Fig 2 (B3) elongation | aniso vs iso | none | **Generate new** from `output/fractal_aniso_vs_iso.json` (`patients` list). |
| Fig 3 (C1) conditioned policy | win rate | `output/rl_kill_conditioned/kill_conditioned_policy.png` | Exists. Check it shows the rule and the day-90 caveat. |
| Fig 3 (C2) probe-then-commit | noise curve | none | **Generate new** from `output/probe_paced_policies/evaluate.json`. |
| Fig 3 (C3) paced TTP | 4 sets | `output/ttp_equal_budget/ttp_equal_budget.png` is script 68, not paced | **Generate new** from `probe_paced_policies/evaluate.json`. Show the real_test +2.1 d bar. |

## Do NOT use
| File | Reason |
|---|---|
| `output/sobol_tornado_plot.png` | Drawn from the retracted reduced-ODE result (S1 = 0.998). Regenerate from `output/sobol_pde_swanson_vs_legacy.json` if a Sobol figure is wanted. |
| `output/ablation_study_figure.png` | Script 60 legacy D (10x too low). Use `ablation_study_figure_swanson.png` or a script 78 figure. |
| `output/figures/fig1..fig8_*.png` | README-era figures (architecture, RL, FNO, radar, roadmap). Numbers not audited. `fig7_cross_track_radar` and `fig2_phase5_rl_results` likely carry stale claims. Fig 1 (architecture) may be reusable as a schematic after a check. |
| `output/65_master_summary_figure.png`, `master_cohort_synthesis.png`, `outputs/validation_summary_dashboard.png` | Summary posters from older runs. Not audited. |
| `output/anisotropic_geometry_summary.png` | Likely shows D_f (retracted metric). Check first. |

## Useful for supplement
`output/mgmt_resistance/km_mgmt_ttp.png` (script 82, KM curves). `output/csgt_mathematical_proof.png`. `output/invasion_dynamics_analysis.png`. `output/fig1_bar_chart.png`, `fig2_scatter.png`, `fig3_delta_boxplot.png` (UCSF tensor study; check they match [[UCSF-PDGM-Tensor-Study]]). `output/batch_digital_twins/` holds about 150 per-patient plots (not needed).

## Missing figures to make (7)
A1 leaderboard, B1 final-mass panel, B2 forecast, B3 elongation, C2 probe-then-commit, C3 paced TTP, plus an optional framework schematic.
