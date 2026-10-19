> **HISTORICAL BACKGROUND (checked 2026-10-19):** Older GBM figures; current RANO paper/poster are separate. Inventory is not publication approval. Patient-fold results supersede random-cell classifier bars; treatment/RL panels are simulations with assumed kill rates. Images below retain historical labels/numbers and must be checked against the corrected ledger before reuse. No images regenerated or edited.

# Paper: Figure Inventory and Plan

## Historical creation status (2026-10-04), inventory checked 2026-10-19
Script: `scripts/make_paper_figures.py`. The generator reads JSON/CSV but also contains handwritten labels/settings; this is not proof every displayed quantity is sourced. Its savefig setting is 300 dpi. Historical visual-review claim is author-reported, not an independent scientific check. Generation command is historical, not current run authorization.
Output in `paper/figures/` (14 PNG; each PNG pHYs chunk reports 11811 pixels/metre, approximately 300 dpi):

| File | Panels | Data source |
|---|---|---|
| `fig1_track_a.png` (+ `fig1a_classification`, `fig1b_spib_saddle`, `fig1c_tcga_survival`) | A zone-classification accuracy (7 methods); B SPIB energy levels + index-1 saddle; C TCGA univariate Cox forest + KM median 432 d | `cgat/gat_metrics.json` (exists at `output/cgat/gat_metrics.json`; random-cell result, superseded for interpretation by `output/cgat_leak_free.json:variants.gat_pca_clean.acc`), `scvi_metrics.json`, `nmf_metrics.json`, `benchmark_comparison.tsv`, `spib_saddle_point_metrics.json`, `survival_stats_summary.json`, `tcga_km_summary.json` |
| `fig2_track_b.png` (+ `fig2a_resistance`, `fig2b_final_mass`, `fig2c_forecast`, `fig2d_elongation`) | A simulated resistant fraction MTD vs adaptive (old 54/61 summary unverified without its exact per-patient field); B final-mass ratio vs rho with 0.02 line; C core-target forecast Dice (all / grew / shrank); D paired elongation | `adaptive_geometry_metrics.json`, `forecast_labels_core/results.json`, `fractal_aniso_vs_iso.json` |
| `fig3_track_c.png` (+ `fig3a_conditioned_policy`, `fig3b_probe_commit`, `fig3c_paced_ttp`) | A win rate by set, blind vs conditioned vs rule; B probe-then-commit vs noise; C paced TTP gain, 4 sets | `rl_kill_conditioned/evaluate.json`, `probe_paced_policies/evaluate.json` |

Required corrections before any reuse (not certification of existing raster text):
- Fig 1A: historical random-cell score with patient/label leakage and test-epoch selection; use patient-fold clean-PCA mean 63.8% (SD 13.3), not the old C-GAT headline. Fig 1B: SPIB found 2 metastable states; the three attractor energies are near zero. Fig 1C: age is a known factor.
- Fig 2B: adaptive ends with more tumour mass in most patients (ratio up to about 4.7), and ratio is 1.0 for the 10 patients above rho 0.02, who all progress earlier (orange).
- Fig 2C: absolute Dice is low; DTI orientation adds nothing (iso is as good); the model loses on shrinking cores. Fig 2D: simulation vs simulation; tract alignment is lower under aniso (0.45 vs 0.50).
- Fig 3A: win is on day-90 volume only, which does not translate to TTP. Fig 3C: orange = real growth inputs, assumed kills (+2.1 simulated days); blue = synthetic parameter sets, also simulations.

Framework schematic made: `paper/figures/fig0_framework.png` (script `scripts/make_framework_figure.py`). Text is hand-written from [[PAPER_FINDINGS_LEDGER]]. Visual inspection on 2026-10-19 confirms the image says "no link" and Track C says "Assumed kill rates". The old 78.7% C-GAT and 54/61 headlines remain in the raster; they are historical, caveated/unsourced claims, not current evidence. The exact per-patient inflammation artifact/provenance is missing here; ledger A24 records the constant 1.0 score. The "Do NOT use" and supplement notes below still apply.

---
(SUPERSEDED older plan below, retained as history; no current generation instructions)

Written 2026-10-04. I listed files and read the JSON that feeds them. **I have not opened any PNG.** Check each image by eye before use.
Original empty-directory statement is superseded: 14 PNG files now exist.

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

## Original missing-figure proposal (superseded)
A1 leaderboard, B1 final-mass panel, B2 forecast, B3 elongation, C2 probe-then-commit, C3 paced TTP, plus an optional framework schematic.

## Intended/reported edits 2026-10-04 (not all independently visually checked)
- Fig 3 (A, C): `real_test` intended label is "real growth, assumed kill".
- Fig 0 (framework): Track A to Track B arrow removed, labelled "no link"; footer states the inflammation score is 1.0 for every patient. Track C data box says assumed kill rates.
- Fig 1 (Track A): C-GAT accuracy bars must carry the script 86 / 88 caveat, or be dropped. Decide after script 88.

Raster-use limitation: other than the framework image, panel labels were not visually rechecked in this correction. The existing warning list remains applicable; source keys and captions must be checked image by image before any paper use.
