# MONTH 10 AUDIT LOG

Generated: 2026-09-03T19:14:59

## Environment
- Python interpreter: `C:\Users\vihan\AppData\Local\Python\pythoncore-3.14-64\python.exe`
- Python version: 3.14.4
- numpy 2.4.4
- scipy 1.18.0
- matplotlib 3.11.1

## Metric JSON inputs (read-only)
| File | Size |
|---|---|
| `output\3d_extension_summary.json` | 5,600 bytes |
| `output\adaptive_geometry_metrics.json` | 6,961 bytes |
| `output\anisotropic_geometry_metrics.json` | 5,076 bytes |
| `output\checkpoint.json` | 33,272 bytes |
| `output\checkpoint_atlas.json` | 18,682 bytes |
| `output\checkpoint_dki.json` | 10,975 bytes |
| `output\checkpoint_standard.json` | 11,324 bytes |
| `output\cv_split.json` | 1,541 bytes |
| `output\dual_drug_comparison.json` | 3,970 bytes |
| `output\final_executive_summary.json` | 2,051 bytes |
| `output\hybrid_cohort_metrics.json` | 8,867 bytes |
| `output\inverse_est_PatientID_0003.json` | 197 bytes |
| `output\inverse_est_PatientID_0005.json` | 197 bytes |
| `output\inverse_est_PatientID_0006.json` | 198 bytes |
| `output\inverse_est_PatientID_0007.json` | 228 bytes |
| `output\inverse_est_PatientID_0008.json` | 196 bytes |
| `output\inverse_est_PatientID_0010.json` | 198 bytes |
| `output\isotropic_baseline_metrics.json` | 1,406 bytes |
| `output\master_cohort_summary.json` | 40,875 bytes |
| `output\mu_glioma_cohort.json` | 228,690 bytes |
| `output\mu_glioma_cohort_smoke.json` | 226,776 bytes |
| `output\phase5_adaptive_metrics.json` | 200 bytes |
| `output\real_patient_validation_summary.json` | 459 bytes |
| `output\rho_calibration.json` | 13,023 bytes |
| `output\rho_calibration_positive.json` | 285 bytes |
| `output\sobol_sensitivity_results.json` | 894 bytes |
| `output\stromal_feedback_metrics.json` | 6,260 bytes |

## Generated PNG deliverables
| File | Size | Pixel dims |
|---|---|---|
| `output\65_master_summary_figure.png` | 405,276 bytes | 4770x3543 |
| `output\adaptive_initial_clones.png` | 151,157 bytes | 3390x1735 |
| `output\adaptive_tensor_validation.png` | 221,678 bytes | 2534x2166 |
| `output\adaptive_therapy_comparison.png` | 955,223 bytes | 3837x3201 |
| `output\adaptive_therapy_dynamics.png` | 451,356 bytes | 3175x1185 |
| `output\adaptive_therapy_metrics_summary.png` | 320,800 bytes | 2971x1773 |
| `output\anisotropic_geometry_summary.png` | 227,731 bytes | 2580x1773 |
| `output\anisotropic_recurrence_maps.png` | 313,378 bytes | 3596x1735 |
| `output\anisotropic_solver_mass_test.png` | 95,454 bytes | 2579x990 |
| `output\anisotropic_tensor_validation.png` | 249,463 bytes | 2514x2166 |
| `output\bbb_permeability_map.png` | 59,622 bytes | 1162x995 |
| `output\fig1_bar_chart.png` | 122,244 bytes | 2970x1770 |
| `output\fig2_scatter.png` | 344,712 bytes | 2887x2369 |
| `output\fig3_delta_boxplot.png` | 144,034 bytes | 2969x1768 |
| `output\master_cohort_synthesis.png` | 760,412 bytes | 4081x3644 |
| `output\phase5_adaptive_steering.png` | 277,924 bytes | 3178x2362 |
| `output\real_patient_evolution.png` | 20,512 bytes | 3000x1000 |
| `output\real_patient_timed_infusion.png` | 3,534,595 bytes | 7539x1552 |
| `output\real_patient_timed_infusion_timeseries.png` | 274,294 bytes | 2969x2368 |
| `output\sobol_tornado_plot.png` | 52,161 bytes | 1737x967 |
| `output\stromal_feedback_init.png` | 30,132 bytes | 2015x990 |
| `output\stromal_feedback_metrics_summary.png` | 310,718 bytes | 2973x1773 |
| `output\stromal_feedback_recurrence_maps.png` | 300,203 bytes | 4282x2221 |
| `output\stromal_tensor_validation.png` | 460,901 bytes | 3940x2389 |

## Validation summary
- Overall validation: **PASS**
- Phase 1 (anisotropic) Df bounds (1.0, 2.0): 8/8 pass
- Phase 2 (stromal) Df bounds (0.4, 1.0), front_corr floor 0.9: 8/8 pass
  (realized front_corr range 0.9355-0.9523)
- Phase 3 (adaptive) drug_reduction in (0,1), TTP>0: 8/8 pass

## Spherical baseline baseline (D3)
- Isotropic baseline cache: `output\isotropic_baseline_metrics.json` (present)

## Spatial-render fallback usage (Panel A)
- anisotropic fallback used: False
- stromal fallback used: False

## `__pycache__` cleanup
- `2290` `__pycache__` directories removed under src/ and venv/ (D5 idempotency — bytecode regenerates on next run).

## Notes
- Bash `run_all.sh` requires Git Bash / WSL on Windows. PowerShell equivalent can be substituted.
- `output/*.npz` are git-ignored per D6 (heavy per-patient arrays not pushed; JSON + PNG evidence trail is tracked).
- Untagged `.npz` arrays remain on disk for reproducibility but are excluded from git tracking.
