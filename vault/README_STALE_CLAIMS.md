# README.md — claims that disagree with files on disk

Checked 2026-10-04. **Do not cite the README in the paper.** Cite the vault notes and output files. The README also still has merge-conflict markers in places.

| README claim | What the file says | Source |
|---|---|---|
| Classical ML reaches ">94%" zone accuracy | LR 69.8%, RF 72.6%, C-GAT 78.7% | `output/benchmark_comparison.tsv`, `output/cgat/gat_metrics.json` |
| "Track A outputs are currently unverifiable" | Outputs exist and were verified (see [[TRACK_A_RESULTS]]) | vault audit |
| GRN has 373 edges; S100B out-degree 45 | 320 edges; S100B 42 | `output/grn_metrics.json` |
| Best single KO is SDE2; best pair S100A11+ZNF106 | TMLHE (collapse 0.0408) by collapse; MT-CO2 (TI +4.82) by therapeutic index; best pair S100A11+S100A6 | `output/single_ko_results.json`, [[TRACK_A_RESULTS]] |
| Anisotropic D_f 1.20-1.55 vs about 0 (d = 10.43) | D_f is lower in aniso (0.925 vs 0.943) and is not a valid metric. Elongation is the real effect. | [[Script-74-Fractal-Elongation]] |
| Sobol: rho_s S1 = 0.607 or 0.998 | PDE Sobol: alpha_sens 0.43, rho_s 0.20, D_w 0.000 | [[Script-80-Sobol-PDE]] |
| Atlas-DTI aniso wins about 35/62 patients | 12/62 | `output/results_atlas.csv`, [[UCSF-PDGM-Tensor-Study]] |
| RL 13.94 mm3 vs Stupp 11.01 mm3 at day 90 | Script 66: PPO 12.67 vs Stupp 14.47 mm3 (different run and setup). Do not mix. | [[Script-66-RL-Equal-Budget]] |
| "Tracks B and C are evaluated on synthetic cohorts" | Track B uses real MU-Glioma-Post parameters. Track C uses a mix (see [[PAPER_DATA_SOURCES]]). | vault |
| Robust MPC: "30% lower cost variance" | Variance reduction is **-31.8% (worse)**; 68.9% dose sparing is the supported number | [[TRACK_C_RESULTS]] |
| Biomarker rho > 0.024 /day (CI 0.0202-0.0249) | Script 62 gives rho* about 0.075 /day, CI [0.0746, 0.0820]. Script 62 is the current analysis (64 real growing patients, equal budget; early start wins for 60/64; only 4 patients lie above rho*). The README's 0.024 comes from an older virtual-cohort analysis. Use script 62. | [[TRACK_C_RESULTS]] |
| 3D adaptive dose sparing 59.8% (one patient) | 8 patients: 81.4% +- 30.4%, collapses from 89-99% to 10.9% as rho rises | [[TRACK_B_RESULTS]] |
| "8-patient synthetic cohort" for Track B | 61 real patients for adaptive therapy; 8 real patients for MPC and 3D | [[TRACK_B_RESULTS]] |
| Adaptive: TTP "non-inferior", dose sparing | Also: final tumour mass is higher under adaptive (315.9 vs 306.4; VR 2.39), and 10/61 progressed earlier | `output/adaptive_cohort_summary.json` |
| Track A "C-GAT 78.7%" as a clean accuracy | Random cell-level split. A patient-ID lookup gets 82.8%; patient-level LR/RF/scVI-LR get 60-66% | `output/patient_level_cv.json`, [[LIMITS_CURRENT]] L0 |
| Track A feeds Track B through an inflammation score | Score is 1.0 for all 61 patients; there is no link | [[PAPER_FINDINGS_LEDGER]] |
