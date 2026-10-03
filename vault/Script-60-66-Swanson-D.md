# Scripts 60/66 re-run at Swanson-type diffusivities: PRE-SPECIFICATION

**Written 2026-10-03, before any re-run output was read.** Results are appended below the line `## RESULTS
All values read back from the output files (listed at the end). Runs: 2026-10-03.

### Headline
1. **Volume endpoints stay flat; spatial endpoints are real.** At Swanson-range D (0.1-0.8 mm²/d) the day-90 total-volume ablation is still small (PPO: no-DTI +0.27% vs +0.007% before; at most 3.6% in any scenario), but the day-90 mask changes materially: Dice(no-DTI vs full) = 0.886 mean (was 0.994), min 0.08. The old null was a parameter artifact for spatial endpoints, not for volume (mass conservation).
2. **The 0.998 does not survive a PDE.** In a real PDE Sobol, rho is not dominant for V_total, and D_w contributes nothing to volume at either D range.
3. **Script 66's study is essentially unchanged**: win rates and mean log-ratios at equal budget barely move.

### Script 60 (primary endpoint: day-90 mean final volume, mm³; 30 LHS scenarios, equal budget)
| Arm | legacy D, full | Swanson D, full | legacy no-DTI | Swanson no-DTI |
|---|---|---|---|---|
| Stupp | 26.452 | 26.659 | 26.455 | 26.810 |
| heuristic_budgeted | 25.388 | 25.486 | 25.390 | 25.552 |
| ppo | 24.085 | 24.182 | 24.087 | 24.247 |
| dagger_oracle | 21.682 | 23.976 | 21.691 | **26.131** |

- Ablation impact on PPO mean final volume: no DTI +0.007% -> **+0.268%**; no mechanics 0.000% (mechanics never enters `pde_step`); pure RD = no DTI.
- Win rate vs Stupp: ppo 100% and heuristic 43% in every ablation. `dagger_oracle`: 100% (full) -> **83% (no DTI)**; its mean volume rises 9% (23.976 -> 26.131) without DTI. The late-block oracle schedule is the arm whose result depends on where the cells are.
- Caveat: the "No DTI" arm is script 60's uniform x-aligned tensor (D_xx = D_white, D_yy = D_zz = D_gray), not a matched isotropic tensor, so it also changes D magnitude. Script 78's `iso_matched` arm isolates orientation (Dice 0.923 at D = 0.1, 90 d, natural growth).
- The script's `best_arm_full_model` field changes from `dagger_oracle` to `ppo`. By mean final volume `dagger_oracle` is still lower (23.976 vs 24.182), but the field uses the script's own criterion. The BC policy was trained without diffusion and is no longer the optimum at the new D (script 66 LHS-30 log ratio -0.103 vs oracle_final -0.198).

### Script 60-style spatial endpoints at Swanson D (script 77 re-run, `ablation_spatial_endpoints_swanson.json`)
| | legacy D | Swanson D |
|---|---|---|
| diffusion length sqrt(2Dt), 90 d | 0.46-1.18 mm | **4.6-11.8 mm** |
| Fisher front 2 sqrt(D rho) t | 0.66-2.85 mm | **6.6-28.5 mm** |
| Stupp: mask Dice no-DTI vs full (mean / min) | 0.994 / 0.904 | **0.886 / 0.080** |
| Stupp: mask Dice no-diffusion vs full (mean) | 0.991 | **0.719** |
| Stupp extent p90 (full / no-DTI / no-diff), mm | 2.28 / 2.27 / 2.30 | **1.79 / 1.46 / 2.30** |
| combo days 56-90: Dice no-DTI vs full (mean) | 0.996 | **0.866** |
| oracle beats Stupp on V_total | 100% | 100% |
| oracle beats Stupp on V_vis(0.01 K) | 43% | **20%** |
The diffusion length is now several voxels (2 mm), so the sub-voxel explanation from script 77 no longer applies.

### Script 66 (primary: equal-budget win rate and mean log(final/Stupp final); evaluate stage, full PDE)
| LHS-30 | win % old -> new | log ratio old -> new |
|---|---|---|
| heuristic59 | 10 -> 13 | +0.071 -> +0.069 |
| combo35_late = oracle_final | 100 -> 100 | -0.210 -> -0.198 |
| ppo_final_s0 | 100 -> 100 | -0.139 -> -0.141 |
| bc_dagger_final | 100 -> 100 | -0.210 -> **-0.103** |
Real-test (21 patients): all numbers identical to 3 decimals except trivial changes; real-test patients use their own estimated `D_mm2_per_day`, unchanged by this re-run.
- Max effect of diffusion on Stupp's final volume (validate stage): 0.029% -> **0.83%**.
- Reused unchanged from the legacy directory: `oracle.json`, `train.json`, PPO/BC weights (D-independent, reaction-only stages).
- Conclusion of script 66 holds: pacing/late-block structure at equal budget is not a D artifact.

### Sobol (script 80, PDE, N = 128 base, 640 runs per regime, Stupp schedule, day 90)
| Output | Parameter | legacy S1 | legacy ST | Swanson S1 | Swanson ST |
|---|---|---|---|---|---|
| **V_total (primary)** | rho | 0.203 ± 0.178 | 0.392 | 0.203 ± 0.179 | 0.393 |
| | D_w | 0.000 | 0.000 | 0.000 | 0.000 |
| | alpha_sens | 0.429 ± 0.140 | 0.760 | 0.426 ± 0.140 | 0.761 |
| extent_p90 | rho | 0.193 ± 0.136 | 0.361 | 0.217 ± 0.159 | 0.453 |
| | D_w | 0.000 | 0.000 | 0.004 ± 0.029 | **0.020 ± 0.012** |
| | alpha_sens | 0.617 ± 0.174 | 0.821 | 0.514 ± 0.190 | 0.782 |
| V_vis(0.04) | all | CIs 0.05-2.1 (many zeros): not interpretable at N = 128 | | | |

- **Does rho_s still dominate?** No. In the PDE with treatment, alpha_sens (kill scale) has the largest S1 for V_total (0.43), rho second (0.20), and ST sums above 1 (interactions). The archived script-46 value (rho_s S1 = 0.998) comes from a toy ODE with k_diff = 15 and D_white held within ±20% of 0.013, so rho dominates by construction. **The Track B headline "rho_s S1 = 0.998" should not be cited as a property of the PDE model.**
- D_w is irrelevant to volume at both ranges (mass conservation). At Swanson D it enters extent only weakly (ST = 0.020 ± 0.012).
- Sobol differences from the earlier "ablation" result: Sobol varies D within 0.1-0.8 on a treated tumour (small by day 90, mean extent ~2 mm), whereas script 78 grew an untreated tumour; so the Sobol D-sensitivity is small while the mask Dice in the ablation is not.

### Verdict
- The ablation null is a parameter artifact for spatial endpoints (Dice 0.99 -> 0.89 with the diffusivity fixed) and still a true null for total volume. The scripts 60/66 conclusions about RL vs Stupp at equal budget hold.
- Caveats: N = 128 Sobol has wide CIs; the unit-corrected range (0.1-0.8) was chosen by converting the script's declared cm²/day, not fitted to data; real-test D estimates were not audited; 30 scenarios, one seed.

### Files
`output/ablation_and_baselines_metrics_swanson.json`, `output/rl_equal_budget_swanson/{validate,theory,evaluate}.json` (+ png, logs), `output/ablation_spatial_endpoints_swanson.json`, `output/sobol_pde_{legacy,swanson}.json`, `output/sobol_pde_swanson_vs_legacy.json`. Legacy outputs are untouched; `GBM_D_REGIME=legacy` reproduces them.
