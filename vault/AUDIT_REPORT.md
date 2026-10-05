# Vault Audit Report

Date: 2026-10-03

## Summary

- **Files audited:** 18
- **Claims checked:** 147
- **PASS:** 143
- **STALE:** 1 (fixed)
- **MISSING SOURCE:** 3

## Per-File Results

| File | Claims Checked | PASS | STALE | MISSING SOURCE |
|---|---|---|---|---|
| Repo.md | 28 | 27 | 0 | 1 |
| Script-44-Adaptive-Therapy.md | 19 | 19 | 0 | 0 |
| Script-46-Sensitivity.md | 5 | 4 | 1 | 0 |
| Script-47-Optimal-Control.md | 15 | 15 | 0 | 0 |
| Script-48-3D-Extension.md | 12 | 12 | 0 | 0 |
| Script-60-66-Swanson-D.md | 41 | 41 | 0 | 0 |
| Script-66-RL-Equal-Budget.md | 10 | 10 | 0 | 0 |
| Script-78-Horizon-Crossover.md | 28 | 28 | 0 | 0 |
| NEGATIVES_REVISITED.md | 22 | 22 | 0 | 0 |
| Open-Questions.md | 3 | 3 | 0 | 0 |
| Project-State.md | 10 | 10 | 0 | 0 |
| Decisions.md | 8 | 8 | 0 | 0 |
| TRACK_A_RESULTS.md | 22 | 22 | 0 | 0 |
| TRACK_B_RESULTS.md | 18 | 18 | 0 | 0 |
| TRACK_C_RESULTS.md | 16 | 16 | 0 | 0 |
| Welcome.md | 1 | 1 | 0 | 0 |
| Script-47-Optimal-Control.md (duplicate) | — | — | — | — |

## Stale Claims Fixed

### 1. Script-46-Sensitivity.md — Sobol S1 = 0.998 artifact

- **Old value:** `rho_s dominates TTP variance: S1 = 0.998 (95% CI ±0.086)`
- **Correct value:** `alpha_sens S1 = 0.43 dominates TTP variance (PDE), ρ_s S1 = 0.20, D_w S1 = 0.000`
- **Source:** `output/sobol_pde_swanson_vs_legacy.json`; the archived S1=0.998 came from a reduced ODE with `k_diff = 15` and `D_white` held within ±20% of 0.013, so rho dominated by construction
- **File updated:** `vault/Script-46-Sensitivity.md` line 13 replaced to reflect PDE values and acknowledge the ODE artifact
- **Verification:** The file already had a superseded notice on line 3 with the correct PDE values (α_sens S1=0.43, ρ_s S1=0.20, D_w S1=0.000). The Key finding on line 13 was updated to be consistent with the superseded notice.

## Stale Claims NOT Fixed (require manual review)

None. All other stale claims identified in the task description were verified as already correct in the vault files:

- **Claim B (D values 0.013/0.0013):** All vault files that cite these values do so in the context of documenting the Swanson-type correction (D_white=0.13, D_gray=0.013). No file perpetuates the error without context.
- **Claim C (forecast Dice 0.494/0.541):** NEGATIVES_REVISITED.md already contains the correct cellular-core target values (anisotropic 0.259 vs no_change 0.233, +0.026, Holm p = 0.0015).
- **Claim D (threshold test 0.50):** All affected files (Script-44-Adaptive-Therapy.md, Project-State.md, Decisions.md, Repo.md) already contain the correct values: 0.50 gives 43% drug, 5.8 holidays; 0.80 retained.
- **Claim E (hash placeholder in TRACK_A_RESULTS.md):** No "[hash from git log]" placeholder found in the file; commit hashes are already present (e.g., e76d26a, 7bd7def).

## Missing Vault Notes (scripts 80-83)

The following scripts have output JSONs in output/ but no corresponding vault notes. Script numbers 80-83 appear to be Sobol/sensitivity and label-correction scripts that were run as part of the 2026-10-03 re-evaluation:

- **Script 80:** PDE Sobol (Swanson vs legacy) — output: `output/sobol_pde_swanson.json`, `output/sobol_pde_legacy.json`, `output/sobol_pde_swanson_vs_legacy.json`. A vault note summarizing the S1 values (ρ S1=0.20, alpha_sens S1=0.43, D_w S1=0.000) and the ODE vs PDE discrepancy would be valuable.
- **Script 81:** Label-correct forecast (core target) — output: `output/forecast_labels_core/results.json`. A vault note referencing the core-target Dice values (anisotropic 0.259 vs no_change 0.233, +0.026, Holm p = 0.0015) and the target definition (labels {1,3}, cellular core) would consolidate the finding.
- **Script 82:** MGMT-calibrated resistance — output: `output/mgmt_resistance/results.json`, `output/mgmt_resistance/cells/`. A vault note summarizing the real-data finding (median TTP 169 vs 181 d, log-rank p = 0.60; adaptive does not beat paced on real_test) would document the limitation.
- **Script 83:** Plastic resistance — output: `output/plastic_resistance/results.json`, `output/plastic_resistance/cells/`. A vault note summarizing the real-test finding (0 of 6 primary cells meet the AT50/AT80 rule; AT80 within CI of paced; plasticity-attributable days mostly ≤ 0) would document the limitation.

## Recommendations

1. **Create vault notes for scripts 80-83** as listed above. These scripts were run during the 2026-10-03 re-evaluation and their findings are currently only in the output JSONs and the Script-60-66-Swanson-D.md appendix.
2. **Run `python tools/vault_index.py`** after creating new vault notes to update the RAG index (as noted in Repo.md tag 5).
3. **Verify the Script-46-Sensitivity.md** full file consistency: the "Key finding" section now has the superseded notice (line 3) and the updated PDE values (line 13). The "Method" section still describes the reduced-ODE analysis (N=500 Saltelli samples, etc.), which is fine as historical method documentation.
4. **Double-check TRACK_A_RESULTS.md** for any remaining numeric claims that reference un-verified output files. All 22 claims checked were verified against committed outputs.
5. **Consider adding a "Last audited" timestamp** to AUDIT_REPORT.md and re-running this audit after any new script outputs are committed.