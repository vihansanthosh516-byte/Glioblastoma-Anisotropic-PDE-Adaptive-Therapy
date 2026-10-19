> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Script 47 — MPC Optimal Control & Dual-Drug

**Status:** Committed [hash]

## Cohort
- Eight selected fitted MU-Glioma inputs for assumed-treatment simulation; not observed treatment outcomes, rho spread from 0.0001 to 0.11 /day
- Selection: positive rho, R² > 0.7, spread across rho range

## Key findings

### Finding 1: Adding a second assumed drug changes model outcomes
- Single-agent progresses earlier than MTD for **1/8** patients (PatientID_0188: single 25d vs MTD 40d). Corrected 2026-10-04 from `output/dual_drug_comparison.json`; the earlier "4/8" did not match.
- Dual-adaptive reaches the 360-day horizon cap in six of eight modeled cases; PatientID_0188 has TTP263 days, PatientID_0123 22 days. Horizon capping is not demonstrated clinical control.

### Finding 2: Simulated model-defined resistance
- MTD: Rf ≈ 1.00 for all patients
- Single-agent: Rf 0.05–0.99 depending on rho (fails at high rho)
- Dual-agent: Rf ≤ 0.05 for 6/8 patients

### Finding 3: Drug cost trade-off
- Dual-agent drug AUC is 944-1280, against 46-144 for single-agent (7-20x) and 385 for MTD (2.5-3.3x). Corrected 2026-10-04: the earlier "2.5x single-agent" was the ratio to MTD.
- **Dual-agent is not dose-sparing.** It uses more drug than MTD. Do not call it sparing in the paper.
- Added assumed drug, exposure and policy all differ; the fair full-dose dual baseline below does not support an adaptivity benefit.

## Limitations
- PatientID_0123 (rho=0.11/day, V0=3.4 mm³) has modeled TTP7/11/22 days for MTD/single/dual; an edge case, not observed treatment failure
- PatientID_0188 (rho=0.037/day) extends TTP from 40 to 263 days but does not reach 360

## Files
- src/47_optimal_control.py
- output/dual_drug_comparison.json
## Fair-baseline check (script 89, 2026-10-04)
MTD plus the same secondary drug at full dose, n=8: mean horizon-capped simulated TTP 358.9 d, mean modeled resistant fraction 0.019741 (range1.22e-13–0.157920; not a cohort1e-12), AUC 1152.5. Dual-adaptive: 305.6 d, 0.13, AUC 1015.5. So the "dual-agent rescue" comes from the second drug, not from adaptive control. Source: `output/dual_mtd_baseline.json`. See [[PAPER_FINDINGS_LEDGER]] 2i.

Source keys: dual_drug_comparison.json:patients[].{rho_per_day,r_squared,mtd,single_adaptive,dual_adaptive}; source89 patients[].mtd_plus_secondary and from_script47. AUC units follow the assumed model. Fitted V0=3.4005 mm3 need not equal the discretized initial seed; inspect the initialization floor before interpreting case0123. No clinical rescue or non-inferiority is established.
