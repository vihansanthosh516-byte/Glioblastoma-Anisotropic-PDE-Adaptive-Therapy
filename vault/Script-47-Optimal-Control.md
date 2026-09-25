# Script 47 — MPC Optimal Control & Dual-Drug

**Status:** Committed [hash]

## Cohort
- 8 real MU-Glioma patients, rho spread from 0.0001 to 0.11 /day
- Selection: positive rho, R² > 0.7, spread across rho range

## Key findings

### Finding 1: Dual-agent MPC rescues single-agent failures
- Single-agent progresses earlier than MTD for 4/8 patients (e.g., PAT_0188: single 25d vs MTD 40d)
- Dual-agent brings all 8 patients to 360d control except the two extreme outliers

### Finding 2: Resistance suppression
- MTD: Rf ≈ 1.00 for all patients
- Single-agent: Rf 0.05–0.99 depending on rho (fails at high rho)
- Dual-agent: Rf ≤ 0.05 for 6/8 patients

### Finding 3: Drug cost trade-off
- Dual-agent AUC ~2.5× single-agent (~1000 vs ~50)
- Second drug pays for itself in TTP extension

## Limitations
- PatientID_0123 (rho=0.11/day, V0=3.4 mm³) fails on all arms — extreme outlier
- PatientID_0188 (rho=0.037/day) extends TTP from 40 to 263 days but does not reach 360

## Files
- src/47_optimal_control.py
- output/dual_drug_comparison.json