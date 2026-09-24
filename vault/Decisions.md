# Decisions Log

## 2026-09-23 — Adaptive therapy high-rho stratification

**Finding:** 23 of 61 patients (rho > 0.02/day) show adaptive ≈ MTD. 
10/23 progress earlier under adaptive.

**Hypothesis:** The 80% controller setpoint is too tight for fast-growing 
tumors — the tumor rebounds from 80% to 100% in days, forcing near-continuous 
dosing, so adaptive degenerates to MTD.

**Test:** Rerun with THRESHOLD_OFF=0.50. If high-rho patients now show real 
holidays and preserved sensitive clone, artifact. If not, real biology.

**Decision criteria:**
- Artifact → rerun primary with 0.50, update committed result
- Real biology → report stratified finding: "adaptive works for rho < 0.02"

## 2026-09-22 — E_MAX_RATIO = 1000

**Problem:** At 200x, MTD shrinks 0-7% (too weak). At 2000x, saturates 
(too strong).

**Test:** Sweep 150x, 500x, 1000x, 2000x on 3 patients.
- 150x → 0-7% shrinkage
- 500x → 25-58%
- 1000x → 40-62% ← matches clinical
- 2000x → 57-62% (saturated)

**Decision:** E_MAX_RATIO = 1000. Same value scales to all patients 
via per-patient rho.

## 2026-09-22 — Per-patient E_max from rho

**Problem:** Global E_max produced 30x spread in effective drug:growth ratio 
across patients — same bug class as the rho scalar.

**Fix:** E_max_patient = patient_rho × E_MAX_RATIO. Every patient gets the 
same nominal drug:growth ratio.