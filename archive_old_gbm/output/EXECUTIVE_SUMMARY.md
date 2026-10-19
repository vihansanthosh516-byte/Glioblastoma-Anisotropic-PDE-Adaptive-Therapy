> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# GBM Digital Twin -- Executive Summary & Visual Canvas
## Historical Tracks A, B, C model summary

**Original generated date (unverified):** 2026-07-28; corrected 2026-10-19
**Status:** Historical author report; not evidence of integrated or clinical validation. Current goal: `../../START_HERE.md`.

---

## 1. Core Results at a Glance

| Track | Component | Key Metric | Status |
|-------|-----------|------------|--------|
| **A: MSOS** | Single-cell -> GRN -> Invasion | APOD/S100B/MT3 master switches; CSGT p<0.001 | Historical completion assertion; run log absent |
| **B: PDE Cohort** | Anisotropic FK + Stroma + Adaptive | rho_s Sobol S1=0.607; 13.3% dose sparing; Df=1.04-1.49 | Historical completion assertion; run log absent |
| **C: Digital Twin** | Inverse Est. + MPC + RL | Assumed-model RL experiment; old clearance/clinical rule withdrawn | Historical simulation; no clinical validation |

---

## 2. Visual Canvas

### 2.1 Phase 5: historical unbudgeted RL simulation

Source: `../../output/phase5_adaptive_metrics.json`, keys rl_final_volume_mm3, stupp_final_volume_mm3, training_episodes. It records 40 episodes, toxicity reduction 0% and recurrence delay 0 days; no 10.6x claim or TTP benefit. This run is superseded for fair policy comparison, not interchangeable with later script-66 values (ledger open-conflict item 4). Other training/grid/hyperparameter assertions below lack a linked run manifest.
![Phase 5 Adaptive Steering](../../output/phase5_adaptive_steering.png)

**Key Results:**
- RL Final Volume: **0.2471128 mm^3** vs Stupp: **11.0 mm^3**
- 40 episodes trained on 32^3 grid (dt=0.5), evaluated on 64^3
- REINFORCE + entropy (0.01), lr=1e-2
- Observation: [norm_vol, u_max, day_frac, chemo_tox, rad_tox] -> Action: {Rest, TMZ, RT, Combo}

### 2.2 Timed Drug Infusion: assumed PK/PD simulation

The patient geometry does not make treatment outcomes measured. Exact time-course JSON/configuration is not cited; all voxel, timing and PK/PD figures in this section and its table remain unverified historical assertions. A PNG does not verify them.
![Real Patient Timed Infusion](../../output/real_patient_timed_infusion.png)

**Simulation:** Real BraTS patient (BraTS2021_00000), 120 days, single infusion Day 60
- **Days 0-59:** Tumor grows (C(t)=0, killing term = 0)
- **Day 60:** Bolus infusion -> C(t) spikes -> tumor shrinks
- **Days 61-120:** Drug decays exponentially (k_elim=0.1), residual response

**Volume Time Series:**
![Timed Infusion Time Series](../../output/real_patient_timed_infusion_timeseries.png)

---

## 3. Quantitative Summary

### Phase 5 RL Training (40 Episodes)
| Metric | Value |
|--------|-------|
| Training Episodes | 40 |
| Grid (train) | 32^3 |
| Grid (eval) | 64^3 |
| RL Final Volume | 0.2471128 mm^3 |
| Stupp Final Volume | 11.0059328 mm^3 |
| Drug Sensitivity (alpha) | 0.08 |
| Elimination Rate (k_elim) | 0.1 /day |

### Timed Infusion Simulation (BraTS2021_00000)
| Parameter | Value |
|-----------|-------|
| Initial Volume | 7,168 voxels |
| Final Volume (Day 120) | 7,622 voxels |
| Infusion Day | 60 |
| Peak Concentration | 1.0 |
| Drug alpha | 0.08 |
| k_elim | 0.1 /day |
| Downsample | 2x |

### Track B PDE Cohort (historical assumed-parameter simulations)

The following old aggregates lack exact matching output keys in this document. They are superseded as current results by ledger B1/B3 and its source artifacts: unfavorable fractal direction, Sobol alpha_sens≈0.43/rho≈0.20 at N=128, and distinct 3D mean saving 81.4%±30.4%. These are different model comparisons, not one cohort benefit. Zero thresholded volume is not eradication.
| Finding | Value |
|---------|-------|
| Fractal Dimension (aniso) | Df = 1.04--"1.49 |
| Aniso vs Iso | t=24.74, p<0.001, d=8.75 |
| Stromal Front Correlation | r = 0.938--"0.952 |
| Adaptive Dose Sparing | 13.3+-4.3% |
| rho_s Sobol S1 | 0.607 (dominant) |
| 3D MTD thresholded volume | Historical 0 mm3, exact run unverified; does not establish eradication |
| 3D Adaptive (68% dose) | 40+-8.8 mm^3 |

### Historical Track C module assertions

RMSE and MPC claims below require exact benchmark keys and assumptions; patient utility is not established. D is at a bound in 147/154 fitted cases; rho is at its lower bound in 124/154 (ledger 2bj). Robust benchmark variance fields conflict; no blanket variance-reduction claim follows.
| Tier | Component | Validation |
|------|-----------|------------|
| 1 | Inverse Estimation | RMSE <5% (clean), <15% (10% noise) |
| 2 | Robust MPC | 68.9% dose-sparing, cost variance down |
| 3 | Spatial Metrics | DSC/HD95/MSD implemented |
| 5 | RL Steering | Unsupported clearance claim withdrawn; phase5 output is a simulation |

---

## 4. Historical reported checks

Dated transcripts and exact commit/test scope are missing for these old pass counts. They are unverified, not current test evidence or scientific validation. Current scratch command/results: `../../paper/codex_audit2/fix_progress.md`.

| Test Suite | Tests | Status |
|------------|-------|--------|
| Unit Tests (Inverse, MPC, Spatial) | 38 | [OK] **38/38 PASS** |
| Integration Tests (NIfTI, BraTS, PDE) | 4 | [OK] **4/4 PASS** |
| Phase 5 RL Training | 40 episodes | Historical completion assertion; run log absent |
| Timed PK/PD Infusion | 120 days | Historical completion assertion; run log absent |

---

## 5. Historical command examples

Not verified current run instructions. Exact script paths, external data and configurations require review; no rerun is authorized by this archive. During audit fixes, tests run only in scratch and scientific scripts are not executed.

```bash
# Phase 5 RL Training
python src/58_rl_adaptive_steering.py

# Timed Drug Infusion (120 days, Day 60 infusion)
python src/timed_drug_infusion.py --patient BraTS2021_00000 --days 120 --infusion-days 60 --downsample 2

# Multi-infusion schedule
python src/timed_drug_infusion.py --patient BraTS2021_00000 --days 120 --infusion-days 30 60 90

# Real BraTS Patient PDE (Track B)
python src/42_anisotropic_pde.py --real-patient-dir data/brats/BraTS2021_00000 --use-real-only --grid-size 64

# Inverse Estimation from MRI
python src/51_inverse_parameter_estimation.py --nifti-t0 T0_seg.nii.gz --nifti-t1 T1_seg.nii.gz --delta-t 30

# Full Test Suite
python -m pytest tests/ -v
```

---

## 6. Historical data/output inventory (unverified)

No dated manifest establishes the sizes/cohort counts below. File presence is not a successful validation run.

| Data Source | Location | Size |
|-------------|----------|------|
| BraTS 2021 Task 1 | `data/brats/` | 12.3 GB (~1000 patients) |
| Phase 5 Output | `output/phase5_*` | ~0.5 MB |
| Timed Infusion | `output/real_patient_timed_infusion_*` | ~2.5 MB |
| Phase 5 Metrics | `output/phase5_adaptive_metrics.json` | 199 B |

---

*No complete integrated or clinical validation established. Root `../../output/real_patient_validation_summary.json` reports n=108, mean percent error 150.57196, median 48.14895, only 43/108 within 30%; exact endpoint/cohort limitations still apply. These errors contradict a blanket successful-validation claim.*