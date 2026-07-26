Here is your completely transformed `README.md`.

It applies the **33 Style Patterns** across the board: it strips out AI buzzwords (*interlocking*, *delve*, *testament*, *game-changer*), cuts artificial drama (*cathedral*, *wake with a hum*), flattens forced triplets and overused em dashes, converts headers to sentence case, removes decorative emojis, and uses clean, direct technical prose.

---

```markdown
# Biophysical Glioblastoma Digital Twin & Reinforcement Learning Adaptive Therapy Framework

[![Python](https://img.shields.io/badge/Python-3.14+-blue.svg?logo=python)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()
[![Status](https://img.shields.io/badge/Status-Research%20Prototype-orange.svg)]()
[![Pipeline](https://img.shields.io/badge/Pipeline-10--Month%20%7C%20MSOS%20%7C%20Digital%20Twin-purple.svg)]()

> **Research prototype.** Not clinically validated. Not for patient care.
> See [Disclaimer](#disclaimer).

---

## Table of contents

- [Overview](#overview)
- [Three parallel research tracks](#three-parallel-research-tracks)
- [Track A: MSOS — Multi-Scale Spatial Oncology Suite](#track-a-msos--multi-scale-spatial-oncology-suite)
- [Track B: 10-Month PDE cohort](#track-b-10-month-pde-cohort)
- [Track C: Digital Twin reactor](#track-c-digital-twin-reactor)
- [Installation & quickstart](#installation--quickstart)
- [Repository structure](#repository-structure)
- [Key findings across tracks](#key-findings-across-tracks)
- [Validation & reproducibility](#validation--reproducibility)
- [Limitations](#limitations)
- [Future work](#future-work)
- [Citation](#citation)
- [Disclaimer](#disclaimer)

---

## Overview

Glioblastoma (GBM) is an aggressive primary brain cancer, with a median survival of approximately 15 months under standard care. Current treatment follows a fixed schedule (Stupp protocol) that does not adjust to patient-specific tumor changes or dynamic treatment responses.

This repository provides a mechanistic computational platform containing three standalone research tracks. Together, they cover single-cell multi-omics, spatial biophysics, virtual clinical trials, and reinforcement learning for adaptive therapy optimization. The framework is designed for hypothesis generation and in-silico trials rather than clinical deployment.

| Track | Scope | Scripts | Execution |
|-------|-------|---------|-----------|
| **A: MSOS** | Single-cell omics, causal GRN, invasion PDE, drug discovery, and clinical validation | 01–41, 50a, 52a, 53b, 54, 55a, 55b, 56a | `src/run_pipeline.py` (partial) |
| **B: 10-Month PDE Cohort** | Anisotropic PDE, stromal coupling, adaptive therapy, Sobol sensitivity, 3D extension, and synthesis | 42–48 | `bash run_all.sh` |
| **C: Digital Twin Reactor** | Inverse estimation, robust MPC, 3D DTI, RL adaptive steering, global SA, virtual cohort, and reporting | 49, 50b, 51, 52b, 53a, 56b, 57–65 | Individual scripts |

---

## Three parallel research tracks

### Track A: MSOS — Multi-Scale Spatial Oncology Suite

**Scope:** Processes single-cell data from the UCSC `multiomic-gbm` dataset (223k cells, 25 patients) through spatial transcriptomics, causal gene regulatory networks, field reconstruction, invasion dynamics, and virtual drug screening, ending with validation on Ivy GAP and TCGA cohorts.

| Month | Scripts | Deliverable |
|-------|---------|-------------|
| **1–3** | 01–09 | Multi-omic ingest, UMAP/DE, baseline models, and unified benchmarking (LR/RF >94% acc) |
| **4** | 10–11 | Contrastive VAE (32-dim latent space) scaled to 140k cells |
| **4–6** | 12–18 | C-GAT (graph attention on cVAE latents, 2.1M edges) with scVI/NMF baselines and CSGT proof (continuous Core → Periphery → Healthy gradient, KW p<0.001) |
| **MSOS M1** | 19–22 | Waddington landscape (Healthy=0.56, Core=0.00; Periphery saddle E=0.865), NEB saddle verification, and drift-diffusion tensors |
| **MSOS M2** | 23–26 | Transfer entropy (100×100 matrix, 373 edges at 95th percentile), PID analysis, causal GRN (master switches: APOD 46, S100B 45, MT3 40), and bootstrap validation |
| **MSOS M3** | 27–30 | ABA lattice (512², wave speed 2.42 µm/hr), FK PDE (ETDRK4 + Strang splitting, analytical 20 µm/hr), and invasion kinetics |
| **MSOS M4** | 31–34 | Virtual single KO (top: SDE2 C=0.0198), dual KO (S100A11+ZNF106 C=0.0143), calibrated TI, and drug gating report |
| **Clinical** | 35–41 | Ivy GAP cohort (120 patients across 3 zones), penalized Elastic Net Cox per zone, FK-PDE spatial recurrence mapping, Hill+Bliss dose optimization, and Clinical Gating Matrix |

**Key outputs:** `output/01_filtered_three_class.h5ad`, `output/cgat/` (C-GAT models + full graph), `output/csgt_*`, `output/te_matrix.npy`, `output/master_switches.tsv`, `output/aba_*`, `output/single_ko_results.json`, `output/dual_ko_ti.json`, `output/real_cohort_*.csv`, `output/penalized_survival_metrics.json`, `output/spatial_recurrence_*`, `output/final_dose_response_matrix.csv`, `output/clinical_actionability_report.md`

---

### Track B: 10-Month PDE cohort

**Scope:** An 8-patient synthetic cohort pipeline executed sequentially via `run_all.sh`. This track adds anisotropic diffusion tensors, tumor–stroma coupling, adaptive therapy with drug holidays, and global sensitivity analysis.

```text
run_all.sh chain:
  42_anisotropic_pde.py           (Month 7)  → anisotropic_geometry_metrics.json
  43_stromal_feedback.py          (Month 8)  → stromal_feedback_metrics.json
  44_adaptive_therapy.py          (Month 9)  → adaptive_geometry_metrics.json
  46_sensitivity_analysis.py       (Phase 2b) → sobol_sensitivity_results.json + tornado plot
  47_optimal_control.py             (Phase 3)  → dual_drug_comparison.json
  48_3d_extension.py               (Phase 3D) → 3d_tumor_volume_patient.npz + 3d_extension_summary.json
  45_validation_synthesis.py       (Month 10) → master_cohort_summary.json + master_cohort_synthesis.png + POSTER_KEY_FINDINGS.md + MONTH10_AUDIT.md

```

| Script | Month/Phase | Core PDE / Algorithm | Result |
| --- | --- | --- | --- |
| **42** | Month 7 | 2D FK with anisotropic tensor (D_∥/D_⊥=10×), gene-driven ρ/D scaling, and box-counting fractal dimension | **Df = 1.04–1.49** (anisotropic) vs ~0 (isotropic); t=24.74, p<0.001, d=8.75 |
| **43** | Month 8 | Coupled tumor-stroma PDE (Michaelis-Menten ρ(G)), D_G=0.13 mm²/day, 200 days | **Front correlation r = 0.938–0.952** across cohort |
| **44** | Month 9 | Dual-clone PDE + TMZ PK, comparing MTD to adaptive dosing (holidays triggered at <80% peak) | **Dose reduction of 9–21% (mean 13.3±4.3%)** with non-inferior TTP ratio (0.50–0.82) |
| **46** | Phase 2b | Reduced ODE + SALib Sobol (5 parameters: ρ_s, aniso_ratio, μ_r, EC50, D_white) | **ρ_s dominates TTP variance** (S1=0.607, ST=0.633) |
| **47** | Phase 3 | 3-arm MPC (MTD / Single-agent / Dual-agent) | **Dual drug control suppresses resistance** (R-frac 0.038 vs 0.99 MTD), TTP=360d for all 8 cases |
| **48** | Phase 3D | 3D FK (50³, 10× anisotropy), MTD vs adaptive over 180 days | **MTD clears local volume**; **Adaptive yields 40±8.8 mm³ at 68% lower dose** |
| **45** | Month 10 | Master synthesis, isotropic baseline, and spatial metrics (DSC/HD/MSD) | **DSC 0.21±0.02, HD 26.3±3.3 mm** (anisotropic vs isotropic structural comparison) |

**Notes on interpretation:** Adaptive therapy achieves equivalent time-to-progression at lower cumulative drug exposure. It does not extend total TTP or reverse selection pressure in this high-selection regime. The main advantage is dynamic dose-sparing with equivalent tumor control, and individual responses scale with background inflammatory markers (Pearson r=0.89, p=0.003).

**Key outputs:** `output/anisotropic_*`, `output/stromal_*`, `output/adaptive_*`, `output/sobol_*`, `output/dual_drug_comparison.json`, `output/3d_*`, `output/master_cohort_summary.json`, `output/master_cohort_synthesis.png`, `output/POSTER_KEY_FINDINGS.md`, `output/MONTH10_AUDIT.md`

---

### Track C: Digital Twin reactor

**Scope:** A 3D model that ingests DTI tensor fields, accounts for poroelastic mechanics, simulates the full Stupp protocol (surgery, radiotherapy, and temozolomide), trains a Gymnasium RL agent for dynamic dose scheduling, and evaluates outcomes across virtual patient cohorts.

| Phase | Scripts | Core Functionality |
| --- | --- | --- |
| **1** | 49, 50b, 51 | 3D visual dashboard and Clinical DSS (inverse L-BFGS-B estimation from 2 timepoints, RMSE <5% noise-free) |
| **2** | 52b | **Robust MPC** (adaptive horizon over 7–21 days; yields 68.9% dose reduction in 50 MC runs) |
| **3** | 53a, 56b | **Spatial metrics** (DSC/HD95/MSD clinical thresholds) and a stable 3D DTI solver (CFL-bounded with dynamic thresholds) |
| **4** | 57 | 90-day virtual Stupp schedule (surgery day 15, RT days 20–50, TMZ days 20–80) |
| **5** | 58 | **Gymnasium RL** (tracks volume, u_max, day, chemo/rad tox; REINFORCE policy evaluated at 64³) |
| **6** | 59 | **Global SA** (LHS N=30); **Biomarker rule: ρ > 0.024 → RL preferred** (win rate >80% for high ρ) |
| **7** | 60–63 | Standard baselines (Stupp, threshold, RL), physical ablations, 5-seed stability tests, and reward sensitivity grids |
| **8** | 64 | **N=20 virtual cohort evaluation** (paired t-test, Wilcoxon signed-rank, and Kaplan-Meier progression curves) |
| **9** | 65 | Summary report and 4-panel dashboard synthesis figure |

**Reinforcement learning results (Phase 5, 64³ evaluation):**

| Metric | RL Adaptive | Standard Stupp | Difference |
| --- | --- | --- | --- |
| **Final Tumor Volume (Day 90)** | **1.04 mm³** | 11.01 mm³ | ~10x lower |
| **Peak Cellularity (u_max)** | 0.02 | 0.15 | 7.5x lower |
| **Time-to-Progression** | >90 days | 42 days | >2x delay |
| **Cumulative Drug Exposure** | 87% | 100% | 13% dose reduction |
| **Cohort Variance (CV)** | Low (CV=0.3) | High (CV=0.8) | More consistent response |

**Biomarker rule (Phases 6–7):**

* **RL Adaptive is preferred when ρ > 0.024 day⁻¹** (faster proliferation rate)
* Standard Stupp performs adequately when ρ ≤ 0.024 day⁻¹ (slower growth)
* RL overall win rate: 36.7% across full range, rising to **>80% for ρ > 0.025**

**Key outputs:** `output/3d_interactive_tumor_dashboard.html`, `output/clinical_reports/`, `output/robust_mpc_benchmark.json`, `output/phase3_3d_dti_metrics.json`, `output/phase4_therapy_metrics.json`, `output/phase5_adaptive_metrics.json`, `output/phase6_sensitivity_metrics.json`, `output/ablation_and_baselines_metrics.json`, `output/rl_convergence_metrics.json`, `output/biomarker_stability_metrics.json`, `output/reward_sensitivity_metrics.json`, `output/phase8_cohort_metrics.json`, `output/final_executive_summary.json`, `output/65_master_summary_figure.png`

---

## Installation & quickstart

```bash
# Clone repository
git clone [https://github.com/vihansanthosh516-byte/Glioblastoma-Anisotropic-PDE-Adaptive-Therapy.git](https://github.com/vihansanthosh516-byte/Glioblastoma-Anisotropic-PDE-Adaptive-Therapy.git)
cd Glioblastoma-Anisotropic-PDE-Adaptive-Therapy

# Create virtual environment (Python 3.14+)
python -m venv venv
source venv/bin/activate          # Linux / macOS
.\venv\Scripts\Activate.ps1       # Windows PowerShell

# Install dependencies
pip install -r Requirements.txt   # Core pinned: SALib>=1.4.7
pip install numpy scipy matplotlib pillow pandas gymnasium torch

# Track B: Run 10-Month PDE Cohort (Months 7 → 10)
bash run_all.sh                   # Linux / macOS / Git Bash / WSL
bash run_all.sh --month10         # Month 10 synthesis only

# Track C: Run Digital Twin phases individually
python src/51_inverse_parameter_estimation.py --test
python src/52_robust_mpc_controller.py --benchmark --n-mc 50
python src/53_spatial_metrics.py --validate
python src/58_rl_adaptive_steering.py
python src/59_sensitivity_analysis.py
python src/60_baselines_and_ablation.py
python src/61_rl_convergence_diagnostics.py
python src/62_biomarker_bootstrap_stability.py
python src/63_reward_sensitivity.py
python src/64_virtual_cohort_simulation.py
python src/65_generate_final_report.py

```

**Platform note:** `run_all.sh` requires a POSIX bash shell. On Windows, run via Git Bash, WSL, or execution through PowerShell:

```powershell
foreach ($s in 42,43,44,45) { & venv\Scripts\python.exe "src\${s}_*.py" }

```

---

## Repository structure

```text
.
├── run_all.sh                      # Track B: Sequential runner script (bash)
├── Requirements.txt                # Pinned dependencies (SALib>=1.4.7)
├── README.md
├── LICENSE                         # MIT License
├── .gitignore                      # Excludes heavy binary output (*.npz, __pycache__/, venv/)
├── .kilo/                          # Configuration directory
├── docs/
│   ├── methodology_upgrade_summary.md    # Inverse estimation, robust MPC, and spatial metrics documentation
│   └── PROJECT_REVIEW_SUMMARY.md         # Review summary for Tracks B and C
├── src/                            # Pipeline scripts (01 through 65)
│   ├── 01–09   Multi-omic ingest, DE analysis, classical ML, and deep learning baselines
│   ├── 10–11   cVAE contrastive pretraining and latent extraction
│   ├── 12–18   C-GAT graph attention models, scVI/NMF baselines, and CSGT proofs
│   ├── 19–22   Waddington energy landscape, NEB saddle verification, and Fokker-Planck solver
│   ├── 23–26   Transfer entropy, causal gene network inference, and bootstrap validation
│   ├── 27–30   ABA tissue lattice mapping, reaction-diffusion PDEs, and invasion kinetics
│   ├── 31–34   Virtual single/dual gene knockouts and drug gating reports
│   ├── 35–41   Clinical cohort mapping (Ivy GAP/TCGA), Cox models, and dose-response profiles
│   ├── 42–48   Track B: Anisotropic PDE, stromal feedback, adaptive dosing, and 3D extension
│   ├── 49      Track C: Interactive 3D visualization dashboard
│   ├── 50a     Track A: Spatial genomics deconvolution via Bayesian ADVI
│   ├── 50b     Track C: Clinical DSS engine
│   ├── 51      Track C: Inverse parameter estimation routines
│   ├── 52a     Track A: 3D DTI solver with mechanical coupling
│   ├── 52b     Track C: Robust MPC solver
│   ├── 53a     Track C: Spatial evaluation metrics (DSC, HD95, MSD)
│   ├── 53b     Track A: 60-day virtual therapy simulations
│   ├── 54      Track A: 2D Neftel state to PDE parameter mappings
│   ├── 55a     Track A: Poroelastic mechanics solver
│   ├── 55b     Track A: Tensor construction for 8-patient cohort
│   ├── 56a     Track A: Mechanical coupling module for 3D DTI
│   ├── 56b     Track C: STABLE 3D DTI solver implementation
│   ├── 57      Track C: 90-day Stupp protocol engine
│   ├── 58      Track C: Gymnasium RL environment for adaptive dosing
│   ├── 59      Track C: Global sensitivity analysis and biomarker rule generation
│   ├── 60      Track C: Baseline comparisons and physical ablation scripts
│   ├── 61      Track C: Multi-seed convergence diagnostics
│   ├── 62      Track C: Bootstrap confidence interval calculations for biomarker cutoffs
│   ├── 63      Track C: Reward function sensitivity experiments
│   ├── 64      Track C: N=20 virtual cohort trial execution
│   ├── 65      Track C: Executive reporting and master figure generation
│   └── run_pipeline.py
├── tests/                          # Test suite (inverse estimation, MPC, spatial metrics)
└── output/                         # Processed data artifacts (JSON, PNG, MD, TSV, CSV)

```

---

## Key findings across tracks

### Track A (MSOS) — Single-cell to spatial dynamics

| Stage | Finding |
| --- | --- |
| **DE (02/03)** | Distinct expression profiles separate Core, Periphery, and Healthy zones |
| **Classical ML (06)** | Logistic Regression and Random Forest reach >94% accuracy on 2.5k highly variable genes |
| **Deep Baselines (07/08)** | Tabular Transformers underperform tree-based models on this dataset size |
| **cVAE (10/11)** | Contrastive pretraining builds a 32-dimensional latent representation for 140k cells |
| **C-GAT (12–14, 16)** | Graph attention network retains accuracy while encoding spatial topology |
| **NMF (16)** | Identifies 4 expression modules corresponding to AC, MES, NPC, and OPC cell states |
| **CSGT (18)** | Confirms a continuous expression gradient extending from Healthy tissue to the tumor core |
| **Waddington (19–22)** | Identifies dual attractor states (Healthy=0.56, Core=0.00) with a saddle point at the periphery |
| **Causal GRN (23–26)** | Highlights APOD, S100B, and MT3 as key upstream regulatory switches |
| **Invasion (27–30)** | ABA model yields front expansion speed of 2.42 µm/hr; FK PDE yields 20 µm/hr |
| **Drug Discovery (31–34)** | Identifies SDE2 (single) and S100A11+ZNF106 (dual) as candidate targets |
| **Clinical (35–41)** | Zone-stratified Cox models map regional recurrence risk across patient cohorts |

### Track B (10-Month PDE cohort) — Spatial PDE & dosing

| Metric | Result | Statistical Notes |
| --- | --- | --- |
| **Anisotropic fractal dimension (Df)** | 1.04–1.49 | Higher complexity vs isotropic (~0); t=24.74, p<0.001 |
| **Stromal front correlation** | 0.938–0.952 | Maintained across all 8 synthetic cases |
| **Adaptive dose-sparing** | 9–21% (mean 13.3±4.3%) | Significant dose reduction vs continuous dosing (p=5.2×10⁻⁵) |
| **TTP ratio (adaptive/MTD)** | 0.50–0.82 | Equivalent time-to-progression maintained with lower exposure |
| **Inflammation vs TTP (MTD)** | Pearson r = -0.98 | Strong inverse relationship (p<0.001) |
| **Sobol Index (ρ_s)** | 0.607 | Proliferation rate is the primary driver of TTP variation |
| **Dual-drug TTP** | 360 days | Outperforms MTD (281–296 days) and single-agent setups |
| **Dual-drug resistant fraction** | 0.038 | Suppresses resistant clone selection compared to MTD (~1.0) |
| **3D MTD final volume** | 0.0 mm³ | Complete local clearance in idealized 3D grid |
| **3D Adaptive final volume** | 40.0 ± 8.8 mm³ | Retains low tumor burden while saving ~67% cumulative dose |

### Track C (Digital Twin reactor) — 3D DTI & reinforcement learning

| Phase | Finding |
| --- | --- |
| **Inverse Estimation** | Parameter reconstruction yields <5% error on noise-free data, <15% under 10% gaussian noise |
| **Robust MPC** | Delivers 68.9% dose-sparing across 50 Monte Carlo trials under parameter uncertainty |
| **Spatial Metrics** | Establishes clinical evaluation thresholds (DSC ≥ 0.70, HD ≤ 5 mm) |
| **RL Adaptive (Phase 5)** | Controls tumor volume down to 1.04 mm³ vs 11.01 mm³ under Stupp protocol, using 13% less drug |
| **Biomarker Rule (Phase 6)** | **ρ > 0.024 day⁻¹** serves as a threshold where RL policy outperforms fixed dosing |
| **Ablations (Phase 7)** | Omitting DTI tensors increases simulated volume errors by 15–20% |
| **Convergence (Phase 7)** | Policy training stabilizes within a 5-seed CV of <5% on final tumor volume |

---

## Validation & reproducibility

* **Idempotency:** Running `45_validation_synthesis.py` outputs deterministic metrics matching stored SHA-256 signatures. Isotropic baselines are cached to `output/isotropic_baseline_metrics.json`.
* **Repository hygiene:** Binary `.npz` arrays are excluded via `.gitignore` to keep the repository lightweight while preserving JSON, PNG, and CSV outputs.
* **Physics checks:** Tensors are verified for symmetric positive-definiteness ($0.0$ residual), and mass conservation errors remain near machine precision ($1.77 \times 10^{-16}$).
* **Uncertainty quantification:** Parameter estimates use 100-sample bootstrap runs; biomarker cutoffs use 1000-sample bootstrap iterations.
* **Sensitivity analysis:** Global variance decomposition uses SALib Sobol indices in Track B and Latin Hypercube Sampling in Track C.
* **Validation tests:** Core routines are tested under `tests/` (`test_inverse_estimation.py`, `test_robust_mpc.py`, `test_spatial_metrics.py`).

---

## Limitations

This repository is an in-silico research framework built using synthetic patient cohorts. Known constraints include:

* Lack of prospective clinical validation.
* Cohorts are computationally generated (8 patients in Track B, 20 in Track C, 120 in Track A).
* PK/PD relationships use simplified one-compartment models.
* Resistance mechanisms and toxicity profiles rely on stylized assumptions.
* Spatial agreement between anisotropic and isotropic models shows low Dice overlap (DSC ~0.21), reflecting fundamental differences in transport mechanics rather than solver error.
* Scripts `56a` and `56b` represent two distinct 3D DTI solvers; `56b` is the active version used in later phases.
* Designed strictly for computational research, not clinical decision-making.

---

## Future work

* Ingestion of real clinical 3x3 DTI patient scans.
* External validation against public datasets (BraTS, TCGA-GBM, Ivy GAP).
* Enforcing brain boundary constraints (zero-flux masks for skull and ventricles).
* Adding multi-clone dynamics to track complex clonal competition.
* Benchmarking RL policies directly against multi-horizon Model Predictive Control.
* Containerization (Docker) for environment setup.

---

## Citation

If you use this repository in your research, please cite:

```bibtex
@software{gbm_pde_rl_cohort_2026,
  title  = {Biophysical Glioblastoma Digital Twin and Reinforcement Learning Adaptive Therapy Framework},
  author = {Vihan},
  year   = {2026},
  note   = {A three-track computational framework spanning single-cell dynamics, spatial PDEs, and RL-driven adaptive therapy optimization.}
}

```

---

## Disclaimer

This software is provided purely for **computational research and educational purposes**. It is **not** clinically validated and must **not** be used for diagnostic or treatment decisions.

```

```
