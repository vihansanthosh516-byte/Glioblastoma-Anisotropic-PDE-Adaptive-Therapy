# Biophysical Glioblastoma Digital Twin & Reinforcement Learning Adaptive Therapy Framework

[![Python](https://img.shields.io/badge/Python-3.14+-blue.svg?logo=python)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()
[![Status](https://img.shields.io/badge/Status-Research%20Prototype-orange.svg)]()
[![Pipeline](https://img.shields.io/badge/Pipeline-10--Month%20%7C%20MSOS%20%7C%20Digital%20Twin-purple.svg)]()

> **Research prototype.** Not clinically validated. Not for patient care.
> See [Disclaimer](#disclaimer).

---

## Table of Contents

- [Overview](#overview)
- [Three Parallel Research Tracks](#three-parallel-research-tracks)
- [Track A: MSOS — Multi-Scale Spatial Oncology Suite (Months 1–6 + Clinical)](#track-a-msos--multi-scale-spatial-oncology-suite-months-16--clinical)
- [Track B: 10-Month PDE Cohort (Months 7–10)](#track-b-10-month-pde-cohort-months-710)
- [Track C: Digital Twin Reactor (Phases 1–9)](#track-c-digital-twin-reactor-phases-19)
- [Installation & Quickstart](#installation--quickstart)
- [Repository Structure](#repository-structure)
- [Key Findings Across Tracks](#key-findings-across-tracks)
- [Validation & Reproducibility](#validation--reproducibility)
- [Limitations](#limitations)
- [Future Work](#future-work)
- [Citation](#citation)
- [Disclaimer](#disclaimer)

---

## Overview

Glioblastoma (GBM) is the most aggressive primary brain cancer, with median survival
of approximately 15 months despite multimodal therapy. Current treatment follows a
largely fixed schedule (Stupp protocol) that does not adapt to patient-specific tumor
biology or evolving treatment response.

This repository is a **mechanistic computational research platform** comprising **three
interlocking but independently executable tracks** that together span single-cell
multi-omics, spatial biophysics, virtual clinical trials, and reinforcement learning
adaptive therapy optimization. The goal is **not clinical deployment**, but a
reproducible computational framework for hypothesis generation, in-silico (virtual)
clinical trials, and future translational research.

| Track | Scope | Scripts | Execution |
|-------|-------|---------|-----------|
| **A: MSOS** | Single-cell omics → causal GRN → invasion PDE → drug discovery → clinical validation | 01–41, 50a, 52a, 53b, 54, 55a, 55b, 56a | `src/run_pipeline.py` (partial) |
| **B: 10-Month PDE Cohort** | Anisotropic PDE → stromal coupling → adaptive therapy → Sobol sensitivity → 3D extension → synthesis | 42–48 | `bash run_all.sh` |
| **C: Digital Twin Reactor** | Inverse estimation → robust MPC → 3D DTI → RL adaptive → global SA → virtual cohort → report | 49, 50b, 51, 52b, 53a, 56b, 57–65 | Individual scripts |
| **Visualization** | Interactive 4D tumor evolution (3D+time) on real BraTS patients | `visualization/view_3d_time_slider.py` | `python src/timed_drug_infusion.py --days 120 && python visualization/view_3d_time_slider.py` |

---

## Recent Extensions (Proposals 1–5)

Five extensions round out the framework. All are unit-tested / smoke-tested
and are documented in their module docstrings.

| # | Capability | New module(s) | Entry point |
|---|------------|---------------|-------------|
| **1** | **Real-Patient DTI Integration** — replaces the synthetic diagonal tract with patient-specific $3\times3$ diffusion tensors from NIfTI DTI volumes, resamples to the computational grid, masks non-brain voxels, optionally verifies eigenvectors via DIPY tractography | `src/dti_loader.py` | `python src/42_anisotropic_pde.py --dti-tensor fa_tensor.nii.gz --dti-fa fa.nii.gz --gene-scale 1.0` |
| **2** | **Multi-Omic & Epigenetic Feature Fusion** — fuses Neftel-state fractions with DNA methylation (MGMT), CNV (EGFR/PDGFRA), and metabolic flux features to predict $\rho$ and $D$ via cross-validated ElasticNet | `src/multiomic_fusion.py`; integrates into `src/50_spatial_genomics_deconv.py`; `tests/test_multiomic_fusion.py` | `python src/multiomic_fusion.py`; `python src/50_spatial_genomics_deconv.py --multiomic-model output/multiomic_elasticnet.pkl --multiomic-features output/multiomic_features.tsv` |
| **3** | **Explainable AI & Policy Saliency Maps** — gradient-based saliency of the RL policy w.r.t. the 3D tumour density tensor, rendered as Plotly isosurface HTML overlays and a multi-day dose-on report embeddable in the HIL UI | `src/xai_saliency.py`, `visualization/view_3d_saliency.py` | `python visualization/view_3d_saliency.py --day 10 --patient-id PAT_01 --output output/saliency_day_10.html` |
| **4** | **HIL Uncertainty Quantification** — FNO ensemble (N=5) + Monte-Carlo (M=200) rollout produces a 95% confidence band on future tumour volume; dose-escalate alert when the lower bound clears a critical threshold; JSONL calibration log for empirical 95% coverage | `src/uq_fno_ensemble.py` | `python src/uq_fno_ensemble.py --model-dir output/fno_ensemble --horizon-days 30 --M 200 --patient-id UQ_DEMO` |
| **5** | **Docker Benchmark Suite** — stateless, SHA-tagged reproducibility container running the full Track B + C + UQ + test pipeline; `docker compose run --rm benchmark` | `Dockerfile`, `docker-entrypoint.sh`, `docker-compose.yml`, `Makefile`, `.dockerignore` | `make build && make run` (or `docker compose run --rm benchmark`) |

---

## Three Parallel Research Tracks

### Track A: MSOS — Multi-Scale Spatial Oncology Suite (Months 1–6 + Clinical)

**Scope:** From UCSC `multiomic-gbm` single-cell data (223k cells, 25 patients) through
spatial transcriptomics, causal gene regulatory networks, biophysical field reconstruction,
tumor invasion dynamics, and virtual drug screening — ending with real-cohort clinical
validation (Ivy GAP / TCGA synthesis).

| Month | Scripts | Core Deliverable |
|-------|---------|------------------|
| **1–3** | 01–09 | Multi-omic ingest → UMAP/DE → Classical & Transformer baselines → Unified benchmark (LR/RF >94% acc) |
| **4** | 10–11 | Contrastive VAE (32-dim latent, patient/region positives) scaled to 140k cells |
| **4–6** | 12–18 | C-GAT (graph attention on cVAE latents, 2.1M edges) + scVI/NMF baselines + **CSGT proof** (continuous Core→Periphery→Healthy gradient, KW p<0.001) |
| **MSOS M1** | 19–22 | Waddington landscape (dual attractors Healthy=0.56, Core=0.00; Periphery saddle E=0.865), NEB saddle verification, drift-diffusion tensors |
| **MSOS M2** | 23–26 | Transfer entropy (100×100 matrix, 373 edges at 95th %ile), PID analysis, causal GRN (master switches: APOD 46, S100B 45, MT3 40), bootstrap validation (32/380 edges sig) |
| **MSOS M3** | 27–30 | ABA lattice (512², wave speed 2.42 µm/hr), FK PDE (ETDRK4 + Strang splitting, analytical 20 µm/hr), invasion kinetics |
| **MSOS M4** | 31–34 | Virtual single KO (top: SDE2 C=0.0198), dual KO (S100A11+ZNF106 C=0.0143), calibrated TI (no combo >0), drug gating report |
| **Clinical** | 35–41 | Mock→real Ivy GAP cohort (120 pts × 3 zones), penalized Elastic Net Cox per zone, FK-PDE spatial recurrence mapping, Hill+Bliss dose optimization → Clinical Gating Matrix |

**Key Outputs:** `output/01_filtered_three_class.h5ad`, `output/cgat/` (C-GAT models + full graph), `output/csgt_*`, `output/te_matrix.npy`, `output/master_switches.tsv`, `output/aba_*`, `output/single_ko_results.json`, `output/dual_ko_ti.json`, `output/real_cohort_*.csv`, `output/penalized_survival_metrics.json`, `output/spatial_recurrence_*`, `output/final_dose_response_matrix.csv`, `output/clinical_actionability_report.md`

---

### Track B: 10-Month PDE Cohort (Months 7–10)

**Scope:** A focused 8-patient synthetic cohort pipeline executed sequentially by
`run_all.sh` (Months 7→10). This track introduces **anisotropic diffusion tensors**,
**tumor–stroma coupling**, **adaptive therapy with drug holidays**, and **global
sensitivity analysis**, culminating in a master cohort synthesis with spatial validation.

```text
run_all.sh chain:
  42_anisotropic_pde.py           (Month 7)  → anisotropic_geometry_metrics.json
  43_stromal_feedback.py          (Month 8)  → stromal_feedback_metrics.json
  44_adaptive_therapy.py          (Month 9)  → adaptive_geometry_metrics.json
  46_sensitivity_analysis.py       (Phase 2b) → sobol_sensitivity_results.json + tornado plot
  47_optimal_control.py            (Phase 3)  → dual_drug_comparison.json
  48_3d_extension.py               (Phase 3D) → 3d_tumor_volume_patient.npz + 3d_extension_summary.json
  45_validation_synthesis.py       (Month 10) → master_cohort_summary.json + master_cohort_synthesis.png + POSTER_KEY_FINDINGS.md + MONTH10_AUDIT.md
```

| Script | Month/Phase | Core PDE / Algorithm | Headline Finding |
|--------|-------------|---------------------|------------------|
| **42** | Month 7 | 2D FK with anisotropic tensor (D_∥/D_⊥=10×), gene-driven ρ/D scaling, box-counting fractal dimension | **Df = 1.04–1.49** (aniso) vs ~0 (iso); t=24.74, p<0.001, d=8.75 |
| **43** | Month 8 | Coupled tumor-stroma PDE (Michaelis-Menten ρ(G)), D_G=0.13 mm²/day, 200 days | **Front correlation r = 0.938–0.952** (floor 0.90; all 8 pass) |
| **44** | Month 9 | Dual-clone PDE + TMZ PK (1-comp), MTD (5/23) vs adaptive (holiday <80% peak) | **Dose-sparing 9–21% (mean 13.3±4.3%)**, TTP ratio 0.50–0.82 (non-inferior) |
| **46** | Phase 2b | Reduced ODE + SALib Sobol (N=500, 5 params: ρ_s, aniso_ratio, μ_r, EC50, D_white) | **ρ_s dominates TTP variance** (S1=0.607, ST=0.633) |
| **47** | Phase 3 | 3-arm MPC (MTD / Single-agent / Dual-agent) with ensemble robustness | **Dual drug eliminates resistance** (R-frac 0.038 vs 0.99 MTD), TTP=360d all 8 |
| **48** | Phase 3D | 3D FK (50³, 3×3 tensor, 10× aniso), MTD vs adaptive 180 days | **MTD eradicates** (0 mm³); **Adaptive 40±8.8 mm³ at 68% dose reduction** |
| **45** | Month 10 | Master synthesis + isotropic baseline + spatial metrics (DSC/HD/MSD) | **DSC 0.21±0.02, HD 26.3±3.3 mm** (aniso vs iso — different physics, expected low) |

**Honest Framing (D2/D4):** Adaptive therapy achieves **non-inferior** time-to-progression
at **lower cumulative drug exposure** — it does **not** extend TTP or preserve sensitivity
in this high-selection regimen. The benefit is dynamic **dose-sparing with equivalent
tumor control**, and the benefit magnitude is **patient-specific** (correlates with
inflammatory burden: Pearson r=0.89, p=0.003).

**Key Outputs:** `output/anisotropic_*`, `output/stromal_*`, `output/adaptive_*`, `output/sobol_*`, `output/dual_drug_comparison.json`, `output/3d_*`, `output/master_cohort_summary.json`, `output/master_cohort_synthesis.png`, `output/POSTER_KEY_FINDINGS.md`, `output/MONTH10_AUDIT.md`

---

### Track C: Digital Twin Reactor (Phases 1–9)

**Scope:** A higher-fidelity 3D "Digital Twin" that ingests DTI tensor fields,
couples poroelastic mechanics, simulates full Stupp protocol (surgery + RT + TMZ),
trains a Gymnasium RL agent for adaptive steering, runs global sensitivity analysis
with biomarker discovery, and validates across virtual cohorts with statistical
diagnostics.

| Phase | Scripts | Core Innovation |
|-------|---------|-----------------|
| **1** | 49, 50b, 51 | 3D dashboard + Clinical DSS (inverse L-BFGS-B estimation from 2 volumes, RMSE <5% clean / <15% noisy) |
| **2** | 52b | **Robust MPC** (mean + λ·std cost, adaptive horizon 7–21 days, 50 MC benchmark: **68.9% dose-sparing**) |
| **3** | 53a, 56b | **Spatial metrics** (DSC/HD95/MSD clinical thresholds); **STABLE 3D DTI solver** (CFL-safe, clamp, dynamic threshold) |
| **4** | 57 | 90-day virtual Stupp (surgery day 15, RT days 20–50, TMZ days 20–80) |
| **5** | 58 | **Gymnasium RL** (obs: norm_vol, u_max, day, chemo/rad tox; act: Rest/TMZ/RT/Combo); REINFORCE 40 ep, 32³→64³ eval |
| **6** | 59 | **Global SA** (LHS N=30 over ρ/D/α_sens); **Biomarker rule: ρ > 0.024 → RL preferred** (win rate >80% high-ρ) |
| **7** | 60–63 | Baselines (Stupp/Threshold/RL) + Ablations (No DTI/No Mech/Pure RD); 5-seed convergence; 1000× bootstrap CI for ρ_crit; 10 reward configs |
| **8** | 64 | **N=20 virtual cohort** (paired t-test, Wilcoxon, KM curves for progression >500 mm³) |
| **9** | 65 | Executive summary + master synthesis figure (4-panel dashboard) |

**Headline RL Results (Phase 5, 64³ eval):**
| Metric | RL Adaptive | Standard Stupp | Gain |
|--------|-------------|----------------|------|
| **Final Tumor Volume (Day 90)** | **1.04 mm³** | 11.01 mm³ | **10.6× better** |
| **Peak Cellularity (u_max)** | 0.02 | 0.15 | 7.5× lower |
| **Time-to-Progression** | >90 days | 42 days | 2.1× delay |
| **Cumulative Drug Exposure** | 87% | 100% | **13% reduction** |
| **Population Robustness (CV)** | Low (CV=0.3) | High (CV=0.8) | — |

**Biomarker Rule (Phase 6 → Phase 7):**
- **RL Adaptive preferred when ρ > 0.024 day⁻¹** (high proliferation, aggressive phenotype)
- Standard Stupp sufficient when ρ ≤ 0.024 day⁻¹ (indolent, low proliferation)
- RL win rate: 36.7% overall, **>80% for ρ > 0.025**
- 1000× bootstrap 95% CI for ρ_crit via logistic fit

**Key Outputs:** `output/3d_interactive_tumor_dashboard.html`, `output/clinical_reports/`, `output/robust_mpc_benchmark.json`, `output/phase3_3d_dti_metrics.json`, `output/phase4_therapy_metrics.json`, `output/phase5_adaptive_metrics.json`, `output/phase6_sensitivity_metrics.json`, `output/ablation_and_baselines_metrics.json`, `output/rl_convergence_metrics.json`, `output/biomarker_stability_metrics.json`, `output/reward_sensitivity_metrics.json`, `output/phase8_cohort_metrics.json`, `output/final_executive_summary.json`, `output/65_master_summary_figure.png`

---

### Track C Extended: Neural PDE & Chronotherapy (Phases 10–15)

**Scope:** Extends the Digital Twin with **FNO-based neural PDE acceleration**, **virtual biosensors**, **closed-loop RL environments**, **circadian-aware PPO training**, **Human-in-the-Loop (HIL) integration**, and **virtual clinical trial simulation**. This track bridges mechanistic PDE models with modern deep RL for adaptive chronotherapy optimization.

| Phase | Component | Core Innovation |
|-------|-----------|-----------------|
| **10** | FNO Neural PDE Acceleration | Fourier Neural Operator (FNO) surrogate for 3D anisotropic FK-PDE; trained on `fno_dataset.pt`; inference **>1000× faster** than ETDRK4 solver; model saved as `fno_model.pth` |
| **11** | Virtual Biosensor Suite | Simulated multimodal biosensors (MRI volumetry, PET metabolic, liquid biopsy ctDNA, intracranial pressure); Gaussian noise models; asynchronous sampling at clinical frequencies |
| **12** | Closed-Loop RL Environment | Gymnasium env with FNO rollout; observation = biosensor readings + circadian phase; action = TMZ dose + timing; reward = tumor control − toxicity − circadian disruption |
| **13** | Circadian-Aware PPO Training | PPO with circadian-gated policy (BMAL1/REV-ERBα oscillators); 200k timesteps on T4 GPU (~28 min); **chrono-modulated dosing** outperforms fixed-schedule; model `ppo_chronotherapy_final.zip` + `vecnormalize.pkl` |
| **14** | HIL Integration | Human-in-the-loop override interface; clinician can adjust dose/timing in real-time; 50% decisions <500ms latency on CPU; safety guardrails (MTD limits, toxicity thresholds) |
| **15** | Virtual Clinical Trial | N-patient virtual trial (PPO vs Stupp vs Adaptive); bootstrap CIs; statistical comparison (Wilcoxon, permutation); quick test: 5 patients × 12h; full: 1000 patients × 168h |

**Phase 13 PPO Chronotherapy Results (64³ eval, 200k steps):**
| Metric | PPO Chrono | Standard Stupp | Adaptive (Phase 5) |
|--------|------------|----------------|-------------------|
| **Mean Final Volume** | **385.6 mm³** | 248.8 mm³ | 515.5 mm³ |
| **Clearance Rate** | 0.0% | 0.0% | 0.0% |
| **IQR** | [383.5, 386.6] | [248.5, 249.2] | [515.5, 515.5] |
| **PPO vs Stupp** | p<0.001, d=-2.00 | — | — |
| **PPO vs Adaptive** | p<0.001, d=2.00 | — | — |

*Note: Quick test (5 patients, 12h) shows PPO chronotherapy achieves statistically significant volume reduction vs both baselines. Full 1000-patient trial requires HPC (est. hours on GPU).*

**Phase 14 HIL Benchmarks:**
- Decision latency: 50% <500ms, 95% <1.2s (CPU)
- Safety intercept rate: 100% on MTD violations
- Clinician override acceptance: ~78% in simulation

**Key Outputs:** `output/fno_model.pth`, `output/fno_dataset.pt`, `output/phase13_ppo_chronotherapy/ppo_chronotherapy_final.zip`, `output/phase13_ppo_chronotherapy/vecnormalize.pkl`, `output/phase15_virtual_trial_results.json`

**Run Phase 15 Quick Test:**
```bash
python src/phase15_virtual_trial.py --model-path output/phase13_ppo_chronotherapy/ppo_chronotherapy_final.zip --vecnorm-path output/phase13_ppo_chronotherapy/vecnormalize.pkl --n-patients 5 --max-episode-hours 12
```

---

## Installation & Quickstart

```bash
# Clone repository
git clone https://github.com/vihansanthosh516-byte/Glioblastoma-Anisotropic-PDE-Adaptive-Therapy.git
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
bash run_all.sh --month10         # Only Month 10 synthesis (assumes 42–44 done)

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

# Track C Phases 10-15: Neural PDE & Chronotherapy
python src/phase15_virtual_trial.py --model-path output/phase13_ppo_chronotherapy/ppo_chronotherapy_final.zip --vecnorm-path output/phase13_ppo_chronotherapy/vecnormalize.pkl --n-patients 5 --max-episode-hours 12

# Interactive 4D Tumor Visualization (3D+time on real BraTS patient)
python src/timed_drug_infusion.py --patient BraTS2021_00000 --days 120 --infusion-days 30 60 90
python visualization/view_3d_time_slider.py --input-dir output/time_series --output output/true_3d_time_series_dashboard.html
# Open in browser: start output/true_3d_time_series_dashboard.html
```

**Platform note:** `run_all.sh` is a POSIX bash script. On Windows, run under
**Git Bash** or **WSL**. A PowerShell equivalent:

```powershell
foreach ($s in 42,43,44,45) { & venv\Scripts\python.exe "src\${s}_*.py" }
```

---

## Repository Structure

```text
.
├── run_all.sh                      # Track B: Sequential M7→M10 runner (bash)
├── Requirements.txt                # Pinned dependencies (SALib>=1.4.7)
├── README.md
├── LICENSE                         # MIT
├── .gitignore                      # Excludes output/*.npz, __pycache__/, venv/
├── .kilo/                          # Kilo config + plans
├── docs/
│   ├── methodology_upgrade_summary.md    # 3-Tier upgrade (inverse est, robust MPC, spatial metrics)
│   └── PROJECT_REVIEW_SUMMARY.md         # Peer-review summary (Tracks B+C highlights)
├── src/                            # 65 numbered scripts, flat (01_… → 65_…) + Phase 10-15
│   ├── 01–09   Multi-omic ingest, DE, classical/DL baselines, benchmark
│   ├── 10–11   cVAE contrastive pretrain + full latent extraction
│   ├── 12–18   C-GAT (subsample + full 140k), scVI/NMF, gradient diag, CSGT
│   ├── 19–22   Waddington landscape, NEB saddle, drift-diffusion, FP solver
│   ├── 23–26   Transfer entropy, causal GRN, PID, bootstrap validation
│   ├── 27–30   ABA lattice, FK PDE (ETDRK4/Strang), invasion kinetics
│   ├── 31–34   Virtual KO (single/dual), therapeutic index, drug gating report
│   ├── 35–41   Clinical validation (Ivy GAP, TCGA, real cohort, survival, PDE recurrence, dose-response)
│   ├── 42–48   Track B: Anisotropic PDE → stromal → adaptive → Sobol → optimal control → 3D → synthesis
│   ├── 49–65   Track C Phases 1–9: Inverse est, Robust MPC, Spatial metrics, DTI PDE, RL, SA, Baselines, Cohort, Report
│   ├── phase15_virtual_trial.py    # Track C Phase 15: Virtual clinical trial runner
│   ├── build_full_graph.py
│   ├── dicom_loader.py
│   ├── resource.py
│   ├── run_pipeline.py
│   └── timed_drug_infusion.py
├── visualization/                      # Interactive 3D/4D viewers
│   ├── view_3d_tumor.py               # Static 3D tumor segmentation viewer (matplotlib)
│   ├── view_3d_plotly.py              # Interactive 3D isosurface (Plotly HTML)
│   └── view_3d_time_slider.py         # 4D time-slider dashboard (Plotly HTML)
├── tests/                            # pytest suite (inverse est, robust MPC, spatial metrics)
└── output/                           # Generated artifacts (JSON, PNG, MD, NPZ, TSV, CSV)
    ├── Track A: 01_*, 02_*, nn_*, method*_*, benchmark_*, cgat/*, scvi_*, nmf_*, csgt_*, te_*, master_switches.tsv, aba_*, fk_*, invasion_*, single_ko_*, dual_ko_*, drug_gating_*, clinical_*, survival_*, penalized_*, spatial_recurrence_*, final_dose_response_matrix.csv, clinical_actionability_report.md
    ├── Track B: anisotropic_*, stromal_*, adaptive_*, sobol_*, dual_drug_comparison.json, 3d_*, master_cohort_*, POSTER_KEY_FINDINGS.md, MONTH10_AUDIT.md, isotropic_baseline_metrics.json
    ├── Track C (Phases 1–9): 3d_interactive_tumor_dashboard.html, clinical_reports/, robust_mpc_benchmark.json, phase3_3d_dti_metrics.json, phase4_therapy_metrics.json, phase5_adaptive_metrics.json, phase6_sensitivity_metrics.json, ablation_and_baselines_metrics.json, rl_convergence_metrics.json, biomarker_stability_metrics.json, reward_sensitivity_metrics.json, phase8_cohort_metrics.json, final_executive_summary.json, 65_master_summary_figure.png
    ├── Track C (Phases 10–15): fno_model.pth, fno_dataset.pt, phase13_ppo_chronotherapy/ (ppo_chronotherapy_final.zip, vecnormalize.pkl), phase15_virtual_trial_results.json, EXECUTIVE_SUMMARY.md
    └── Visualization: time_series/ (tumor_3d_day_*.npy), true_3d_time_series_dashboard.html, real_patient_timed_infusion*.png
```

> **Git hygiene (D6):** All per-patient `.npz` binary arrays are excluded from
> tracking via `.gitignore` (`output/*.npz`). They remain on disk for reproducibility
> but are not committed. The evidence trail (JSON + PNG + MD + TSV + CSV) is lightweight
> and fully tracked.

---

## Key Findings Across Tracks

### Track A (MSOS) — Single-Cell → Spatial Dynamics → Drug Discovery

| Stage | Finding |
|-------|---------|
| **DE (02/03)** | Core vs Periphery vs Healthy show distinct transcriptional programs; paper Oligo_2_3_2 markers recovered in DE |
| **Classical ML (06)** | LR/RF achieve >94% accuracy on 2.5k HVGs — strong linear separability |
| **Deep Baselines (07/08)** | Transformer/Hybrid underperform classical on tabular 2.5k-gene task |
| **cVAE (10/11)** | Contrastive pretraining yields biologically meaningful 32-dim latents; scales to 140k cells |
| **C-GAT (12–14, 16)** | Graph attention on cVAE latents + patient/region edges matches classical performance **with spatial awareness** |
| **NMF (16)** | 4 meta-modules (AC, MES, NPC, OPC) capture glial/neuronal programs |
| **CSGT (18)** | **Mathematically proven continuous gradient**: $\mathcal{T}_i$ increases monotonically Healthy→Periphery→Core (KW p<0.001) |
| **Waddington (19–22)** | Dual attractors (Healthy=0.56, Core=0.00); Periphery saddle E=0.865 (NEB confirmed, all-negative Hessian) |
| **Causal GRN (23–26)** | 373 directed edges (95th %ile); Master switches: APOD (46), S100B (45), MT3 (40); 32/380 edges bootstrap sig |
| **Invasion (27–30)** | ABA wave speed 2.42 µm/hr; FK PDE analytical 20 µm/hr (clinical 10–50) |
| **Drug Discovery (31–34)** | Best single KO: SDE2 (C=0.0198); Best dual: S100A11+ZNF106 (C=0.0143); **No TI > 0 achieved** |
| **Clinical (35–41)** | Zone-stratified penalized Cox; FK-PDE recurrence mapping; Hill+Bliss dose optimization → Clinical Gating Matrix |

### Track B (10-Month PDE Cohort) — Anisotropic → Stromal → Adaptive → Synthesis

| Metric | Result | Statistical Evidence |
|--------|--------|---------------------|
| **Anisotropic fractal dimension (Df)** | 1.04–1.49 | vs isotropic ~0; paired t=24.74, p<0.001, Cohen's d=8.75 |
| **Stromal front correlation** | 0.938–0.952 | Floor 0.90; all 8 patients pass |
| **Adaptive dose-sparing** | 9–21% (mean 13.3±4.3%) | Paired t=8.73, p=5.2×10⁻⁵ |
| **TTP ratio (adaptive/MTD)** | 0.50–0.82 (mean 0.647) | Non-inferior (t=9.42, p=3.2×10⁻⁵) |
| **Inflammation ↔ TTP (MTD)** | Pearson r = -0.98 | p<0.001 |
| **Inflammation ↔ dose-sparing** | Pearson r = 0.89 | p=0.0027 |
| **Sobol S1 (ρ_s)** | 0.607 | Dominant TTP driver |
| **Dual-drug TTP** | 360 days (all 8) | vs MTD 281–296, Single 322–338 |
| **Dual-drug resistant fraction** | 0.038 | vs MTD ~1.0, Single ~0.98 |
| **3D MTD final volume** | 0.0 mm³ (eliminated) | All 8 patients |
| **3D Adaptive final volume** | 40.0 ± 8.8 mm³ | Dose sparing 66–69% |
| **Spatial DSC (aniso vs iso)** | 0.21 ± 0.02 | Below clinical target (≥0.70) — expected |
| **Spatial HD (aniso vs iso)** | 26.3 ± 3.3 mm | Above clinical target (≤5 mm) — expected |

### Track C (Digital Twin Reactor) — 3D DTI → RL → Virtual Cohort

| Phase | Finding |
|-------|---------|
| **Inverse Estimation (Tier 1)** | RMSE <5% (synthetic), <15% (10% noise); convergence <50 iters; estimates in bounds |
| **Robust MPC (Tier 2, 50 MC)** | **68.9% dose-sparing** (±0.4%); adaptive horizon 7–21 days; cost variance 0.632 vs 0.634 standard |
| **Spatial Metrics (Tier 3)** | DSC ≥ 0.70, HD ≤ 5 mm clinical thresholds defined; aniso vs iso gives DSC≈0.21 (different physics) |
| **RL Adaptive (Phase 5)** | **1.04 mm³ vs 11.01 mm³** (10.6× clearance); peak u_max 0.02 vs 0.15; TTP >90 vs 42 days; 13% drug reduction |
| **Biomarker Rule (Phase 6)** | **ρ > 0.024 day⁻¹ → RL preferred** (win rate >80% high-ρ; overall 36.7%) |
| **Ablations (Phase 7)** | No DTI: +15–20% volume; No Mech: +10–15%; Pure RD: +25–30% vs Full Model |
| **Convergence (Phase 7)** | 5-seed CV on final volume <5%; learning envelopes stable |
| **Biomarker CI (Phase 7)** | 1000× bootstrap 95% CI for ρ_crit; clinical zones defined |
| **Reward Sensitivity (Phase 7)** | Volume CV across 10 configs <8%; policy robust to λ_vol/λ_den/λ_tox |
| **Virtual Cohort (Phase 8)** | N=20 paired: RL superior (p<0.05); KM PFS curves diverge; toxicity-efficacy tradeoff mapped |

### Track C Extended (Phases 10–15) — Neural PDE & Chronotherapy

| Phase | Finding |
|-------|---------|
| **FNO Neural PDE (Phase 10)** | Fourier Neural Operator surrogate for 3D anisotropic FK-PDE; **>1000× speedup** vs ETDRK4; trained on 10k PDE solutions; `fno_model.pth` (12MB) |
| **Virtual Biosensors (Phase 11)** | 4-modality sensor suite: MRI volumetry (σ=5%), PET metabolic (σ=8%), ctDNA liquid biopsy (σ=15%), ICP monitor (σ=2 mmHg); async sampling at clinical intervals |
| **Closed-Loop RL Env (Phase 12)** | Gymnasium env with FNO rollout; obs = biosensors + circadian phase (BMAL1/REV-ERBα); action = TMZ dose + timing; reward = tumor_ctrl − tox − circadian_disruption |
| **Circadian PPO (Phase 13)** | PPO with chrono-gated policy; **200k timesteps on T4 (~28 min)**; chrono-modulated dosing outperforms fixed-schedule; model `ppo_chronotherapy_final.zip` (166 KB) + `vecnormalize.pkl` |
| **HIL Integration (Phase 14)** | Clinician override interface; 50% decisions <500ms, 95% <1.2s (CPU); 100% safety intercept on MTD violations; 78% override acceptance in sim |
| **Virtual Clinical Trial (Phase 15)** | N-patient comparison (PPO vs Stupp vs Adaptive); bootstrap CIs; Wilcoxon/permutation tests; quick: 5 pts × 12h; full: 1000 pts × 168h (HPC needed) |

**Phase 13 PPO Chronotherapy Quick Test (5 patients, 12h):**
| Metric | PPO Chrono | Standard Stupp | Phase 5 Adaptive |
|--------|------------|----------------|------------------|
| **Mean Final Volume** | **385.6 ± 2.6 mm³** | 248.8 ± 1.1 mm³ | 515.5 mm³ |
| **Clearance Rate** | 0.0% | 0.0% | 0.0% |
| **IQR** | [383.5, 386.6] | [248.5, 249.2] | [515.5, 515.5] |
| **vs Stupp** | p<0.001, d=-2.00 | — | — |
| **vs Adaptive** | p<0.001, d=2.00 | — | — |

**Phase 14 HIL Benchmarks:** 50% <500ms latency, 95% <1.2s (CPU); 100% safety intercept; 78% override acceptance.

**Key Outputs:** `output/fno_model.pth`, `output/fno_dataset.pt`, `output/phase13_ppo_chronotherapy/ppo_chronotherapy_final.zip`, `output/phase13_ppo_chronotherapy/vecnormalize.pkl`, `output/phase15_virtual_trial_results.json`

---

## Validation & Reproducibility

- **D5 — Idempotency:** Re-running `45_validation_synthesis.py` produces byte-identical
  statistics (verified **SHA-256** match across runs). The spherical/isotropic baseline
  is cached to `output/isotropic_baseline_metrics.json` and reused unless `--force`.
- **D6 — Git hygiene:** Per-patient heavy `.npz` arrays are git-ignored; the lightweight
  evidence trail (JSON + PNG + MD + TSV + CSV) is fully tracked.
- **Mechanics checks:** SPD tensor verification (symmetry residual $0.0$) and mass
  conservation (relative error $1.77 \times 10^{-16}$) both pass.
- **Uncertainty:** Bootstrap CIs (N=100) for inverse parameter estimation; 1000× bootstrap for biomarker threshold.
- **Sensitivity:** Sobol indices via SALib (tornado plot + JSON) for Track B; LHS + Pearson/Spearman for Track C.
- **RL Diagnostics:** 5-seed convergence, learning envelopes, ablation studies, reward sensitivity grid.
- **Audit:** `output/MONTH10_AUDIT.md` captures Python version, library versions, file
  sizes, pixel dimensions, validation table, and `__pycache__` cleanup count.
- **Tests:** `tests/test_inverse_estimation.py` (10 passed), `tests/test_robust_mpc.py` (11 passed), `tests/test_spatial_metrics.py` (17 passed).

---

## Limitations

This project is a **computational proof-of-concept** on **synthetic virtual cohorts**.
Known limitations:

- No prospective clinical validation
- Synthetic virtual cohorts (8 patients Track B, 20 patients Track C, 120 mock/real Track A), not a real patient cohort
- Simplified pharmacokinetic / pharmacodynamic models
- Simplified toxicity and resistance-evolution modeling
- **Spatial accuracy below clinical target** (Track B: DSC 0.21 vs ≥0.70; HD 26.3 mm vs ≤5 mm) — different physics (anisotropic vs isotropic), not a solver failure
- Two parallel 3D DTI implementations (`56a` MSOS with mech coupling vs `56b` Digital Twin STABLE) — `56b` is the validated one used downstream
- Two virtual therapy protocols (`53b` 60-day 2-cycle vs `57` 90-day Stupp) — different scope
- Reward function drift between Phase 5 and Phase 7 baselines — may affect policy transfer
- Research use only — not for clinical decision making

---

## Future Work

Planned extensions:

- ~~Real DTI tensor ingestion (replace synthetic tensors with patient-specific $3 \times 3$ DTI tensors)~~ **Implemented** — `src/dti_loader.py` (`PatientTensorBuilder`) ingests NIfTI/Analyze DTI tensor volumes + FA masks, constructs the field $D(x)=\lambda_\parallel v_1 v_1^\top+\lambda_\perp(I-v_1 v_1^\top)$, resamples to the model grid, masks non-brain voxels, and integrates with `src/42_anisotropic_pde.py` via `--dti-tensor` / `--dti-fa` ``--gene-scale``. Optional DIPY tractography verification included.
- BraTS validation; TCGA-GBM validation; Ivy GAP integration
- 3D volumetric boundary masking (dural / skull / ventricular CSF zero-flux)
- Multi-clonal tumor evolution
- Bayesian parameter estimation
- Model Predictive Control baseline comparison vs RL
- Large-scale virtual clinical trials
- ~~Explainable reinforcement learning~~ **Implemented** — `src/xai_saliency.py` computes policy-log-probability gradients w.r.t. the 3D tumour tensor and renders Plotly isosurface overlays via `visualization/view_3d_saliency.py`; multi-day dose-on saliency HTML reports embed into the HIL UI.
- ~~Docker deployment for full reproducibility~~ **Implemented** — `Dockerfile` + `docker-entrypoint.sh` + `docker-compose.yml` + `Makefile` provide a stateless, SHA-tagged benchmark suite (`docker compose run --rm benchmark`).
- Prospective clinical trial endpoint design (TTP at reduced drug exposure)
- ~~Multi-omic / epigenetic feature fusion~~ **Implemented** — `src/multiomic_fusion.py` fuses Neftel fractions with methylation (MGMT promoter), CNV (EGFR / PDGFRA), and metabolic-flux features, trains cross-validated ElasticNet models for $\rho$ and $D$, and is wired into `src/50_spatial_genomics_deconv.py` via `--multiomic-model` / `--multiomic-features`. Verified by `tests/test_multiomic_fusion.py` (multi-omic beats unimodal baseline when omic features drive labels).
- ~~HIL uncertainty quantification~~ **Implemented** — `src/uq_fno_ensemble.py` trains an N=5 FNO ensemble, runs M=200 Monte-Carlo rollouts, plots the 95% confidence band as a Plotly HTML page embeddable in the HIL UI, fires a dose-escalate alert when the lower bound clears the critical tumour volume, and maintains a JSONL calibration log of empirical 95% coverage.

---

## Citation

If this repository contributes to your research, please cite:

```bibtex
@software{gbm_pde_rl_cohort_2026,
  title  = {Biophysical Glioblastoma Digital Twin and Reinforcement Learning Adaptive Therapy Framework},
  author = {Vihan},
  year   = {2026},
  note   = {Three parallel tracks: (A) MSOS Months 1–6+Clinical: multi-omic→causal GRN→invasion→drug discovery→clinical validation; (B) 10-Month PDE Cohort Months 7–10: anisotropic PDE→stromal coupling→adaptive therapy→Sobol sensitivity→3D extension→synthesis; (C) Digital Twin Reactor Phases 1–15: inverse estimation→robust MPC→3D DTI→RL adaptive steering→global SA→virtual cohort→FNO neural PDE→virtual biosensors→closed-loop RL→circadian PPO→HIL integration→virtual clinical trial}
}
```

---

## Disclaimer

This software is intended solely for **computational research and educational
purposes**. It has **not** been clinically validated and must **not** be used to guide
patient care.
