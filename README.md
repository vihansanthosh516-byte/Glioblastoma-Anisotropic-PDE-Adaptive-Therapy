# Glioblastoma Digital Twin & Multi-Modal Anisotropic Diffusion Validation Platform

[![Python](https://img.shields.io/badge/Python-3.12+-blue.svg?logo=python)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()
[![ISEF](https://img.shields.io/badge/ISEF-2027-purple.svg)]()

> **Research prototype for a science-fair computational biology project.** Not clinically validated. Not for patient care.

---

## Table of Contents

- [Overview](#overview)
- [Track A — MSOS: Single-Cell to Systems Biology](#track-a--msos-single-cell-to-systems-biology)
- [Track B — 10-Month Anisotropic PDE Cohort](#track-b--10-month-anisotropic-pde-cohort)
- [Track C — Digital Twin Reactor](#track-c--digital-twin-reactor)
- [Track C-Extended — Phases 10–15](#track-c-extended--phases-1015)
- [Validation Study: Three-Way Tensor Comparison on Real Patients](#validation-study-three-way-tensor-comparison-on-real-patients)
- [Clinical Decision Rules](#clinical-decision-rules)
- [Installation](#installation)
- [Repository Structure](#repository-structure)
- [Reproducing the Pipeline](#reproducing-the-pipeline)
- [Testing](#testing)
- [Citation](#citation)
- [Disclaimer](#disclaimer)

---

## Overview

Glioblastoma (GBM) is the most aggressive primary brain cancer in adults, with median survival of approximately 15 months despite maximal surgical resection, radiotherapy, and temozolomide (the Stupp protocol). The Stupp protocol is **open-loop** — it applies the same fixed schedule to every patient regardless of tumor genotype, invasion geometry, or how the tumor is actually responding.

This repository is a **multi-track computational research platform** built around that gap. It spans four pieces of work:

| Track | Scripts | What it does |
|---|---|---|
| **A — MSOS** (Multi-Scale Spatial Oncology Suite) | `src/01`–`41` | Single-cell multi-omics → causal gene-regulatory network → continuous tumor-state landscape → invasion PDE → virtual drug screening → clinical survival validation |
| **B — 10-Month PDE Cohort** | `src/42`–`48` | Anisotropic diffusion PDE (tumor spreads faster along white-matter tracts), tumor–stroma coupling, adaptive drug-holiday dosing, global sensitivity analysis, dual-drug optimal control, 3D extension |
| **C — Digital Twin Reactor** | `src/49`–`65` | Inverse estimation of patient-specific parameters from imaging, uncertainty-aware model-predictive control, a reinforcement-learning agent that steers therapy day-by-day, biomarker discovery, 20-patient virtual cohort validation |
| **C-Extended** | `neural_pde/`, `rl/`, `sensing/`, `hil/` | Neural-PDE surrogate (>1000× faster inference), virtual biosensors, circadian-aware PPO policy, hardware-in-the-loop pump interface, virtual clinical trial harness |

In addition, a **real-patient validation study** — the Standard vs. Atlas vs. DKI tensor comparison — tests the core physical assumption behind Track B and C (that anisotropic, tract-guided diffusion better predicts tumor shape than isotropic diffusion) against 62 real patients from the public UCSF-PDGM glioma imaging cohort. That study is described in its own section below.

---

## Track A — MSOS: Single-Cell to Systems Biology

Built from single-cell RNA-seq of ~223,000 cells across 25 patients (three tissue zones: Core, Periphery, Healthy).

| Stage | Scripts | Result |
|---|---|---|
| Classification benchmarks | 05–16 | Classical ML (Logistic Regression / Random Forest) reaches >94% zone-classification accuracy on 2,500 HVGs; a contrastive VAE + graph attention network ("C-GAT") matches this while adding spatial structure at 140k-cell scale |
| Continuous tumor-state proof (CSGT) | 17–18 | Healthy → Periphery → Core forms a statistically continuous transition, not three discrete clusters (Kruskal-Wallis p < 0.001) |
| Waddington landscape | 19–22 | Energy landscape reconstruction via Fokker-Planck density estimation; a Climbing-Image NEB solver formally proves the Periphery state is a saddle point between the Healthy and Core attractors |
| Causal gene-regulatory network | 23–26 | Directed transfer-entropy network (373 edges); master switch genes ranked by out-degree: **APOD (46), S100B (45), MT3 (40)** |
| Invasion modeling | 27–30 | Agent-based automaton + Fisher-Kolmogorov PDE; analytical wave speed 20 µm/hr, within the clinically observed 10–50 µm/hr range |
| Virtual drug screening | 31–34 | Single- and dual-gene knockout screens via latent-space perturbation; best single knockout SDE2, best pair S100A11+ZNF106 — no combination reached a positive therapeutic index under the calibrated model |
| Clinical validation | 35–41 | Zone-stratified penalized (Elastic Net) Cox survival model; FK-PDE spatial recurrence mapping; Hill+Bliss dose-response optimization feeding a Clinical Gating Matrix |

---

## Track B — 10-Month Anisotropic PDE Cohort

An 8-patient synthetic cohort modeled with a 2D anisotropic reaction-diffusion PDE (parallel diffusivity 10× perpendicular, approximating white-matter-tract-guided spread).

| Result | Value | Statistical support |
|---|---|---|
| Anisotropic vs. isotropic invasion geometry | Fractal dimension 1.04–1.49 vs. ≈0 | paired t = 24.74, p < 0.001, Cohen's d = 8.75 |
| Tumor–stroma front correlation | r = 0.938–0.952 across all 8 patients | above the 0.90 validation floor |
| Adaptive dose-sparing at equivalent control | 9–21% dose reduction (mean 13.3 ± 4.3%), non-inferior time-to-progression (ratio 0.50–0.82) | t = 8.73/9.42, p < 10⁻⁴ |
| Dominant driver of outcome variance | Proliferation rate ρ (Sobol S1 = 0.607) | SALib, N=500 base samples |
| Dual-drug resistance suppression | Resistant fraction 0.038 vs. ≈1.0 under MTD dosing | 3-arm MPC comparison |
| 3D extension (180-day horizon) | MTD eradicates (0 mm³); adaptive reaches 40.0 ± 8.8 mm³ at ~66–69% dose sparing | 50³ grid, full 3×3 diffusion tensor |

---

## Track C — Digital Twin Reactor

A higher-fidelity 3D digital twin wrapped in a full control loop: imaging → inverse-estimated parameters → PDE simulation → treatment decision.

| Component | Result |
|---|---|
| Inverse parameter estimation | RMSE < 5% (clean data), < 15% (10% noise); converges in < 50 iterations |
| Robust (uncertainty-aware) MPC | 68.9 ± 0.4% dose-sparing across 50 Monte-Carlo trials, no increase in outcome variance |
| RL adaptive therapy vs. Stupp protocol (day-90 volume, 64³ grid) | RL: 13.94 mm³ vs. Stupp: 11.01 mm³ — statistically non-inferior tumor control at ~13% lower cumulative drug exposure |
| Clinically actionable biomarker rule | ρ > 0.024 day⁻¹ → RL adaptive preferred (win rate > 80% in that regime, 36.7% overall) |
| Model-component ablations | Removing DTI anisotropy: +15–20% final volume; removing mechanics: +10–15%; pure reaction-diffusion: +25–30% |
| 20-patient virtual cohort | RL vs. Stupp non-inferiority confirmed (paired t p = 0.00067, Wilcoxon p = 0.00026, Cohen's d = −0.93; 73.8 mm³ vs. 137.8 mm³, 46.5% reduction in this cohort configuration) |

---

## Track C-Extended — Phases 10–15

Infrastructure aimed at a real-time, sensor-driven closed loop:

| Phase | Module | Result |
|---|---|---|
| 10 | `neural_pde/` | Fourier Neural Operator surrogate reproduces the 3D PDE solver's output at >1000× inference speed |
| 11 | `sensing/virtual_sensor.py` | 4-modality virtual biosensor suite (MRI, PET, ctDNA, ICP) with realistic per-modality noise models |
| 12–13 | `rl/chronotherapy_env.py`, `rl/train_chronotherapy.py` | Closed-loop Gymnasium environment + PPO training with a circadian-gated policy |
| 14 | `hil/pump_interface.py` | Hardware-in-the-loop pump interface; 50% of decisions <500ms, 95% <1.2s (CPU); safety watchdog intercepts simulated MTD violations |
| 15 | `phase15_virtual_trial.py` | Virtual clinical trial harness, scalable to large simulated cohorts |

---

## Validation Study: Three-Way Tensor Comparison on Real Patients

This study tests the core physical assumption behind Tracks B and C — that anisotropic, tract-guided diffusion predicts tumor shape better than isotropic diffusion — directly against real tumor segmentations, using the public **UCSF-PDGM** glioma imaging cohort (WHO Grade 2–3, IDH-mutant astrocytoma), 62 patients, 42/20 train/test split (seed=42, stratified).

Three ways of constructing the diffusion tensor that drives the PDE are compared:

| Mode | Tensor construction |
|---|---|
| **Patient-Specific DTI** | The patient's own fitted diffusion-tensor eigenvalues, used directly |
| **Anisotropy-Enhanced DTI** | Patient eigenvalues with anisotropy scaled up (λ₁×1.5, λ₂/λ₃×0.5) toward population-atlas-level values |
| **Kurtosis-Adjusted DTI** | Patient eigenvalues with a different scaling (λ₁×1.3, λ₂/λ₃×0.8) approximating a diffusion-kurtosis correction |

**Results** (Dice Similarity Coefficient between predicted and ground-truth segmentation, Wilcoxon signed-rank test on paired per-patient aniso-vs-iso DSC):

| Mode | DSC (aniso) | DSC (iso) | Δ (aniso − iso) | Patients where aniso wins | p-value |
|---|---|---|---|---|---|
| Patient-Specific DTI | 0.7843 | 0.9464 | −0.1621 | 0 / 62 | < 0.001 |
| Anisotropy-Enhanced DTI | 0.8229 | 0.7987 | +0.0242 | ~35 / 62 | < 0.05 |
| Kurtosis-Adjusted DTI | 0.1065 | 0.0275 | +0.0790 | 60 / 62 | < 0.001 |

**Train/test consistency (42/20 split):**

| Mode | Δ, train (n=42) | Δ, test (n=20) |
|---|---|---|
| Patient-Specific DTI | −0.162 | −0.162 |
| Anisotropy-Enhanced DTI | +0.018 | +0.036 |
| Kurtosis-Adjusted DTI | +0.067 | +0.103 |

**Reading these results:** using a patient's own raw DTI eigenvalues directly, anisotropic modeling *underperforms* a simple isotropic model on every single patient. Once anisotropy is heuristically enhanced — either toward atlas-typical values or via a kurtosis-style adjustment — the anisotropic model overtakes the isotropic one, most decisively under the kurtosis-style adjustment (60/62 patients). The absolute DSC values for the Kurtosis-Adjusted mode are low in both arms because the larger eigenvalue perturbation shifts both predictions further from the ground-truth mask; the relative comparison (which one tracks the true shape better) is what the win-rate and Wilcoxon statistics test, and both point the same direction. Reproducing this comparison end-to-end requires `run_improved_aniso.py` followed by `isef_figures.py`; the driver scripts currently have machine-specific data paths that need to be pointed at a local UCSF-PDGM download before rerunning.

**Figures:**

| Figure | Description |
|---|---|
| `fig1_bar_chart.png` | Anisotropic vs. isotropic DSC by mode, with error bars |
| `fig2_scatter.png` | Patient-level comparison, aniso DSC vs. iso DSC |
| `fig3_delta_boxplot.png` | Distribution of per-patient performance difference (aniso − iso) |

---

## Clinical Decision Rules

Two distinct, independently derived decision rules come out of this platform:

**1. Adaptive dosing rule (Track C, virtual cohort):**
> If a patient's estimated tumor proliferation rate ρ exceeds 0.024 day⁻¹ (estimable from two longitudinal MRI volumes via the platform's inverse-estimation module), the RL-adaptive dosing policy is preferred over standard fixed-schedule Stupp therapy (win rate > 80% in this regime, vs. 36.7% overall). Below this threshold, standard Stupp is sufficient. Bootstrapped (1000×) 95% CI on the threshold: ρ_crit ∈ [0.0202, 0.0249] day⁻¹.

**2. Tensor-construction rule (real-patient validation study):**
> Do not drive the anisotropic PDE from raw patient-specific DTI eigenvalues alone — on this cohort it underperformed isotropic modeling in every patient. Anisotropy-enhanced or kurtosis-adjusted tensor construction should be used instead when anisotropic modeling is desired, with the kurtosis-adjusted variant showing the most consistent aniso-over-iso advantage (60/62 patients).

Both are **research hypotheses generated from model output** (a virtual cohort in the first case, a single real cohort with no external replication in the second), not clinically validated recommendations.

---

## Installation

```bash
git clone https://github.com/vihansanthosh516-byte/Glioblastoma-Anisotropic-PDE-Adaptive-Therapy.git
cd Glioblastoma-Anisotropic-PDE-Adaptive-Therapy

python -m venv venv
source venv/bin/activate          # Linux/macOS
# .\venv\Scripts\Activate.ps1     # Windows PowerShell

pip install -r Requirements.txt
```

If you hit missing imports beyond what `Requirements.txt` pins, install them individually (`scanpy`, `anndata`, `umap-learn`, `torch`, `torch-geometric`, `gymnasium`, `stable-baselines3`, `networkx`, `seaborn`, `pymc`, `pytensor`, `nibabel`, `dipy`, `SALib`, and `plotly` cover the full pipeline across all tracks).

**Docker (Track B + C):**
```bash
make build && make run
```

---

## Repository Structure

```text
.
├── src/                          105 Python files: numbered pipeline (01–65)
│   ├── 01–41                     Track A: single-cell multi-omics → causal GRN → virtual drug screening
│   ├── 42–48                     Track B: anisotropic PDE cohort → adaptive therapy → sensitivity analysis
│   ├── 49–65                     Track C: digital twin → RL → virtual cohort → final report
│   └── neural_pde/, rl/, sensing/, hil/    Track C-Extended (Phases 10–15)
├── run_improved_aniso.py         Three-way tensor validation driver (Standard / Atlas / DKI)
├── cross_validation.py           42/20 train-test split + Wilcoxon testing for the validation study
├── isef_figures.py               Figure generation for the validation study
├── tests/                        9 pytest files, 56 tests
├── visualization/                5 interactive Plotly/matplotlib viewers (4D time-slider, saliency maps)
├── docs/                         Project review, methodology notes, dataset info
├── manuscript/                   Draft paper (PAPER.md, manuscript.v3.{md,docx,pdf}, paper.tex, references.bib)
├── scripts/                      BraTS downloader, figure generator, real-data validator
├── data/                         Clinical metadata; UCSF-PDGM imaging must be obtained separately
├── output/                       Generated results, figures, and CSVs (created when the pipeline runs)
├── run_all.sh, Dockerfile, docker-compose.yml, Makefile
└── REPO_SUMMARY.md, ISEF_FINAL_REPORT.md   Companion documents: full script-by-script inventory and submission narrative
```

---

## Reproducing the Pipeline

**Validation study (Standard vs. Atlas vs. DKI):**
```bash
python run_improved_aniso.py
python isef_figures.py
```

**Track B (10-Month PDE Cohort):**
```bash
bash run_all.sh          # runs scripts 42 → 48 → 45 in sequence
```

**Track C (Digital Twin Reactor):**
```bash
python src/51_inverse_parameter_estimation.py --test
python src/52_robust_mpc_controller.py --benchmark --n-mc 50
python src/58_rl_adaptive_steering.py
python src/64_virtual_cohort_simulation.py
python src/65_generate_final_report.py
```

**Track C-Extended (quick test, Phase 15):**
```bash
python src/phase15_virtual_trial.py \
    --model-path output/phase13_ppo_chronotherapy/ppo_chronotherapy_final.zip \
    --n-patients 5 --max-episode-hours 12
```

**Track A** requires the external UCSC Cell Browser `multiomic-gbm` dataset (~1.4 GB, https://cells.ucsc.edu/?ds=multiomic-gbm) and data-path configuration in `src/01`–`16`.

---

## Testing

```bash
pytest tests/
```
56 tests across 9 files, covering inverse estimation, joint inverse estimation, robust MPC, spatial metrics, multi-omic fusion, hybrid controller, MU-Glioma data loading, retrospective validation, and treatment models.

---

## Citation

```bibtex
@software{gbm_platform_2026,
  title  = {Glioblastoma Digital Twin and Multi-Modal Anisotropic Diffusion Validation Platform},
  author = {Vihan Santhosh},
  year   = {2026},
  note   = {Three research tracks (single-cell multi-omics, anisotropic PDE cohort,
            digital-twin reinforcement learning) plus a real-patient validation of
            Standard, Atlas, and DKI tensor construction on 62 UCSF-PDGM patients.}
}
```

---

## Disclaimer

This software is intended solely for **computational research and educational purposes** as part of an independent science-fair research project. It has **not** been clinically validated in any prospective or retrospective patient study and must **not** be used to guide patient care. Tracks B and C are evaluated on synthetic virtual cohorts generated by the models themselves. The tensor-validation study uses real, de-identified imaging data from a public research cohort (UCSF-PDGM); its results have not been independently replicated on any other cohort. All quantitative results in this README are as produced by the corresponding scripts in this repository — see `REPO_SUMMARY.md` for the full script-by-script reference.
