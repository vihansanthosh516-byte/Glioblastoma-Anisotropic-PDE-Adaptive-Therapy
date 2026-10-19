> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

> **Audit corrections, 2026-10-19:** Original layout retained. Struck claims are withdrawn or superseded. The evidence appendix gives corrected values, sources and specific missing evidence. Unstruck historical claims retain their original scope.

# Multi-Omic Anisotropic PDE Modeling & Adaptive Therapy Dynamics in Glioblastoma Cohorts

![Python](https://img.shields.io/badge/Python-3.14+-blue?logo=python)
![License](https://img.shields.io/badge/License-MIT-green)
~~![Month](https://img.shields.io/badge/Status-10--Month%20End--to--End%20Complete-brightgreen)~~
~~![Pipeline](https://img.shields.io/badge/Pipeline-Multi--Omic%20%7C%20Anisotropic%20PDE%20%7C%20Adaptive%20Therapy-purple)~~

---

## Abstract

~~We developed a patient-specific computational oncology platform for glioblastoma treatment planning that integrates:~~

~~1. **3D Anisotropic DTI Tract Modeling**: Patient-specific diffusion tensor fields aligned to white matter tracts derived from DTI-MRI (validated across 8 patients, Months 7–8)~~

~~2. **Inverse Biophysical Parameter Estimation**: Bounded optimization framework solving for growth rate (ρ) and diffusion coefficient (D) directly from longitudinal imaging volumes (T₀, T₁) (Tier 1: scripts 51, RMSE < 15% with noise)~~

~~3. **Uncertainty-Aware Adaptive-Horizon MPC**: Robust model predictive control optimizing dosing over confidence intervals (ρ±15%, D±15%) with dynamic horizon adjustment (7–21 days) based on tumor growth dynamics (Tier 2: scripts 52-53, dose-sparing ≥60%, non-inferior TTP)~~

~~4. **Spatial Validation Framework**: Quantitative assessment using Dice Similarity Coefficient (DSC) and Hausdorff Distance (HD) between simulated and observed tumor masks (Tier 3: script 53, cohort mean DSC = 0.21 ± 0.02 vs isotropic baseline)~~

~~5. **Deterministic Reproducibility**: SHA-256 mathematical provenance certification for offline execution audit~~

~~> **Note:** Track A (scripts 01–41, Months 1–6 omics & biomarkers) outputs are currently unverifiable. The presentation board and abstract focus on Tracks B, C (MPC adaptive dosing) and the DTI validation study where real, reproducible evidence exists.~~

~~The platform transforms multi-omic risk stratification (S100A8/S100A11/LST1 inflammatory signature) into a clinical decision support tool for neuro-oncology treatment planning.~~

---

## Key Findings & Clinical Summary

| Discovery / Metric | Result / Finding | Statistical Evidence |
| :--- | :--- | :--- |
| ~~**Inflammatory Stratification**~~ | ~~$S100A8/S100A11/LST1$ tri-gene signature stratifies patient risk tiers~~ | ~~Pearson $r = -0.99$, $p < 0.001$ vs TTP~~ |
| ~~**Adaptive TTP Non-Inferiority**~~ | ~~Equal progression-free survival vs MTD~~ | ~~Paired $t = -1.00$, $p = 0.35$ (non-inferior)~~ |
| ~~**Drug Toxicity Reduction**~~ | ~~**63.4% – 73.0%** cumulative dose-sparing (mean $68.1\% \pm 3.5\%$)~~ | ~~Paired $t = 55.0$, $p < 0.001$; 95% CI $[65.2\%, 71.1\%]$~~ |
| ~~**Anisotropic Fractal Fronts**~~ | ~~$D_f = 1.20 – 1.55$ (vs isotropic $\approx 0.0$)~~ | ~~Paired $t = 29.5$, $p < 0.001$, Cohen's $d = 10.43$~~ |
| ~~**Stromal Front Alignment**~~ | ~~$r \in [0.93, 0.96]$ tract correlation~~ | ~~Exceeds $r \ge 0.90$ floor across all 8 patients~~ |
| ~~**Parameter Estimation Accuracy**~~ | ~~RMSE < 5% (synthetic), < 15% (noisy)~~ | ~~N=100 bootstrap, 95% CI~~ |
| ~~**Robust MPC Cost Variance**~~ | ~~30% lower vs standard MPC~~ | ~~F-test p < 0.01~~ |
| ~~**Spatial DSC**~~ | ~~0.85 ± 0.08 (cohort mean)~~ | ~~Paired t-test vs isotropic baseline~~ |
| ~~**Spatial HD**~~ | ~~3.2 ± 1.1 mm (cohort mean)~~ | ~~Below 5mm clinical threshold~~ |

~~> **Honest Framing:** Adaptive therapy achieves **non-inferior** time-to-progression at substantially lower drug exposure — it does **not** extend TTP or preserve sensitivity in this high-selection regimen. The benefit is dynamic dose-sparing with equivalent tumor control.~~

---

## Full 10-Month Project Architecture

> Historical example only; not certified as current instructions or verified output. See evidence appendix.

```text
[ Months 1–3: Omics & Biomarkers ] ──> [ Months 4–6: ML Risk Stratification ]
                                                        │
                                                        ▼
[ Months 7–8: Anisotropic PDE & Stroma ] <── [ Patient Cohort Parameterization ]
            │
            ▼
[ Month 9: Evolutionary Adaptive Dosing ] ──> [ Month 10: Cohort Synthesis & Stats ]
---

## Installation & Quickstart

```bash
# Clone repository
~~git clone https://github.com/vihansanthosh516-byte/Glioblastoma-Anisotropic-PDE-Adaptive-Therapy.git~~
~~cd Glioblastoma-Anisotropic-PDE-Adaptive-Therapy~~

# ~~Create virtual environment (Python 3.14+)~~
~~python -m venv venv~~
~~source venv/bin/activate          # Linux / macOS~~
~~.\venv\Scripts\Activate.ps1       # Windows PowerShell~~

# ~~Install dependencies~~
~~pip install numpy scipy matplotlib pillow pandas~~

> Historical example only; not certified as current instructions or verified output. See evidence appendix.

```

**Platform note:** `run_all.sh` is a POSIX bash script. On Windows, run under **Git Bash** or **WSL**.

---

## Key Deliverables & Artifacts

All outputs are written to `output/` by the pipeline:

| Artifact | Description |
|----------|-------------|
| `master_cohort_synthesis.png` | Publication-grade **3821 × 3503 px** 4-panel master poster canvas (Panels A–D). |
| `master_cohort_summary.json` | Unified 8-patient JSON dataset with per-phase metrics, spherical baseline comparison, and full statistical test results ($t$-stats, $p$-values, effect sizes, CIs). |
| `POSTER_KEY_FINDINGS.md` | Dynamically generated, mathematically verified bullet points for presentation boards (no placeholder text). |
| `MONTH10_AUDIT.md` | Clean-environment execution audit log: Python environment, file inventory, validation summary, pycache cleanup record. |
| `isotropic_baseline_metrics.json` | Cached isotropic Fisher–Kolmogorov baseline (D3 idempotency). |
| `adaptive_*.npz`, `stromal_*.npz`, `anisotropic_*.npz` | Per-patient heavy binary arrays (excluded from git via `.gitignore`). |

---

## 3D Volumetric Extension (Phase 3D)

While the primary validation framework operates on high-resolution 2D patient slices ($100 \times 100\text{ mm}$ grid), the mathematical solver naturally extends to full 3D voxel meshes ($\mathbb{R}^3$). The 3D extension module (`src/48_3d_extension.py`) demonstrates this architectural upgrade:

| Feature | 2D Implementation | 3D Extension |
|---------|-------------------|--------------|
| **Grid** | $100 \times 100$ voxels (1.0 mm) | $50 \times 50 \times 50$ voxels (1.0 mm) |
| **Tensor** | $2 \times 2$ symmetric $\mathbf{D}_{2\times2}$ | $3 \times 3$ symmetric $\mathbf{D}_{3\times3}$ |
| **PDE** | $\partial_t u = \nabla_{2D} \cdot (\mathbf{D} \nabla_{2D} u) + \rho u (1 - u/K)$ | $\partial_t u = \nabla_{3D} \cdot (\mathbf{D} \nabla_{3D} u) + \rho u (1 - u/K) - \gamma C(t) u$ |
| **BCs** | Neumann zero-flux (4 faces) | Neumann zero-flux (6 faces) |
| **CFL** | $dt \leq dx^2 / (2 \cdot \max(D))$ | $dt \leq dx^2 / (2 \cdot 3 \cdot \max(D))$ |

### 3D Mathematical Formulation

The 3D anisotropic Fisher-Kolmogorov PDE for tumor density $u(\mathbf{x}, t)$ on a 3D voxel mesh $\Omega \subset \mathbb{R}^3$:

$$\frac{\partial u}{\partial t} = \nabla \cdot \left( \mathbf{D}(\mathbf{x}) \nabla u \right) + \rho u \left( 1 - \frac{u}{K} \right) - \gamma C(t) u$$

where the diffusion tensor expands to:

$$\mathbf{D}(\mathbf{x}) = \begin{bmatrix} D_{xx} & D_{xy} & D_{xz} \\ D_{yx} & D_{yy} & D_{yz} \\ D_{zx} & D_{zy} & D_{zz} \end{bmatrix}$$

### 3D Execution Results

| Metric | MTD | Adaptive Therapy |
|--------|-----|------------------|
| **Final Volume (180 days)** | 0.0 mm³ (eliminated) | 203.0 mm³ (controlled) |
| **Drug Exposure** | 100% (continuous) | 40.2% (59.8% sparing) |
| **Sphericity Index** | N/A | 1.000 (compact sphere) |

**Interpretation:** In 3D, MTD eliminates the tumor within 180 days, while adaptive therapy maintains stable disease at ~200 mm³ with **59.8% dose reduction** — consistent with the 2D cohort findings.

### Future Clinical Roadmap

1. **3D DTI Tensor Ingestion:** Replace synthetic tensor fields with patient-specific $3 \times 3$ diffusion tensors derived from clinical DTI sequences.
2. **Volumetric Boundary Masking:** Enforce 3D zero-flux conditions along dural/skull boundaries and ventricular CSF spaces.
3. **Digital Twin Integration:** Scale the 14-day MPC controller to 3D patient digital twins for prospective neuro-oncology treatment planning.

---

## Repository Tree & Git Exclusions

```text
.
~~├── .gitignore                    # Excludes output/*.npz, __pycache__/, bio_env/~~
~~├── .kilo/~~
~~│   ├── kilo.jsonc               # Kilo config (allow-all permissions)~~
~~│   └── plans/~~
~~│       └── 1784578649972-month-10-cohort-synthesis.md~~
~~├── run_all.sh                    # Sequential M7→M10 runner (bash)~~
~~├── src/~~
~~│~~
~~└── output/                       # Generated artifacts (JSON, PNG, MD)~~
~~    ├── master_cohort_summary.json~~
~~    ├── master_cohort_synthesis.png~~
~~    ├── POSTER_KEY_FINDINGS.md~~
~~    ├── MONTH10_AUDIT.md~~
~~    ├── isotropic_baseline_metrics.json~~
~~    └── *.npz                     # (git-ignored; kept on disk for reproducibility)~~
> Historical example only; not certified as current instructions or verified output. See evidence appendix.

```

> **Git hygiene (D6):** All 32 per-patient `.npz` binary arrays are excluded from tracking via `.gitignore` (`output/*.npz`). They remain on disk for reproducibility but are not committed. The evidence trail (JSON + PNG + MD) is lightweight and fully tracked.

---

## Reproducibility & Idempotency

* **D5 Idempotency:** Re-running `45_validation_synthesis.py` produces byte-identical statistics (verified SHA-256 hash match across runs). The spherical baseline is cached to `output/isotropic_baseline_metrics.json` and reused unless `--force` is passed.
* **Validation:** All three phase metric JSONs pass schema verification and range checks (8/8 patients PASS per phase).
* **Audit:** Full environment capture (`MONTH10_AUDIT.md`) includes Python version, library versions, file sizes, pixel dimensions, validation table, and `__pycache__` cleanup count.
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
~~git clone https://github.com/vihansanthosh516-byte/Glioblastoma-Anisotropic-PDE-Adaptive-Therapy.git~~
~~cd Glioblastoma-Anisotropic-PDE-Adaptive-Therapy~~

~~python -m venv venv~~
~~source venv/bin/activate          # Linux/macOS~~
# ~~.\venv\Scripts\Activate.ps1     # Windows PowerShell~~

~~pip install -r Requirements.txt~~
> Historical example only; not certified as current instructions or verified output. See evidence appendix.

```

If you hit missing imports beyond what `Requirements.txt` pins, install them individually (`scanpy`, `anndata`, `umap-learn`, `torch`, `torch-geometric`, `gymnasium`, `stable-baselines3`, `networkx`, `seaborn`, `pymc`, `pytensor`, `nibabel`, `dipy`, `SALib`, and `plotly` cover the full pipeline across all tracks).

**Docker (Track B + C):**
```bash
~~make build && make run~~
> Historical example only; not certified as current instructions or verified output. See evidence appendix.

```

---

## Repository Structure

```text
.
~~├── src/                          105 Python files: numbered pipeline (01–65)~~
~~│   ├── 01–41                     Track A: single-cell multi-omics → causal GRN → virtual drug screening~~
~~│   ├── 42–48                     Track B: anisotropic PDE cohort → adaptive therapy → sensitivity analysis~~
~~│   ├── 49–65                     Track C: digital twin → RL → virtual cohort → final report~~
~~│   └── neural_pde/, rl/, sensing/, hil/    Track C-Extended (Phases 10–15)~~
~~├── run_improved_aniso.py         Three-way tensor validation driver (Standard / Atlas / DKI)~~
~~├── cross_validation.py           42/20 train-test split + Wilcoxon testing for the validation study~~
~~├── isef_figures.py               Figure generation for the validation study~~
~~├── tests/                        9 pytest files, 56 tests~~
~~├── visualization/                5 interactive Plotly/matplotlib viewers (4D time-slider, saliency maps)~~
~~├── docs/                         Project review, methodology notes, dataset info~~
~~├── manuscript/                   Draft paper (PAPER.md, manuscript.v3.{md,docx,pdf}, paper.tex, references.bib)~~
~~├── scripts/                      BraTS downloader, figure generator, real-data validator~~
~~├── data/                         Clinical metadata; UCSF-PDGM imaging must be obtained separately~~
~~├── output/                       Generated results, figures, and CSVs (created when the pipeline runs)~~
~~├── run_all.sh, Dockerfile, docker-compose.yml, Makefile~~
~~└── REPO_SUMMARY.md, ISEF_FINAL_REPORT.md   Companion documents: full script-by-script inventory and submission narrative~~
> Historical example only; not certified as current instructions or verified output. See evidence appendix.

```

---

## Reproducing the Pipeline

**Validation study (Standard vs. Atlas vs. DKI):**
```bash
~~python run_improved_aniso.py~~
~~python isef_figures.py~~
> Historical example only; not certified as current instructions or verified output. See evidence appendix.

```

**Track B (10-Month PDE Cohort):**
```bash
~~bash run_all.sh          # runs scripts 42 → 48 → 45 in sequence~~
> Historical example only; not certified as current instructions or verified output. See evidence appendix.

```

**Track C (Digital Twin Reactor):**
```bash
~~python src/51_inverse_parameter_estimation.py --test~~
~~python src/52_robust_mpc_controller.py --benchmark --n-mc 50~~
~~python src/58_rl_adaptive_steering.py~~
~~python src/64_virtual_cohort_simulation.py~~
~~python src/65_generate_final_report.py~~
> Historical example only; not certified as current instructions or verified output. See evidence appendix.

```

**Track C-Extended (quick test, Phase 15):**
```bash
~~python src/phase15_virtual_trial.py \~~
~~    --model-path output/phase13_ppo_chronotherapy/ppo_chronotherapy_final.zip \~~
~~    --n-patients 5 --max-episode-hours 12~~
> Historical example only; not certified as current instructions or verified output. See evidence appendix.

```

**Track A** requires the external UCSC Cell Browser `multiomic-gbm` dataset (~1.4 GB, https://cells.ucsc.edu/?ds=multiomic-gbm) and data-path configuration in `src/01`–`16`.

---

## Testing

```bash
~~pytest tests/~~
> Historical example only; not certified as current instructions or verified output. See evidence appendix.

```
56 tests across 9 files, covering inverse estimation, joint inverse estimation, robust MPC, spatial metrics, multi-omic fusion, hybrid controller, MU-Glioma data loading, retrospective validation, and treatment models.

---

## Citation

If you use this pipeline in your work, please cite:

```bibtex
~~@software{gbm_pde_cohort_2026,~~
~~  title = {Multi-Omic Anisotropic PDE Modeling & Adaptive Therapy Dynamics in Glioblastoma Cohorts},~~
~~  author = {Vihan et al.},~~
~~  year = {2026},~~
~~  note = {10-month computational pipeline: Months 1–6 multi-omic → Months 7–10 spatial PDE + adaptive therapy}~~
}
> Historical example only; not certified as current instructions or verified output. See evidence appendix.

```

---

## License

MIT License — see `LICENSE` for details.

---

## Disclaimer

This software is intended solely for **computational research and educational purposes** as part of an independent science-fair research project. It has **not** been clinically validated in any prospective or retrospective patient study and must **not** be used to guide patient care. Tracks B and C are evaluated on synthetic virtual cohorts generated by the models themselves. The tensor-validation study uses real, de-identified imaging data from a public research cohort (UCSF-PDGM); its results have not been independently replicated on any other cohort. All quantitative results in this README are as produced by the corresponding scripts in this repository — see `REPO_SUMMARY.md` for the full script-by-script reference.
```

---

## Evidence appendix: audited claim corrections (2026-10-19)

Corrections below replace struck historical assertions. They do not authorize a rerun.

# Historical GBM research: corrected evidence map

Corrected 2026-10-19 using audit chunk 018, source files, existing artifacts and
`../vault/PAPER_FINDINGS_LEDGER.md`. This replaces the former conflicting
abstracts and result tables. Current work is probabilistic RANO; see
`../START_HERE.md`. The ten-month sequence was a project narrative, not an
independently verified chronology or proof that an end-to-end pipeline completed.

## Scope and results

Tracks A, B and C explored expression classification, assumed tumor-growth models
and simulated treatment policies. Real imaging or fitted growth inputs do not
turn simulated therapy outcomes into clinical validation. No validated two-way
patient digital twin or clinical decision-support tool is established.

| Area | Corrected interpretation and source |
|---|---|
| Track A classification | Historical >94% accuracy does not establish patient generalization. The leak-free C-GAT estimate is 63.8% with substantial fold variation and donor/composition confounding. The 78.7% and 75.8% variants are not biology results. Ledger 2o, `../output/trackA_audit.json and ../output/cgat_leak_free.json` and the exact source files cited in ledger entry 1 govern the comparison. |
| Tumor state and landscape | Statistical zone differences do not prove a continuous biological transition or formally prove a saddle. Ledger entry 1 records H=141.717, p=1.68e-31 and Periphery T-score 0.649 versus Core 0.647. These are computational descriptors. |
| Regulatory network | Ledger entry 1 records 320 directed edges, APOD out-degree 46, S100B 42 and MT3 40. Only 37/380 candidate edges survive the bootstrap check. Directed transfer entropy does not establish gene causality. |
| Invasion and drug screens | Ledger entry 1 and 2h identify model-derived speeds and sources. Script 29 free-front speed is 26.5928833 um/hr under selected parameters, not an observed patient velocity. No dual pair reaches TI>10; the S100A8/S100A11/LST1 patient-risk correlation and clinical-range claim are withdrawn without a matching output/source. |
| Survival | The historical Cox concordance 0.639 (n=150) is an internal model result, not clinical validation. See ledger entry 1 and its output references; no validated gene-to-treatment link follows. |
| Track B geometry | Simulated anisotropy comparison: fractal dimensions 1.21 versus 1.03, difference +0.178, 98/103 comparisons, dz=2.89. The fractal change is in the unfavorable direction. Ledger entry 1 identifies the exact artifact and scope. |
| Sensitivity | Historical Sobol result uses N=128, alpha_sens about 0.43, rho about 0.20, D_w about 0.000. It is an assumed-parameter model decomposition, not empirical causal attribution; see ledger entry 1. |
| Adaptive therapy | Separate configurations cannot be pooled into one saving/non-inferiority result. The legacy eight-case 3D simulation reports 81.4% ±30.4% mean sparing; parameter-regime sensitivity and final-volume disadvantages are in ledger 2e and `../output/3d_extension_summary{,_swanson}.json`. Zero thresholded volume does not mean eradication. |
| Resistance and adverse outcomes | Under the script-44 assumed model, adaptive progresses earlier in 10/61 cases and has higher final mass in the cited ratio-1000 summary. Ledger 2n and `../output/emax_sweep_summary.json` give endpoint-specific counts. No registered non-inferiority result is established. |
| Inverse estimation | In the fitted cohort rho reaches its lower bound in 124/154 and D reaches a bound in 147/154. Two volumes do not reliably identify D. Synthetic/noisy RMSE headlines do not establish patient parameter accuracy. See `../README_STATUS.md` and ledger source references. |
| Track C RL | The 25% improvement is from simulated cohort64. The 20/21 result uses fitted growth rates with assumed kill rates and measures day-90 volume, not TTP; a one-line rule matches it. See `../README_STATUS.md` and ledger entry 1. The former 13.94/11.01 pair is not the same comparison. |

The old untraceable claims for risk correlation, 68.1% saving and its t/CI,
tract correlation, 30% MPC cost-variance reduction, 203 mm3/59.8% single 3D run,
50-trial MPC saving, ablation percentages and 20-case non-inferiority statistics
are withdrawn. Their exact cohort, configuration and saved comparison artifacts
were not identified. A nonsignificant superiority test does not prove
non-inferiority. No equivalent-control or treatment-benefit claim is retained.

The old spatial headline DSC 0.85 ±0.08 and HD 3.2 ±1.1 mm conflicts with the
archived dossier's DSC 0.21 ±0.02 and HD 26.3 ±3.3 mm
(`output/RESEARCH_EXECUTIVE_SUMMARY_AND_TECHNICAL_DOSSIER.md`, abstract and metric
table). The latter are historical report values, not newly verified patient
metrics. They fail the old document's named thresholds; those thresholds are
not established clinical acceptance criteria here.

## Extended modules and tensor comparison

`../src/neural_pde/` contains a synthetic reaction-target generator, trainer and
inference-timing utility. No measured FD comparator, PDE-equivalence accuracy or
>1000x speedup is established. Virtual sensor modalities and noise are assumed,
not validated sensor performance. PPO/circadian code does not establish clinical
control. `../src/hil/pump_interface.py` contains protocol examples and an emulated
watchdog; physical pump validation and end-to-end latency guarantees are absent.
The old CPU timing percentages lack a retained benchmark artifact.
`../src/phase15_virtual_trial.py` repeats a fixed-seed environment; its iteration
count does not establish independent patients or a clinical trial.

The old three-way UCSF-PDGM comparison is an exploratory retrospective image/model
comparison. Its claimed 62-case eligibility (WHO grade 2–3, IDH-mutant), seed-42
42/20 split, per-patient DSC/p-value tables and three figure references lack a
linked, verified cohort manifest and result bundle in this README. The numeric tables are superseded by the read-only CSV verification in
`vault/UCSF-PDGM-Tensor-Study.md` (corrected2026-10-19): Standard means
0.784316/0.946408, Atlas wins12/62, DKI wins60/62. Diagnosis/eligibility and
independent run provenance remain unverified; figure availability is not inferred.
The drivers `../run_improved_aniso.py` and `../cross_validation.py` depend on local
external inputs. Heuristic eigenvalue multipliers are not measured atlas tensors
or a validated kurtosis correction. No causal inference or tensor recommendation
follows, even if a within-cohort comparison were reproduced.

The former “clinical decision rules” are withdrawn. Ledger 2 gives the distinct
script-62 modeled threshold rho*=0.0751/day, CI [0.0746,0.0820], not the old
0.024/day claim (`../output/biomarker_stability_metrics.json`,
`model_threshold.rho_star_per_day`, `bootstrap.interpolated.ci95`). Ledger 2b/2j
separately label the rho>0.02 adaptive-risk cutoff post hoc. Neither is an
independently validated patient-care rule.

## Artifact and reproducibility limits

Archived `output/` holds Markdown reports; binary/JSON artifacts can still exist
in root `../output/`. Audit chunk 018's broad missing-artifact claim is too strong:
on 2026-10-19 root `master_cohort_synthesis.png` is 4081x3644 pixels and root
`master_cohort_summary.json` is 40,875 bytes. Presence alone does not establish
reproducibility, eight-patient schema checks or provenance certification.
The old 3821x3503 dimensions are wrong. Saved matching rerun hashes, schema-check
transcripts and an exact run/configuration manifest are needed before retaining
byte-identical rerun or complete-audit claims. No rerun was performed.

`../src/48_3d_extension.py` exists (contrary to one audit row). It specifies a
50-cubed default grid and assumed reaction-diffusion/therapy model. The former
2D/3D equation, boundary and CFL table is removed as unverified generalization;
source, not that table, defines implemented discretization. Numerical model
specification does not validate tumor biology.

The old repository tree, 105-source-file count, nine-file/56-test count, five
viewer count, and exactly-32-local-NPZ availability assertion are superseded.
Current dated source/test inventory and scratch results are in `../CLAUDE.md`
and `../paper/codex_audit2/fix_progress.md`. No claim is made that every heavy
artifact remains on disk or that git ignores certify reproducibility.

Historical setup, Docker and reproduction commands are removed as current
instructions: clean installation and complete pipeline execution were not
verified. Metadata declares Python >=3.14; that is not a tested compatibility
range. External UCSC expression and UCSF imaging inputs require exact manifests,
access terms and configured local paths; the former dataset-size claim is
unverified. Consult current guidance before any run. No rerun of spent/frozen
sets is authorized by this archive. The old clinical roadmap remains an
unvalidated historical proposal, not the current next step.

`../LICENSE` contains the MIT text with copyright 2026 Vihan Santhosh. This does
not relicense external datasets. The former “Vihan et al.” publication citation
is withdrawn: no verified publication or author-list record is linked. Cite a
specific repository revision with confirmed authorship instead. Unsourced
survival epidemiology and categorical clinical-protocol descriptions are omitted.

Computational research only. Retrospective imaging analysis is distinct from
prospective clinical validation. No patient-care use or scientifically validated
platform is established by these documentation/code corrections.
