# Multi-Scale Spatial Oncology Suite (MSOS) — 5-Month Master Plan

**Project:** Science Fair 2026 / ISEF Grand Award Track  
**Dataset:** 15,000 Single-Cell Spatial Matrix (Glioblastoma: Healthy, Periphery, Core)  
**Validated Foundation:** VAE (32D latent), NMF Module State Mapper, C-GAT, Trajectory Framework $\mathcal{T}_i$ ($p = 1.17 \times 10^{-31}$)  
**Timeline:** 5 Months | **Start Date:** July 15, 2026

---

## 🎯 Executive Mission

Build a **production-grade, self-contained predictive software ecosystem** that moves from descriptive modeling → autonomous clinical forecasting. Every module must be mathematically rigorous, computationally efficient (vectorized NumPy/SciPy/PyTorch), and publication-ready.

---

## 🗓️ Month-by-Month Architecture

### MONTH 1: Biophysical Velocity Fields & Waddington Landscape Physics
**Goal:** Map the explicit trajectory vector fields of cellular fate decisions; derive the Waddington energy landscape from the validated VAE latent space.

#### Mathematical Engine
- Treat 32D VAE latent space as continuous fluid dynamical system
- Compute localized **Drift $A(z)$** and **Diffusion $B(z)$** coefficients per cell
- Solve **Fokker-Planck Equation** → derive **Energy Potential $U(z) = -\ln P_{\text{steady}}(z)$**
- Prove Periphery = unstable thermodynamic saddle point (high energy, high entropy)

#### Deliverables
| Script | Output | Validation |
|--------|--------|------------|
| `src/19_phenotypic_velocity.py` | `output/tumor_phenotypic_flux.png` (quiver plot) | Visual vector field over UMAP |
| `src/20_fokker_planck_solver.py` | `output/waddington_landscape.npy`, `output/energy_potential.png` | 2D energy contour map |
| `src/21_drift_diffusion_analysis.py` | `output/drift_vectors.npy`, `output/diffusion_tensors.npy` | Per-cell $A(z), B(z)$ |
| `src/22_saddle_point_proof.py` | `output/saddle_point_metrics.json` | Hessian eigenvalues at Periphery centroid |

#### Key Equations
```
Drift: A(z) = lim_{Δt→0} E[Δz | z] / Δt
Diffusion: B(z) = lim_{Δt→0} E[Δz Δz^T | z] / Δt
Fokker-Planck: ∂ₜP = -∇·(A P) + ½ ∇²:(B P)
Steady-state: U(z) = -ln P_ss(z)
```

#### Success Criteria
- Vector field shows coherent flow Healthy → Periphery → Core
- Energy landscape reveals Periphery as saddle (1 positive, 31 negative Hessian eigenvalues)
- All computations vectorized (no Python loops over 15k cells)

---

### MONTH 2: Information-Theoretic Gene Regulatory Networks (GRNs)
**Goal:** Replace correlation with causality. Compute Directed Transfer Entropy to identify master regulator TFs driving the Periphery transition.

#### Mathematical Engine
- **Transfer Entropy** $\mathcal{T}_{X \to Y}$ between TFs ($X$) and functional programs ($Y$) across spatial gradient
- **Directed Information** for causal wiring diagram
- **Partial Information Decomposition (PID)** to separate unique/synergistic/redundant information

#### Deliverables
| Script | Output | Validation |
|--------|--------|------------|
| `src/23_transfer_entropy_engine.py` | `output/te_matrix.npy` (TF × Target) | $k$-NN Kraskov estimator |
| `src/24_causal_grn_builder.py` | `output/causal_grn.graphml`, `output/master_switches.tsv` | Top 10 drivers by out-degree |
| `src/25_pid_analysis.py` | `output/pid_decomposition.npy` | Unique info > Redundant info |
| `src/26_grn_validation.py` | `output/grn_bootstrap_ci.json` | 1000 bootstrap resamples |

#### Key Equations
```
Transfer Entropy: T_{X→Y} = Σ p(y_{t+1}, y_t, x_t) log₂ [p(y_{t+1}|y_t,x_t) / p(y_{t+1}|y_t)]
Directed Info: I(Xⁿ→Yⁿ) = Σ I(Xᵢ; Yᵢ | Y^{i-1})
PID: I(X₁,X₂;Y) = Unq₁ + Unq₂ + Red + Syn
```

#### Success Criteria
- Master switches: ≤ 5 TFs controlling >80% of Periphery program variance
- Bootstrap 95% CI excludes 0 for top edges
- Causal direction Periphery-specific (not in Healthy/Core)

---

### MONTH 3: Stochastic Agent-Based Tissue Invasion Engine
**Goal:** Convert 15,000 cells → active agents on 2D/3D lattice; simulate tumor invasion over time with morphogen PDEs.

#### Mathematical Engine
- **Cellular Automaton** on 500×500 (2D) or 200×200×50 (3D) lattice
- **Agent Rules** derived from empirical graph metrics (C-GAT edge weights, transition scores)
- **Fisher-Kolmogorov PDE** for morphogen diffusion: $\frac{\partial \rho}{\partial t} = D \nabla^2 \rho + r \rho (1 - \rho)$
- **Stochastic Gillespie Algorithm** for reaction-diffusion events

#### Deliverables
| Script | Output | Validation |
|--------|--------|------------|
| `src/27_aba_lattice.py` | `output/aba_state_*.npy` (time series) | Cell counts match analytical FK solution |
| `src/28_fisher_kolmogorov_pde.py` | `output/morphogen_field_*.npy` | Wave speed $c = 2\sqrt{Dr}$ verified |
| `src/29_invasion_simulator.py` | `output/invasion_movie.mp4`, `output/invasion_metrics.json` | Front velocity vs. empirical data |
| `src/30_aba_analysis.py` | `output/invasion_clinical_correlation.tsv` | Correlation with patient survival |

#### Agent State Machine
```
HEALTHY: P(proliferate) = 0.01, P(differentiate) = 0.05
PERIPHERY: P(invade) = 0.15, P(transition→CORE) = 0.08, P(secrete_morphogen) = 0.3
CORE: P(proliferate) = 0.25, P(necrose) = 0.02
```

#### Success Criteria
- Simulated invasion front velocity matches analytical $2\sqrt{Dr}$ within 10%
- Spatial heterogeneity emerges (fingering, clustering) matching histology
- Runtime: < 5 min for 30 simulated days on 500×500 lattice

---

### MONTH 4: In Silico Combinatorial Drug Gating (Therapeutic Discovery)
**Goal:** Algorithmic screening platform for virtual single/dual knockouts; compute Network Collapse Score $\mathcal{C}$.

#### Mathematical Engine
- Virtual knockouts: zero target columns in expression matrix → re-run CSGT/C-GAT
- **Network Collapse Score:** $\mathcal{C} = 1 - \frac{\text{Tr}(\Sigma_{\text{perturbed}})}{\text{Tr}(\Sigma_{\text{baseline}})}$
- **Synergy Metric:** Bliss independence / Loewe additivity for dual combos
- **Therapeutic Index:** $\mathcal{TI} = \frac{\mathcal{C}_{\text{tumor}}}{\mathcal{C}_{\text{healthy}}}$

#### Deliverables
| Script | Output | Validation |
|--------|--------|------------|
| `src/31_virtual_knockout_engine.py` | `output/single_ko_scores.tsv` (all genes) | Top 20: $\mathcal{C} > 0.7$ |
| `src/32_combinatorial_screen.py` | `output/dual_ko_synergy_matrix.npy` | 500×500 matrix, FDR < 0.05 |
| `src/33_therapeutic_index.py` | `output/therapeutic_index_ranking.tsv` | TI > 5 for lead combos |
| `src/34_drug_gating_report.py` | `output/drug_gating_report.pdf` | Publication-ready figures |

#### Screening Space
- **Single KO:** ~2,000 TFs + signaling genes (from GRN Month 2)
- **Dual KO:** Top 100 singles → 4,950 pairs (feasible)
- **Readout:** $\mathcal{T}_i$ shift, velocity field collapse, energy landscape flattening

#### Success Criteria
- Identify ≥ 3 dual-target combos with $\mathcal{C} > 0.85$ and $\mathcal{TI} > 10$
- Healthy centroid structural integrity preserved (PCA Procrustes distance < 0.1)
- Results reproducible across 10 bootstrap resamples

---

### MONTH 5: Clinical Generalization & Manuscript Compilation
**Goal:** Validate MSOS on independent human cohorts (Ivy GAP, TCGA-GBM); produce publication-ready manuscript + public repo.

#### Mathematical Engine
- **Domain Adaptation:** MMD / CORAL alignment of latent spaces across datasets
- **Survival Analysis:** Cox PH model with MSOS-derived features (velocity magnitude, $\mathcal{C}$, master switch expression)
- **Cross-Cohort GRN Conservation:** Jaccard similarity of top regulatory edges

#### Deliverables
| Script | Output | Validation |
|--------|--------|------------|
| `src/35_ivygap_ingest.py` | `output/ivygap_latent.npy`, `output/ivygap_metrics.json` | 500+ spatial samples |
| `src/36_tcga_validation.py` | `output/tcga_survival_analysis.tsv` | HR > 2, p < 0.001 for MSOS risk score |
| `src/37_cross_cohort_grn.py` | `output/grn_conservation.json` | Jaccard > 0.6 for top 50 edges |
| `src/38_manuscript_compiler.py` | `manuscript/main.tex`, `manuscript/figures/`, `manuscript/supplementary/` | Compiles to PDF |
| `src/39_repo_packager.py` | `MSOS/` (public repo structure) | `pip install -e .` works |

#### Manuscript Structure (Nature Methods / Cell Systems format)
1. **Title:** Multi-Scale Spatial Oncology Suite: From Biophysical Velocity Fields to In Silico Combinatorial Therapy
2. **Abstract:** 250 words, 4 key results
3. **Main Figures:** 6 (Velocity, Landscape, GRN, Invasion, Drug Screen, Clinical Validation)
4. **Supplementary:** 12 figures, 3 tables, mathematical derivations
5. **Code Availability:** GitHub + Zenodo DOI

#### Success Criteria
- MSOS risk score stratifies TCGA-GBM (n=150+) into high/low risk (log-rank p < 10⁻⁵)
- Ivy GAP spatial patterns recapitulate Periphery saddle point (energy correlation > 0.7)
- Manuscript compiles without errors; repo passes `pytest`, `ruff`, `mypy`
- All 38 scripts documented with docstrings, type hints, unit tests

---

## 🏗️ Technical Infrastructure Standards

### Code Quality (Enforced via CI)
```bash
# Required in every script
#!/usr/bin/env python3
"""Module docstring with mathematical basis."""
from __future__ import annotations
import numpy as np
import numpy.typing as npt
from typing import Tuple, Dict, List, Optional
```

### Data Contracts
| Artifact | Format | Schema |
|----------|--------|--------|
| Latent spaces | `.npy` | `(N, 32)` float32 |
| Transition scores | `.npy` | `(N,)` float32 |
| Graphs | `.graphml` | NetworkX compatible |
| Metrics | `.json` | Strict schema (pydantic) |
| Figures | `.png` (300 DPI) + `.svg` | Publication ready |

### Computational Requirements
- **Vectorization:** No Python loops over cells (use `scipy.spatial.cKDTree`, `sklearn.neighbors`, broadcasting)
- **Memory:** < 8 GB peak (stream large matrices via `mmap_mode='r'`)
- **Parallelism:** `joblib.Parallel` for bootstrap/knockout loops
- **GPU:** PyTorch for C-GAT retraining (Month 4); CPU-only for Months 1-3,5

### Version Control
- **Branching:** `main` (stable), `month-{1..5}-*` (feature), `manuscript-*` (writing)
- **Commits:** Conventional Commits (`feat:`, `fix:`, `docs:`, `refactor:`)
- **Tags:** `v1.0-month1`, `v2.0-month2`, ..., `v5.0-final`

---

## 📊 Milestone Tracking Dashboard

| Month | Week | Primary Deliverable | Status | Blockers |
|-------|------|---------------------|--------|----------|
| 1 | 1-2 | `19_phenotypic_velocity.py` → flux plot | 🟢 Ready to run | — |
| 1 | 2-3 | `20_fokker_planck_solver.py` | ⬜ | Need drift/diffusion validation |
| 1 | 3-4 | `21_drift_diffusion_analysis.py` + `22_saddle_point_proof.py` | ⬜ | Hessian computation at scale |
| 2 | 1-2 | `23_transfer_entropy_engine.py` (Kraskov) | ⬜ | k-NN MI estimator choice |
| 2 | 2-3 | `24_causal_grn_builder.py` + `25_pid_analysis.py` | ⬜ | PID library (dit / idtxl) |
| 2 | 3-4 | `26_grn_validation.py` (bootstrap) | ⬜ | Compute time (1000×15k) |
| 3 | 1-2 | `27_aba_lattice.py` + `28_fisher_kolmogorov_pde.py` | ⬜ | Lattice size vs. resolution |
| 3 | 2-3 | `29_invasion_simulator.py` (video) | ⬜ | Matplotlib animation perf |
| 3 | 3-4 | `30_aba_analysis.py` (clinical corr) | ⬜ | Survival data access |
| 4 | 1-2 | `31_virtual_knockout_engine.py` | ⬜ | CSGT re-run speed |
| 4 | 2-3 | `32_combinatorial_screen.py` (4950 pairs) | ⬜ | Parallelization strategy |
| 4 | 3-4 | `33_therapeutic_index.py` + `34_drug_gating_report.py` | ⬜ | Healthy baseline definition |
| 5 | 1-2 | `35_ivygap_ingest.py` + `36_tcga_validation.py` | ⬜ | Data download / MTA |
| 5 | 2-3 | `37_cross_cohort_grn.py` | ⬜ | Orthology mapping |
| 5 | 3-4 | `38_manuscript_compiler.py` + `39_repo_packager.py` | ⬜ | LaTeX template, CI setup |

---

## 🚨 Risk Mitigation

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Fokker-Planck numerical instability | Medium | High | Implicit Euler scheme; validate against known potentials |
| Transfer Entropy compute time (O(N²)) | High | Medium | Subsample 5k cells for TE; use GPU k-NN (RAPIDS/cuML) |
| ABA lattice memory (500³ × 30 days) | Medium | High | Sparse representation; chunked time series |
| TCGA/Ivy GAP data access delays | Medium | High | Start MTA Week 1 Month 5; use public subsets first |
| Manuscript scope creep | High | Medium | Freeze figure list Week 2 Month 5 |

---

## 🔬 Immediate Next Step (Month 1, Week 1)

**Execute:** `python src/19_phenotypic_velocity.py`

**Expected Output:** `output/tumor_phenotypic_flux.png` — quiver plot showing directional velocity vectors over UMAP, colored by spatial zone (Healthy/Periphery/Core).

**Validation Checklist:**
- [ ] Arrows point Healthy → Periphery → Core
- [ ] Periphery shows highest vector magnitude (transition zone)
- [ ] No NaN/inf in velocity vectors
- [ ] Runtime < 30 seconds (vectorized KDTree)

---

## 📝 Plan Governance

- **Review Cadence:** Weekly (Sunday 20:00) — update dashboard, adjust next week
- **Decision Log:** Append to `DECISIONS.md` for all architectural choices
- **Retrospective:** End of each month — what worked, what didn't, velocity adjustment
- **Plan Exit:** Call `plan_exit` when Month 1 deliverables pass validation

---

**Author:** Kilo (Computational Oncology Research Lab)  
**Version:** 1.0 | **Date:** 2026-07-15  
**Status:** **ACTIVE** — Month 1 execution begins now.