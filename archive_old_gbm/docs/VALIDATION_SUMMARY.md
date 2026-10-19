> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

> **Audit corrections, 2026-10-19:** Original layout retained. Struck claims are withdrawn or superseded. The evidence appendix gives corrected values, sources and specific missing evidence. Unstruck historical claims retain their original scope.

# ~~MSOS Validation Summary — Final Report (Months 1-4)~~

~~**Project:** Multi-Scale Spatial Oncology Suite (MSOS)  ~~
~~**Dataset:** 15,000 Single-Cell Spatial Matrix (Glioblastoma: Healthy, Periphery, Core)  ~~
~~**Hardware:** NVIDIA RTX 4050 Laptop GPU (6.4 GB VRAM), CUDA 13.2  ~~
~~**Period:** July 15, 2026 — Completed in single session~~

---

## Execution Status Matrix

| ~~Script~~ | ~~Status~~ | ~~Runtime~~ | ~~Key Output~~ |
|--------|--------|---------|------------|
| ~~`19_phenotypic_velocity.py`~~ | ~~Complete~~ | ~~2.2s~~ | ~~`tumor_phenotypic_flux.png`~~ |
| ~~`20_fokker_planck_solver.py`~~ | ~~Complete~~ | ~~1.2s~~ | ~~`energy_potential.png`, `waddington_landscape.npy`~~ |
| ~~`21_drift_diffusion_analysis.py`~~ | ~~Complete~~ | ~~0.3s~~ | ~~`drift_vectors.npy`, `diffusion_tensors.npy`~~ |
| ~~`22_saddle_point_proof.py`~~ | ~~Complete~~ | ~~0.5s~~ | ~~`saddle_point_metrics.json` (NEB on Periphery)~~ |
| ~~`23_transfer_entropy_engine.py`~~ | ~~Complete~~ | ~~15 min~~ | ~~`te_matrix.npy` (100×100), 685 non-zero TE edges~~ |
| ~~`24_causal_grn_builder.py`~~ | ~~Complete~~ | ~~0.01s~~ | ~~`causal_grn.graphml`, `master_switches.tsv`~~ |
| ~~`25_pid_analysis.py`~~ | ~~Complete~~ | ~~<0.01s~~ | ~~`pid_decomposition.tsv` (simplified)~~ |
| ~~`26_grn_validation.py`~~ | ~~Complete~~ | ~~20s~~ | ~~`grn_bootstrap_ci.json` (32/380 edges sig)~~ |
| ~~`27_aba_lattice.py`~~ | ~~Complete~~ | ~~2.3s~~ | ~~`aba_grid_history.npy`, `aba_metrics.json`~~ |
| ~~`28_fisher_kolmogorov_pde.py`~~ | ~~Complete~~ | ~~1.0s~~ | ~~`fk_field_history.npy`, **31.6% error** (needs dt=0.005)~~ |
| ~~`29_invasion_simulator.py`~~ | ~~Complete~~ | ~~21s~~ | ~~`invasion_metrics.json` (wave speed 4.3 µm/hr)~~ |
| ~~`30_aba_analysis.py`~~ | ~~Complete~~ | ~~5s~~ | ~~`invasion_dynamics_analysis.png`~~ |
| ~~`31_virtual_knockout_engine.py`~~ | ~~Complete~~ | ~~8s~~ | ~~`single_ko_results.json` (200 genes)~~ |
| ~~`32_combinatorial_screen.py`~~ | ~~Complete~~ | ~~0.2s~~ | ~~`dual_ko_results.json` (15 pairs)~~ |
| ~~`33_therapeutic_index.py`~~ | ~~Complete~~ | ~~19s~~ | ~~`single_ko_ti.json` (4 genes), `dual_ko_ti.json` (15 pairs)~~ |
| ~~`34_drug_gating_report.py`~~ | ~~Complete~~ | ~~1.2s~~ | ~~`drug_gating_report.md`, `optimization_matrix.png`~~ |

~~**Total Scripts Executed:** 16/39 (Months 1-4 complete, Month 5 pending)~~

---

## Key Scientific Results

### Month 1: Biophysical Fields & Waddington Landscape
~~- **Phenotypic Velocity Field:** 15,000 cells × 32D latent space → quiver plot~~
~~- **Waddington Energy Landscape:** Dual attractors (Healthy=0.56, Core=0.00), **Periphery saddle at E=5.74** (NEB confirmed)~~
~~- **Drift/Diffusion Tensors:** 15,000 × 32×32 drift, 15,000 × 32×32 diffusion~~
~~- **Saddle Point Proof:** NEB on Periphery confirms mixed Hessian (±λ)~~

### ~~Month 2: Causal GRN & Master Switches~~
~~- **Transfer Entropy Matrix:** 100×100 genes, 685 significant directed edges (p<0.01)~~
~~- **Top Master Switches (Out-Degree):** APOD(46), S100B(45), MT3(40), S100A8(26), S100A9(23)~~
~~- **Bootstrap Validation:** 32/380 edges significant (95% CI > 0)~~

### ~~Month 3: Invasion Dynamics~~
~~- **FK PDE (ETDRK2, dt=0.0001):** Wave speed 4.3 µm/hr (clinical: 10-50 µm/hr) — **15.6% error**~~
~~- **Integrated CA+PDE (512×512):** 400 PDE steps × 10 CA sub-steps~~
~~- **Necrotic Fraction:** 6.7% (clinical: 10-40%) — needs `core_necrose=0.005`~~

### ~~Month 4: Drug Discovery & Therapeutic Index~~
| ~~Metric~~ | ~~Best Single KO~~ | ~~Best Dual KO~~ |
|--------|----------------|--------------|
| ~~**Tumor Collapse (C)**~~ | ~~TMLHE: 0.0378~~ | ~~MT-ATP6+S100A6: 1.0000~~ |
| ~~**Healthy Collapse (C)**~~ | ~~All: 1.0000~~ | ~~All: 1.0000~~ |
| ~~**Therapeutic Index (TI)**~~ | ~~MT-ATP6: 0.02~~ | ~~All dual: 1.00~~ |
| ~~**Bliss Synergy**~~ | ~~—~~ | ~~LST1+CCL3L1: -0.4359~~ |

~~**Clinical Reality Check:** No combination achieved TI > 10 with tumor_collapse > 0.05~~

---

## ~~Remaining Calibration Gaps (P0)~~

| ~~Issue~~ | ~~Current~~ | ~~Target~~ | ~~Fix Strategy~~ |
|-------|---------|--------|--------------|
| ~~**FK Wave Speed Error**~~ | ~~31.6%~~ | ~~<10%~~ | ~~Use dt=0.0001 + RK4 or implicit diffusion~~ |
| ~~**Necrotic Fraction**~~ | ~~6.7%~~ | ~~10-40%~~ | ~~Increase `core_necrose` to 0.005~~ |
| ~~**TI Clinical Threshold**~~ | ~~TI_max=1.0~~ | ~~TI > 10~~ | ~~Increase transition rates, add dose-response~~ |
| ~~**Healthy Zone Collapse**~~ | ~~C=1.0 (max)~~ | ~~C < 0.2~~ | ~~Fix encoder saturation, add homeostatic term~~ |

---

## ~~Final Artifact Inventory~~

> Historical example only; not certified as current instructions or verified output. See evidence appendix.

```
output/
├── Month 1: Biophysical Fields
│   ├── tumor_phenotypic_flux.png
│   ├── energy_potential.png
│   ├── waddington_landscape.npy
│   ├── drift_vectors.npy (15000, 32)
│   ├── diffusion_tensors.npy (15000, 32, 32)
│   └── saddle_point_metrics.json
├── Month 2: Causal GRN
│   ├── te_matrix.npy (100, 100)
│   ├── causal_grn.graphml
│   ├── master_switches.tsv
│   ├── pid_decomposition.tsv
│   ├── grn_bootstrap_ci.json
│   └── grn_metrics.json
├── Month 3: Invasion Engine
│   ├── fk_field_history.npy (11, 512, 512)
│   ├── fk_metrics.json (error: 31.6%)
│   ├── invasion_metrics.json
│   ├── invasion_dynamics_analysis.png
│   └── aba_analysis_results.json
├── Month 4: Drug Discovery
│   ├── single_ko_results.json (200 genes)
│   ├── single_ko_ti.json (4 genes)
│   ├── dual_ko_results.json (15 pairs)
│   ├── dual_ko_ti.json (15 pairs)
│   ├── optimization_matrix.npy (100×6)
│   ├── optimization_matrix.png
│   └── drug_gating_report.md
└── Pipeline Metadata
    ├── VALIDATION_SUMMARY.md (this file)
    └── BENCHMARK_SUMMARY.md
```

---

## ~~Next Phase: Month 5 — Clinical Validation & Manuscript~~

| ~~Script~~ | ~~Purpose~~ | ~~Priority~~ |
|--------|---------|----------|
| ~~`35_ivygap_ingest.py`~~ | ~~Load Ivy GAP spatial transcriptomics~~ | ~~P0~~ |
| ~~`36_tcga_validation.py`~~ | ~~TCGA-GBM survival analysis with MSOS risk score~~ | ~~P0~~ |
| ~~`37_cross_cohort_grn.py`~~ | ~~Validate master switches across cohorts~~ | ~~P1~~ |
| ~~`38_manuscript_compiler.py`~~ | ~~Auto-generate LaTeX manuscript~~ | ~~P1~~ |
| ~~`39_repo_packager.py`~~ | ~~GitHub release + Zenodo DOI + CI/CD~~ | ~~P1~~ |

~~**Estimated Timeline:** 1 week for Month 5 execution~~

---

~~**Pipeline Status:** **Months 1-4 Complete** | **Month 5 Ready to Launch**~~

---

## Evidence appendix: audited claim corrections (2026-10-19)

Corrections below replace struck historical assertions. They do not authorize a rerun.

# Historical model-output summary: evidence corrections

Corrected 2026-10-19. The old “Months 1–4 complete” matrix lists scripts 19–34
(16 entries), not a verified 16/39 execution manifest. Its July-15 single-session,
GPU/VRAM/CUDA details and runtime values lack saved environment/command logs and
are withdrawn. The data used are a 15000-cell subsample, not a complete spatial
matrix; full cohort dimensions need input metadata. These are legacy simulations,
not the current RANO goal (`../../START_HERE.md`).

The audit searched archived output only. Several corresponding root artifacts
do exist, so blanket “missing” claims would be wrong. Read-only root array
headers: phenotypic_velocity.npy (15000,32), drift_vectors.npy (15000,32),
diffusion_tensors.npy (15000,32,32), te_matrix.npy (100,100). Shapes do not prove
execution timing, biological validity or a complete reproduction chain.

| Root artifact (under ../../output/) | Corrected evidence and scope |
|---|---|
| `grn_metrics.json` | n_edges=320; inferred out-degrees APOD 46, S100B 42, MT3 40. “685 significant edges” and the old gene-degree list are not supported by this graph summary. No causal biology or unreported multiplicity control is inferred. |
| `grn_bootstrap_ci.json` | Present in root output. The old 32/380 count is superseded by ledger's historical 37/380; exact interval selection and resampling method must accompany any inferential use. |
| `fk_metrics.json` | Selected simulation: numerical speed 3.623743 versus analytical 4.0, error 9.40642%; configured 20000 steps, completed 20010. Old 15.6% and 31.6% estimates are unsupported for this saved artifact and not interchangeable. |
| `invasion_summary.json` | Selected coupled simulation: free speed 26.592883 versus analytical 28.284271 um/hr, 5.98% relative error. This is not a clinical wave-speed measurement. |
| `aba_metrics.json`, `aba_analysis_results.json` | Tumor-cell snapshot necrosis 22.2/30.4/38.1% at 90/120/150; separate whole-grid final fraction 17.19%. Old 6.7% is not a verified current calibration result. Denominators and timepoints matter. |
| `single_ko_results.json`, `dual_ko_results.json` | 200 single-gene rows, 21 dual pairs. Highest single collapse is TMLHE 0.0407984; no positive Bliss score in the dual output. These are latent-model perturbation scores, not drug response. |
| `single_ko_ti.json`, `dual_ko_ti.json` | MT-CO2 model score 4.819325; dual range -0.978599 to +0.442764. Floor-based score field named therapeutic_index is not a pharmacologic TI; old all-dual=1 and healthy-collapse=1 table is withdrawn. |

The old NEB energy/attractor assertions and “saddle proof” are not established by
this report. The distinct SPIB output has internal model-consistency checks,
not independent biological validation. See `../TRACK_A_RESULTS.md` for corrected
CSGT, latent/saddle and missing-provenance details.

The old final artifact tree is removed: its mix of exact files, intended outputs,
obsolete filenames and shape assertions did not establish completeness. Markdown
`drug_gating_report.md` exists under this archive's output directory; root output
contains the named machine-readable artifacts above. `invasion_metrics.json`
and old `saddle_point_metrics.json` must be individually checked before reuse;
neither is replaced by another run's similarly named summary. Artifact presence
alone is not a successful-execution or reproducibility certificate.

P0 proposals to tune wave speed, necrosis, collapse or a TI threshold were
historical modeling ideas. Their old numerical targets were not verified
measurements or clinical acceptance rules. No threshold tuning or rerun is
authorized here. The Month-5 one-week schedule, “Ready to Launch” and clinical
validation/manuscript plan are superseded as active next steps. Current ledger
and START_HERE govern; missing execution/provenance evidence is not regenerated.
