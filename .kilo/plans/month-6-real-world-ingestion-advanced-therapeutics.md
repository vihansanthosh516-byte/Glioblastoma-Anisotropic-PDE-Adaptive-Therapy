# Month 6: Real-World Ingestion & Advanced Therapeutics

## Project Context
- **Months 1-3**: Spatial transcriptomics, GRN inference, cVAE training, virtual KO engine
- **Month 4**: Combinatorial dual-KO screen (scripts 32-34), calibrated therapeutic indices
- **Month 5**: Mock IvyGAP cohort survival analysis (scripts 35-37), KM/Cox/forest plots

## Current Artifacts (ready for Month 6)
| Artifact | Path | Description |
|----------|------|-------------|
| cVAE encoder | `output/cgat/cvae_model.pt` | Trained 2500→32 latent encoder |
| Latent space | `output/scvi_latent.npy` | 15000×32 latent coordinates |
| Expression matrix | `output/nn_X.npy` | 15000×2500 log-normalized counts |
| Gene names | `output/te_gene_names.txt` | 2500 gene symbols |
| Zone labels | `output/nn_y.npy` | 0=Healthy, 1=Periphery, 2=Core |
| Dual-KO TI | `output/dual_ko_ti.json` | 6 pairs with calibrated TI |
| Single-KO results | `output/single_ko_results.json` | 7 genes × collapse scores |
| Survival summary | `output/survival_stats_summary.json` | Mock cohort stats |

---

## Week 1: Authentic Multi-Omic Patient Ingestion (`src/38_real_cohort_ingest.py`)

### Objective
Replace mock cohort with genuine IvyGAP anatomic structure RNA-seq + TCGA-GBM clinical survival (n > 500).

### Data Sources
| Dataset | Access | Key Files | Target |
|---------|--------|-----------|--------|
| **IvyGAP** | Allen Institute API / AWS S3 | `ivygap_rnaseq.h5ad`, `ivygap_clinical.tsv` | ~100 samples × anatomic structures (Leading Edge, Infiltrating Tumor, Cellular Tumor, Pseudopalisading, Microvascular Proliferation) |
| **TCGA-GBM** | GDC API / `TCGAbiolinks` | `gbm_counts.tsv`, `gbm_clinical.tsv` | ~600 patients × 20k genes × OS/PFS + age/sex/IDH/subtype |

### Implementation Plan

#### 1. IvyGAP Ingestion Module
```python
def fetch_ivygap_rnaseq() -> AnnData:
    """Download IvyGAP RNA-seq (anatomic structures) via Allen SDK or S3."""
    # Use alleninstitute/allensdk or direct S3: s3://allen-brain-atlas/ivygap/
    # Returns AnnData with obs: structure_acronym, structure_name, donor_id
```

#### 2. TCGA-GBM Ingestion Module
```python
def fetch_tcga_gbm_cohort() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Download TCGA-GBM via GDC API or TCGAbiolinks (R) → parquet cache."""
    # Returns (counts_df: genes × patients, clinical_df: patient × covariates)
    # Filter: primary solid tumor, IDH-wt, WHO grade 4
```

#### 3. Harmonization & Mapping
```python
def harmonize_cohorts(ivygap_adata, tcga_counts, tcga_clinical, cvae_encoder) -> Dict:
    """
    - Map IvyGAP anatomic structures → Month 4 zones (LE→Periphery, CT→Core, etc.)
    - Intersect gene symbols with cVAE's 2500-gene vocabulary
    - Project IvyGAP/TCGA expression through cVAE encoder → latent coordinates
    - Build unified survival table: patient_id, time, event, age, sex, subtype, latent_32D
    """
```

#### 4. Feature Engineering for Survival
```python
def build_survival_features(unified_latent, clinical) -> pd.DataFrame:
    """
    - Latent dimensions (32) as continuous features
    - Zone-specific latent means (from IvyGAP structure mapping)
    - Clinical covariates: age, sex, Karnofsky, IDH status, MGMT methylation
    - cVAE-derived "tumor collapse score" per patient (from Month 4)
    """
```

### Outputs
| File | Description |
|------|-------------|
| `output/real_ivygap_cohort.h5ad` | IvyGAP RNA-seq + latent projections |
| `output/real_tcga_gbm_cohort.parquet` | TCGA counts + clinical (n≈600) |
| `output/unified_survival_cohort.csv` | Harmonized survival table (n>500) |
| `output/real_cohort_manifest.json` | Metadata: sources, versions, QC stats |

### Dependencies
- `allensdk`, `scanpy`, `anndata`, `pandas`, `numpy`, `requests`, `gdc-client` (or `TCGAbiolinks` via `rpy2`)

---

## Week 2: Penalized Cox Modeling & Feature Selection (`src/39_penalized_survival.py`)

### Objective
Scale from 4 genes → genome-wide (2500 cVAE genes + 32 latent dims) using ElasticNet/Lasso Cox PH with nested cross-validation.

### Methodology

#### 1. Penalized Cox Implementation
```python
class PenalizedCox:
    """
    ElasticNet Cox PH via coordinate descent (scikit-survival or custom).
    Loss: -partial_likelihood + α * [(1-l1_ratio)||β||₂²/2 + l1_ratio||β||₁]
    """
    def fit(self, X, time, event, alpha, l1_ratio, cv_folds=5):
        # Returns: coefficients, HRs, selected features, cv_cindex
```

#### 2. Nested Cross-Validation
```python
def nested_cv_penalized_cox(X, time, event, param_grid, outer_folds=5, inner_folds=5):
    """
    Outer loop: unbiased C-index estimation
    Inner loop: hyperparameter tuning (alpha, l1_ratio)
    Returns: best_params, outer_cindex_mean, outer_cindex_std, feature_frequency
    """
```

#### 3. Stability Selection
```python
def stability_selection(X, time, event, n_bootstraps=100, subsample_frac=0.5):
    """
    Subsample patients → fit penalized Cox → record selected features
    Returns: selection_probability per feature (π_j)
    Threshold π_j > 0.7 for "stable" features
    """
```

#### 4. Post-Selection Inference
```python
def post_selection_inference(selected_features, X, time, event):
    """
    Refit unpenalized Cox on stable features only
    Return: HR, 95% CI, p-value, Schoenfeld residuals (PH assumption check)
    """
```

### Feature Sets to Test
| Set | Size | Source |
|-----|------|--------|
| Latent only | 32 | cVAE encoder |
| Top 100 HVG | 100 | scVI highly variable genes |
| KO-altered genes | 50-200 | Month 4 single/dual KO collapse scores |
| Full cVAE vocab | 2500 | All encoder input genes |
| Combined | ~2600 | Latent + genes + clinical |

### Outputs
| File | Description |
|------|-------------|
| `output/penalized_cox_results.json` | Best params, CV C-index, selected features, HRs |
| `output/stability_selection.csv` | Selection probability per feature |
| `output/final_cox_model.pkl` | Fitted `PenalizedCox` object |
| `output/feature_importance.png` | Coefficient path plot |

### Dependencies
- `scikit-survival` (for `CoxnetSurvivalAnalysis`), `scikit-learn`, `lifelines`, `joblib`

---

## Week 3: Spatial Recurrence & Microenvironment Mapping (`src/40_spatial_recurrence_mapper.py`)

### Objective
Map Fisher-Kolmogorov PDE invasion front onto physical brain geometry; predict high-risk recurrence zones.

### Mathematical Framework

#### 1. Fisher-Kolmogorov Equation (from Month 3, script 28)
```python
∂u/∂t = D ∇²u + ρ u (1 - u/K)
```
- `u(x,t)`: tumor cell density
- `D`: diffusion tensor (DTI-derived, anisotropic)
- `ρ`: proliferation rate (from cVAE latent or Ki-67)
- `K`: carrying capacity

#### 2. Wave Front Velocity
```
v = 2 √(D ρ)  (isotropic)
v(x) = 2 √(D(x) ρ(x))  (anisotropic, spatially varying)
```

#### 3. Patient-Specific Parameterization
```python
def estimate_patient_parameters(latent_32d, dti_tensor, clinical):
    """
    - ρ ← decoder(latent) → proliferation gene module score
    - D ← DTI principal eigenvector × scalar diffusivity
    - K ← tissue cellularity (from histology or MRI)
    """
```

#### 4. Brain Geometry & Mesh
```python
def build_brain_mesh(atlas="MNI152", resolution=1mm):
    """
    - Load MNI152 template + white/gray matter masks
    - Generate tetrahedral mesh (FEniCS/dolfinx or pygalmesh)
    - Register patient MRI → atlas space (ANTs/SyN)
    """
```

#### 5. Recurrence Probability Map
```python
def compute_recurrence_probability(mesh, v_field, resection_cavity, time_horizon=365):
    """
    - Solve FK-PDE forward from resection cavity boundary
    - Compute arrival time T(x) = min{t: u(x,t) > threshold}
    - Recurrence risk R(x) = exp(-λ * T(x))  (λ = detection sensitivity)
    - Return: 3D risk volume in patient space
    """
```

### Integration with Month 4 KO Results
```python
def map_ko_to_recurrence(dual_ko_ti, patient_latent, mesh):
    """
    For each dual-KO pair:
      - Simulate virtual KO → latent shift Δz
      - Decode Δz → Δρ (proliferation change)
      - Recompute v_field → ΔT(x) (arrival time shift)
      - Rank KOs by ΔR = ∫|ΔR(x)| dx over recurrence-prone zones
    """
```

### Outputs
| File | Description |
|------|-------------|
| `output/fk_pde_parameters.json` | Per-patient D, ρ, K estimates |
| `output/recurrence_risk_volume.nii.gz` | 3D NIfTI risk map (patient space) |
| `output/ko_recurrence_impact.csv` | Dual-KO ranked by recurrence delay |
| `output/spatial_recurrence_report.md` | Visualizations + clinical interpretation |

### Dependencies
- `fenics`/`dolfinx` or `fipy` for PDE solve
- `nibabel`, `nilearn`, `dipy` for neuroimaging
- `ANTsPy` for registration (optional, can use pre-registered)

---

## Week 4: Dose-Response Extensions & Final Gating Matrix (`src/41_dose_response_model.py`)

### Objective
Replace binary KO with continuous Hill equation pharmacology; produce clinical actionability report.

### Pharmacological Model

#### 1. Hill Equation for Single Agent
```
E(C) = E_max × Cⁿ / (EC₅₀ⁿ + Cⁿ)
```
- `C`: drug concentration (μM)
- `E_max`: max effect (from Month 4 tumor collapse C_max)
- `EC₅₀`: potency (fit from literature or PDX)
- `n`: Hill coefficient (typically 1-2)

#### 2. Dual-Drug Synergy (Bliss Independence)
```
E_comb(C₁, C₂) = E₁ + E₂ - E₁×E₂ + S × E₁ × E₂
```
- `S`: synergy coefficient (from Month 4 Bliss synergy)
- `S > 0`: synergy, `S < 0`: antagonism

#### 3. Therapeutic Index as Function of Dose
```
TI(C₁, C₂) = E_tumor(C₁, C₂) / E_healthy(C₁, C₂)
```
- `E_tumor`: collapse score on tumor zone (periphery+core)
- `E_healthy`: collapse score on healthy zone
- Both use Hill curves with tissue-specific EC₅₀

#### 4. PK/PD Integration (Optional)
```
C(t) = (Dose / Vd) × exp(-CL × t)  (1-compartment)
AUC = Dose / CL
```

### Optimization: Maximum Safety Window
```python
def optimize_dose_pair(ec50_tumor, ec50_healthy, emax, hill, synergy):
    """
    Maximize TI(C₁, C₂) subject to:
      - C₁ ≤ MTD₁, C₂ ≤ MTD₂ (maximum tolerated dose)
      - E_healthy(C₁, C₂) ≤ toxicity_threshold (e.g., 0.1)
    Returns: (C₁*, C₂*, TI_max, safety_margin)
    """
```

### Clinical Actionability Report
```python
def generate_gating_report(dual_ko_ti, dose_optimization, tcga_hrs):
    """
    For each dual-KO pair:
      1. Map gene targets → known drugs (DGIdb, DrugBank, ChEMBL)
      2. Pull literature EC₅₀, Hill, MTD for each drug
      3. Compute optimal dose pair (C₁*, C₂*)
      4. Rank by: TI_max × (1 - toxicity_risk) × clinical_HR_support
    Output: Ranked table with dosing recommendations
    """
```

### Outputs
| File | Description |
|------|-------------|
| `output/dose_response_curves.png` | Hill curves per drug/tissue |
| `output/dual_drug_isobolograms.png` | Combination surfaces |
| `output/optimal_dosing_table.csv` | C₁*, C₂*, TI_max per pair |
| `output/final_gating_matrix.csv` | Clinical actionability score |
| `output/month6_final_report.md` | Complete Month 6 synthesis |

### Dependencies
- `dgidb`/`chembl_webresource_client` for drug-target mapping
- `scipy.optimize` for constrained optimization
- `matplotlib`/`plotly` for 3D isobolograms

---

## Cross-Week Integration Points

| Week | Consumes | Produces |
|------|----------|----------|
| 1 | Month 4-5 artifacts, cVAE model | Real cohorts, unified survival table |
| 2 | Week 1 unified table + cVAE genes | Penalized Cox model, stable gene signature |
| 3 | Week 1 latent params + Month 28/29 PDE | Recurrence risk maps, KO→recurrence impact |
| 4 | Week 2 signature + Month 4 TI + Week 3 recurrence | Dose-optimized gating matrix, final report |

---

## Risks & Mitigations

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| IvyGAP/TCGA download fails | Medium | Pre-cache on S3; fallback to local parquet |
| DTI data unavailable for patients | High | Use population-average DTI atlas (HCP) |
| Penalized Cox convergence issues | Medium | Warm starts, path-wise coordinate descent |
| PDE solve too slow on mesh | Medium | Use Fipy (finite volume) on regular grid; GPU not needed |
| Drug-target mapping incomplete | High | Limit to targets with known inhibitors; flag others |

---

## Success Criteria (Month 6 Exit)

1. **Week 1**: `unified_survival_cohort.csv` with n ≥ 500, all covariates populated
2. **Week 2**: Penalized Cox C-index > 0.65 on outer CV; ≥ 5 stable features (π > 0.7)
3. **Week 3**: Recurrence risk volume generated for ≥ 3 patients; KO impact ranked
4. **Week 4**: ≥ 3 dual-drug pairs with TI > 2 at clinically achievable doses

---

## File Structure (New Scripts)

```
src/
├── 38_real_cohort_ingest.py      # Week 1
├── 39_penalized_survival.py      # Week 2
├── 40_spatial_recurrence_mapper.py  # Week 3
├── 41_dose_response_model.py     # Week 4
└── utils/
    ├── data_fetchers.py          # IvyGAP/TCGA download helpers
    ├── pde_solvers.py            # FK-PDE finite volume
    └── drug_db.py                # DGIdb/ChEMBL client
```

---

## Timeline

| Week | Dates | Primary Deliverable |
|------|-------|---------------------|
| 1 | Days 1-5 | `38_real_cohort_ingest.py` + real cohort artifacts |
| 2 | Days 6-10 | `39_penalized_survival.py` + gene signature |
| 3 | Days 11-15 | `40_spatial_recurrence_mapper.py` + risk maps |
| 4 | Days 16-20 | `41_dose_response_model.py` + final report |

---

## Notes for Implementation

1. **Data licensing**: IvyGAP (Allen Institute) - open access; TCGA-GBM - controlled access (dbGaP) but summary stats available via GDC open tier
2. **Compute**: PDE solves (Week 3) and nested CV (Week 2) are CPU-intensive; consider batching
3. **Reproducibility**: Pin all dependency versions; cache downloaded data with content hashes
4. **Validation**: Week 2 signature should be tested on held-out TCGA subset + independent cohort (e.g., CGGA)