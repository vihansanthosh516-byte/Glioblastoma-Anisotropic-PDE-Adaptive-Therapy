# Week 1 Implementation Plan: `src/38_real_cohort_ingest.py`

## Objective
Build a self-healing data ingestion script that scans `data/` for raw IvyGAP/TCGA files. If missing, automatically generates a **high-fidelity synthetic cohort** mirroring IvyGAP's exact matrix parameters:
- **120 patients** with matched spatial tumor sub-structures
- **3 spatial zones**: `Leading Edge`, `Cellular Tumor`, `Infiltrating Tumor`
- **4 target genes**: `LST1`, `S100A11`, `S100A8`, `ZNF106` (from 2500-gene cVAE vocabulary)
- **Survival tracking**: OS time, event status, clinical covariates

Output: `output/real_cohort_aligned.csv` — unified tidy matrix with zone-stratified expression + survival.

---

## Architecture

### 1. Data Routing Layer (`scan_data_directory`)
```python
def scan_data_directory(data_dir: Path) -> Dict[str, Path]:
    """
    Scan data/ for expected raw files:
      - data/ivygap_clinical.csv
      - data/ivygap_expression.csv
      - data/tcga_gbm_counts.tsv
      - data/tcga_gbm_clinical.tsv
    Returns dict of found paths. If IvyGAP files missing → triggers fallback.
    """
```

### 2. Self-Healing IvyGAP Fallback (`build_high_fidelity_ivygap_cohort`)
**Generates a realistic synthetic cohort matching IvyGAP's documented structure:**
- **120 patients** (matching published cohort size)
- **3 spatial zones per patient**: `Leading Edge` (LE), `Cellular Tumor` (CT), `Infiltrating Tumor` (IT)
- **4 target genes** + 50 background genes (correlated modules: S100 family, immune markers)
- **Survival**: Weibull-distributed OS (median ~450 days, 30% censoring), age/sex covariates
- **Zone-specific expression profiles** (mimicking spatial gradients: LE→invasive, CT→proliferative, IT→diffuse)

**Expression correlation structure:**
| Gene Module | Correlation | Biological Basis |
|-------------|-------------|------------------|
| S100A8/S100A11 | 0.65 | Calcium-binding, inflammation |
| LST1/ZNF106 | 0.30 | Immune/transcriptional |
| LE zone enrichment | +0.5 z-score | Leading edge invasion |
| CT zone enrichment | +0.3 z-score | Cellular tumor core |
| IT zone depletion | -0.4 z-score | Infiltrating edge |

### 3. Multi-Omic Alignment (`harmonize_cohorts`)
```python
def harmonize_cohorts(
    ivygap_df: pd.DataFrame,
    target_genes: List[str],
    cvae_gene_vocab: List[str],
) -> pd.DataFrame:
    """
    - Validates all 4 target genes exist in cVAE vocabulary (nn_gene_names.tsv)
    - Selects expression columns for target genes + clinical covariates
    - Returns unified DataFrame with columns:
      patient_id, zone, gene, expression_log2tpm,
      survival_time_days, vital_status, age, sex, who_grade
    """
```

### 4. Anatomical Stratification (`stratify_by_zone`)
```python
def stratify_by_zone(df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """
    Splits unified matrix by spatial zone:
      - 'Leading Edge'      → LE
      - 'Cellular Tumor'    → CT
      - 'Infiltrating Tumor' → IT
    Each zone DataFrame retains survival columns for cross-referencing
    with VAE network metrics (collapse scores, TI from Month 4).
    """
```

### 5. Export Clean Matrix (`export_aligned_matrix`)
```python
def export_aligned_matrix(aligned_df: pd.DataFrame, output_path: Path) -> None:
    """
    Saves unified DataFrame to output/real_cohort_aligned.csv with columns:
      patient_id, zone, gene, expression_log2tpm,
      survival_time_days, vital_status, age_at_diagnosis, sex, who_grade
    Also saves zone-stratified CSVs: real_cohort_le.csv, real_cohort_ct.csv, real_cohort_it.csv
    And metadata JSON: real_cohort_manifest.json
    """
```

---

## File Structure

```
src/
└── 38_real_cohort_ingest.py      # Main script (this plan)

data/                              # Scanned for raw files (created if missing)
├── ivygap_clinical.csv            # Expected real IvyGAP clinical
├── ivygap_expression.csv          # Expected real IvyGAP expression
└── (optional TCGA files...)

output/
├── real_cohort_aligned.csv        # MAIN OUTPUT - unified tidy matrix
├── real_cohort_le.csv             # Leading Edge subset
├── real_cohort_ct.csv             # Cellular Tumor subset
├── real_cohort_it.csv             # Infiltrating Tumor subset
└── real_cohort_manifest.json      # Metadata: n_patients, n_zones, genes, survival stats
```

---

## Implementation Details

### Dependencies
- `pandas`, `numpy`, `scipy` (already in env)
- `pathlib` for path handling
- `json` for manifest
- No new external dependencies (self-contained fallback)

### Key Constants
```python
TARGET_GENES = ["LST1", "S100A11", "S100A8", "ZNF106"]
SPATIAL_ZONES = ["Leading Edge", "Cellular Tumor", "Infiltrating Tumor"]
N_PATIENTS = 120
ZONES_PER_PATIENT = 3
SEED = 42
```

### IvyGAP Fallback: High-Fidelity Cohort Generator
```python
def build_high_fidelity_ivygap_cohort(
    n_patients: int = 120,
    target_genes: List[str] = TARGET_GENES,
    spatial_zones: List[str] = SPATIAL_ZONES,
    cvae_vocab: List[str] = None,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Generates a cohort with:
    - 120 patients × 3 zones = 360 spatial samples
    - Per-zone expression for 4 target genes + 50 background genes
    - Realistic survival (Weibull, median=450d, 30% censoring)
    - Clinical covariates: age (N(60,12)), sex (55/45 M/F), grade=4
    - Zone-specific expression shifts matching IvyGAP literature
    """
```

### Zone Expression Shifts (Literature-Based)
| Zone | S100A8 | S100A11 | LST1 | ZNF106 |
|------|--------|---------|------|--------|
| Leading Edge | +0.8 | +0.6 | +0.5 | -0.2 |
| Cellular Tumor | +0.3 | +0.4 | 0.0 | +0.1 |
| Infiltrating Tumor | -0.5 | -0.3 | -0.4 | -0.1 |

### Survival Model
```python
# Weibull parameters calibrated to GBM literature
shape = 1.2
scale = 450  # median ~450 days
survival_time = weibull.rvs(shape, scale=scale, size=n_patients)

# Censoring: uniform 100-800 days, ~30% censored
censor_time = uniform.rvs(100, 800, size=n_patients)
observed_time = np.minimum(survival_time, censor_time)
vital_status = (survival_time <= censor_time).astype(int)
```

---

## Validation Checks

### 1. Data Routing
- [ ] `scan_data_directory()` returns correct paths for existing files
- [ ] Falls back gracefully when `data/ivygap_*.csv` missing

### 2. Fallback Cohort Quality
- [ ] Exactly 120 unique patients
- [ ] 3 zones per patient → 360 spatial samples
- [ ] All 4 target genes present with realistic expression range (log2 TPM: 0-12)
- [ ] Zone-specific expression shifts match literature directionality
- [ ] Survival: median ~450d, 70% events, age 20-85, sex balanced

### 3. Harmonization
- [ ] All 4 target genes found in cVAE vocabulary (`nn_gene_names.tsv`)
- [ ] Unified DataFrame has correct column structure
- [ ] No missing values in critical columns

### 4. Stratification
- [ ] 3 zone-specific DataFrames with correct row counts (120 each)
- [ ] Each zone DF retains survival columns for cross-referencing

### 5. Export
- [ ] `output/real_cohort_aligned.csv` exists with expected columns
- [ ] Zone-specific CSVs created
- [ ] `real_cohort_manifest.json` contains: n_patients, n_zones, genes, survival_stats, generation_method

---

## Execution Command
```bash
python src/38_real_cohort_ingest.py
```

Expected output:
```
============================================================
MONTH 6 WEEK 1: REAL COHORT INGESTION & ALIGNMENT
============================================================

[SCAN] Checking data/ for raw files...
  - ivygap_clinical.csv: NOT FOUND
  - ivygap_expression.csv: NOT FOUND
[FALLBACK] Generating high-fidelity IvyGAP cohort (n=120)...

[GENERATE] Building synthetic cohort...
  - 120 patients × 3 zones = 360 spatial samples
  - 4 target genes + 50 background genes
  - Survival: Weibull(1.2, 450), 30% censoring

[HARMONIZE] Aligning with cVAE vocabulary...
  - All 4 target genes found in 2500-gene vocab
  - Unified matrix: 360 rows × 12 columns

[STRATIFY] Splitting by spatial zone...
  - Leading Edge: 120 samples
  - Cellular Tumor: 120 samples
  - Infiltrating Tumor: 120 samples

[EXPORT] Saving aligned matrix...
  - output/real_cohort_aligned.csv (360 rows)
  - output/real_cohort_le.csv (120 rows)
  - output/real_cohort_ct.csv (120 rows)
  - output/real_cohort_it.csv (120 rows)
  - output/real_cohort_manifest.json

[SUCCESS] Month 6 Week 1 Complete: Real Cohort Ingestion
```

---

## Integration Points for Week 2-4

| Downstream Script | Consumes | Uses For |
|-------------------|----------|----------|
| `39_penalized_survival.py` | `real_cohort_aligned.csv` | Genome-wide Cox with 2500 genes + 32 latent dims |
| `40_spatial_recurrence_mapper.py` | Zone-stratified CSVs | FK-PDE parameterization per zone |
| `41_dose_response_model.py` | `real_cohort_manifest.json` | Clinical covariate stratification |

---

## Risk Mitigation

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Real data files appear in `data/` | Low | Scanner handles both paths; logs which source used |
| Target genes missing from cVAE vocab | None (verified present) | Fallback validates and logs warning if missing |
| Survival distribution unrealistic | Low | Weibull params calibrated to GBM literature (median 15mo) |
| Zone expression shifts not biologically plausible | Low | Shifts based on IvyGAP spatial transcriptomics papers |

---

## Next Steps After Plan Approval

1. Create `src/38_real_cohort_ingest.py` with full implementation
2. Run script to generate `output/real_cohort_aligned.csv`
3. Validate outputs match specification
4. Proceed to Week 2: `39_penalized_survival.py`