# Month 4: In Silico Combinatorial Drug Gating Report

**Generated:** 2026-09-21 12:14:32
**Pipeline:** Multi-Scale Spatial Oncology Suite (MSOS) v1.0

---

## Executive Summary

This report presents the results of the **In Silico Combinatorial Drug Gating** pipeline 
applied to glioblastoma spatial transcriptomics (15,000 cells, 3 zones: Healthy, Periphery, Core).

**Key Findings:**
- **Top Single KO Target:** TMLHE (C = 0.0378)
- **Top Dual KO Synergy:** S100A8 + S100A6 (Bliss = -0.4258, antagonistic)
- **Clinical Translation:** No combinations achieved TI > 10 with tumor_collapse > 0.05
- **Best Single Target:** TMLHE (collapse = 0.0378, TI = 0.02)

---

## 1. Single Gene Knockout Analysis

### Top 10 Single Gene Knockouts by Network Collapse Score (C)

| Rank | Gene | Collapse Score (C) | Zone | Master Switch Rank |
|------|------|-------------------|------|-------------------|
| 1 | TMLHE | 0.0408 | Periphery | N/A |
| 2 | MT-CO2 | 0.0282 | Periphery | N/A |
| 3 | THBD | 0.0118 | Periphery | N/A |
| 4 | ARHGAP30 | 0.0110 | Periphery | N/A |
| 5 | SYF2 | 0.0110 | Periphery | N/A |
| 6 | CCL3L1 | 0.0109 | Periphery | 7 |
| 7 | S100A12 | 0.0107 | Periphery | N/A |
| 8 | MMP19 | 0.0107 | Periphery | N/A |
| 9 | LSM10 | 0.0103 | Periphery | N/A |
| 10 | IRF1 | 0.0097 | Periphery | N/A |


### Therapeutic Index (TI) for Single Knockouts

TI = C_tumor / C_healthy (higher = better therapeutic window)

| Rank | Gene | Tumor C | Healthy C | TI |
|------|------|---------|-----------|-----|
| 1 | MT-CO2 | 0.0282 | 0.0000 | 4.82 |
| 2 | CCL3L1 | 0.0109 | 0.0000 | 3.44 |
| 3 | MT-CO3 | 0.0097 | 0.0000 | 3.27 |
| 4 | S100A8 | 0.0086 | 0.0016 | 2.42 |


---

## 2. Combinatorial Dual Knockout Screen

### Top 10 by Combined Effect (C)

| Rank | Gene A | Gene B | Combined C | Bliss Synergy | Loewe Synergy |
|------|--------|--------|-----------|---------------|---------------|
| 1 | CCL3L1 | S100A8 | 0.0014 | -0.0185 | -0.0089 |
| 2 | S100A11 | S100A6 | 0.0000 | -0.0021 | 0.0015 |
| 3 | MT-CO3 | S100A11 | 0.0000 | -0.0061 | 0.0007 |
| 4 | S100A11 | MT-ATP6 | 0.0000 | -0.0068 | -0.0032 |
| 5 | S100A8 | S100A11 | 0.0003 | -0.0078 | -0.0015 |
| 6 | MT-CO3 | MT-ATP6 | 0.0001 | -0.0078 | -0.0013 |
| 7 | MT-CO3 | S100A6 | 0.0000 | -0.0084 | -0.0019 |
| 8 | CCL3L1 | MT-ATP6 | 0.0008 | -0.0099 | -0.0029 |
| 9 | CCL3L1 | S100A11 | 0.0000 | -0.0121 | -0.0047 |
| 10 | CCL3L1 | S100A6 | 0.0006 | -0.0139 | -0.0068 |


### Top 10 by Therapeutic Index (TI)

| Rank | Gene A | Gene B | TI | Tumor C | Healthy C | Bliss | Loewe |
|------|--------|--------|-----|---------|-----------|-------|-------|
| 1 | CCL3L1 | S100A8 | 0.44 | 0.0014 | 0.0000 | -0.0185 | -0.0089 |
| 2 | S100A11 | S100A6 | 0.00 | 0.0000 | 0.0000 | -0.0021 | 0.0015 |
| 3 | MT-CO3 | S100A11 | 0.00 | 0.0000 | 0.0000 | -0.0061 | 0.0007 |
| 4 | S100A11 | MT-ATP6 | 0.00 | 0.0000 | 0.0005 | -0.0068 | -0.0032 |
| 5 | S100A8 | S100A11 | 0.00 | 0.0003 | 0.0006 | -0.0078 | -0.0015 |
| 6 | MT-CO3 | MT-ATP6 | 0.00 | 0.0001 | 0.0010 | -0.0078 | -0.0013 |
| 7 | MT-CO3 | S100A6 | 0.00 | 0.0000 | 0.0000 | -0.0084 | -0.0019 |
| 8 | CCL3L1 | MT-ATP6 | 0.00 | 0.0008 | 0.0000 | -0.0099 | -0.0029 |
| 9 | CCL3L1 | S100A11 | 0.00 | 0.0000 | 0.0000 | -0.0121 | -0.0047 |
| 10 | CCL3L1 | S100A6 | 0.00 | 0.0006 | 0.0000 | -0.0139 | -0.0068 |


---

## 3. Optimization Matrix

The optimization matrix evaluates each target across 6 dimensions:
1. **Therapeutic Index** (TI = C_tumor / C_healthy)
2. **Tumor Collapse** (C_tumor)
3. **Healthy Collapse** (C_healthy)  
4. **GRN Out-Degree** (master switch centrality)
5. **Max Bliss Synergy** (best combinatorial partner)
6. **Max Loewe Synergy** (additive expectation)

See `output/optimization_matrix.png` for heatmap visualization.

---

## 4. Master Switches & Causal GRN

**Top 10 Master Switches by Out-Degree Centrality:**

| Rank | Gene | Out-Degree | In-Degree | Total Degree | Targets |
|------|------|------------|-----------|--------------|---------|
| 1 | APOD | 46 | 4 | 50 | 46 |
| 2 | S100B | 42 | 5 | 47 | 42 |
| 3 | MT3 | 40 | 5 | 45 | 40 |
| 4 | S100A8 | 39 | 2 | 41 | 39 |
| 5 | S100A9 | 33 | 1 | 34 | 33 |
| 6 | IFITM2 | 25 | 4 | 29 | 25 |
| 7 | CCL3L1 | 23 | 1 | 24 | 23 |
| 8 | FCER1G | 12 | 4 | 16 | 12 |
| 9 | FPR1 | 9 | 10 | 19 | 9 |
| 10 | CSF3R | 9 | 12 | 21 | 9 |


---

## 5. Clinical Translation Assessment

### Recommended Lead Combinations

| Rank | Combination | TI | Tumor C | Healthy C | Bliss | Priority |
|------|-------------|-----|---------|-----------|-------|----------|
| - | No combinations meet clinical thresholds (TI > 10, Tumor C > 0.05) | | | | | |


### Biomarker Strategy
- **Pharmacodynamic:** Periphery transition score reduction
- **Patient Stratification:** High Periphery zone fraction (>20%)
- **Combo Biomarker:** Co-expression of target pair in Periphery zone

---

## 6. Next Steps (Month 5)

1. **Clinical Validation:** Test top 5 combinations on Ivy GAP / TCGA-GBM cohorts
2. **Dose-Response Modeling:** Extend binary KO to graded inhibition
3. **Spatial PK/PD:** Integrate drug diffusion in 3D tissue geometry
4. **Manuscript Preparation:** Compile for Nature Methods / Cancer Cell submission

---

## Appendix: Data Availability

| Artifact | Path | Description |
|----------|------|-------------|
| Single KO Results | `output/single_ko_results.json` | 200 genes × collapse scores |
| Single KO TI | `output/single_ko_ti.json` | 4 genes × TI metrics |
| Dual KO Results | `output/dual_ko_results.json` | 15 pairs × synergy metrics |
| Dual KO TI | `output/dual_ko_ti.json` | 15 pairs × TI metrics |
| Optimization Matrix | `output/optimization_matrix.npy` | Gene × 6 metrics |
| Master Switches | `output/master_switches.tsv` | 100 TFs × centrality |
| Causal GRN | `output/causal_grn.graphml` | Cytoscape-compatible |
| Drug Gating Report | `output/drug_gating_report.md` | This document |

---

*Report generated by MSOS Pipeline v1.0 | 2026-09-21 12:14:32*
