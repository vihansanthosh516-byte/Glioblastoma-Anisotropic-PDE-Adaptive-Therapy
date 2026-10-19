> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Month 4: In Silico Combinatorial Drug Gating Report

**Generated:** 2026-09-21 12:14:32
**Pipeline:** Multi-Scale Spatial Oncology Suite (MSOS) v1.0

---

## Executive Summary

This report presents the results of the **In Silico Combinatorial Drug Gating** pipeline 
applied to glioblastoma spatial transcriptomics (15,000 cells, 3 zones: Healthy, Periphery, Core).

**Key Findings:**
- **Top Single KO Target:** TMLHE (model C=0.0407984; `../../output/single_ko_results.json`)
- **Top Dual KO Synergy:** Highest listed Bliss: S100A11+S100A6, -0.00212064 (`../../output/dual_ko_results.json`); all Bliss values negative
- **Model-only interpretation:** No combinations achieved TI > 10 with tumor_collapse > 0.05
- **Separate rankings:** TMLHE highest collapse; MT-CO2 highest single model score 4.819325. No TMLHE TI=0.02 source exists.

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

therapeutic_index is a calibrated model score, not literal division by displayed healthy collapse; raw zeros and the 0.001 floor affect it. No pharmacologic therapeutic window measured.

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
1. **Therapeutic Index** (TI is a floor/calibration-dependent model score)
2. **Tumor Collapse** (C_tumor)
3. **Healthy Collapse** (C_healthy)  
4. **GRN Out-Degree** (master switch centrality)
5. **Max Bliss Synergy** (best combinatorial partner)
6. **Max Loewe Synergy** (additive expectation)

See `output/optimization_matrix.png` for heatmap visualization.

---

## 4. Inferred Network Centrality (not causal evidence)

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

## 5. Model-only interpretation Assessment

### Recommended Lead Combinations

| Rank | Combination | TI | Tumor C | Healthy C | Bliss | Priority |
|------|-------------|-----|---------|-----------|-------|----------|
| - | No combinations meet chosen model cutoffs (TI > 10, Tumor C > 0.05) | | | | | |


### Unvalidated model hypotheses
- **Pharmacodynamic:** Periphery transition score reduction
- **Patient Stratification:** High Periphery zone fraction (>20%)
- **Combo Biomarker:** Co-expression of target pair in Periphery zone

---

## 6. Historical proposed experiments (not clinical recommendations)

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
| Dual KO Results | `output/dual_ko_results.json` | 21 pairs × synergy metrics |
| Dual KO TI | `output/dual_ko_ti.json` | 21 pairs × TI metrics |
| Optimization Matrix | `output/optimization_matrix.npy` | Gene × 6 metrics |
| Master Switches | `output/master_switches.tsv` | 10 supplied centrality rows × centrality |
| Causal GRN | `output/causal_grn.graphml` | Cytoscape-compatible |
| Drug Gating Report | `output/drug_gating_report.md` | This document |

---

*Report generated by MSOS Pipeline v1.0 | 2026-09-21 12:14:32*

Audit correction 2026-10-19: Virtual perturbation scores are not drug efficacy. Root arrays: 200 single-KO rows, four single-score rows, 21 dual pairs (audit 15-pair count incorrect). Highest dual score: CCL3L1+S100A8 0.442764, Bliss -0.0185325, Loewe -0.00887558; distinct from combined-effect/Bliss ranking. Root optimization_matrix.npy and causal_grn.graphml exist, contrary to missing-file assertions. No causal validation, clinical threshold, patient treatment recommendation or journal suitability follows.
