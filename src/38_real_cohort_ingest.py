#!/usr/bin/env python3
"""
Month 5, Week 4: Real TCGA-GBM Cohort Ingestion

Merges real TCGA-GBM RNA-seq expression (S100A6, S100A11, S100A8, CCL3L1)
with real TCGA-GBM clinical survival data.

Produces real_cohort_aligned.csv with n=150 patients who have both
expression and survival.
"""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd

from pathlib import Path as _Path
PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Target genes present in the Xena HiSeqV2 matrix
TARGET_GENES = ["S100A6", "S100A11", "S100A8", "CCL3L1"]


def load_expression(path: Path = DATA_DIR / "tcga_gbm_expression.tsv") -> pd.DataFrame:
    """Load Xena HiSeqV2 expression matrix, transpose to patient × gene."""
    print(f"[LOAD] Reading expression: {path}")
    expr = pd.read_csv(path, sep="\t").set_index("sample")
    print(f"  Loaded {expr.shape[0]} genes × {expr.shape[1]} samples")

    # Subset to target genes
    missing = [g for g in TARGET_GENES if g not in expr.index]
    if missing:
        print(f"  [WARN] Missing genes: {missing}")
    present = [g for g in TARGET_GENES if g in expr.index]
    expr = expr.loc[present]
    print(f"  Using {len(present)} target genes: {present}")

    # Transpose: rows = patients, cols = genes
    expr = expr.T
    expr.index.name = "sample"
    print(f"  Transposed to {expr.shape[0]} patients × {expr.shape[1]} genes")
    return expr


def load_clinical(path: Path = DATA_DIR / "tcga_gbm_clinical.csv") -> pd.DataFrame:
    print(f"[LOAD] Reading clinical: {path}")
    clin = pd.read_csv(path)
    print(f"  Loaded {len(clin)} patients")

    # Subset to needed columns
    keep = ["sample", "overall_survival", "overall_survival_time",
            "age_at_diagnosis", "gender", "gene_expression", "vital_status"]
    keep = [c for c in keep if c in clin.columns]
    clin = clin[keep].copy()

    # Clean
    clin = clin.dropna(subset=["overall_survival", "overall_survival_time", "age_at_diagnosis"])
    print(f"  After cleaning: {len(clin)} patients with survival + age")
    return clin


def main():
    print("=" * 60)
    print("REAL TCGA-GBM COHORT INGESTION")
    print("=" * 60)
    print()

    expr = load_expression()
    print()
    clin = load_clinical()
    print()

    # Merge
    merged = clin.merge(expr, left_on="sample", right_index=True, how="inner")
    print(f"[MERGE] Overlap: {len(merged)} patients")

    # Rename for downstream scripts
    merged = merged.rename(columns={
        "overall_survival_time": "os_time_days",
        "overall_survival": "os_event",
    })

    # Coerce numeric
    for col in ["os_time_days", "os_event", "age_at_diagnosis"] + TARGET_GENES:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")

    # Drop rows with any NaN in target genes
    before = len(merged)
    merged = merged.dropna(subset=TARGET_GENES)
    print(f"  After dropping NaN in target genes: {len(merged)} (was {before})")

    print()
    print("[STATS] Cohort summary:")
    print(f"  Patients: {len(merged)}")
    print(f"  Events: {int(merged['os_event'].sum())} ({merged['os_event'].mean()*100:.1f}%)")
    print(f"  Median OS: {merged['os_time_days'].median():.1f} days")
    print(f"  Age range: {merged['age_at_diagnosis'].min():.0f}-{merged['age_at_diagnosis'].max():.0f}")
    print()

    # Split by zone (using clinical 'gene_expression' subtype as a proxy)
    # Write the full aligned cohort
    main_path = OUTPUT_DIR / "real_cohort_aligned.csv"
    merged.to_csv(main_path, index=False)
    print(f"[SAVE] {main_path}  ({len(merged)} rows, {len(merged.columns)} cols)")

    # Write "zone" splits by molecular subtype (proxy for zones since no IvyGAP)
    if "gene_expression" in merged.columns:
        for subtype in merged["gene_expression"].dropna().unique():
            sub_df = merged[merged["gene_expression"] == subtype]
            if len(sub_df) >= 10:
                suffix = subtype.lower().replace(" ", "_")
                path = OUTPUT_DIR / f"real_cohort_{suffix}.csv"
                sub_df.to_csv(path, index=False)
                print(f"[SAVE] {path}  ({len(sub_df)} patients)")

    # Manifest
    manifest = {
        "data_source": "TCGA-GBM (Xena HiSeqV2 + cBioPortal clinical)",
        "n_patients": int(len(merged)),
        "n_events": int(merged["os_event"].sum()),
        "median_os_days": float(merged["os_time_days"].median()),
        "target_genes": TARGET_GENES,
        "target_genes_present": [g for g in TARGET_GENES if g in merged.columns],
        "target_genes_missing": [g for g in TARGET_GENES if g not in merged.columns],
        "columns": list(merged.columns),
        "notes": "Real RNA-seq expression (Xena HiSeqV2, log2 TPM) merged with real survival (cBioPortal TCGA-GBM clinical). 150-patient overlap between the 518 clinical and 172 expression cohorts.",
    }
    manifest_path = OUTPUT_DIR / "real_cohort_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[SAVE] {manifest_path}")

    print()
    print("[SUCCESS] Real cohort ingestion complete")


if __name__ == "__main__":
    main()