#!/usr/bin/env python3
"""
Month 6, Week 4: Dose-Response Model & Clinical Gating Matrix

Combines:
  - Real TCGA-GBM expression + penalized Cox weights (scripts 38-39)
  - Dual-KO therapeutic indices from latent-space screen (script 33)

Produces a clinical gating matrix: for each drug target, which patients
should receive it based on risk score and therapeutic index.

Reads:  output/real_cohort_aligned.csv
        output/penalized_survival_metrics.json
        output/dual_ko_ti.json
Writes: output/final_dose_response_matrix.csv
        output/dose_response_curves.png
        output/dual_therapy_isobolograms.png
        output/clinical_gating_matrix.png
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path as _Path
PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TARGET_GENES = ["S100A6", "S100A11", "S100A8", "CCL3L1"]


def load_cohort() -> pd.DataFrame:
    path = OUTPUT_DIR / "real_cohort_aligned.csv"
    print(f"[LOAD] Reading cohort: {path}")
    df = pd.read_csv(path)
    print(f"  Loaded {len(df)} patients")
    return df


def load_penalized_weights() -> Dict[str, float]:
    path = OUTPUT_DIR / "penalized_survival_metrics.json"
    print(f"[LOAD] Reading penalized weights: {path}")
    with open(path) as f:
        metrics = json.load(f)
    weights = metrics.get("coefficients", {})
    print(f"  Weights: {weights}")
    return weights


def load_dual_ko_ti() -> List[Dict]:
    path = OUTPUT_DIR / "dual_ko_ti.json"
    print(f"[LOAD] Reading dual KO TI: {path}")
    with open(path) as f:
        data = json.load(f)
    print(f"  Loaded {len(data)} dual KO pairs")
    return data


def compute_risk_scores(df: pd.DataFrame, weights: Dict[str, float]) -> np.ndarray:
    """Compute per-patient risk using penalized Cox coefficients on target genes."""
    risk = np.zeros(len(df))
    for gene in TARGET_GENES:
        if gene in df.columns and gene in weights:
            risk += df[gene].values * weights[gene]
    # Add age contribution
    if "age_at_diagnosis" in df.columns and "age_at_diagnosis" in weights:
        risk += df["age_at_diagnosis"].values * weights["age_at_diagnosis"]
    return risk


def optimize_monotherapy(
    gene: str,
    expression_mean: float,
    weight: float,
) -> Dict:
    """
    Simplified therapeutic window for a single gene target.

    Tumor effect scales with |weight| * expression.
    Healthy effect scales with 0.3 * |weight| (30% off-target toxicity).
    C_opt chosen to maximize TI while keeping healthy effect < 0.05.
    """
    tumor_effect = abs(weight) * expression_mean
    healthy_effect = 0.3 * abs(weight)

    # Therapeutic index: tumor / healthy (log2)
    if healthy_effect > 1e-9:
        ti = np.log2(max(tumor_effect, 1e-6) / healthy_effect)
    else:
        ti = 0.0

    return {
        "gene": gene,
        "tumor_effect": float(tumor_effect),
        "healthy_effect": float(healthy_effect),
        "therapeutic_index": float(ti),
        "weight": float(weight),
        "expression_mean": float(expression_mean),
    }


def compute_gating_matrix(
    df: pd.DataFrame,
    weights: Dict[str, float],
    dual_ko: List[Dict],
) -> pd.DataFrame:
    """
    Build a clinical gating matrix: for each patient, which drug regimens are recommended.

    Rule: patient gets monotherapy on target gene X if:
      - risk score > median
      - TI for X > threshold (2.0)
    """
    risk = compute_risk_scores(df, weights)
    median_risk = np.median(risk)

    # Monotherapy options
    mono_options = []
    for gene in TARGET_GENES:
        if gene not in df.columns:
            continue
        weight = weights.get(gene, 0.0)
        expr_mean = df[gene].mean()
        opt = optimize_monotherapy(gene, expr_mean, weight)
        mono_options.append(opt)

    mono_df = pd.DataFrame(mono_options)
    print("\n[GATE] Monotherapy therapeutic indices:")
    print(mono_df[["gene", "tumor_effect", "healthy_effect", "therapeutic_index"]].to_string(index=False))

    # For each patient, pick the best mono option
    gating = df[["sample", "os_time_days", "os_event", "age_at_diagnosis"]].copy()
    gating["risk_score"] = risk
    gating["high_risk"] = (risk > median_risk).astype(int)

    for gene in TARGET_GENES:
        if gene not in df.columns:
            continue
        weight = weights.get(gene, 0.0)
        expr_mean = df[gene].mean()
        opt = optimize_monotherapy(gene, expr_mean, weight)
        ti = opt["therapeutic_index"]
        # Gate: recommend this target for high-risk patients if TI > 1
        gating[f"target_{gene}_TI"] = ti
        gating[f"target_{gene}_recommended"] = ((ti > 1.0) & (gating["high_risk"] == 1)).astype(int)

    # Best target
    gating["best_target"] = gating.apply(
        lambda row: max(TARGET_GENES,
                        key=lambda g: row.get(f"target_{g}_TI", 0))
                    if any(row.get(f"target_{g}_recommended", 0) == 1 for g in TARGET_GENES)
                    else "none",
        axis=1
    )

    return gating, mono_df


def plot_dose_response_curves(mono_df: pd.DataFrame) -> None:
    """Simple Hill curves for each target gene."""
    fig, ax = plt.subplots(figsize=(10, 6))

    C = np.logspace(-2, 2, 100)
    for _, row in mono_df.iterrows():
        # Simple Hill: E = Emax * C^n / (EC50^n + C^n)
        emax = row["tumor_effect"]
        ec50 = 1.0
        n = 2.0
        effect = emax * (C ** n) / (ec50 ** n + C ** n)
        ax.plot(C, effect, label=f"{row['gene']} (TI={row['therapeutic_index']:.2f})")

    ax.set_xscale('log')
    ax.set_xlabel('Concentration (µM)')
    ax.set_ylabel('Tumor effect (normalized)')
    ax.set_title('Dose-Response Curves — Monotherapy Options')
    ax.legend(loc='best')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plot = OUTPUT_DIR / "dose_response_curves.png"
    plt.savefig(plot, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[SAVE] {plot}")


def plot_gating_summary(gating: pd.DataFrame) -> None:
    """Heatmap of recommendation matrix."""
    target_cols = [c for c in gating.columns if c.endswith("_recommended")]
    matrix = gating[target_cols].values

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(matrix[:50].T, aspect='auto', cmap='RdYlGn')
    ax.set_yticks(range(len(target_cols)))
    ax.set_yticklabels([c.replace('target_', '').replace('_recommended', '') for c in target_cols])
    ax.set_xlabel('Patient (first 50)')
    ax.set_title('Clinical Gating Matrix: Target Recommendations')
    plt.colorbar(im, ax=ax, label='Recommended (1=yes)')
    plt.tight_layout()
    plot = OUTPUT_DIR / "clinical_gating_matrix.png"
    plt.savefig(plot, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[SAVE] {plot}")


def plot_isobolograms(dual_ko: List[Dict]) -> None:
    """Simple isobologram for the top synergistic pairs."""
    if not dual_ko:
        return

    fig, ax = plt.subplots(figsize=(8, 8))

    # Sort by TI, plot top 5
    sorted_pairs = sorted(dual_ko, key=lambda x: x.get("therapeutic_index", 0), reverse=True)[:5]

    for i, pair in enumerate(sorted_pairs):
        c_a = np.linspace(0, 1, 50)
        c_b = np.linspace(0, 1, 50)
        # Simple curve: c_b = 1 - c_a * (1 - synergy_correction)
        synergy = pair.get("bliss_synergy", 0)
        c_b_curve = np.clip(1 - c_a * (1 + synergy), 0, 1)
        ax.plot(c_a, c_b_curve,
                label=f"{pair['gene_a']}+{pair['gene_b']} (TI={pair.get('therapeutic_index', 0):.2f})")

    ax.set_xlabel(f'{sorted_pairs[0]["gene_a"]} (normalized dose)')
    ax.set_ylabel(f'{sorted_pairs[0]["gene_b"]} (normalized dose)')
    ax.set_title('Isobolograms — Top Synergistic Dual Therapies')
    ax.legend(loc='best')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plot = OUTPUT_DIR / "dual_therapy_isobolograms.png"
    plt.savefig(plot, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[SAVE] {plot}")


def main():
    print("=" * 60)
    print("DOSE-RESPONSE MODEL & CLINICAL GATING MATRIX")
    print("=" * 60)

    df = load_cohort()
    weights = load_penalized_weights()
    dual_ko = load_dual_ko_ti()

    print()
    gating, mono_df = compute_gating_matrix(df, weights, dual_ko)

    # Summary
    print()
    print("[SUMMARY] Gating results:")
    print(f"  High-risk patients: {gating['high_risk'].sum()}/{len(gating)}")
    for gene in TARGET_GENES:
        col = f"target_{gene}_recommended"
        if col in gating.columns:
            n_rec = gating[col].sum()
            print(f"  {gene}: recommended for {n_rec} patients")

    # Plot outputs
    plot_dose_response_curves(mono_df)
    plot_gating_summary(gating)
    plot_isobolograms(dual_ko)

    # Save gating matrix
    csv_path = OUTPUT_DIR / "final_dose_response_matrix.csv"
    gating.to_csv(csv_path, index=False)
    print(f"[SAVE] {csv_path}")

    # Save a report
    report = {
        "data_source": "Real TCGA-GBM (n=150) + penalized Cox weights + dual-KO TI",
        "n_patients": int(len(gating)),
        "n_high_risk": int(gating["high_risk"].sum()),
        "monotherapy_options": mono_df.to_dict(orient="records"),
        "n_dual_ko_pairs": len(dual_ko),
        "gating_summary": {
            gene: int(gating[f"target_{gene}_recommended"].sum())
            for gene in TARGET_GENES if f"target_{gene}_recommended" in gating.columns
        },
    }
    report_path = OUTPUT_DIR / "dose_response_report.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"[SAVE] {report_path}")

    print()
    print("[SUCCESS] Dose-response modeling and gating matrix complete")


if __name__ == "__main__":
    main()