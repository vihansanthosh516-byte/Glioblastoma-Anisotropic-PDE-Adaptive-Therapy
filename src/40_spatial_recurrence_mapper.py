#!/usr/bin/env python3
"""
Month 6, Week 3: Recurrence Risk Mapper (adapted for TCGA-GBM molecular subtypes)

Original design: IvyGAP spatial zones (Leading Edge / Cellular Tumor / Infiltrating Tumor).
Adapted design: molecular subtypes (Proneural / Classical / Mesenchymal / Neural) as
pseudo-zones, since TCGA-GBM is bulk tumor (one sample per patient, no spatial dimensions).

Reads:  output/real_cohort_proneural.csv, _classical.csv, _mesenchymal.csv, _neural.csv
        output/penalized_survival_metrics.json (from script 39, for risk weights)
Writes: output/subtype_recurrence_summary.json
        output/subtype_recurrence_risk.png
        output/subtype_risk_profile.png
        output/subtype_invasion_scores.png
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path as _Path
PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Molecular subtypes act as pseudo-zones
SUBTYPES = ["Proneural", "Classical", "Mesenchymal", "Neural"]
SUBTYPE_FILES = {
    "Proneural": "real_cohort_proneural.csv",
    "Classical": "real_cohort_classical.csv",
    "Mesenchymal": "real_cohort_mesenchymal.csv",
    "Neural": "real_cohort_neural.csv",
}

TARGET_GENES = ["S100A6", "S100A11", "S100A8", "CCL3L1"]


def load_subtype_data(subtype: str) -> pd.DataFrame:
    """Load subtype-stratified cohort."""
    path = OUTPUT_DIR / SUBTYPE_FILES[subtype]
    if not path.exists():
        print(f"  [WARN] Missing {path}")
        return pd.DataFrame()
    df = pd.read_csv(path)
    print(f"  Loaded {subtype}: {len(df)} patients")
    return df


def compute_invasion_score(df: pd.DataFrame, weights: Dict[str, float]) -> float:
    """
    Compute invasion/risk score from expression + penalized Cox weights.

    Uses the coefficients from script 39's Elastic Net model.
    """
    score = 0.0
    for gene in TARGET_GENES:
        if gene in df.columns:
            mean_expr = df[gene].mean()
            weight = weights.get(gene, 0.0)
            score += weight * mean_expr
    return score


def load_penalized_weights(path: Path = OUTPUT_DIR / "penalized_survival_metrics.json") -> Dict[str, float]:
    """Load penalized Cox coefficients from script 39."""
    if not path.exists():
        print(f"  [WARN] {path} missing — using default weights")
        return {g: 0.0 for g in TARGET_GENES}
    with open(path) as f:
        metrics = json.load(f)
    return metrics.get("coefficients", {g: 0.0 for g in TARGET_GENES})


def main():
    print("=" * 60)
    print("RECURRENCE RISK MAPPER — TCGA-GBM MOLECULAR SUBTYPES")
    print("=" * 60)
    print()

    # Load penalized Cox weights from script 39
    print("[LOAD] Reading penalized survival weights...")
    weights = load_penalized_weights()
    print(f"  Weights: {weights}")

    print()
    print("[LOAD] Reading subtype-stratified cohorts...")
    subtype_dfs = {}
    for subtype in SUBTYPES:
        subtype_dfs[subtype] = load_subtype_data(subtype)

    # Compute per-subtype invasion score
    print()
    print("[COMPUTE] Invasion scores per subtype")
    subtype_scores = {}
    for subtype in SUBTYPES:
        df = subtype_dfs[subtype]
        if df.empty:
            subtype_scores[subtype] = 0.0
            continue
        score = compute_invasion_score(df, weights)
        subtype_scores[subtype] = score
        n_events = int(df["os_event"].sum()) if "os_event" in df.columns else 0
        median_os = df["os_time_days"].median() if "os_time_days" in df.columns else float("nan")
        print(f"  {subtype}: n={len(df)}, events={n_events}, "
              f"median_OS={median_os:.0f}d, score={score:.4f}")

    # Plot 1: Invasion score by subtype
    fig, ax = plt.subplots(figsize=(8, 6))
    subtypes_list = list(subtype_scores.keys())
    scores_list = list(subtype_scores.values())
    colors = ['green', 'blue', 'red', 'orange']
    ax.bar(subtypes_list, scores_list, color=colors)
    ax.set_ylabel('Invasion Score (weighted expression)')
    ax.set_title('Invasion Score by Molecular Subtype (TCGA-GBM, real data)')
    plt.xticks(rotation=15)
    plt.tight_layout()
    score_plot = OUTPUT_DIR / "subtype_invasion_scores.png"
    plt.savefig(score_plot, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[SAVE] {score_plot}")

    # Plot 2: Risk profile — expression heatmap by subtype
    fig, axes = plt.subplots(1, 4, figsize=(16, 5))
    for i, subtype in enumerate(SUBTYPES):
        df = subtype_dfs[subtype]
        if df.empty:
            continue
        means = [df[g].mean() for g in TARGET_GENES if g in df.columns]
        stds = [df[g].std() for g in TARGET_GENES if g in df.columns]
        present_genes = [g for g in TARGET_GENES if g in df.columns]
        axes[i].bar(present_genes, means, yerr=stds, color=colors[i], capsize=5)
        axes[i].set_title(f'{subtype} (n={len(df)})')
        axes[i].set_ylabel('Expression (log2 TPM)')
        axes[i].tick_params(axis='x', rotation=30)
    plt.suptitle('Target Gene Expression by Molecular Subtype (Real TCGA-GBM)', fontsize=14)
    plt.tight_layout()
    profile_plot = OUTPUT_DIR / "subtype_risk_profile.png"
    plt.savefig(profile_plot, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[SAVE] {profile_plot}")

    # Plot 3: Survival curves by subtype
    fig, ax = plt.subplots(figsize=(8, 6))
    for i, subtype in enumerate(SUBTYPES):
        df = subtype_dfs[subtype]
        if df.empty or "os_time_days" not in df.columns:
            continue
        t_g = df["os_time_days"].values
        e_g = df["os_event"].values
        order = np.argsort(t_g)
        t_sorted = t_g[order]
        e_sorted = e_g[order]
        survival = 1.0
        times = [0]
        survs = [1.0]
        for t in np.unique(t_sorted):
            at_risk = (t_sorted >= t).sum()
            events = ((t_sorted == t) & (e_sorted == 1)).sum()
            if at_risk > 0 and events > 0:
                survival *= (1 - events / at_risk)
            times.append(t)
            survs.append(survival)
        ax.step(times, survs, where='post', label=f'{subtype} (n={len(df)})', color=colors[i])

    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Survival probability')
    ax.set_title('Kaplan-Meier by Molecular Subtype (TCGA-GBM, real data)')
    ax.legend(loc='best')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    risk_plot = OUTPUT_DIR / "subtype_recurrence_risk.png"
    plt.savefig(risk_plot, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[SAVE] {risk_plot}")

    # Save summary
    summary = {
        "data_source": "Real TCGA-GBM (Xena expression + cBioPortal clinical)",
        "n_subtypes": len(SUBTYPES),
        "subtypes": SUBTYPES,
        "target_genes": TARGET_GENES,
        "penalized_weights": weights,
        "subtype_stats": {
            subtype: {
                "n_patients": int(len(subtype_dfs[subtype])),
                "n_events": int(subtype_dfs[subtype]["os_event"].sum())
                if not subtype_dfs[subtype].empty else 0,
                "median_os_days": float(subtype_dfs[subtype]["os_time_days"].median())
                if not subtype_dfs[subtype].empty else None,
                "invasion_score": float(subtype_scores[subtype]),
            }
            for subtype in SUBTYPES
        },
        "notes": "Spatial recurrence mapping adapted to molecular subtype stratification, "
                 "since TCGA-GBM is bulk tumor (no IvyGAP spatial zones available).",
    }
    summary_path = OUTPUT_DIR / "subtype_recurrence_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[SAVE] {summary_path}")

    print()
    print("[STATS] Summary by subtype:")
    for subtype in SUBTYPES:
        s = summary["subtype_stats"][subtype]
        print(f"  {subtype}: n={s['n_patients']}, events={s['n_events']}, "
              f"median_OS={s['median_os_days']}, score={s['invasion_score']:.4f}")

    print()
    print("[SUCCESS] Subtype-stratified recurrence mapping complete")


if __name__ == "__main__":
    main()