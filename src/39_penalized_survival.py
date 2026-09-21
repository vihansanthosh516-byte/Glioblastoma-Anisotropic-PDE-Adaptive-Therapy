#!/usr/bin/env python3
"""
Month 6, Week 2: Penalized Survival Modeling

Elastic Net Cox regression on real TCGA-GBM expression + clinical covariates.
Uses lifelines (industry-standard survival analysis library).

Reads:  output/real_cohort_aligned.csv  (150 patients, real RNA-seq + survival)
Writes: output/penalized_survival_metrics.json
        output/penalized_coefficients.png
        output/penalized_regularization_paths.png
        output/penalized_survival_curves.png
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

from pathlib import Path as _Path
PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


EXPRESSION_GENES = ["S100A6", "S100A11", "S100A8", "CCL3L1"]
PENALIZER = 0.1
L1_RATIO = 0.5


def load_cohort(path: Path = OUTPUT_DIR / "real_cohort_aligned.csv") -> pd.DataFrame:
    print(f"[LOAD] Reading {path}")
    df = pd.read_csv(path)
    print(f"  Loaded {len(df)} patients, {len(df.columns)} columns")
    return df


def build_feature_frame(df: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
    """
    Build a pandas DataFrame with survival columns + features.
    CoxPHFitter wants everything in a single frame with named columns.
    """
    df = df.dropna(subset=["os_time_days", "os_event", "age_at_diagnosis"]).copy()

    features = df[EXPRESSION_GENES + ["age_at_diagnosis"]].copy()

    if "gender" in df.columns:
        features["gender_male"] = (df["gender"].astype(str).str.upper() == "MALE").astype(float)

    if "gene_expression" in df.columns:
        subtype_dummies = pd.get_dummies(df["gene_expression"], prefix="subtype")
        if subtype_dummies.shape[1] > 1:
            subtype_dummies = subtype_dummies.iloc[:, 1:]
        features = pd.concat([features, subtype_dummies], axis=1)

    # Survival columns
    features["T"] = df["os_time_days"].values.astype(float)
    features["E"] = df["os_event"].values.astype(int)

    feature_names = [c for c in features.columns if c not in ("T", "E")]

    print(f"[DESIGN] Frame: {len(features)} patients × {len(feature_names)} features")
    print(f"  Features: {feature_names}")
    print(f"  Events: {int(features['E'].sum())}/{len(features)}")
    return features, feature_names


def main():
    print("=" * 60)
    print("PENALIZED SURVIVAL MODELING — REAL TCGA-GBM COHORT")
    print("=" * 60)
    print()

    df = load_cohort()
    features, feature_names = build_feature_frame(df)

    print()
    print(f"[FIT] Elastic Net Cox (lifelines, penalizer={PENALIZER}, l1_ratio={L1_RATIO})")

    cph = CoxPHFitter(penalizer=PENALIZER, l1_ratio=L1_RATIO)
    cph.fit(features, duration_col="T", event_col="E")

    cph.print_summary()

    # Extract coefficients — note: .values is a property, not a method
    coefs = cph.params_
    coef_values = coefs.values  # numpy array
    coef_names = list(coefs.index)
    c_index = cph.concordance_index_

    print()
    print(f"[BEST] C-index: {c_index:.3f}")
    print(f"[BEST] Nonzero coefficients (|coef| > 1e-4):")
    for name, coef in zip(coef_names, coef_values):
        if abs(coef) > 1e-4:
            direction = "risk" if coef > 0 else "protective"
            print(f"  {name}: {coef:+.4f}  ({direction})")

    # Plot 1: Coefficient bar chart
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ['red' if c > 0 else 'blue' for c in coef_values]
    ax.barh(coef_names, coef_values, color=colors)
    ax.axvline(0, color='black', linewidth=0.5)
    ax.set_xlabel('Coefficient (log hazard ratio)')
    ax.set_title(f'Elastic Net Cox Coefficients (C-index={c_index:.3f})')
    plt.tight_layout()
    coef_plot = OUTPUT_DIR / "penalized_coefficients.png"
    plt.savefig(coef_plot, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[SAVE] {coef_plot}")

    # Plot 2: Regularization path
    print()
    print("[FIT] Regularization path (multiple penalizer values)")
    penalizers = np.logspace(-2, 1, 15)
    coef_path = []
    for p in penalizers:
        try:
            cph_p = CoxPHFitter(penalizer=p, l1_ratio=L1_RATIO)
            cph_p.fit(features, duration_col="T", event_col="E")
            coef_path.append(cph_p.params_.values)
        except Exception as e:
            print(f"  penalizer={p:.4f} failed: {e}")
            coef_path.append(np.zeros(len(feature_names)))
    coef_path = np.array(coef_path)

    fig, ax = plt.subplots(figsize=(10, 6))
    for i, name in enumerate(feature_names):
        ax.plot(penalizers, coef_path[:, i], label=name, marker='o', markersize=3)
    ax.axvline(PENALIZER, color='red', linestyle='--', label=f'Chosen penalizer = {PENALIZER}')
    ax.set_xscale('log')
    ax.set_xlabel('Penalizer (log scale)')
    ax.set_ylabel('Coefficient')
    ax.set_title('Elastic Net Cox Regularization Path')
    ax.legend(loc='best', fontsize=8)
    ax.axhline(0, color='black', linewidth=0.5)
    plt.tight_layout()
    reg_plot = OUTPUT_DIR / "penalized_regularization_paths.png"
    plt.savefig(reg_plot, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[SAVE] {reg_plot}")

    # Plot 3: Kaplan-Meier by risk quartile
    risk = cph.predict_partial_hazard(features)
    quartiles = np.percentile(risk, [25, 50, 75])
    group = np.digitize(risk, quartiles)

    fig, ax = plt.subplots(figsize=(8, 6))
    for g, color, label in zip([0, 1, 2, 3],
                                ['green', 'blue', 'orange', 'red'],
                                ['Q1 (lowest risk)', 'Q2', 'Q3', 'Q4 (highest risk)']):
        mask = group == g
        if mask.sum() == 0:
            continue
        t_g = features["T"].values[mask]
        e_g = features["E"].values[mask]
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
        ax.step(times, survs, where='post', label=f'{label} (n={mask.sum()})', color=color)

    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Survival probability')
    ax.set_title(f'Kaplan-Meier by Risk Quartile (C-index={c_index:.3f})')
    ax.legend(loc='best')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    km_plot = OUTPUT_DIR / "penalized_survival_curves.png"
    plt.savefig(km_plot, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[SAVE] {km_plot}")

    # Save metrics
    metrics = {
        "data_source": "Real TCGA-GBM (Xena expression + cBioPortal clinical)",
        "n_patients": int(len(features)),
        "n_events": int(features["E"].sum()),
        "n_features": len(feature_names),
        "feature_names": feature_names,
        "penalizer": float(PENALIZER),
        "l1_ratio": float(L1_RATIO),
        "c_index": float(c_index),
        "coefficients": {name: float(c) for name, c in zip(coef_names, coef_values)},
        "nonzero_features": [
            {"name": name, "coefficient": float(c),
             "direction": "risk" if c > 0 else "protective"}
            for name, c in zip(coef_names, coef_values) if abs(c) > 1e-4
        ],
        "concordance_index_lifelines": float(c_index),
    }
    metrics_path = OUTPUT_DIR / "penalized_survival_metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2, default=str)
    print(f"[SAVE] {metrics_path}")

    print()
    print("[SUCCESS] Penalized survival modeling complete")
    print(f"  C-index: {c_index:.3f}")
    print(f"  Nonzero features: {sum(1 for c in coef_values if abs(c) > 1e-4)}")


if __name__ == "__main__":
    main()