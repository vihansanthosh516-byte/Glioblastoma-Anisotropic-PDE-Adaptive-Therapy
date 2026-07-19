#!/usr/bin/env python3
"""
Month 6, Week 2: Penalized Survival Modeling & Feature Selection

Loads the aligned cohort from Week 1, pivots to wide format per zone,
implements Elastic Net penalized Cox PH via coordinate descent,
runs cross-validation for alpha selection, and exports zone-stratified
hazard ratios for the 4 target genes.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from sklearn.model_selection import KFold

warnings.filterwarnings("ignore")


# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #
TARGET_GENES = ["LST1", "S100A11", "S100A8", "ZNF106"]
SPATIAL_ZONES = ["Leading Edge", "Cellular Tumor", "Infiltrating Tumor"]
ALPHAS = np.logspace(-4, 1, 20)  # 20 alphas from 1e-4 to 10
L1_RATIO = 0.5  # Elastic Net mixing (0.5 = equal L1/L2)
N_FOLDS = 5
SEED = 42


# --------------------------------------------------------------------------- #
# Data Loading & Pivoting
# --------------------------------------------------------------------------- #
def load_aligned_cohort(path: Path = Path("output/real_cohort_aligned.csv")) -> pd.DataFrame:
    """Load the unified long-format cohort from Week 1."""
    df = pd.read_csv(path)
    print(f"  Loaded {df.shape[0]} rows, {df.shape[1]} columns from {path}")
    print(f"  Patients: {df['patient_id'].nunique()}, Zones: {df['zone'].nunique()}, Genes: {df['gene'].nunique()}")
    return df


def pivot_to_wide(
    long_df: pd.DataFrame,
    zone: str,
    target_genes: List[str] = TARGET_GENES,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]:
    """
    Pivot long-format data to wide matrix for a specific zone.

    Returns:
        X: (n_patients, n_features) design matrix
        time: (n_patients,) survival time
        event: (n_patients,) event indicator (1=deceased, 0=censored)
        feature_names: list of column names
    """
    zone_df = long_df[long_df["zone"] == zone].copy()

    # Pivot expression: rows=patient, cols=gene
    expr_wide = zone_df.pivot_table(
        index="patient_id",
        columns="gene",
        values="expression_log2tpm",
        aggfunc="first",
    )

    # Ensure target genes are present
    for g in target_genes:
        if g not in expr_wide.columns:
            raise ValueError(f"Target gene {g} not found in zone {zone}")

    # Clinical covariates (take first row per patient)
    clinical = zone_df.drop_duplicates("patient_id").set_index("patient_id")[
        ["survival_time_days", "vital_status", "age_at_diagnosis", "sex"]
    ]

    # Align
    common_patients = expr_wide.index.intersection(clinical.index)
    expr_wide = expr_wide.loc[common_patients]
    clinical = clinical.loc[common_patients]

    # Build design matrix: target genes + age + sex (encoded)
    X_genes = expr_wide[target_genes].values.astype(float)
    age = clinical["age_at_diagnosis"].values.astype(float).reshape(-1, 1)
    sex = (clinical["sex"].values == "M").astype(float).reshape(-1, 1)

    X = np.hstack([X_genes, age, sex])
    feature_names = target_genes + ["age_at_diagnosis", "sex_M"]

    time = clinical["survival_time_days"].values.astype(float)
    event = clinical["vital_status"].values.astype(int)

    # Center/scale features for penalized regression
    X_mean = X.mean(axis=0)
    X_std = X.std(axis=0)
    X_std[X_std == 0] = 1.0
    X_scaled = (X - X_mean) / X_std

    print(f"  Zone '{zone}': {X.shape[0]} patients, {X.shape[1]} features, {event.sum()} events")
    return X_scaled, time, event, feature_names, X_mean, X_std


# --------------------------------------------------------------------------- #
# Penalized Cox PH (Elastic Net via Coordinate Descent)
# --------------------------------------------------------------------------- #
class PenalizedCox:
    """
    Elastic Net Cox Proportional Hazards model via coordinate descent.

    Loss: -partial_log_likelihood + alpha * [(1-l1_ratio)/2 * ||beta||_2^2 + l1_ratio * ||beta||_1]
    """

    def __init__(
        self,
        alpha: float = 1.0,
        l1_ratio: float = 0.5,
        max_iter: int = 100,
        tol: float = 1e-4,
    ):
        self.alpha = alpha
        self.l1_ratio = l1_ratio
        self.max_iter = max_iter
        self.tol = tol
        self.coef_ = None
        self.baseline_hazard_ = None

    def _partial_log_likelihood(self, beta: np.ndarray, X: np.ndarray, time: np.ndarray, event: np.ndarray) -> float:
        """Negative partial log-likelihood (Breslow)."""
        n = X.shape[0]
        if n == 0:
            return np.inf

        eta = X @ beta
        # Sort by time
        order = np.argsort(time)
        X_sorted = X[order]
        eta_sorted = eta[order]
        event_sorted = event[order]
        time_sorted = time[order]

        # Unique event times
        event_times = np.unique(time_sorted[event_sorted == 1])
        if len(event_times) == 0:
            return np.inf

        loglik = 0.0
        for t in event_times:
            # Risk set: all with time >= t
            risk_mask = time_sorted >= t
            if not risk_mask.any():
                continue

            X_risk = X_sorted[risk_mask]
            eta_risk = eta_sorted[risk_mask]

            # Events at time t
            events_at_t = (time_sorted == t) & (event_sorted == 1)
            X_events = X_sorted[events_at_t]
            n_events = events_at_t.sum()

            if n_events == 0:
                continue

            # Breslow approximation
            exp_eta_risk = np.exp(eta_risk - eta_risk.max())  # numerical stability
            s0 = exp_eta_risk.sum()
            s1 = (X_risk.T @ exp_eta_risk) / s0

            loglik += (X_events @ beta).sum() - n_events * (np.log(s0) + eta_risk.max())

        return -loglik / n  # average negative log-likelihood

    def _gradient(self, beta: np.ndarray, X: np.ndarray, time: np.ndarray, event: np.ndarray) -> np.ndarray:
        """Gradient of negative partial log-likelihood."""
        n, p = X.shape
        eta = X @ beta
        order = np.argsort(time)
        X_sorted = X[order]
        eta_sorted = eta[order]
        event_sorted = event[order]
        time_sorted = time[order]

        event_times = np.unique(time_sorted[event_sorted == 1])
        if len(event_times) == 0:
            return np.zeros(p)

        grad = np.zeros(p)
        for t in event_times:
            risk_mask = time_sorted >= t
            if not risk_mask.any():
                continue

            X_risk = X_sorted[risk_mask]
            eta_risk = eta_sorted[risk_mask]
            events_at_t = (time_sorted == t) & (event_sorted == 1)
            X_events = X_sorted[events_at_t]
            n_events = events_at_t.sum()

            if n_events == 0:
                continue

            exp_eta_risk = np.exp(eta_risk - eta_risk.max())
            s0 = exp_eta_risk.sum()
            s1 = (X_risk.T @ exp_eta_risk) / s0

            grad += X_events.sum(axis=0) - n_events * s1

        return -grad / n

    def _coordinate_descent_step(
        self,
        beta: np.ndarray,
        X: np.ndarray,
        time: np.ndarray,
        event: np.ndarray,
        j: int,
    ) -> float:
        """Single coordinate descent update for coordinate j."""
        n, p = X.shape

        # Compute gradient at current beta
        grad = self._gradient(beta, X, time, event)

        # Coordinate-wise update for Elastic Net
        # beta_j <- S(beta_j - grad_j / H_jj, alpha * l1_ratio) / (1 + alpha * (1 - l1_ratio))
        # where S is soft-thresholding, H_jj is diagonal of Hessian

        # Approximate diagonal Hessian (Fisher information)
        eta = X @ beta
        order = np.argsort(time)
        X_sorted = X[order]
        eta_sorted = eta[order]
        event_sorted = event[order]
        time_sorted = time[order]

        event_times = np.unique(time_sorted[event_sorted == 1])
        h_jj = 0.0

        for t in event_times:
            risk_mask = time_sorted >= t
            if not risk_mask.any():
                continue

            X_risk = X_sorted[risk_mask]
            eta_risk = eta_sorted[risk_mask]
            events_at_t = (time_sorted == t) & (event_sorted == 1)
            n_events = events_at_t.sum()

            if n_events == 0:
                continue

            exp_eta_risk = np.exp(eta_risk - eta_risk.max())
            s0 = exp_eta_risk.sum()
            s1_j = (X_risk[:, j] * exp_eta_risk).sum() / s0
            s2_jj = (X_risk[:, j] ** 2 * exp_eta_risk).sum() / s0

            var_j = s2_jj - s1_j ** 2
            h_jj += n_events * var_j

        h_jj = max(h_jj / n, 1e-8)

        # Current coefficient
        beta_j = beta[j]
        grad_j = grad[j]

        # Soft-thresholding update
        z = beta_j - grad_j / h_jj
        lambda_l1 = self.alpha * self.l1_ratio
        lambda_l2 = self.alpha * (1 - self.l1_ratio)

        if z > lambda_l1 / h_jj:
            beta_new = (z - lambda_l1 / h_jj) / (1 + lambda_l2 / h_jj)
        elif z < -lambda_l1 / h_jj:
            beta_new = (z + lambda_l1 / h_jj) / (1 + lambda_l2 / h_jj)
        else:
            beta_new = 0.0

        return beta_new

    def fit(
        self,
        X: np.ndarray,
        time: np.ndarray,
        event: np.ndarray,
        warm_start: Optional[np.ndarray] = None,
    ) -> "PenalizedCox":
        n, p = X.shape

        if warm_start is not None:
            beta = warm_start.copy()
        else:
            beta = np.zeros(p)

        for iteration in range(self.max_iter):
            beta_old = beta.copy()

            # Coordinate descent cycle
            for j in range(p):
                beta[j] = self._coordinate_descent_step(beta, X, time, event, j)

            # Check convergence
            if np.max(np.abs(beta - beta_old)) < self.tol:
                break

        self.coef_ = beta
        self._compute_baseline_hazard(X, time, event)
        return self

    def _compute_baseline_hazard(self, X: np.ndarray, time: np.ndarray, event: np.ndarray):
        """Compute Breslow baseline hazard."""
        if self.coef_ is None:
            return

        eta = X @ self.coef_
        order = np.argsort(time)
        time_sorted = time[order]
        event_sorted = event[order]
        eta_sorted = eta[order]

        event_times = np.unique(time_sorted[event_sorted == 1])
        baseline_hazard = []

        for t in event_times:
            risk_mask = time_sorted >= t
            if not risk_mask.any():
                continue

            exp_eta_risk = np.exp(eta_sorted[risk_mask] - eta_sorted[risk_mask].max())
            s0 = exp_eta_risk.sum()
            d_t = ((time_sorted == t) & (event_sorted == 1)).sum()

            if s0 > 0:
                baseline_hazard.append((t, d_t / s0))

        self.baseline_hazard_ = np.array(baseline_hazard) if baseline_hazard else np.array([])

    def predict_risk(self, X: np.ndarray) -> np.ndarray:
        """Predict risk scores (linear predictor)."""
        if self.coef_ is None:
            raise ValueError("Model not fitted")
        return X @ self.coef_

    def score(self, X: np.ndarray, time: np.ndarray, event: np.ndarray) -> float:
        """Concordance index (C-index)."""
        risk = self.predict_risk(X)
        return concordance_index(risk, time, event)


def concordance_index(risk: np.ndarray, time: np.ndarray, event: np.ndarray) -> float:
    """Harrell's C-index for survival data."""
    n = len(time)
    if n < 2:
        return 0.5

    concordant = 0
    total = 0

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            # Comparable pairs: i has event before j, or i has event and j is censored after i's time
            if event[i] == 1 and time[i] < time[j]:
                total += 1
                if risk[i] > risk[j]:
                    concordant += 1
                elif risk[i] == risk[j]:
                    concordant += 0.5
            elif event[i] == 1 and event[j] == 0 and time[i] <= time[j]:
                total += 1
                if risk[i] > risk[j]:
                    concordant += 1
                elif risk[i] == risk[j]:
                    concordant += 0.5

    return concordant / total if total > 0 else 0.5


# --------------------------------------------------------------------------- #
# Cross-Validation for Alpha Selection
# --------------------------------------------------------------------------- #
def cross_validate_alpha(
    X: np.ndarray,
    time: np.ndarray,
    event: np.ndarray,
    alphas: np.ndarray,
    l1_ratio: float = L1_RATIO,
    n_folds: int = N_FOLDS,
    seed: int = SEED,
) -> Dict:
    """K-fold CV to select optimal alpha."""
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    cv_results = {alpha: [] for alpha in alphas}

    for fold, (train_idx, val_idx) in enumerate(kf.split(X)):
        X_train, X_val = X[train_idx], X[val_idx]
        time_train, time_val = time[train_idx], time[val_idx]
        event_train, event_val = event[train_idx], event[val_idx]

        # Warm start: fit from largest to smallest alpha
        warm_start = None
        for alpha in sorted(alphas, reverse=True):
            model = PenalizedCox(alpha=alpha, l1_ratio=l1_ratio)
            model.fit(X_train, time_train, event_train, warm_start=warm_start)
            warm_start = model.coef_

            cidx = model.score(X_val, time_val, event_val)
            cv_results[alpha].append(cidx)

    # Average C-index per alpha
    mean_cidx = {a: np.mean(cv_results[a]) for a in alphas}
    std_cidx = {a: np.std(cv_results[a]) for a in alphas}

    # Best alpha: highest mean C-index
    best_alpha = max(alphas, key=lambda a: mean_cidx[a])

    return {
        "alphas": alphas.tolist(),
        "mean_cindex": [mean_cidx[a] for a in alphas],
        "std_cindex": [std_cidx[a] for a in alphas],
        "best_alpha": float(best_alpha),
        "best_cindex": float(mean_cidx[best_alpha]),
        "fold_results": cv_results,
    }


# --------------------------------------------------------------------------- #
# Main Pipeline
# --------------------------------------------------------------------------- #
def main():
    print("=" * 60)
    print("MONTH 6 WEEK 2: PENALIZED SURVIVAL MODELING")
    print("=" * 60)

    # 1. Load data
    print("\n[LOAD] Reading aligned cohort...")
    cohort_df = load_aligned_cohort()

    # 2. Zone-stratified modeling
    print("\n[MODEL] Fitting penalized Cox per zone...")
    zone_results = {}

    for zone in SPATIAL_ZONES:
        print(f"\n  --- Zone: {zone} ---")
        X, time, event, feature_names, X_mean, X_std = pivot_to_wide(cohort_df, zone)

        # Cross-validation for alpha
        print(f"    [CV] Searching alpha over {len(ALPHAS)} values...")
        cv_result = cross_validate_alpha(X, time, event, ALPHAS)
        print(f"    [CV] Best alpha = {cv_result['best_alpha']:.4f}, C-index = {cv_result['best_cindex']:.4f}")

        # Fit final model on full data with best alpha
        best_alpha = cv_result["best_alpha"]
        final_model = PenalizedCox(alpha=best_alpha, l1_ratio=L1_RATIO)
        final_model.fit(X, time, event)

        # Extract coefficients for target genes
        coef_dict = dict(zip(feature_names, final_model.coef_))
        target_coefs = {g: coef_dict.get(g, 0.0) for g in TARGET_GENES}

        # Hazard ratios
        target_hrs = {g: np.exp(coef_dict.get(g, 0.0)) for g in TARGET_GENES}

        # C-index on full data
        cindex_full = final_model.score(X, time, event)

        zone_results[zone] = {
            "best_alpha": float(best_alpha),
            "cindex_cv": float(cv_result["best_cindex"]),
            "cindex_full": float(cindex_full),
            "coefficients": coef_dict,
            "target_coefficients": target_coefs,
            "target_hazard_ratios": target_hrs,
            "feature_names": feature_names,
            "n_patients": int(len(time)),
            "n_events": int(event.sum()),
        }

        print(f"    [RESULT] C-index (full) = {cindex_full:.4f}")
        for g in TARGET_GENES:
            print(f"      {g}: beta={target_coefs[g]:.4f}, HR={target_hrs[g]:.4f}")

    # 3. Leading Edge vs Cellular Tumor comparison
    print("\n[COMPARE] Leading Edge vs Cellular Tumor:")
    le_coefs = zone_results["Leading Edge"]["target_coefficients"]
    ct_coefs = zone_results["Cellular Tumor"]["target_coefficients"]
    for g in TARGET_GENES:
        diff = le_coefs[g] - ct_coefs[g]
        print(f"  {g}: LE beta={le_coefs[g]:.4f}, CT beta={ct_coefs[g]:.4f}, diff={diff:.4f}")

    # 4. Export metrics JSON
    print("\n[EXPORT] Saving metrics...")
    output_dir = Path("output")
    output_dir.mkdir(exist_ok=True)

    metrics = {
        "target_genes": TARGET_GENES,
        "zones": SPATIAL_ZONES,
        "l1_ratio": L1_RATIO,
        "cv_folds": N_FOLDS,
        "alphas_tested": ALPHAS.tolist(),
        "zone_results": zone_results,
        "comparison": {
            "leading_edge_vs_cellular_tumor": {
                g: {
                    "le_beta": float(le_coefs[g]),
                    "ct_beta": float(ct_coefs[g]),
                    "difference": float(le_coefs[g] - ct_coefs[g]),
                    "le_hr": float(zone_results["Leading Edge"]["target_hazard_ratios"][g]),
                    "ct_hr": float(zone_results["Cellular Tumor"]["target_hazard_ratios"][g]),
                }
                for g in TARGET_GENES
            }
        },
    }

    metrics_path = output_dir / "penalized_survival_metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2, default=str)
    print(f"  Saved: {metrics_path}")

    # 5. Comparative coefficient footprint plot
    print("\n[PLOT] Generating coefficient footprint plot...")
    fig, axes = plt.subplots(1, len(SPATIAL_ZONES), figsize=(5 * len(SPATIAL_ZONES), 5), sharey=True)

    for idx, zone in enumerate(SPATIAL_ZONES):
        ax = axes[idx]
        coefs = zone_results[zone]["target_coefficients"]
        genes = list(coefs.keys())
        values = list(coefs.values())

        colors = ['#e74c3c' if v > 0 else '#3498db' for v in values]
        bars = ax.barh(genes, values, color=colors, alpha=0.7, edgecolor='k', linewidth=0.5)
        ax.axvline(x=0, color='k', linestyle='-', linewidth=0.5)
        ax.set_xlabel('Coefficient (beta)', fontsize=10)
        ax.set_title(f'{zone}\n(α={zone_results[zone]["best_alpha"]:.3f}, C={zone_results[zone]["cindex_full"]:.3f})', fontsize=11)
        ax.grid(True, axis='x', alpha=0.3)

        # Add value labels
        for bar, val in zip(bars, values):
            ax.text(val + (0.01 if val >= 0 else -0.01), bar.get_y() + bar.get_height()/2,
                    f'{val:.3f}', va='center', ha='left' if val >= 0 else 'right', fontsize=9)

    axes[0].set_ylabel('Target Gene', fontsize=11)
    plt.suptitle('Penalized Cox Coefficients: Target Genes by Spatial Zone\n(Elastic Net, α=CV-selected, L1_ratio=0.5)',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()

    plot_path = output_dir / "penalized_coefficients.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {plot_path}")

    # 6. Also save coefficient path plot (alpha vs coefficient)
    print("\n[PLOT] Generating regularization path...")
    fig, axes = plt.subplots(1, len(SPATIAL_ZONES), figsize=(5 * len(SPATIAL_ZONES), 5), sharey=True)

    for idx, zone in enumerate(SPATIAL_ZONES):
        ax = axes[idx]
        X, time, event, feature_names, _, _ = pivot_to_wide(cohort_df, zone)

        # Fit path
        coef_paths = {g: [] for g in TARGET_GENES}
        for alpha in sorted(ALPHAS, reverse=True):
            model = PenalizedCox(alpha=alpha, l1_ratio=L1_RATIO)
            model.fit(X, time, event)
            for g in TARGET_GENES:
                g_idx = feature_names.index(g) if g in feature_names else -1
                coef_paths[g].append(model.coef_[g_idx] if g_idx >= 0 else 0.0)

        alphas_sorted = sorted(ALPHAS, reverse=True)
        for g in TARGET_GENES:
            ax.plot(alphas_sorted, coef_paths[g], 'o-', label=g, linewidth=2, markersize=4)

        ax.set_xscale('log')
        ax.set_xlabel('Alpha (regularization strength)', fontsize=10)
        ax.set_title(f'{zone}', fontsize=11)
        ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    axes[0].set_ylabel('Coefficient (beta)', fontsize=11)
    plt.suptitle('Regularization Paths: Target Gene Coefficients vs Alpha\n(Elastic Net, L1_ratio=0.5)',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()

    path_plot_path = output_dir / "penalized_regularization_paths.png"
    plt.savefig(path_plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path_plot_path}")

    print("\n" + "=" * 60)
    print("[SUCCESS] Month 6 Week 2 Complete: Penalized Survival Modeling")
    print("=" * 60)
    print(f"  - {metrics_path}")
    print(f"  - {plot_path}")
    print(f"  - {path_plot_path}")


if __name__ == "__main__":
    main()