#!/usr/bin/env python3
"""
Month 5, Week 2: Survival Analysis Suite (Real TCGA-GBM Clinical Data)

Loads the clinical mapped cohort from script 35 (real TCGA-GBM data),
computes Kaplan-Meier curves, log-rank tests, Cox PH models on clinical
covariates (age, gender, molecular subtype), and generates publication-ready
multi-panel survival plots.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

from pathlib import Path as _Path
PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")


# --------------------------------------------------------------------------- #
# Survival Analysis Core (Manual Implementation)
# --------------------------------------------------------------------------- #
def logrank_test(
    time: np.ndarray,
    event: np.ndarray,
    group: np.ndarray,
) -> Tuple[float, float]:
    """Two-sample log-rank test (Mantel-Cox)."""
    unique_times = np.unique(time[event == 1])
    if len(unique_times) == 0:
        return 0.0, 1.0

    g1, g2 = np.unique(group)
    mask1 = group == g1
    mask2 = group == g2

    O1 = O2 = E1 = E2 = 0.0

    for t in unique_times:
        at_risk1 = np.sum((time >= t) & mask1)
        at_risk2 = np.sum((time >= t) & mask2)
        at_risk = at_risk1 + at_risk2

        events1 = np.sum((time == t) & (event == 1) & mask1)
        events2 = np.sum((time == t) & (event == 1) & mask2)
        events = events1 + events2

        if at_risk > 0:
            exp1 = events * at_risk1 / at_risk
            exp2 = events * at_risk2 / at_risk
            O1 += events1
            O2 += events2
            E1 += exp1
            E2 += exp2

    V = 0.0
    for t in unique_times:
        at_risk1 = np.sum((time >= t) & mask1)
        at_risk2 = np.sum((time >= t) & mask2)
        at_risk = at_risk1 + at_risk2
        events = np.sum((time == t) & (event == 1))
        if at_risk > 1 and events > 0:
            v_t = (at_risk1 * at_risk2 * events * (at_risk - events)) / \
                  (at_risk**2 * (at_risk - 1))
            V += v_t

    if V == 0:
        return 0.0, 1.0

    chi2 = (O1 - E1)**2 / V
    p = 1 - stats.chi2.cdf(chi2, df=1)
    return chi2, p


def kaplan_meier(
    time: np.ndarray,
    event: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Kaplan-Meier survival curve."""
    order = np.argsort(time)
    time_sorted = time[order]
    event_sorted = event[order]

    unique_times = np.unique(time_sorted[event_sorted == 1])
    if len(unique_times) == 0:
        return np.array([time_sorted.max()]), np.array([1.0])

    surv = 1.0
    surv_curve = [1.0]
    time_curve = [0.0]

    for t in unique_times:
        at_risk = np.sum(time_sorted >= t)
        events = np.sum((time_sorted == t) & (event_sorted == 1))
        if at_risk > 0:
            surv *= (1 - events / at_risk)
        surv_curve.append(surv)
        time_curve.append(t)

    return np.array(time_curve), np.array(surv_curve)


def cox_ph_fit(
    time: np.ndarray,
    event: np.ndarray,
    X: np.ndarray,
    max_iter: int = 100,
    tol: float = 1e-6,
) -> Dict:
    """Cox Proportional Hazards via Breslow partial likelihood."""
    n, p = X.shape
    if p == 0:
        return {"coefficients": np.array([]), "hr": np.array([]),
                "ci_lower": np.array([]), "ci_upper": np.array([]),
                "p_values": np.array([])}

    X_centered = X - X.mean(axis=0)
    event_times = np.unique(time[event == 1])

    beta = np.zeros(p)
    for _ in range(max_iter):
        U = np.zeros(p)
        I = np.zeros((p, p))

        for t in event_times:
            risk_mask = time >= t
            X_risk = X_centered[risk_mask]
            n_risk = X_risk.shape[0]
            if n_risk == 0:
                continue

            exp_xb = np.exp(X_risk @ beta)
            s0 = exp_xb.sum()
            # Weighted mean
            s1 = (X_risk.T @ exp_xb) / s0
            # Weighted covariance
            s2 = np.zeros((p, p))
            for k in range(n_risk):
                diff = X_risk[k] - s1
                s2 += exp_xb[k] * np.outer(diff, diff)
            s2 /= s0

            events_mask = (time == t) & (event == 1)
            X_events = X_centered[events_mask]
            n_events = X_events.shape[0]
            if n_events == 0:
                continue

            U += X_events.sum(axis=0) - n_events * s1
            I += n_events * s2

        try:
            delta = np.linalg.solve(I + 1e-8 * np.eye(p), U)
        except np.linalg.LinAlgError:
            delta = np.linalg.pinv(I + 1e-8 * np.eye(p)) @ U

        beta += delta
        if np.max(np.abs(delta)) < tol:
            break

    try:
        var = np.diag(np.linalg.inv(I + 1e-8 * np.eye(p)))
    except np.linalg.LinAlgError:
        var = np.diag(np.linalg.pinv(I + 1e-8 * np.eye(p)))

    se = np.sqrt(np.maximum(var, 0))
    hr = np.exp(beta)
    z = beta / se
    p_values = 2 * (1 - stats.norm.cdf(np.abs(z)))

    ci_lower = np.exp(beta - 1.96 * se)
    ci_upper = np.exp(beta + 1.96 * se)

    return {
        "coefficients": beta,
        "hr": hr,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "p_values": p_values,
        "standard_errors": se,
    }


def stratify_by_covariate(
    df: pd.DataFrame,
    covariate: str,
    method: str = "median",
) -> np.ndarray:
    """Stratify patients by clinical covariate (median for continuous, value for categorical)."""
    values = df[covariate].values
    if covariate in ["age_at_diagnosis"]:
        if method == "median":
            cutoff = np.median(values)
        elif method == "tertile":
            cutoff = np.percentile(values, 66.67)
        elif method == "quartile":
            cutoff = np.percentile(values, 75)
        else:
            raise ValueError(f"Unknown method: {method}")
        return (values > cutoff).astype(int)
    elif covariate in ["gender"]:
        return values.astype(int)
    elif covariate in ["molecular_subtype"]:
        unique_vals = np.unique(values)
        if len(unique_vals) <= 2:
            return (values > np.median(values)).astype(int)
        ref = unique_vals[0]
        return (values != ref).astype(int)
    else:
        raise ValueError(f"Unknown covariate: {covariate}")


def run_univariate_survival_clinical(
    patient_df: pd.DataFrame,
    covariates: List[str],
    time: np.ndarray,
    event: np.ndarray,
) -> pd.DataFrame:
    """Run log-rank test and Cox PH for each clinical covariate."""
    results = []

    for cov in covariates:
        if cov not in patient_df.columns:
            continue

        # Stratify by covariate
        group = stratify_by_covariate(patient_df, cov)

        # Log-rank test
        chi2, p_lr = logrank_test(time, event, group)

        # Cox PH
        if cov == "molecular_subtype":
            # Dummy encode: one-hot with reference
            dummies = pd.get_dummies(patient_df[cov], prefix=cov, drop_first=True)
            X = dummies.values.astype(float)
        else:
            X = patient_df[cov].values.reshape(-1, 1).astype(float)
        cox = cox_ph_fit(time, event, X)

        results.append({
            "covariate": cov,
            "logrank_chi2": float(chi2),
            "logrank_p": float(p_lr),
            "cox_hr": float(cox["hr"][0]) if len(cox["hr"]) > 0 else 1.0,
            "cox_ci_lower": float(cox["ci_lower"][0]) if len(cox["ci_lower"]) > 0 else 1.0,
            "cox_ci_upper": float(cox["ci_upper"][0]) if len(cox["ci_upper"]) > 0 else 1.0,
            "cox_p": float(cox["p_values"][0]) if len(cox["p_values"]) > 0 else 1.0,
            "n_high": int(group.sum()),
            "n_low": int((1 - group).sum()),
        })

    return pd.DataFrame(results)


def run_multivariate_survival_clinical(
    patient_df: pd.DataFrame,
    covariates: List[str],
    time: np.ndarray,
    event: np.ndarray,
) -> Dict:
    """Run multivariate Cox PH with clinical covariates."""
    # Build design matrix - handle categorical
    X_df = pd.get_dummies(patient_df[covariates], columns=["molecular_subtype"], drop_first=True)
    feature_names = X_df.columns.tolist()
    X = X_df.values.astype(float)
    # Handle missing
    if np.isnan(X).any():
        X = np.nan_to_num(X, nan=np.nanmean(X))

    cox = cox_ph_fit(time, event, X)

    return {
        "features": feature_names,
        "coefficients": dict(zip(feature_names, cox["coefficients"])),
        "hr": dict(zip(feature_names, cox["hr"])),
        "ci_lower": dict(zip(feature_names, cox["ci_lower"])),
        "ci_upper": dict(zip(feature_names, cox["ci_upper"])),
        "p_values": dict(zip(feature_names, cox["p_values"])),
    }


# --------------------------------------------------------------------------- #
# Visualization
# --------------------------------------------------------------------------- #
def plot_km_panels(
    patient_df: pd.DataFrame,
    covariates: List[str],
    time: np.ndarray,
    event: np.ndarray,
    output_path: str = str(OUTPUT_DIR / "km_survival_curves.png"),
) -> None:
    """Generate multi-panel Kaplan-Meier survival curves for clinical covariates."""
    n_cov = len(covariates)
    n_cols = min(2, n_cov)
    n_rows = (n_cov + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
    if n_cov == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for idx, cov in enumerate(covariates):
        ax = axes[idx]
        
        group = stratify_by_covariate(patient_df, cov)

        # High group
        mask_high = group == 1
        t_high, s_high = kaplan_meier(time[mask_high], event[mask_high])
        ax.step(t_high, s_high, where='post', label=f'High (n={mask_high.sum()})',
                color='#e74c3c', linewidth=2)

        # Low group
        mask_low = group == 0
        t_low, s_low = kaplan_meier(time[mask_low], event[mask_low])
        ax.step(t_low, s_low, where='post', label=f'Low (n={mask_low.sum()})',
                color='#3498db', linewidth=2)

        # Log-rank test
        chi2, p = logrank_test(time, event, group)
        ax.text(0.02, 0.02, f'Log-rank p = {p:.4f}',
                transform=ax.transAxes, fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        ax.set_xlabel('Survival Time (days)', fontsize=11)
        ax.set_ylabel('Survival Probability', fontsize=11)
        ax.set_title(f'{cov} (High vs Low)', fontsize=12, fontweight='bold')
        ax.legend(loc='lower left', fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 1.05)

    # Hide unused subplots
    for idx in range(len(covariates), len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle('Kaplan-Meier Survival Curves by Clinical Covariate\n'
                 '(High vs Low, median split)',
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[PLOT] Kaplan-Meier panels saved to {output_path}")


def plot_forest(
    multivariate_results: Dict,
    output_path: str = str(OUTPUT_DIR / "forest_plot.png"),
) -> None:
    """Generate forest plot of multivariate Cox HRs."""
    features = multivariate_results["features"]
    hr = multivariate_results["hr"]
    ci_lower = multivariate_results["ci_lower"]
    ci_upper = multivariate_results["ci_upper"]
    p_values = multivariate_results["p_values"]

    # Convert dicts to ordered arrays
    hr_arr = np.array([hr[f] for f in features])
    lo_arr = np.array([ci_lower[f] for f in features])
    hi_arr = np.array([ci_upper[f] for f in features])
    p_arr = np.array([p_values[f] for f in features])

    n = len(features)
    y_pos = np.arange(n)

    fig, ax = plt.subplots(figsize=(8, max(4, n * 0.4)))

    colors = ['#e74c3c' if p < 0.05 else '#95a5a6' for p in p_arr]

    for i, (feat, h, lo, hi, p, c) in enumerate(zip(features, hr_arr, lo_arr, hi_arr, p_arr, colors)):
        ax.plot([lo, hi], [i, i], color=c, linewidth=2)
        ax.plot(h, i, 'o', color=c, markersize=6)
        ax.text(hi + 0.05, i, f'{h:.2f} ({lo:.2f}-{hi:.2f})',
                va='center', fontsize=9)
        ax.text(0.98, i, f'{feat}\np={p:.3f}',
                ha='right', va='center', fontsize=9,
                transform=ax.get_yaxis_transform())

    ax.axvline(x=1.0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('Hazard Ratio (95% CI)', fontsize=11)
    ax.set_title('Multivariate Cox PH: Clinical Covariates (Age, Gender, Subtype)',
                 fontsize=12, fontweight='bold')
    ax.set_yticks([])
    ax.set_xlim(0.5, max(hi_arr) * 1.3)
    ax.grid(True, axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[PLOT] Forest plot saved to {output_path}")


# --------------------------------------------------------------------------- #
# Main Pipeline
# --------------------------------------------------------------------------- #
def main():
    print("=" * 60)
    print("MONTH 5 WEEK 2: SURVIVAL ANALYSIS SUITE (Real TCGA-GBM Data)")
    print("=" * 60)

    # Load clinical cohort (real TCGA-GBM data from script 35)
    print("\n[LOAD] Reading clinical_mapped_cohort.csv...")
    patient_df = pd.read_csv(str(OUTPUT_DIR / "clinical_mapped_cohort.csv"))
    print(f"  Cohort: {patient_df.shape[0]} patients, {patient_df.shape[1]} features")

    # Clinical covariates for survival analysis
    clinical_covariates = ["age_at_diagnosis", "gender", "molecular_subtype"]
    print(f"\n[COVARIATES] Clinical covariates: {clinical_covariates}")

    # Fix column names
    time = patient_df["os_time_days"].values
    event = patient_df["os_event"].values

    print(f"\n  Survival time range: {time.min():.0f} - {time.max():.0f} days")
    print(f"  Events: {event.sum()}/{len(event)} ({100*event.mean():.1f}%)")
    print(f"  Median survival: {np.median(time[event == 1]):.1f} days")

    # Univariate survival on clinical covariates
    print("\n[SURVIVAL] Univariate analysis (log-rank + Cox)...")
    univariate_results = []
    for cov in clinical_covariates:
        if cov not in patient_df.columns:
            print(f"  {cov}: NOT IN DATA")
            continue

        # Stratify by covariate
        group = stratify_by_covariate(patient_df, cov)
        
        # Log-rank test
        chi2, p_lr = logrank_test(time, event, group)

        # Cox PH
        if cov == "molecular_subtype":
            # Dummy encode: one-hot with reference
            dummies = pd.get_dummies(patient_df[cov], prefix=cov, drop_first=True)
            X = dummies.values.astype(float)
        else:
            X = patient_df[cov].values.reshape(-1, 1).astype(float)
        cox = cox_ph_fit(time, event, X)

        univariate_results.append({
            "covariate": cov,
            "logrank_chi2": float(chi2),
            "logrank_p": float(p_lr),
            "cox_hr": float(cox["hr"][0]) if len(cox["hr"]) > 0 else 1.0,
            "cox_ci_lower": float(cox["ci_lower"][0]) if len(cox["ci_lower"]) > 0 else 1.0,
            "cox_ci_upper": float(cox["ci_upper"][0]) if len(cox["ci_upper"]) > 0 else 1.0,
            "cox_p": float(cox["p_values"][0]) if len(cox["p_values"]) > 0 else 1.0,
            "n_high": int(group.sum()),
            "n_low": int((1 - group).sum()),
        })

        print(f"  {cov}: HR={cox['hr'][0]:.3f} "
              f"({cox['ci_lower'][0]:.3f}-{cox['ci_upper'][0]:.3f}), "
              f"log-rank p={p_lr:.4f}, Cox p={cox['p_values'][0]:.4f}")

    uni_df = pd.DataFrame(univariate_results)
    print(uni_df.to_string(index=False))

    # Multivariate Cox
    print("\n[SURVIVAL] Multivariate Cox PH (age + gender + molecular subtype)...")
    X_df = pd.get_dummies(patient_df[["age_at_diagnosis", "gender", "molecular_subtype"]], 
                          columns=["molecular_subtype"], drop_first=True)
    feature_names = X_df.columns.tolist()
    X = X_df.values.astype(float)
    if np.isnan(X).any():
        X = np.nan_to_num(X, nan=np.nanmean(X))

    cox = cox_ph_fit(time, event, X)
    multivariate_results = {
        "features": feature_names,
        "coefficients": dict(zip(feature_names, cox["coefficients"])),
        "hr": dict(zip(feature_names, cox["hr"])),
        "ci_lower": dict(zip(feature_names, cox["ci_lower"])),
        "ci_upper": dict(zip(feature_names, cox["ci_upper"])),
        "p_values": dict(zip(feature_names, cox["p_values"])),
    }

    for feat in feature_names:
        print(f"  {feat}: HR={cox['hr'][feature_names.index(feat)]:.3f} "
              f"({cox['ci_lower'][feature_names.index(feat)]:.3f}-{cox['ci_upper'][feature_names.index(feat)]:.3f}), "
              f"p={cox['p_values'][feature_names.index(feat)]:.4f}")

    # Generate KM plots for clinical covariates
    print("\n[PLOTS] Generating visual artifacts...")
    OUTPUT_DIR.mkdir(exist_ok=True)
    plot_km_panels(patient_df, clinical_covariates, time, event,
                   str(OUTPUT_DIR / "km_survival_curves.png"))
    plot_forest(multivariate_results, str(OUTPUT_DIR / "forest_plot.png"))

    # Export summaries
    print("\n[EXPORT] Saving statistical summaries...")
    uni_df = pd.DataFrame(univariate_results)
    uni_df.to_csv(str(OUTPUT_DIR / "univariate_survival.csv"), index=False)
    
    with open(str(OUTPUT_DIR / "multivariate_survival.json"), "w") as f:
        json.dump(multivariate_results, f, indent=2)

    # Combined summary JSON
    summary = {
        "univariate": uni_df.to_dict(orient="records"),
        "multivariate": multivariate_results,
        "cohort_size": int(patient_df.shape[0]),
        "n_events": int(event.sum()),
        "median_survival_days": float(np.median(time[event == 1])),
        "covariates_analyzed": clinical_covariates,
    }
    with open(str(OUTPUT_DIR / "survival_stats_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print("\n[SUCCESS] Month 5 Week 2 Complete: Survival Analysis (Real TCGA-GBM)")
    print("  - output/km_survival_curves.png (clinical covariates)")
    print("  - output/forest_plot.png (multivariate Cox)")
    print("  - output/univariate_survival.csv")
    print("  - output/multivariate_survival.json")
    print("  - output/survival_stats_summary.json")


if __name__ == "__main__":
    main()