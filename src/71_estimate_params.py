#!/usr/bin/env python3
"""
Estimate per-patient tumor growth parameters (ρ, D) from longitudinal
MU-Glioma-Post tumor volumes.

Model: dV/dt = ρ * V * (1 - V/K)  (logistic growth)
For early-stage, exponential: V(t) = V0 * exp(ρ*t)
Fit ρ from log-linear regression on log(V) vs t.

Timepoints → approximate days:
  Timepoint_1 → 0 days
  Timepoint_2 → ~90 days
  Timepoint_3 → ~180 days
  Timepoint_4 → ~270 days
  Timepoint_5 → ~360 days
  Timepoint_6 → ~450 days
(assumes ~90-day intervals; adjust if metadata available)
"""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd

from pathlib import Path as _Path
PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"

# Approximate days for each timepoint (90-day spacing)
TIMEPOINT_DAYS = {
    "Timepoint_1": 0,
    "Timepoint_2": 90,
    "Timepoint_3": 180,
    "Timepoint_4": 270,
    "Timepoint_5": 360,
    "Timepoint_6": 450,
}


def fit_exponential_growth(days, volumes):
    """Fit log(V) = log(V0) + ρ*t. Returns (V0, rho, r_squared)."""
    if len(days) < 2:
        return None
    days = np.array(days, dtype=float)
    volumes = np.array(volumes, dtype=float)
    log_v = np.log(volumes)
    coeffs = np.polyfit(days, log_v, 1)
    rho = coeffs[0]
    v0 = np.exp(coeffs[1])
    predicted = np.polyval(coeffs, days)
    ss_res = np.sum((log_v - predicted) ** 2)
    ss_tot = np.sum((log_v - log_v.mean()) ** 2)
    r_squared = 1 - ss_res / max(ss_tot, 1e-10)
    return v0, rho, r_squared


def estimate_D_from_radius(growth_rate_daily, rho_per_day):
    """
    For a spherical tumor, radius grows as r(t) = r0 * exp(rho*t/3).
    FK wave speed c = 2*sqrt(D*rho) → D = c²/(4*rho).
    Approximate c from volume growth: c ≈ (delta_r) / delta_t.
    """
    # Simplified: assume Fisher-Kolmogorov relation
    # c ≈ 2*sqrt(D*rho) → D ≈ c²/(4*rho)
    # For a volumetric growth rate, we approximate c from the rate
    if rho_per_day <= 0:
        return None
    # Use empirical relation: c ≈ 0.1 * rho_volume (calibrated)
    # This is a placeholder — will need real calibration
    return None  # TODO


def main():
    print("=" * 60)
    print("PARAMETER ESTIMATION FROM LONGITUDINAL MU-GLIOMA")
    print("=" * 60)
    print()

    df = pd.read_csv(OUTPUT_DIR / "mu_glioma_volumes.csv")
    print(f"[LOAD] {len(df)} patients")

    results = []
    for _, row in df.iterrows():
        patient_id = row["patient_id"]
        days = []
        volumes = []
        timepoints_used = []
        for tp, day in TIMEPOINT_DAYS.items():
            if tp in df.columns and pd.notna(row[tp]) and row[tp] > 0:
                days.append(day)
                volumes.append(row[tp])
                timepoints_used.append(tp)

        if len(days) < 2:
            continue

        fit = fit_exponential_growth(days, volumes)
        if fit is None:
            continue
        v0, rho, r2 = fit

        results.append({
            "patient_id": patient_id,
            "n_timepoints": len(days),
            "timepoints": ",".join(timepoints_used),
            "V0_mm3": float(v0),
            "rho_per_day": float(rho),
            "r_squared": float(r2),
            "first_volume": float(volumes[0]),
            "last_volume": float(volumes[-1]),
            "trajectory": "growing" if rho > 0 else "shrinking",
        })

    out_df = pd.DataFrame(results)
    out_path = OUTPUT_DIR / "mu_glioma_params.csv"
    out_df.to_csv(out_path, index=False)
    print(f"\n[SAVE] {out_path}")
    print(f"  Patients with fit: {len(out_df)}")
    print()

    # Statistics
    print("[STATS] Growth rate ρ (per day):")
    print(f"  Mean: {out_df['rho_per_day'].mean():.5f}")
    print(f"  Median: {out_df['rho_per_day'].median():.5f}")
    print(f"  Std: {out_df['rho_per_day'].std():.5f}")
    print(f"  Range: [{out_df['rho_per_day'].min():.5f}, {out_df['rho_per_day'].max():.5f}]")
    print()
    print("[STATS] Growing vs shrinking:")
    print(out_df["trajectory"].value_counts())
    print()
    print("[STATS] R² distribution:")
    print(f"  Median R²: {out_df['r_squared'].median():.3f}")
    print(f"  Patients with R² > 0.8: {(out_df['r_squared'] > 0.8).sum()}/{len(out_df)}")
    print()
    print("First 10 rows:")
    print(out_df.head(10).to_string())

    print()
    print("[SUCCESS] Parameter estimation complete")


if __name__ == "__main__":
    main()