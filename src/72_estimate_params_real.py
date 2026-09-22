#!/usr/bin/env python3
"""
Estimate per-patient tumor growth parameters using REAL days from diagnosis.
"""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd

from pathlib import Path as _Path
PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"


def fit_exponential_growth(days, volumes):
    days = np.array(days, dtype=float)
    volumes = np.array(volumes, dtype=float)
    if len(days) < 2:
        return None
    log_v = np.log(volumes)
    coeffs = np.polyfit(days, log_v, 1)
    rho = coeffs[0]
    v0 = np.exp(coeffs[1])
    predicted = np.polyval(coeffs, days)
    ss_res = np.sum((log_v - predicted) ** 2)
    ss_tot = np.sum((log_v - log_v.mean()) ** 2)
    r_squared = 1 - ss_res / max(ss_tot, 1e-10)
    return v0, rho, r_squared


def main():
    print("=" * 60)
    print("REAL PARAMETER ESTIMATION (using actual MRI dates)")
    print("=" * 60)
    print()

    cohort_path = OUTPUT_DIR / "mu_glioma_cohort.json"
    with open(cohort_path) as f:
        cohort = json.load(f)
    print(f"[LOAD] {len(cohort)} patients from {cohort_path}")

    results = []
    n_skipped = 0
    for patient in cohort:
        valid = []
        for tp in patient["timepoints"]:
            if tp["volume_mm3"] is not None and tp["day_from_diagnosis"] is not None:
                valid.append((tp["day_from_diagnosis"], tp["volume_mm3"], tp["number"]))

        valid.sort()

        if len(valid) < 2:
            n_skipped += 1
            continue

        days = [v[0] for v in valid]
        volumes = [v[1] for v in valid]

        fit = fit_exponential_growth(days, volumes)
        if fit is None:
            n_skipped += 1
            continue

        v0, rho, r2 = fit

        results.append({
            "patient_id": patient["patient_id"],
            "n_timepoints": len(days),
            "days_from_dx": ",".join(str(int(d)) for d in days),
            "volumes_mm3": ",".join(str(int(v)) for v in volumes),
            "V0_mm3": float(v0),
            "rho_per_day": float(rho),
            "doubling_time_days": float(np.log(2) / abs(rho)) if abs(rho) > 1e-10 else None,
            "D_mm2_per_day": float(abs(rho) / 10.0),
            "r_squared": float(r2),
            "trajectory": "growing" if rho > 0 else "shrinking",
        })

    out_df = pd.DataFrame(results)
    out_path = OUTPUT_DIR / "mu_glioma_params_real.csv"
    out_df.to_csv(out_path, index=False)
    print(f"\n[SAVE] {out_path}")
    print(f"  Patients with valid fit: {len(out_df)} (skipped {n_skipped})")

    if len(out_df) == 0:
        print("\n[ERROR] No patients with 2+ valid timepoints")
        return

    print()
    print("[STATS] Growth rate rho (per day, from real MRI dates):")
    print(f"  Mean:   {out_df['rho_per_day'].mean():+.5f}")
    print(f"  Median: {out_df['rho_per_day'].median():+.5f}")
    print(f"  Std:    {out_df['rho_per_day'].std():.5f}")
    print(f"  Range:  [{out_df['rho_per_day'].min():+.5f}, {out_df['rho_per_day'].max():+.5f}]")

    print()
    print("[STATS] Doubling time (days):")
    dt = out_df["doubling_time_days"].dropna()
    if len(dt) > 0:
        print(f"  Median: {dt.median():.0f} days")
        print(f"  Range:  [{dt.min():.0f}, {dt.max():.0f}] days")

    print()
    print("[STATS] Trajectory distribution:")
    print(out_df["trajectory"].value_counts().to_string())

    print()
    print("[STATS] R^2 distribution:")
    print(f"  Median R^2:              {out_df['r_squared'].median():.3f}")
    print(f"  Patients with R^2 > 0.8: {(out_df['r_squared'] > 0.8).sum()}/{len(out_df)}")

    print()
    print("Sample of 10 patients:")
    print(out_df.head(10)[["patient_id", "n_timepoints", "rho_per_day",
                            "doubling_time_days", "r_squared", "trajectory"]].to_string(index=False))

    summary = {
        "n_patients": int(len(out_df)),
        "n_skipped": n_skipped,
        "rho_stats": {
            "mean": float(out_df["rho_per_day"].mean()),
            "median": float(out_df["rho_per_day"].median()),
            "std": float(out_df["rho_per_day"].std()),
            "min": float(out_df["rho_per_day"].min()),
            "max": float(out_df["rho_per_day"].max()),
        },
        "doubling_time_stats": {
            "median_days": float(dt.median()) if len(dt) > 0 else None,
            "min_days": float(dt.min()) if len(dt) > 0 else None,
            "max_days": float(dt.max()) if len(dt) > 0 else None,
        },
        "trajectory_counts": out_df["trajectory"].value_counts().to_dict(),
    }
    summary_path = OUTPUT_DIR / "mu_glioma_params_real_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n[SAVE] {summary_path}")
    print()
    print("[SUCCESS] Real parameter estimation complete")


if __name__ == "__main__":
    main()