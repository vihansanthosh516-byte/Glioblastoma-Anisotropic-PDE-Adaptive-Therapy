#!/usr/bin/env python3
"""Cohort-level forward validation using joint treatment-aware inverse estimation.

Loads all MU-Glioma-Post patients with >=3 timepoints, fits rho/D jointly across
all adjacent intervals (with a physiological prior), then validates forward by
refitting on the first N-1 timepoints and predicting the final one.

Outputs:
    output/real_patient_validation.csv  per-patient metrics
    output/real_patient_validation_summary.json  aggregate statistics
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

CLINICAL_EXCEL = ROOT / "data" / "tcia" / "MU-Glioma-Post_ClinicalData-July2025.xlsx"
COHORT_JSON = ROOT / "output" / "mu_glioma_cohort.json"
DATA_ROOT = ROOT / "data" / "tcia" / "MU-Glioma-Post"

TP_COLUMNS = {
    1: "Number of Days from Diagnosis to 1st MRI (Timepoint_1) ",
    2: "Number of Days from Diagnosis to 2nd MRI (Timepoint_2) ",
    3: "Number of Days from Diagnosis to 3rd MRI (Timepoint_3) ",
    4: "Number of Days from Diagnosis to 4th MRI (Timepoint_4) ",
    5: "Number of Days from Diagnosis to 5th MRI (Timepoint_5) ",
    6: "Number of Days from Diagnosis to 6th MRI (Timepoint_6) ",
}

TREATMENT_RANGE_KEYS = [
    (
        "Number of days from Diagnosis to Radiation Therapy Start date",
        "Number of days from Diagnosis to Radiation Therapy end date",
    ),
    (
        " Number of days from Diagnosis to Initial Chemo Therapy Start date",
        " Number of days from Diagnosis to Initial Chemo Therapy end date",
    ),
    (
        "Number of Days from Diagnosis to Starting Additional Therapy ",
        "Number of Days from Diagnosis to Complete Additional Therapy ",
    ),
    (
        "Number of Days from Diagnosis to Starting 2nd_Additional Therapy ",
        "Number of Days from Dagnosis to Complete 2nd_Additional Therapy ",
    ),
    (
        "Number of Days from Diagnosis to Start Immunotherapy ",
        "Number of Days from Diagnosis to Complete Immunotherapy ",
    ),
    (
        "Number of Days from Diagnosis to Start Other Additional Therapy ",
        "Number of Days from Diagnosis to Complete Other Additional Therapy ",
    ),
]


def _load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "src" / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


_inv = _load_module("inv_est", "51_inverse_parameter_estimation.py")
_joint = _load_module("joint_est", "51_joint_inverse_estimation.py")

from treatment_aware_pde import TreatmentSchedule, treatment_aware_ode_model  # noqa: E402
from radiation_model import RadiationSchedule  # noqa: E402


def get_value(row: pd.Series, key: str) -> Optional[float]:
    value = row.get(key)
    if value is None or pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def load_volume_mm3(patient_id: str, timepoint: int) -> Optional[float]:
    path = DATA_ROOT / patient_id / f"Timepoint_{timepoint}" / f"{patient_id}_Timepoint_{timepoint}_tumorMask.nii.gz"
    try:
        import nibabel as nib
        image = nib.load(str(path))
        voxel_vol = abs(float(np.linalg.det(image.affine[:3, :3])))
        return float(np.count_nonzero(image.get_fdata() > 0) * voxel_vol)
    except Exception:
        return None


def build_treatment_ranges(row: pd.Series) -> List[Tuple[float, float]]:
    ranges: List[Tuple[float, float]] = []
    surgery = get_value(row, "Number of days from Diagnosis to First surgery or procedure ")
    if surgery is not None:
        ranges.append((surgery, surgery))
    for start_key, end_key in TREATMENT_RANGE_KEYS:
        start = get_value(row, start_key)
        end = get_value(row, end_key)
        if start is not None and end is not None:
            ranges.append((min(start, end), max(start, end)))
    return ranges


def interval_is_treatment_free(ranges: Sequence[Tuple[float, float]], t0: float, t1: float) -> bool:
    for start, end in ranges:
        if start < t1 and end > t0:
            return False
    return True


def build_tmz_schedule(row: pd.Series) -> Tuple[float]:
    start = get_value(row, " Number of days from Diagnosis to Initial Chemo Therapy Start date")
    end = get_value(row, " Number of days from Diagnosis to Initial Chemo Therapy end date")
    if start is None or end is None:
        return ()
    return tuple(float(day) for day in range(int(start), int(end) + 1))


def build_radiation_schedule(row: pd.Series) -> Optional[RadiationSchedule]:
    start = get_value(row, "Number of days from Diagnosis to Radiation Therapy Start date")
    end = get_value(row, "Number of days from Diagnosis to Radiation Therapy end date")
    if start is None or end is None or end <= start:
        return None
    return RadiationSchedule(start, end)


def build_patient_schedule(row: pd.Series) -> TreatmentSchedule:
    return TreatmentSchedule(
        tmz_bolus_days=build_tmz_schedule(row),
        radiation=build_radiation_schedule(row),
    )


def validate_patient(
    patient_id: str,
    timepoints: Sequence[int],
    row: pd.Series,
) -> Optional[Dict[str, Any]]:
    volumes = [load_volume_mm3(patient_id, tp) for tp in timepoints]
    if any(v is None for v in volumes) or any(v <= 0 for v in volumes):
        return None

    days: List[float] = []
    for tp in timepoints:
        day = get_value(row, TP_COLUMNS[tp])
        if day is None:
            return None
        days.append(day)
    if not all(b > a for a, b in zip(days, days[1:])):
        return None

    ranges = build_treatment_ranges(row)
    schedule = build_patient_schedule(row)

    joint = _joint.estimate_joint_parameters(volumes, days, schedule=schedule, regularization=0.01)
    rho, D = joint["rho"], joint["D"]

    if len(volumes) >= 4:
        refit = _joint.estimate_joint_parameters(volumes[:-1], days[:-1], schedule=schedule, regularization=0.01)
    else:
        refit = joint
    pred = treatment_aware_ode_model(
        refit["rho"], refit["D"], volumes[-2], days[-1] - days[-2],
        schedule=schedule, start_day=days[-2],
    )

    actual = volumes[-1]
    mae = abs(pred - actual)
    pct = mae / max(actual, 1.0) * 100.0

    free_intervals = [
        (days[i], days[i + 1])
        for i in range(len(days) - 1)
        if interval_is_treatment_free(ranges, days[i], days[i + 1])
    ]

    return {
        "patient_id": patient_id,
        "n_timepoints": len(volumes),
        "n_treatment_free_intervals": len(free_intervals),
        "rho": rho,
        "D": D,
        "rho_refit": refit["rho"],
        "D_refit": refit["D"],
        "predicted_volume_mm3": pred,
        "actual_volume_mm3": actual,
        "mae_mm3": mae,
        "percent_error": pct,
        "objective": joint["objective"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Cohort forward validation")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--output", type=Path, default=ROOT / "output" / "real_patient_validation.csv")
    parser.add_argument("--summary", type=Path, default=ROOT / "output" / "real_patient_validation_summary.json")
    args = parser.parse_args()

    with open(COHORT_JSON, "r", encoding="utf-8") as fh:
        cohort = json.load(fh)

    df = pd.read_excel(CLINICAL_EXCEL, sheet_name="MU Glioma Post")
    df["Patient_ID"] = df["Patient_ID"].astype(str)

    patients = [p for p in cohort if p["n_timepoints"] >= 3]
    if args.limit is not None:
        patients = patients[: args.limit]

    results: List[Dict[str, Any]] = []
    for patient in patients:
        pid = patient["patient_id"]
        tps = sorted(tp["number"] for tp in patient["timepoints"])
        matches = df[df["Patient_ID"] == pid]
        if matches.empty:
            continue
        row = matches.iloc[0]
        result = validate_patient(pid, tps, row)
        if result is not None:
            results.append(result)
            flag = "OK " if result["percent_error"] < 30 else "BAD"
            print(f"[{flag}] {pid}: pred={result['predicted_volume_mm3']:.0f} "
                  f"actual={result['actual_volume_mm3']:.0f} "
                  f"({result['percent_error']:.1f}%) rho={result['rho']:.4f} D={result['D']:.4f}")

    if not results:
        print("No patients validated")
        return 1

    df_out = pd.DataFrame(results)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(args.output, index=False)

    errors = np.asarray([r["percent_error"] for r in results])
    summary = {
        "n_patients": len(results),
        "mean_percent_error": float(np.mean(errors)),
        "median_percent_error": float(np.median(errors)),
        "within_30pct": int(np.sum(errors < 30)),
        "within_10pct": int(np.sum(errors < 10)),
        "within_30pct_frac": float(np.mean(errors < 30)),
        "within_10pct_frac": float(np.mean(errors < 10)),
        "mean_mae_mm3": float(np.mean([r["mae_mm3"] for r in results])),
        "median_mae_mm3": float(np.median([r["mae_mm3"] for r in results])),
        "mean_rho": float(np.mean([r["rho"] for r in results])),
        "mean_D": float(np.mean([r["D"] for r in results])),
        "mean_treatment_free_intervals": float(np.mean([r["n_treatment_free_intervals"] for r in results])),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"\nSaved per-patient results to {args.output}")
    print(f"Saved summary to {args.summary}")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())