#!/usr/bin/env python3
"""Historical within-cohort temporal comparison of assumed treatment models.

Real imaging volumes are inputs; this is not an independent locked validation.
Model regularization/treatment effects are assumptions, not clinical evidence.
Output when explicitly run: output/naive_vs_treatment_comparison.json.
Importing the module performs no data loading or fitting.
"""
import importlib.util
import json
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from treatment_aware_pde import TreatmentSchedule, treatment_aware_ode_model  # noqa: E402

TP = {i: f"Number of Days from Diagnosis to {o} MRI (Timepoint_{i}) "
      for i, o in enumerate(["1st", "2nd", "3rd", "4th", "5th", "6th"], 1)}


def gv(r, k):
    v = r.get(k)
    if v is None or pd.isna(v):
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def load_vol(pid, tp):
    try:
        p = ROOT / f"data/tcia/MU-Glioma-Post/{pid}/Timepoint_{tp}/{pid}_Timepoint_{tp}_tumorMask.nii.gz"
        img = nib.load(str(p))
        return np.sum(img.get_fdata() > 0) * abs(np.linalg.det(img.affine[:3, :3]))
    except Exception:
        return None


def build_sched(row):
    ts = gv(row, " Number of days from Diagnosis to Initial Chemo Therapy Start date")
    te = gv(row, " Number of days from Diagnosis to Initial Chemo Therapy end date")
    tmz = tuple(float(d) for d in range(int(ts), int(te) + 1)) if (ts is not None and te is not None) else ()
    return TreatmentSchedule(tmz_bolus_days=tmz)


def main():
    spec = importlib.util.spec_from_file_location("joint", ROOT / "src" / "51_joint_inverse_estimation.py")
    joint = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(joint)

    df = pd.read_excel(ROOT / "data/tcia/MU-Glioma-Post_ClinicalData-July2025.xlsx", sheet_name="MU Glioma Post")
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    with open(ROOT / "output/mu_glioma_cohort.json") as fh:
        cohort = json.load(fh)

    err_naive, err_aware = [], []
    counts = {"n": 0, "both_ok": 0}

    for p in cohort:
        pid = p["patient_id"]
        tps = sorted(tp["number"] for tp in p["timepoints"])
        if len(tps) < 3:
            continue
        vols = [load_vol(pid, tp) for tp in tps]
        if any(v is None for v in vols):
            continue
        row = df[df["Patient_ID"] == pid]
        if row.empty:
            continue
        r = row.iloc[0]
        days = []
        for tp in tps:
            d = gv(r, TP[tp])
            if d is None:
                break
            days.append(d)
        if len(days) < 3 or not all(b > a for a, b in zip(days, days[1:])):
            continue
        counts["n"] += 1
        sched = build_sched(r)
        fit_days = days[:-1]
        fit_vols = vols[:-1]
        target = vols[-1]
        dt_pred = days[-1] - days[-2]

        if len(fit_vols) >= 3:
            est_n = joint.estimate_joint_parameters(fit_vols, fit_days, schedule=None, regularization=0.01)
            pred_n = treatment_aware_ode_model(est_n["rho"], est_n["D"], fit_vols[-1], dt_pred, start_day=fit_days[-1])
            est_a = joint.estimate_joint_parameters(fit_vols, fit_days, schedule=sched, regularization=0.01)
            pred_a = treatment_aware_ode_model(
                est_a["rho"], est_a["D"], fit_vols[-1], dt_pred, schedule=sched, start_day=fit_days[-1]
            )
            counts["both_ok"] += 1
            err_naive.append(abs(pred_n - target) / max(target, 1) * 100)
            err_aware.append(abs(pred_a - target) / max(target, 1) * 100)

    err_naive = np.array(err_naive)
    err_aware = np.array(err_aware)
    print(f"N patients with >=4 timepoints (proper refit): {counts['n']}, both_ok: {counts['both_ok']}")
    if len(err_naive):
        print(f"NAIVE: median {np.median(err_naive):.1f}% | within30 {100*np.mean(err_naive<30):.1f}% | mean {np.mean(err_naive):.1f}%")
        print(f"AWARE: median {np.median(err_aware):.1f}% | within30 {100*np.mean(err_aware<30):.1f}% | mean {np.mean(err_aware):.1f}%")
        print(f"IMPROVEMENT: {np.median(err_naive)-np.median(err_aware):.1f}pp median, {100*(np.mean(err_naive<30)-np.mean(err_aware<30)):.1f}pp within30")

    output_path = ROOT / "output" / "naive_vs_treatment_comparison.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps({
        "label": "Exploratory within-cohort temporal comparison; assumed treatment model; not locked clinical validation",
        "counts": counts,
        "relative_error_percent_naive": err_naive.tolist(),
        "relative_error_percent_treatment_aware": err_aware.tolist(),
        "regularization_assumed": 0.01,
    }, indent=2), encoding="utf-8")
    if output_path.stat().st_size == 0:
        raise RuntimeError(f"Empty comparison artifact: {output_path}")


if __name__ == "__main__":
    main()
