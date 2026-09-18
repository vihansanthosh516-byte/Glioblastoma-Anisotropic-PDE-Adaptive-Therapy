#!/usr/bin/env python3
"""
Calibrate GBM growth rate rho from real MU-Glioma-Post longitudinal volumes.

rho = ln(V1/V0) / (t1 - t0)   [exponential growth phase assumption]

Only uses patients with >=2 timepoints, valid volumes (>0), and actual
clinical timing. Reports the cohort median rho. Honest real data only.
"""
import sys
import json
from pathlib import Path
import numpy as np

sys.path.insert(0, "src")
from src.mu_glioma_loader import load_patient_record, list_patient_ids

DATA_ROOT = Path("data/tcia/MU-Glioma-Post")
CLINICAL = Path("data/tcia/MU-Glioma-Post_ClinicalData-July2025.xlsx")

patient_ids = list_patient_ids(DATA_ROOT)
print(f"Total patients: {len(patient_ids)}")

pairs = []  # (patient_id, t0_day, t1_day, V0, V1, rho)
eligible = 0

for pid in sorted(patient_ids):
    rec = load_patient_record(pid, DATA_ROOT, CLINICAL, include_volumes=True)
    # Need >=2 timepoints with known volume and known timing
    tps = [(tp.day_from_diagnosis, tp.volume_mm3) for tp in rec.timepoints
           if tp.volume_mm3 is not None and tp.volume_mm3 > 0
           and tp.day_from_diagnosis is not None]
    if len(tps) < 2:
        continue
    eligible += 1
    # use earliest and latest timepoints
    tps.sort(key=lambda x: x[0])
    t0, v0 = tps[0]
    t1, v1 = tps[-1]
    dt = t1 - t0
    if dt <= 0 or v0 <= 0 or v1 <= 0:
        continue
    rho = np.log(v1 / v0) / dt
    pairs.append((pid, t0, t1, v0, v1, rho))

print(f"Patients with >=2 timed volume pairs: {eligible}")
print(f"Pairs with valid dt>0, V>0: {len(pairs)}")
print(f"\n{'Patient':14s} {'t0':>6s} {'t1':>6s} {'V0(m3)':>9s} {'V1(m3)':>9s} {'rho(/day)':>10s}")
for pid, t0, t1, v0, v1, rho in pairs:
    print(f"{pid:14s} {t0:6.0f} {t1:6.0f} {v0:9.1f} {v1:9.1f} {rho:10.5f}")

if pairs:
    rhos = np.array([p[5] for p in pairs])
    print("\n=== CALIBRATED RHO ===")
    print(f"median rho = {np.median(rhos):.5f} /day")
    print(f"mean   rho = {np.mean(rhos):.5f} /day")
    print(f"q25    rho = {np.percentile(rhos,25):.5f} /day")
    print(f"q75    rho = {np.percentile(rhos,75):.5f} /day")
    print(f"min/max   = {rhos.min():.5f} / {rhos.max():.5f} /day")
    print(f"negative rho count (shrinking tumors): {(rhos<0).sum()}/{len(rhos)}")
    with open("output/rho_calibration.json", "w") as f:
        json.dump({
            "n_pairs": len(pairs),
            "median_rho": float(np.median(rhos)),
            "mean_rho": float(np.mean(rhos)),
            "q25": float(np.percentile(rhos,25)),
            "q75": float(np.percentile(rhos,75)),
            "min": float(rhos.min()),
            "max": float(rhos.max()),
            "negative_count": int((rhos<0).sum()),
            "per_patient": [{"pid": p[0], "rho": p[5]} for p in pairs],
        }, f, indent=2)
    print("\nSaved output/rho_calibration.json")