#!/usr/bin/env python3
"""Script 106: measured boundary disagreement between two segmenters (removes the ASSUMED 2 mm in script 103). Exploratory.

Declared before the run (this file is committed before it is executed):
  Data: LUMIERE, every scan with both automated label maps in native CT1 space
        (HD-GLIO-AUTO segmentation_CT1_origspace.nii.gz; DeepBraTumIA ct1_seg_mask.nii.gz).
  Label codes (found by matching voxel counts to the LUMIERE pyradiomics CSVs, 3 scans, within 1-3%):
        HD-GLIO-AUTO 1 = non-enhancing, 2 = contrast-enhancing; DeepBraTumIA 1 = necrosis, 2 = contrast-enhancing, 3 = edema.
  Regions compared: enhancing (HD 2 vs DBT 2); whole abnormality (HD 1|2 vs DBT 1|2|3).
  Per scan, both native and after nearest-neighbour resampling to 2 mm isotropic (the MU analysis grid):
        Dice; average symmetric surface distance (ASSD, mm); 95th percentile Hausdorff (HD95, mm).
  Scans where either mask is empty for a region are counted and excluded for that region.
  Summary: median and IQR over scans, patient-bootstrap CI of the median (patient_stats.cluster_bootstrap, 10,000 draws),
        for all scans and for RANO-rated follow-ups (PD/SD/PR/CR, as script 104).
  Comparison: script 103's assumed noise is a 1-voxel (2 mm) shift plus boundary jitter; its ceilings are read from
        output/segmentation_perturbation.json and reported next to the measured 2 mm Dice.
Weakest points (stated before the run): two automated tools are not a scan-rescan pair; LUMIERE enhancing label is not
  MU core (labels 1+3); native voxel sizes differ between scans (resampling to 2 mm makes the Dice comparable to MU).
Input : data/external/lumiere/Imaging-v202211.zip (read in place), LUMIERE-ExpertRating-v202211.csv
Output: output/boundary_noise_measured.json, output/boundary_noise_measured_scans.csv (+ manifest)
"""
from __future__ import annotations

import gzip
import json
import sys
import time
import zipfile
from pathlib import Path as _Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
LUM = PROJECT_ROOT / "data" / "external" / "lumiere"
sys.path.insert(0, str(PROJECT_ROOT))

from src import patient_stats as ps  # noqa: E402
from src.run_manifest import write_run_manifest  # noqa: E402

ZIP = LUM / "Imaging-v202211.zip"
HD = "HD-GLIO-AUTO-segmentation/native/segmentation_CT1_origspace.nii.gz"
DBT = "DeepBraTumIA-segmentation/native/segmentation/ct1_seg_mask.nii.gz"
REGIONS = {"enhancing": ((2,), (2,)), "whole": ((1, 2), (1, 2, 3))}
GRID_MM = 2.0


def _load(z, name):
    img = nib.Nifti1Image.from_bytes(gzip.decompress(z.read(name)))
    return np.asarray(img.dataobj).astype(np.int16), np.asarray(img.header.get_zooms()[:3], float)


def _to_grid(mask, zooms):
    return ndimage.zoom(mask.astype(np.uint8), zooms / GRID_MM, order=0).astype(bool)


def _surface_dist(a, b, spacing):
    """Distances (mm) from the surface voxels of a to the surface of b, and back."""
    pad = 3
    idx = np.argwhere(a | b)
    lo = np.maximum(idx.min(0) - pad, 0)
    hi = np.minimum(idx.max(0) + pad + 1, a.shape)
    sl = tuple(slice(int(l), int(h)) for l, h in zip(lo, hi))
    a, b = a[sl], b[sl]
    sa = a & ~ndimage.binary_erosion(a)
    sb = b & ~ndimage.binary_erosion(b)
    da = ndimage.distance_transform_edt(~sb, sampling=spacing)[sa]
    db = ndimage.distance_transform_edt(~sa, sampling=spacing)[sb]
    return np.concatenate([da, db])


def _metrics(a, b, spacing):
    inter = int((a & b).sum())
    d = _surface_dist(a, b, spacing)
    return {"dice": 2 * inter / (a.sum() + b.sum()), "assd_mm": float(d.mean()), "hd95_mm": float(np.quantile(d, 0.95)),
            "vol_a_mm3": float(a.sum() * np.prod(spacing)), "vol_b_mm3": float(b.sum() * np.prod(spacing))}


def main():
    t0 = time.time()
    z = zipfile.ZipFile(ZIP)
    names = set(z.namelist())
    scans = sorted({n[: -len(HD)] for n in names if n.endswith(HD)})
    rows, n_missing = [], 0
    for k, base in enumerate(scans):
        if base + DBT not in names:
            n_missing += 1
            continue
        _, pid, tp, _ = base.split("/")
        hd, zh = _load(z, base + HD)
        dbt, zd = _load(z, base + DBT)
        if hd.shape != dbt.shape or not np.allclose(zh, zd, atol=1e-3):
            rows.append({"Patient": pid, "tp": tp, "skip": "grid mismatch"})
            continue
        for reg, (lh, ld) in REGIONS.items():
            a, b = np.isin(hd, lh), np.isin(dbt, ld)
            r = {"Patient": pid, "tp": tp, "region": reg}
            if not a.any() or not b.any():
                r["skip"] = "empty " + ("hd" if not a.any() else "dbt")
                rows.append(r)
                continue
            r.update({f"native_{m}": v for m, v in _metrics(a, b, zh).items()})
            a2, b2 = _to_grid(a, zh), _to_grid(b, zh)
            if a2.any() and b2.any():
                r.update({f"g2_{m}": v for m, v in _metrics(a2, b2, np.full(3, GRID_MM)).items()})
            r["voxel_mm"] = "x".join(f"{v:.2f}" for v in zh)
            rows.append(r)
        if (k + 1) % 50 == 0:
            print(f"[{time.time() - t0:5.0f}s] {k + 1}/{len(scans)}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "boundary_noise_measured_scans.csv", index=False)
    rat = pd.read_csv(LUM / "LUMIERE-ExpertRating-v202211.csv")
    rat.columns = ["Patient", "tp", "lt3", "nonmeas", "rating", "rationale"]
    rated = set(map(tuple, rat[rat["rating"].isin({"PD", "SD", "PR", "CR"})][["Patient", "tp"]].to_numpy()))
    ok = df[df["skip"].isna()] if "skip" in df else df
    res = {"script": "106_boundary_noise_measured", "status": "exploratory (declared before run)",
           "n_scans_listed": len(scans), "n_scans_without_dbt": n_missing,
           "skips": df["skip"].value_counts().to_dict() if "skip" in df else {},
           "label_codes": {"hd_glio_auto": {"1": "non-enhancing", "2": "contrast-enhancing"},
                           "deepbratumia": {"1": "necrosis", "2": "contrast-enhancing", "3": "edema"}},
           "n_boot": ps.DEFAULT_N_BOOT, "summary": {}}
    for reg in REGIONS:
        for subset, d in (("all_scans", ok[ok["region"] == reg]),
                          ("rated_followups", ok[(ok["region"] == reg) & [(p, t) in rated for p, t in zip(ok["Patient"], ok["tp"])]])):
            block = {"n_scans": int(len(d)), "n_patients": int(d["Patient"].nunique())}
            for col in ("native_dice", "native_assd_mm", "native_hd95_mm", "g2_dice", "g2_assd_mm", "g2_hd95_mm"):
                x = d[col].dropna()
                if len(x) == 0:
                    continue
                lo, hi = ps.cluster_bootstrap(d.dropna(subset=[col]), lambda g, c=col: g[c].median(), patient="Patient")
                block[col] = {"median": float(x.median()), "q25": float(x.quantile(0.25)), "q75": float(x.quantile(0.75)),
                              "median_ci95_patient_bootstrap": [float(lo), float(hi)]}
            res["summary"].setdefault(reg, {})[subset] = block
    pert = json.loads((OUT / "segmentation_perturbation.json").read_text())
    res["script103_assumed_ceilings_gbm_only"] = pert["A_ceiling_patient_mean_dice"]["gbm_only"]
    res["weakest_points"] = [
        "Two automated tools are not a scan-rescan pair: this measures method disagreement, which includes systematic bias.",
        "LUMIERE enhancing label is not MU core (labels 1+3, includes necrosis); HD-GLIO-AUTO has no necrosis label, so a core match is not possible.",
        "Native voxel sizes differ between scans; the 2 mm resampled metrics are the ones comparable to MU."]
    (OUT / "boundary_noise_measured.json").write_text(json.dumps(res, indent=1, default=float))
    write_run_manifest("boundary_noise_106", OUT / "boundary_noise_measured.manifest.json",
                       script="src/106_boundary_noise_measured.py", seed=ps.DEFAULT_SEED,
                       config={"grid_mm": GRID_MM, "regions": {k: [list(v[0]), list(v[1])] for k, v in REGIONS.items()},
                               "n_boot": ps.DEFAULT_N_BOOT},
                       inputs=[LUM / "LUMIERE-ExpertRating-v202211.csv", OUT / "segmentation_perturbation.json"],
                       dataset="LUMIERE", patient_split="none (descriptive)",
                       primary_endpoint="median inter-tool ASSD (mm) and Dice at 2 mm, enhancing region")
    for reg in REGIONS:
        s = res["summary"][reg]["all_scans"]
        print(reg, {k: round(v["median"], 3) for k, v in s.items() if isinstance(v, dict)}, "n", s["n_scans"])


if __name__ == "__main__":
    main()
