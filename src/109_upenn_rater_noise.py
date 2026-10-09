#!/usr/bin/env python3
"""Script 109: expert-correction noise on UPENN-GBM (GRAND_PLAN L5/L6 substitute). Exploratory, declared before the run.

Why: no open glioma scan-rescan data exists (RIDER Neuro MRI and QIN-GBM are NIH controlled access). UPENN-GBM
(TCIA, CC BY 4.0, via NCI Imaging Data Commons) has the AIMI/BAMF AI segmentation of each scan and, for 220 scans
(55 patients), a version corrected by a board-certified radiologist (radiologist 3: 100 scans, radiologist 11: 120;
no scan was corrected by both). So this measures how much an expert changes a tool's mask on the SAME scan.
Labels: 1 necrosis, 2 edema, 3 enhancing lesion. Regions: core = necrosis + enhancing (same definition as MU labels
1 + 3), enhancing, whole = all three.
Per scan and region (both masks non-empty): Dice, ASSD (mm), HD95 (mm), |ln(V_ai / V_rad)|; also Dice after
nearest-neighbour resampling to 2 mm isotropic (the MU grid). Volume floor analogue to A6.5: median |ln ratio| / 2.
Summary: median, IQR, patient-bootstrap CI of the median (10,000 draws), overall and per radiologist.
Pairing: the AI and corrected SEG of a study must reference the same source series and have the same geometry
(orientation, pixel spacing, rows, columns); frames are matched by their ImagePositionPatient. Others are counted
and excluded.
Weakest points (stated before the run): expert correction of an AI mask is anchored on that mask, so it understates
independent rater disagreement; same scan, so scanner/registration noise between visits is not included; the AI
tool (AIMI/BAMF) differs from the MU and LUMIERE pipelines; 55 patients, 2 raters, no rater overlap.
Input : data/external/upenn_gbm/seg/ (DICOM-SEG from IDC), data/external/upenn_gbm/seg_series.csv
Output: output/upenn_rater_noise.json, output/upenn_rater_noise_scans.csv (+ manifest)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path as _Path

import numpy as np
import pandas as pd
import pydicom
from scipy import ndimage

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
UP = PROJECT_ROOT / "data" / "external" / "upenn_gbm"
sys.path.insert(0, str(PROJECT_ROOT))

from src import patient_stats as ps  # noqa: E402
from src.run_manifest import write_run_manifest  # noqa: E402

REGIONS = {"core": (1, 3), "enhancing": (3,), "whole": (1, 2, 3)}
GRID_MM = 2.0


def read_seg(path):
    d = pydicom.dcmread(str(path))
    sh = d.SharedFunctionalGroupsSequence[0]
    px = [float(x) for x in sh.PixelMeasuresSequence[0].PixelSpacing]
    orient = tuple(round(float(x), 4) for x in sh.PlaneOrientationSequence[0].ImageOrientationPatient)
    arr = d.pixel_array.reshape(int(d.NumberOfFrames), d.Rows, d.Columns).astype(bool)
    frames = []
    for i, fg in enumerate(d.PerFrameFunctionalGroupsSequence):
        seg = int(fg.SegmentIdentificationSequence[0].ReferencedSegmentNumber)
        pos = tuple(round(float(x), 2) for x in fg.PlanePositionSequence[0].ImagePositionPatient)
        frames.append((seg, pos, i))
    ref = d.ReferencedSeriesSequence[0].SeriesInstanceUID if "ReferencedSeriesSequence" in d else None
    return {"geom": (orient, tuple(round(p, 4) for p in px), d.Rows, d.Columns), "px": px, "orient": orient,
            "frames": frames, "arr": arr, "ref": ref}


def to_volumes(s, positions):
    """{segment: (Z, R, C) bool} on a shared ordered list of slice positions."""
    zi = {p: k for k, p in enumerate(positions)}
    vols = {}
    for seg, pos, i in s["frames"]:
        v = vols.setdefault(seg, np.zeros((len(positions), s["arr"].shape[1], s["arr"].shape[2]), bool))
        v[zi[pos]] |= s["arr"][i]
    return vols


def slice_order(positions, orient):
    n = np.cross(np.array(orient[:3]), np.array(orient[3:]))
    dist = {p: float(np.dot(n, p)) for p in positions}
    ordered = sorted(positions, key=dist.get)
    d = np.diff([dist[p] for p in ordered])
    return ordered, float(np.median(d)) if len(d) else 1.0


def surface_dist(a, b, spacing):
    idx = np.argwhere(a | b)
    lo, hi = np.maximum(idx.min(0) - 3, 0), np.minimum(idx.max(0) + 4, a.shape)
    sl = tuple(slice(int(x), int(y)) for x, y in zip(lo, hi))
    a, b = a[sl], b[sl]
    sa, sb = a & ~ndimage.binary_erosion(a), b & ~ndimage.binary_erosion(b)
    return np.concatenate([ndimage.distance_transform_edt(~sb, sampling=spacing)[sa],
                           ndimage.distance_transform_edt(~sa, sampling=spacing)[sb]])


def metrics(a, b, spacing):
    d = surface_dist(a, b, spacing)
    vv = float(np.prod(spacing))
    return {"dice": float(2 * (a & b).sum() / (a.sum() + b.sum())), "assd_mm": float(d.mean()),
            "hd95_mm": float(np.quantile(d, 0.95)), "v_ai_mm3": float(a.sum() * vv), "v_rad_mm3": float(b.sum() * vv),
            "abs_ln_ratio": float(abs(np.log(a.sum() / b.sum())))}


def main():
    ser = pd.read_csv(UP / "seg_series.csv")
    ser["who"] = ser["SeriesDescription"].str.extract(r"(AI segmentation|radiologist 11|radiologist 3)")[0]
    files = {p.parent.name.replace("SEG_", ""): p for p in (UP / "seg").rglob("*.dcm")}
    rows, skips = [], {}
    for study, g in ser.groupby("StudyInstanceUID"):
        rads = g[g["who"].str.startswith("radiologist", na=False)]
        ais = g[g["who"] == "AI segmentation"]
        if rads.empty or ais.empty:
            continue
        rad_row = rads.iloc[0]
        rad = read_seg(files[rad_row["SeriesInstanceUID"]])
        match = [read_seg(files[u]) for u in ais["SeriesInstanceUID"]]
        match = [s for s in match if s["ref"] == rad["ref"] and s["geom"] == rad["geom"]]
        if not match:
            skips["no AI SEG with same source series and geometry"] = skips.get("no AI SEG with same source series and geometry", 0) + 1
            continue
        ai = match[0]
        positions, dz = slice_order(sorted({p for _, p, _ in ai["frames"]} | {p for _, p, _ in rad["frames"]}), rad["orient"])
        va, vr = to_volumes(ai, positions), to_volumes(rad, positions)
        spacing = np.array([abs(dz), rad["px"][0], rad["px"][1]])
        for reg, labs in REGIONS.items():
            shape = (len(positions),) + ai["arr"].shape[1:]
            a = np.zeros(shape, bool)
            b = np.zeros(shape, bool)
            for lab in labs:
                a |= va.get(lab, np.zeros(shape, bool))
                b |= vr.get(lab, np.zeros(shape, bool))
            r = {"PatientID": rad_row["PatientID"], "study": study, "rater": rad_row["who"], "region": reg,
                 "spacing_mm": "x".join(f"{x:.2f}" for x in spacing)}
            if not a.any() or not b.any():
                r["skip"] = "empty " + ("ai" if not a.any() else "rad")
                rows.append(r)
                continue
            r.update(metrics(a, b, spacing))
            a2 = ndimage.zoom(a.astype(np.uint8), spacing / GRID_MM, order=0).astype(bool)
            b2 = ndimage.zoom(b.astype(np.uint8), spacing / GRID_MM, order=0).astype(bool)
            r["g2_dice"] = float(2 * (a2 & b2).sum() / max(a2.sum() + b2.sum(), 1))
            rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "upenn_rater_noise_scans.csv", index=False)
    ok = df[df["skip"].isna()] if "skip" in df else df
    res = {"script": "109_upenn_rater_noise", "status": "exploratory (declared before run)",
           "dataset": "UPENN-GBM (TCIA, CC BY 4.0) via NCI IDC; AIMI/BAMF AI vs radiologist-corrected DICOM-SEG",
           "n_studies_paired": int(ok["study"].nunique()), "n_patients": int(ok["PatientID"].nunique()),
           "skips_studies": skips, "skips_regions": df["skip"].value_counts().to_dict() if "skip" in df else {},
           "n_boot": ps.DEFAULT_N_BOOT, "summary": {}}
    for reg in REGIONS:
        for sub, d in [("all", ok[ok["region"] == reg])] + [(w, ok[(ok["region"] == reg) & (ok["rater"] == w)])
                                                            for w in ("radiologist 3", "radiologist 11")]:
            blk = {"n_scans": int(len(d)), "n_patients": int(d["PatientID"].nunique())}
            for col in ("dice", "g2_dice", "assd_mm", "hd95_mm", "abs_ln_ratio"):
                x = d[col].dropna()
                if len(x) < 2:
                    continue
                lo, hi = ps.cluster_bootstrap(d.dropna(subset=[col]), lambda g, c=col: g[c].median(), patient="PatientID")
                blk[col] = {"median": float(x.median()), "q25": float(x.quantile(0.25)), "q75": float(x.quantile(0.75)),
                            "median_ci95_patient_bootstrap": [float(lo), float(hi)]}
            if "abs_ln_ratio" in blk:
                blk["volume_floor_ln_A6.5_analogue"] = blk["abs_ln_ratio"]["median"] / 2
            res["summary"].setdefault(reg, {})[sub] = blk
    res["weakest_points"] = [
        "Correction of an AI mask is anchored on that mask: understates independent rater disagreement.",
        "Same scan: no between-visit scanner or registration noise.",
        "AIMI/BAMF tool differs from MU and LUMIERE pipelines; 55 patients, 2 radiologists, no overlap between them."]
    (OUT / "upenn_rater_noise.json").write_text(json.dumps(res, indent=1, default=float))
    write_run_manifest("upenn_rater_noise_109", OUT / "upenn_rater_noise.manifest.json", script="src/109_upenn_rater_noise.py",
                       seed=ps.DEFAULT_SEED, config={"regions": {k: list(v) for k, v in REGIONS.items()}, "grid_mm": GRID_MM,
                                                     "n_boot": ps.DEFAULT_N_BOOT},
                       inputs=[UP / "seg_series.csv"], dataset="UPENN-GBM", patient_split="none (descriptive)",
                       primary_endpoint="core region: median Dice, ASSD and |ln volume ratio|, AI vs radiologist-corrected")
    for reg in REGIONS:
        b = res["summary"][reg]["all"]
        print(reg, b["n_scans"], {k: round(v["median"], 3) for k, v in b.items() if isinstance(v, dict)})


if __name__ == "__main__":
    main()
