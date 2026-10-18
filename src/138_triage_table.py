#!/usr/bin/env python3
"""Script 138: triage table (Amendment 14, A14.2-A14.4). Features and outcome per patient. No modelling here.

1 mm grid. core = labels 1|3, edema = label 2, all recurrence = labels 1-3 inside the pre-op brain.
standard plan = predictgbm_io.standard_plan (EDT <= 15 mm from the core, inside the brain).
Output: output/triage/triage_table.csv, output/triage/triage_exclusions.json
"""
from __future__ import annotations

import gzip
import json
import os
import sys
import zipfile
from pathlib import Path as _Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "triage"
H3 = _Path(os.environ.get("UCSF_H3_DIR", r"C:\Users\vihan\data\ucsf_h3"))
ZIP = _Path(os.environ.get("UCSF_ALPTDG_ZIP", r"C:\Users\vihan\Downloads\UCSF_POSTOP_GLIOMA_DATASET_FINAL_v1.0.zip"))
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import predictgbm_io as pg  # noqa: E402

MARGIN = pg.CTV_MARGIN


def ucsf_label(z, sid):
    a = np.asarray(nib.Nifti1Image.from_bytes(gzip.decompress(z.read(f"{sid}/{sid}_time1_seg.nii.gz"))).dataobj)
    return np.ascontiguousarray(np.rint(a).astype(np.uint8)[::-1, ::-1])


def row(pid, centre, seg, brain, rec):
    core = pg.core(seg)
    if not core.any():
        return None, "no core voxel"
    ra = np.isin(rec, (1, 2, 3)) & brain
    if not ra.any():
        return None, "no all-recurrence voxel inside the pre-op brain"
    std = pg.standard_plan(seg, brain)
    dil = (ndimage.distance_transform_edt(~core) <= MARGIN)
    idx = np.argwhere(core)
    cen = idx.mean(0)
    bb = np.argwhere(brain)
    lo, hi = bb.min(0), bb.max(0)
    dsurf = ndimage.distance_transform_edt(brain)[tuple(np.clip(np.rint(cen).astype(int), 0, np.array(brain.shape) - 1))]
    lab, n = ndimage.label(core)
    sizes = ndimage.sum(core, lab, range(1, n + 1)) if n else []
    cov = float((ra & std).sum() / ra.sum())
    return {
        "pid": pid, "centre": centre,
        "f1_log_core_ml": float(np.log(core.sum() / 1000.0)),
        "f2_log1p_edema_ml": float(np.log1p((seg == 2).sum() / 1000.0)),
        "f3_enh_share": float((seg == 3).sum() / core.sum()),
        "f4_lateral_mm": float(abs(cen[0] - (lo[0] + hi[0]) / 2.0)),
        "f5_ap_pos": float((cen[1] - lo[1]) / max(hi[1] - lo[1], 1)),
        "f6_si_pos": float((cen[2] - lo[2]) / max(hi[2] - lo[2], 1)),
        "f7_dist_surface_mm": float(dsurf),
        "f8_margin_truncation": float(1.0 - std.sum() / max(dil.sum(), 1)),
        "f9_max_extent_mm": float((idx.max(0) - idx.min(0) + 1).max()),
        "f10_n_components": int((np.array(sizes) >= 500).sum()),
        "n_rec_vox": int(ra.sum()), "standard_coverage": cov,
        "fail": int(cov < 0.80), "fail90": int(cov < 0.90)}, None


def main():
    rows, excl = [], []

    def add(pid, centre, seg, brain, rec):
        r, why = row(pid, centre, seg, brain, rec)
        if r is None:
            excl.append({"pid": pid, "centre": centre, "reason": why})
        else:
            rows.append(r)

    for pid in pg.dev_ids():
        c = pg.load_case(pid, with_recurrence=True)
        add(pid, "TUM", c["seg"], c["brain"], c["rec"])
    for pid in sorted(pg.test_ids()):
        c = pg.load_case(pid, with_recurrence=True)
        add(pid, "RHUH" if pid.startswith("RHUH") else "LUMIERE", c["seg"], c["brain"], c["rec"])
        print(pid, flush=True) if len(rows) % 20 == 0 else None
    z = zipfile.ZipFile(ZIP)
    for f in sorted((H3 / "inputs").glob("*.npz")):
        x = np.load(f)
        add(f"ucsf_{f.stem}", "UCSF", np.rint(x["seg"]).astype(np.uint8), x["brain"] > 0, ucsf_label(z, f.stem))
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "triage_table.csv", index=False)
    (OUT / "triage_exclusions.json").write_text(json.dumps(excl, indent=1))
    print(df.groupby("centre").agg(n=("pid", "count"), fail=("fail", "mean"), fail90=("fail90", "mean"),
                                   cov=("standard_coverage", "median")).round(3))
    print("total", len(df), "excluded", len(excl))


if __name__ == "__main__":
    main()
