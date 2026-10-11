#!/usr/bin/env python3
"""Script 124: H-3 test INPUTS from UCSF-PDGM pre-op scans, aligned to UCSF-ALPTDG time1 space (Amendment 9).

No ALPTDG tumour label is read here (time1_seg, time2_seg, subtraction segs are never opened). Per grade-4 patient in the
ALPTDG "UCSF-PDGM Mapping" sheet with a local PDGM folder:
  1. rigid (6 dof, Mattes MI) registration of PDGM T1 -> ALPTDG time1 T1 (SimpleITK); pre-op tumour seg (BraTS 1/2/4 ->
     1/2/3), brain mask, T1c, FLAIR resampled (labels nearest, images linear);
  2. tissue: 3-class Gaussian mixture on T1 inside the brain minus tumour -> csf < gm < wm soft maps (declared proxy;
     PREDICT-GBM provides pbmaps from its own pipeline);
  3. arrays put in PREDICT-GBM order (PDGM/ALPTDG are LPS, PREDICT-GBM is RAS: flip x and y).
Exclusion (inputs only): registration brain-mask Dice to ALPTDG time1 brain (T1 > 0) < 0.80, or empty pre-op core.
Output (outside git): C:/Users/vihan/data/ucsf_h3/inputs/<lptdg_id>.npz ; manifest output/predictgbm/h3/inputs_manifest.json
"""
from __future__ import annotations

import gzip
import io
import json
import os
import sys
import zipfile
from pathlib import Path as _Path

import nibabel as nib
import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
MAN = PROJECT_ROOT / "output" / "predictgbm" / "h3"
ZIP = _Path(os.environ.get("UCSF_ALPTDG_ZIP", r"C:\Users\vihan\Downloads\UCSF_POSTOP_GLIOMA_DATASET_FINAL_v1.0.zip"))
PDGM = _Path(os.environ.get("UCSF_PDGM_V5", r"C:\Users\vihan\Downloads\PKG - UCSF-PDGM Version 5\UCSF-PDGM-v5"))
OUTD = _Path(os.environ.get("UCSF_H3_DIR", r"C:\Users\vihan\data\ucsf_h3")) / "inputs"
TABLE = "UCSF_PostopGlioma_Table S1 R1 V5.0_UNBLINDED_FINAL.xlsx"
FORBIDDEN = ("seg", "subtraction")   # ALPTDG files never opened before the H-3 freeze


def zread(z, name):
    assert not any(f in name for f in FORBIDDEN), name
    return nib.Nifti1Image.from_bytes(gzip.decompress(z.read(name)))


def cohort(z):
    x = pd.read_excel(io.BytesIO(z.read(TABLE)), sheet_name=["UCSF-PDGM Mapping", "Clinical Info"])
    m = x["UCSF-PDGM Mapping"].merge(x["Clinical Info"][["SubjectID", "Grade"]], left_on="UCSF-LPTDG ID", right_on="SubjectID")
    m = m[m["Grade"] == 4]
    m["pdgm_dir"] = [PDGM / f"UCSF-PDGM-{int(p.split('-')[-1]):04d}_nifti" for p in m["UCSF-PDGM ID"]]
    return m


def to_pg(a):
    return np.ascontiguousarray(a[::-1, ::-1])


def tissue(t1, brain, tumour):
    m = brain & ~tumour
    v = t1[m].astype(np.float64)
    mu = np.percentile(v, [15, 50, 85])           # 1-D 3-class EM (numpy; sklearn DLLs are blocked on this PC)
    sd, w = np.full(3, v.std() / 3), np.full(3, 1 / 3)
    for _ in range(100):
        lik = w * np.exp(-0.5 * ((v[:, None] - mu) / sd) ** 2) / sd + 1e-300
        r = lik / lik.sum(1, keepdims=True)
        n = r.sum(0)
        mu, w = (r * v[:, None]).sum(0) / n, n / len(v)
        sd = np.sqrt((r * (v[:, None] - mu) ** 2).sum(0) / n) + 1e-6
    order = np.argsort(mu)                        # csf, gm, wm by T1 brightness
    p = np.zeros((3,) + t1.shape, np.float32)
    p[:, m] = r.T[order]
    return p


def main():
    import SimpleITK as sitk
    OUTD.mkdir(parents=True, exist_ok=True)
    MAN.mkdir(parents=True, exist_ok=True)
    z = zipfile.ZipFile(ZIP)
    rows = []
    for _, r in cohort(z).iterrows():
        sid, d = int(r["UCSF-LPTDG ID"]), r["pdgm_dir"]
        rec = {"lptdg_id": sid, "pdgm_id": r["UCSF-PDGM ID"]}
        if not d.exists():
            rows.append({**rec, "excluded": "no PDGM folder"})
            continue
        pre = d.name.replace("_nifti", "")
        fix_nii = zread(z, f"{sid}/{sid}_time1_t1.nii.gz")
        fix = sitk.GetImageFromArray(np.asarray(fix_nii.dataobj, np.float32).T)
        mov = sitk.GetImageFromArray(np.asarray(nib.load(str(d / f"{pre}_T1.nii.gz")).dataobj, np.float32).T)
        reg = sitk.ImageRegistrationMethod()
        reg.SetMetricAsMattesMutualInformation(50)
        reg.SetMetricSamplingStrategy(reg.RANDOM)
        reg.SetMetricSamplingPercentage(0.05, 0)
        reg.SetInterpolator(sitk.sitkLinear)
        reg.SetOptimizerAsRegularStepGradientDescent(2.0, 1e-4, 300)
        reg.SetOptimizerScalesFromPhysicalShift()
        reg.SetShrinkFactorsPerLevel([4, 2, 1])
        reg.SetSmoothingSigmasPerLevel([2, 1, 0])
        init = sitk.CenteredTransformInitializer(fix, mov, sitk.Euler3DTransform(), sitk.CenteredTransformInitializerFilter.MOMENTS)
        reg.SetInitialTransform(init, inPlace=False)
        tx = reg.Execute(fix, mov)

        def warp(name, nearest):
            a = sitk.GetImageFromArray(np.asarray(nib.load(str(d / f"{pre}_{name}.nii.gz")).dataobj, np.float32).T)
            o = sitk.Resample(a, fix, tx, sitk.sitkNearestNeighbor if nearest else sitk.sitkLinear, 0.0)
            return sitk.GetArrayFromImage(o).T

        seg = np.rint(warp("tumor_segmentation", True)).astype(np.uint8)
        seg[seg == 4] = 3
        brain = warp("brain_segmentation", True) > 0
        fb = np.asarray(fix_nii.dataobj) > 0
        dice = float(2 * (brain & fb).sum() / max(brain.sum() + fb.sum(), 1))
        rec.update({"reg_brain_dice": round(dice, 4), "metric": round(reg.GetMetricValue(), 4)})
        if dice < 0.80:
            rows.append({**rec, "excluded": "registration brain Dice < 0.80"})
            continue
        if not np.isin(seg, (1, 3)).any():
            rows.append({**rec, "excluded": "empty pre-op core"})
            continue
        t1 = warp("T1", False)
        csf, gm, wm = tissue(t1, brain, seg > 0)
        out = {"seg": seg, "brain": brain, "wm": wm, "gm": gm, "csf": csf,
               "t1c": warp("T1c", False), "flair": warp("FLAIR", False)}
        np.savez_compressed(OUTD / f"{sid}.npz", **{k: to_pg(v) for k, v in out.items()})
        rows.append({**rec, "excluded": ""})
        print(sid, rec["reg_brain_dice"], flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(MAN / "inputs_cohort.csv", index=False)
    man = {"script": "124_pack_ucsf_h3_inputs", "n_grade4_mapped": int(len(df)),
           "n_included": int((df["excluded"] == "").sum()), "exclusions": df["excluded"].value_counts().to_dict(),
           "alptdg_labels_read": 0, "orientation": "LPS -> RAS (flip x, y)"}
    (MAN / "inputs_manifest.json").write_text(json.dumps(man, indent=1))
    print(json.dumps(man, indent=1))


if __name__ == "__main__":
    sys.exit(main())
