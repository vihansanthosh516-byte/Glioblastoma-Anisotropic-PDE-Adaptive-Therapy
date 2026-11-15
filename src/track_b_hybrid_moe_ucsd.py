"""Frozen hybrid finalist (moe) on UCSD-PTGBM, run once.

Lock: paper/trackB_hybrid_finalist_ucsd_2026-11-14.md. Rules from paper/trackB_v8_ucsd_external_protocol_2026-11-14.md.
"""
from __future__ import annotations

import glob
import json
import os
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage

import track_b_hybrid_harness as h
import track_b_hybrid_moe as moe
import track_b_mu_tcks as t
import track_b_mu_v3 as v3
import track_b_v8_ucsd_external as ext

OUT = h.ROOT / "output" / "track_b_hybrid_moe_ucsd_20261114"
V5C = ["log_edema", "nec_frac", "enh_frac", "log_whole", "compact"]


def extra_features(seg):
    base = seg.replace("_BraTS_tumor_seg.nii.gz", "_")
    a = np.asarray(nib.load(seg).dataobj)
    core, tumor = (a == 1) | (a == 3), np.isin(a, [1, 2, 3])
    v = int(core.sum())
    f = {"surf_ratio": 0.0, "sphericity": 0.0, "edema_reach": 0.0}
    if v:
        area = int((core & ~ndimage.binary_erosion(core)).sum())
        f = {"surf_ratio": area / v, "sphericity": float(np.pi ** (1 / 3) * (6 * v) ** (2 / 3) / max(area, 1)),
             "edema_reach": float(np.log((int((a == 2).sum()) + v) / v))}
    t1c = np.asarray(nib.load(base + "T1post.nii.gz").dataobj, dtype=np.float32)
    fl = np.asarray(nib.load(base + "FLAIR.nii.gz").dataobj, dtype=np.float32)
    norm = lambda im: im / max(float(np.median(im[im > 0])), 1e-6) if (im > 0).any() else im
    t1c, fl = norm(t1c), norm(fl)
    mean = lambda im, m: float(im[m].mean()) if m.any() else 0.0
    shell = ndimage.binary_dilation(tumor, iterations=5) & ~tumor & (fl > 0)
    f.update({"t1c_enh": mean(t1c, a == 3), "t1c_core_sd": float(t1c[core].std()) if v else 0.0,
              "flair_edema": mean(fl, a == 2), "flair_shell": mean(fl, shell)})
    return f


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "result.json").exists():
        raise SystemExit("locked run already done; refusing to overwrite")
    mu, comp_mu = h.load_dev()
    u = ext.ucsd_pairs()
    c = pd.read_excel(ext.UCLIN, sheet_name="Cinical info")
    c.columns = [str(x).strip() for x in c.columns]
    col = next(x for x in c.columns if "date of initial surgery, treatment or diagnosis" in x.lower())
    sess = pd.DataFrame({"session": c["ID"].astype(str).str.strip(), "day_in": -pd.to_numeric(c[col], errors="coerce")})
    sess["patient_id"] = sess.session.str.rsplit("_", n=1).str[0]
    segs = {}
    for p in sorted(glob.glob(str(ext.UCSD / "**" / "*_BraTS_tumor_seg.nii.gz"), recursive=True)):
        segs.setdefault(os.path.basename(p).split("_BraTS")[0], p)
    u = u.merge(sess, on=["patient_id", "day_in"], how="left")
    assert u.session.notna().all() and not u.duplicated(["patient_id", "day_in"]).any()
    extra = pd.DataFrame([extra_features(segs[s]) for s in u.session])
    comp_u = pd.concat([u[V5C].reset_index(drop=True), extra], axis=1)[list(comp_mu.columns)]
    te = u.drop(columns=["y_out", "core_vox_out", "r_obs"])
    pred = np.asarray(moe.fit_predict(mu.reset_index(drop=True), te.reset_index(drop=True), comp_mu, comp_u), float)
    base, _ = ext.fit_predict(mu, u, comp_mu[V5C], u[V5C])
    y = u.y_out.values
    errs = {k: np.abs(v - y) for k, v in base.items()}
    errs["moe"] = np.abs(pred - y)
    pat = pd.DataFrame(errs).assign(pid=u.patient_id.values).groupby("pid").mean()
    diff = (pat.moe - pat.pop_rate).values
    ci = v3.boot_ci(diff, np.random.default_rng(t.SEED), 0.95)
    res = {"script": "src/track_b_hybrid_moe_ucsd.py", "lock": "paper/trackB_hybrid_finalist_ucsd_2026-11-14.md",
           "label": "second exposure of UCSD; finalist failed dev rule",
           "n_pairs": int(len(u)), "n_patients": int(len(pat)), "patient_mae_log_volume": pat.mean().to_dict(),
           "moe_minus_pop_rate": {"mean": float(diff.mean()), "ci95": ci}, "success": bool(ci[1] < 0),
           "moe_minus_v8": float((pat.moe - pat.v8).mean())}
    (OUT / "result.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    pd.DataFrame({"patient_id": u.patient_id, "y_out": y, "moe": pred}).to_csv(OUT / "predictions.csv", index=False)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    run()
