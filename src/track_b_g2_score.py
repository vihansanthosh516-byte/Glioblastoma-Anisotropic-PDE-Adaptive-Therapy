"""Track B gate G2 (protocol section 4): share of tool-core volume that lies inside the human t0 GTV dilated by 5 mm, on CFB t0 scans.
Cases with tool core = 0 are left out (protocol). Pass: median share >= 0.70. The t0 GTV is used here only; it is never a target.
Usage: py -3.14 src/track_b_g2_score.py <folder with pred/ and g2_status.json from the Kaggle kernel>
"""
import csv, glob, json, os, sys
from pathlib import Path
import numpy as np, nibabel as nib
from scipy.ndimage import distance_transform_edt
ROOT = Path(__file__).resolve().parent.parent; OUT = ROOT / "output" / "track_b_vig_cfb_waw_20261115"
CFB = Path(os.environ.get("CFB_DIR", r"C:\Users\vihan\Downloads\PKG - CFB-GBM version 4\CFB-GBM"))

def main(folder):
    folder = Path(folder); st = json.load(open(folder / "g2_status.json")); info = json.load(open(OUT / "g2_cases.json")); rows = []
    for c in info["cases"]:
        n = int(c[3:6]); pf = glob.glob(str(folder / "pred" / f"{c}.nii*"))
        if not pf: rows.append({"case": c, "status": "missing"}); continue
        gt = nib.load(CFB / f"{n:03d}" / "t0" / f"{n}_t0_gtv.nii.gz"); g = np.asarray(gt.dataobj) > 0
        lo, hi = info["boxes"][c]; sl = tuple(slice(a, b) for a, b in zip(lo, hi))
        p = np.asarray(nib.load(pf[0]).dataobj); assert p.shape == tuple(b - a for a, b in zip(lo, hi)), c
        core = np.zeros(g.shape, bool); core[sl] = (p == 1) | (p == 3)
        if core.sum() == 0: rows.append({"case": c, "status": "tool_core_zero", "gtv_empty": int(not g.any())}); continue
        if g.any():
            near = distance_transform_edt(~g, sampling=gt.header.get_zooms()[:3]) <= 5.0
            share = float((core & near).sum() / core.sum())
        else: share = 0.0
        rows.append({"case": c, "status": "scored", "share": share, "tool_core_vox": int(core.sum()), "gtv_vox": int(g.sum()), "gtv_empty": int(not g.any())})
    sc = [r["share"] for r in rows if r["status"] == "scored"]
    res = {"n_cases": len(info["cases"]), "n_scored": len(sc), "n_tool_core_zero": sum(r["status"] == "tool_core_zero" for r in rows), "n_missing": sum(r["status"] == "missing" for r in rows),
           "n_gtv_empty_among_scored": sum(r.get("gtv_empty", 0) for r in rows if r["status"] == "scored"),
           "median_share": float(np.median(sc)), "mean_share": float(np.mean(sc)), "p25_share": float(np.quantile(sc, .25)),
           "pass": bool(len(sc) > 0 and np.median(sc) >= 0.70 and not any(r["status"] == "missing" for r in rows)), "threshold_median": 0.70,
           "checkpoint_sha256": st["checkpoint_sha256"], "nnunetv2": st["nnunetv2"], "torch": st["torch"]}
    json.dump(res, open(OUT / "gate_g2.json", "w"), indent=1)
    with open(OUT / "gate_g2_cases.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["case", "status", "share", "tool_core_vox", "gtv_vox", "gtv_empty"]); w.writeheader(); w.writerows(rows)
    print(json.dumps(res, indent=1))

if __name__ == "__main__":
    main(sys.argv[1])
