"""Track B gate G1 (protocol section 4): tool core vs released core on UCSD (spent cohort, not in training).
d = |log(tool core + 1) - log(released core + 1)|, core = labels 1 + 3, volume = voxels x voxel volume (mm3). Pass: median d <= 0.20 AND mean d <= 0.50.
Usage: py -3.14 src/track_b_g1_score.py <folder with pred/ and g1_status.json from the Kaggle kernel>
"""
import glob, json, os, sys
from pathlib import Path
import numpy as np, nibabel as nib
ROOT = Path(__file__).resolve().parent.parent
UCSD = ROOT / "data" / "external" / "UCSD-PTGBM-v3"; OUT = ROOT / "output" / "track_b_vig_cfb_waw_20261115"

def core_mm3(a, zooms):
    return float(((a == 1) | (a == 3)).sum()) * float(np.prod(zooms))

def main(folder):
    folder = Path(folder); st = json.load(open(folder / "g1_status.json")); cases = json.load(open(OUT / "g1_cases.json"))
    rel = {os.path.basename(p).split("_BraTS")[0]: p for p in glob.glob(str(UCSD / "**" / "*_BraTS_tumor_seg.nii.gz"), recursive=True)}
    rows = []
    for c in cases:
        pf = glob.glob(str(folder / "pred" / f"{c}.nii*"))
        if not pf: rows.append({"case": c, "missing": True}); continue
        p = nib.load(pf[0]); r = nib.load(rel[c])
        assert p.shape == r.shape, c
        t, g = core_mm3(np.asarray(p.dataobj), p.header.get_zooms()), core_mm3(np.asarray(r.dataobj), r.header.get_zooms())
        rows.append({"case": c, "tool_mm3": t, "released_mm3": g, "d": abs(np.log(t + 1) - np.log(g + 1))})
    ok = [x for x in rows if "d" in x]; d = np.array([x["d"] for x in ok])
    res = {"n_cases": len(cases), "n_scored": len(ok), "n_missing": len(rows) - len(ok), "median_d": float(np.median(d)), "mean_d": float(d.mean()),
           "p90_d": float(np.quantile(d, 0.9)), "pass": bool(len(ok) == len(cases) and np.median(d) <= 0.20 and d.mean() <= 0.50),
           "thresholds": {"median": 0.20, "mean": 0.50}, "checkpoint_sha256": st["checkpoint_sha256"], "nnunetv2": st["nnunetv2"], "torch": st["torch"]}
    OUT.mkdir(exist_ok=True); json.dump(res, open(OUT / "gate_g1.json", "w"), indent=1)
    import csv
    with open(OUT / "gate_g1_cases.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["case", "tool_mm3", "released_mm3", "d", "missing"]); w.writeheader(); w.writerows(rows)
    print(json.dumps(res, indent=1))

if __name__ == "__main__":
    main(sys.argv[1])
