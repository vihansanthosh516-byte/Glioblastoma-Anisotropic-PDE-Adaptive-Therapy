"""Track B, lock chain step 5b: tool core volumes of CFB scans -> output/.../cfb_tool_volumes.csv (id_patient, temporality, core_mm3).
core = labels 1 + 3; volume = voxels x voxel volume (mm3). Reads prediction masks only (no image, no GTV). Prints counts only, never volumes:
the file holds the targets of the sealed test. Inputs: folder of the G2 run (reused) and folder of the CFB run, each with pred/ and a status json.
Usage: py -3.14 src/track_b_cfb_volumes.py <g2 folder> <cfb folder>
"""
import glob, json, re, sys
from pathlib import Path
import numpy as np, nibabel as nib, pandas as pd
ROOT = Path(__file__).resolve().parent.parent; OUT = ROOT / "output" / "track_b_vig_cfb_waw_20261115"
SHA = "bf4a8825b1555117504e9ecd61592fabce31674be05589c3cda3a65115addb1f"

def main(g2, cfb):
    rows = []
    for folder, stf in ((g2, "g2_status.json"), (cfb, "cfb_status.json")):
        st = json.load(open(Path(folder) / stf)); assert st["checkpoint_sha256"] == SHA and st["returncode"] == 0 and st["n_predictions"] == st["n_inputs"]
        for p in sorted(glob.glob(str(Path(folder) / "pred" / "CFB*_t?.nii*"))):
            mt = re.match(r"CFB(\d+)_(t\d)\.nii", Path(p).name); im = nib.load(p); a = np.asarray(im.dataobj)
            rows.append({"id_patient": int(mt.group(1)), "temporality": mt.group(2), "core_mm3": float(((a == 1) | (a == 3)).sum()) * float(np.prod(im.header.get_zooms()[:3]))})
    d = pd.DataFrame(rows); assert not d.duplicated(["id_patient", "temporality"]).any()
    d.to_csv(OUT / "cfb_tool_volumes.csv", index=False)
    print("scans", len(d), "patients", d.id_patient.nunique(), "by visit", d.temporality.value_counts().to_dict(), "scans with core 0:", int((d.core_mm3 == 0).sum()))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
