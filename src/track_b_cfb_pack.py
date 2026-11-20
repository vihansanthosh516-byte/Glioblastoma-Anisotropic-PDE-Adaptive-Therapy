"""Track B, lock chain step 5a: pack every other eligible CFB scan (T1Gd + FLAIR + released brain mask, any visit) for the frozen tool.
The 185 t0 scans already run for gate G2 (same method, same inputs) are not packed again; their predictions are reused.
Images only; no GTV voxel and no volume is read here. Same masking and cropping as src/track_b_g2_pack.py. Files get the ending .niigz.
"""
import json, os, shutil, zipfile
from pathlib import Path
import numpy as np, nibabel as nib, pandas as pd
ROOT = Path(__file__).resolve().parent.parent
CFB = Path(os.environ.get("CFB_DIR", r"C:\Users\vihan\Downloads\PKG - CFB-GBM version 4\CFB-GBM"))
MRI = ROOT / "research_notes" / "GBM patient analysis readiness 3" / "source_tables" / "CFB-GBM_mri_availability_v02_20260129.tsv"
RUN = ROOT / "kaggle_run" / "out_ptglioma_train_v4" / "nnUNet_results" / "Dataset501_PTGlioma" / "nnUNetTrainer_VIG250__nnUNetPlans__3d_fullres"
PACK = ROOT / "kaggle_run" / "pack_cfb"; OUT = ROOT / "output" / "track_b_vig_cfb_waw_20261115"

def main():
    m = pd.read_csv(MRI, sep="\t"); m = m[(m.t1gd == 1) & (m.flair == 1)]
    done = set(json.load(open(OUT / "g2_cases.json"))["cases"])
    shutil.rmtree(PACK, ignore_errors=True); (PACK / "files").mkdir(parents=True)
    cases, skipped = [], []
    for r in m.sort_values(["id_patient", "temporality"]).itertuples():
        n, t = int(r.id_patient), r.temporality; c = f"CFB{n:03d}_{t}"
        if c in done: continue
        d = CFB / f"{n:03d}" / t
        f = {k: d / f"{n}_{t}_{k}.nii.gz" for k in ("t1gd", "flair", "brain_mask")}
        if not all(p.exists() for p in f.values()): skipped.append(c); continue
        bm = nib.load(f["brain_mask"]); b = np.asarray(bm.dataobj) > 0
        if not b.any(): skipped.append(c); continue
        idx = np.argwhere(b); lo, hi = idx.min(0), idx.max(0) + 1; sl = tuple(slice(int(a), int(z)) for a, z in zip(lo, hi))
        shift = np.eye(4); shift[:3, 3] = lo
        for ch, k in (("0000", "t1gd"), ("0001", "flair")):
            im = nib.load(f[k]); assert im.shape == bm.shape, c
            nib.save(nib.Nifti1Image((np.asarray(im.dataobj, dtype=np.float32) * b)[sl], im.affine @ shift), PACK / "files" / f"{c}_{ch}.nii.gz")
            os.rename(PACK / "files" / f"{c}_{ch}.nii.gz", PACK / "files" / f"{c}_{ch}.niigz")
        cases.append(c)
    for fn in ("dataset.json", "dataset_fingerprint.json", "plans.json"): shutil.copyfile(RUN / fn, PACK / "files" / fn)
    shutil.copyfile(RUN / "fold_all" / "checkpoint_final.pth", PACK / "files" / "checkpoint_final.pth")
    with zipfile.ZipFile(PACK / "cfb_inputs.zip", "w", zipfile.ZIP_STORED) as z:
        for p in sorted((PACK / "files").iterdir()): z.write(p, p.name)
    json.dump({"title": "trackb-ptglioma-nnunet-raw-private", "id": "vihansanthosh/trackb-ptglioma-nnunet-raw-private", "licenses": [{"name": "other"}]}, open(PACK / "dataset-metadata.json", "w"))
    json.dump({"packed_cases": cases, "skipped": skipped, "reused_from_g2": sorted(done), "n_eligible_by_table": int(len(m))}, open(OUT / "cfb_label_cases.json", "w"), indent=1)
    print(len(cases), "packed;", len(done), "reused;", len(skipped), "skipped;", "table rows", len(m), "; zip MB", round(os.path.getsize(PACK / "cfb_inputs.zip") / 1e6))

if __name__ == "__main__":
    main()
