"""Track B, lock chain step 4 (gate G2): pack CFB t0 scans that have a HUMAN GTV, as input of the frozen segmentation tool.
Reads: CFB t0 T1Gd, FLAIR, released brain mask (images; no GTV voxel is read here), availability tables (metadata).
Images are multiplied by the brain mask and cropped to the mask bounding box (affine shifted; voxel size unchanged). Crop box is saved for the score step.
Files get the ending .niigz so that Kaggle does not unpack them; the kernel links them back to .nii.gz.
"""
import json, os, shutil, zipfile
from pathlib import Path
import numpy as np, nibabel as nib, pandas as pd
ROOT = Path(__file__).resolve().parent.parent
CFB = Path(os.environ.get("CFB_DIR", r"C:\Users\vihan\Downloads\PKG - CFB-GBM version 4\CFB-GBM"))
TAB = ROOT / "research_notes" / "GBM patient analysis readiness 3" / "source_tables" / "CFB-GBM_treatment_imaging_availability_v03_20260619.tsv"
MRI = ROOT / "research_notes" / "GBM patient analysis readiness 3" / "source_tables" / "CFB-GBM_mri_availability_v02_20260129.tsv"
RUN = ROOT / "kaggle_run" / "out_ptglioma_train_v4" / "nnUNet_results" / "Dataset501_PTGlioma" / "nnUNetTrainer_VIG250__nnUNetPlans__3d_fullres"
PACK = ROOT / "kaggle_run" / "pack_g2"; OUT = ROOT / "output" / "track_b_vig_cfb_waw_20261115"

def main():
    g = pd.read_csv(TAB, sep="\t"); m = pd.read_csv(MRI, sep="\t")
    g = g[(g.temporality == "t0") & (g.gtv_type == "humanSegmentation") & (g.gtv == 1)]
    m = m[(m.temporality == "t0") & (m.t1gd == 1) & (m.flair == 1)]
    ids = sorted(set(g.id_patient) & set(m.id_patient))
    shutil.rmtree(PACK, ignore_errors=True); (PACK / "files").mkdir(parents=True)
    cases, boxes, skipped = [], {}, []
    for pid in ids:
        n = int(pid); d = CFB / f"{n:03d}" / "t0"
        f = {k: d / f"{n}_t0_{k}.nii.gz" for k in ("t1gd", "flair", "brain_mask")}
        if not all(p.exists() for p in f.values()): skipped.append(int(pid)); continue
        bm = nib.load(f["brain_mask"]); b = np.asarray(bm.dataobj) > 0
        if not b.any(): skipped.append(int(pid)); continue
        idx = np.argwhere(b); lo, hi = idx.min(0), idx.max(0) + 1; sl = tuple(slice(int(a), int(c)) for a, c in zip(lo, hi))
        c = f"CFB{n:03d}_t0"
        for ch, k in (("0000", "t1gd"), ("0001", "flair")):
            im = nib.load(f[k]); assert im.shape == bm.shape
            a = (np.asarray(im.dataobj, dtype=np.float32) * b)[sl]
            shift = np.eye(4); shift[:3, 3] = lo                      # new voxel (0,0,0) = old voxel lo
            nib.save(nib.Nifti1Image(a, im.affine @ shift), PACK / "files" / f"{c}_{ch}.nii.gz")
            os.rename(PACK / "files" / f"{c}_{ch}.nii.gz", PACK / "files" / f"{c}_{ch}.niigz")
        cases.append(c); boxes[c] = [[int(x) for x in lo], [int(x) for x in hi]]
    for fn in ("dataset.json", "dataset_fingerprint.json", "plans.json"): shutil.copyfile(RUN / fn, PACK / "files" / fn)
    shutil.copyfile(RUN / "fold_all" / "checkpoint_final.pth", PACK / "files" / "checkpoint_final.pth")
    with zipfile.ZipFile(PACK / "g2_inputs.zip", "w", zipfile.ZIP_STORED) as z:
        for p in sorted((PACK / "files").iterdir()): z.write(p, p.name)
    json.dump({"title": "trackb-ptglioma-nnunet-raw-private", "id": "vihansanthosh/trackb-ptglioma-nnunet-raw-private", "licenses": [{"name": "other"}]}, open(PACK / "dataset-metadata.json", "w"))
    json.dump({"cases": cases, "boxes": boxes, "skipped_ids": skipped, "n_human_t0_gtv_in_table": int(len(g)), "n_with_t1gd_flair": len(ids)}, open(OUT / "g2_cases.json", "w"), indent=1)
    print(len(cases), "cases; skipped", skipped, "; zip MB", round(os.path.getsize(PACK / "g2_inputs.zip") / 1e6))

if __name__ == "__main__":
    main()
