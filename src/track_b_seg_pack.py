"""Track B, lock chain step 3: pack the training set of the two-sequence segmentation model in nnU-Net raw format.

Protocol: paper/trackB_vig_cfb_waw_protocol_2026-11-15.md, section 4. Spent cohorts only (MU, UCSF-ALPTDG).
Channels: 0 = T1Gd, 1 = FLAIR. Labels: 1 necrosis / non-enhancing core, 2 FLAIR abnormality, 3 enhancing, 4 resection cavity.
Images are copied byte for byte. Masks are rewritten as uint8 on the grid of the T1Gd image (same shape asserted).
"""
import glob, gzip, json, os, shutil, zipfile
from pathlib import Path
import numpy as np, nibabel as nib

ROOT = Path(__file__).resolve().parent.parent
MU = ROOT / "data" / "tcia" / "MU-Glioma-Post"
UCSF_ZIP = Path(os.environ.get("UCSF_ALPTDG_ZIP", r"C:\Users\vihan\Downloads\UCSF_POSTOP_GLIOMA_DATASET_FINAL_v1.0.zip"))
PACK = ROOT / "kaggle_run" / "pack_ptglioma"
DS = PACK / "nnUNet_raw" / "Dataset501_PTGlioma"


def put_label(arr, ref_img, dst, name):
    a = np.asarray(arr); u = set(np.unique(a).astype(int).tolist())
    assert u <= {0, 1, 2, 3, 4}, f"unexpected label {u} in {name}"
    assert a.shape == ref_img.shape, f"mask and image shape differ in {name}"
    nib.save(nib.Nifti1Image(a.astype(np.uint8), ref_img.affine), dst)


def main():
    (DS / "imagesTr").mkdir(parents=True, exist_ok=True); (DS / "labelsTr").mkdir(exist_ok=True)
    rows, moved = [], 0
    for m in sorted(glob.glob(str(MU / "*" / "Timepoint_*" / "*_tumorMask.nii.gz"))):
        t1c, fl = m.replace("tumorMask", "brain_t1c"), m.replace("tumorMask", "brain_t2f")
        if not (os.path.exists(t1c) and os.path.exists(fl)): continue
        p = Path(m).name.split("_"); case = f"MU_{p[1]}_{p[3]}"
        if not (DS / "labelsTr" / f"{case}.nii.gz").exists():
            ref, mk = nib.load(t1c), nib.load(m)
            moved += float(np.abs(ref.affine - mk.affine).max()) > 1e-3 or float(np.abs(ref.affine - nib.load(fl).affine).max()) > 1e-3
            shutil.copyfile(t1c, DS / "imagesTr" / f"{case}_0000.nii.gz"); shutil.copyfile(fl, DS / "imagesTr" / f"{case}_0001.nii.gz")
            put_label(mk.dataobj, ref, DS / "labelsTr" / f"{case}.nii.gz", case)
        rows.append({"case": case, "cohort": "MU"})
    with zipfile.ZipFile(UCSF_ZIP) as z:
        names = set(z.namelist())
        for s in sorted(x for x in names if x.endswith("_seg.nii.gz") and "subtraction" not in x):
            t1c, fl = s.replace("_seg", "_t1ce"), s.replace("_seg", "_flair")
            if not (t1c in names and fl in names): continue
            b = Path(s).name.split("_"); case = f"UCSF_{b[0]}_{b[1]}"
            if not (DS / "labelsTr" / f"{case}.nii.gz").exists():
                raw = z.read(t1c); ref = nib.Nifti1Image.from_bytes(gzip.decompress(raw)); mk = nib.Nifti1Image.from_bytes(gzip.decompress(z.read(s)))
                moved += float(np.abs(ref.affine - mk.affine).max()) > 1e-3
                (DS / "imagesTr" / f"{case}_0000.nii.gz").write_bytes(raw); (DS / "imagesTr" / f"{case}_0001.nii.gz").write_bytes(z.read(fl))
                put_label(mk.dataobj, ref, DS / "labelsTr" / f"{case}.nii.gz", case)
            rows.append({"case": case, "cohort": "UCSF"})
    json.dump({"channel_names": {"0": "T1Gd", "1": "FLAIR"}, "labels": {"background": 0, "netc": 1, "snfh": 2, "et": 3, "rc": 4},
               "numTraining": len(rows), "file_ending": ".nii.gz"}, open(DS / "dataset.json", "w"), indent=1)
    json.dump({"title": "trackb-ptglioma-nnunet-raw-private", "id": "vihansanthosh/trackb-ptglioma-nnunet-raw-private", "licenses": [{"name": "other"}]},
              open(PACK / "dataset-metadata.json", "w"), indent=1)
    n = {c: sum(r["cohort"] == c for r in rows) for c in ("MU", "UCSF")}
    size = sum(f.stat().st_size for f in DS.rglob("*") if f.is_file())
    out = {"cases": len(rows), "by_cohort": n, "images": len(list((DS / "imagesTr").glob("*.nii.gz"))), "labels": len(list((DS / "labelsTr").glob("*.nii.gz"))),
           "bytes": size, "cases_with_affine_difference_over_1e-3_this_run": int(moved)}
    json.dump({"summary": out, "cases": rows}, open(ROOT / "output" / "track_b_vig_cfb_waw_20261115" / "seg_train_manifest.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
