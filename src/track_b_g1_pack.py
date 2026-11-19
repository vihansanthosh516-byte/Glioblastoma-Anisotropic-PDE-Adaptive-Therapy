"""Track B, lock chain step 3b: pack UCSD images (T1post, FLAIR) and the frozen segmentation weights for gate G1.
UCSD is a spent cohort and was NOT in the training set. Images are copied byte for byte. Released masks stay on this PC.
"""
import glob, json, os, shutil
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
UCSD = ROOT / "data" / "external" / "UCSD-PTGBM-v3"
RUN = ROOT / "kaggle_run" / "out_ptglioma_train_v4" / "nnUNet_results" / "Dataset501_PTGlioma" / "nnUNetTrainer_VIG250__nnUNetPlans__3d_fullres"
PACK = ROOT / "kaggle_run" / "pack_g1"
IM = PACK / "ucsd_images"; MD = PACK / "model"

def main():
    IM.mkdir(parents=True, exist_ok=True); MD.mkdir(parents=True, exist_ok=True)
    cases = []
    for m in sorted(glob.glob(str(UCSD / "**" / "*_BraTS_tumor_seg.nii.gz"), recursive=True)):
        t1, fl = m.replace("BraTS_tumor_seg", "T1post"), m.replace("BraTS_tumor_seg", "FLAIR")
        if not (os.path.exists(t1) and os.path.exists(fl)): continue
        s = os.path.basename(m).split("_BraTS")[0]
        if s in cases: continue
        shutil.copyfile(t1, IM / f"{s}_0000.nii.gz"); shutil.copyfile(fl, IM / f"{s}_0001.nii.gz"); cases.append(s)
    for f in ("dataset.json", "dataset_fingerprint.json", "plans.json"): shutil.copyfile(RUN / f, MD / f)
    shutil.copyfile(RUN / "fold_all" / "checkpoint_final.pth", MD / "checkpoint_final.pth")
    json.dump({"title": "trackb-ucsd-g1-private", "id": "vihansanthosh/trackb-ucsd-g1-private", "licenses": [{"name": "other"}]}, open(PACK / "dataset-metadata.json", "w"), indent=1)
    json.dump(cases, open(ROOT / "output" / "track_b_vig_cfb_waw_20261115" / "g1_cases.json", "w"), indent=1)
    print(len(cases), "cases;", len(list(IM.glob("*.nii.gz"))), "image files")

if __name__ == "__main__":
    main()
