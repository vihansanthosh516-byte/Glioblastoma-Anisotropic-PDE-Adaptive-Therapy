"""Read-only diagnostics for the Track B review. Post hoc, descriptive. Writes nothing."""
import glob, gzip, io, os, sys, zipfile
from pathlib import Path
import numpy as np, pandas as pd, nibabel as nib

ROOT = Path(sys.argv[1])
p = pd.read_csv(ROOT / "data/manifests/forecast_pairs_mu.csv")
s = pd.read_csv(ROOT / "data/manifests/split_mu.csv")
print("PAIRS cols", list(p.columns)); print("SPLIT cols", list(s.columns), "rows", len(s), "unique pid", s.patient_id.nunique())
for c in s.columns:
    if s[c].nunique() < 12: print(" split", c, s[c].value_counts(dropna=False).to_dict())
m = p.merge(s, on="patient_id", how="left")
prim = m[m.in_primary & m.input_core_nonempty]
print("primary pairs", len(prim), "patients", prim.patient_id.nunique(), "fold NaN", int(prim.fold.isna().sum()) if "fold" in prim else "nofold")
print("folds per patient max", prim.groupby("patient_id").fold.nunique().max())
print("pairs/patient by fold", prim.groupby("fold").agg(n=("patient_id", "size"), pat=("patient_id", "nunique")).to_dict())
print("dt primary", prim.dt_days.describe().round(1).to_dict())
print("dt all", p.dt_days.describe().round(1).to_dict(), "n all", len(p))
print("in_primary false reasons: n", int((~m.in_primary).sum()), "; input empty", int((~m.input_core_nonempty).sum()))
yi, yo = np.log(prim.core_vox_in + 1), np.log(prim.core_vox_out + 1)
e = (yo - yi).abs()
print("MU persistence pair-level: mean", round(e.mean(), 3), "median", round(e.median(), 3))
z = prim.core_vox_out == 0
print("MU y_out==0 pairs", int(z.sum()), "share of pair abs err", round(e[z].sum() / e.sum(), 3))
print("MU core_vox_in quantiles", prim.core_vox_in.quantile([0, .05, .1, .25, .5, .75, 1]).round(0).to_dict())
small = prim.core_vox_in < 500
print("MU input<500 vox: n", int(small.sum()), "mean err", round(e[small].mean(), 3), "vs rest", round(e[~small].mean(), 3), "share err", round(e[small].sum() / e.sum(), 3))
pp = prim.groupby("patient_id").size()
print("pairs per patient", pp.describe().round(2).to_dict())

def hdr(img, tag):
    a = np.asarray(img.dataobj)
    print(tag, "shape", img.shape, "zooms", tuple(round(float(x), 3) for x in img.header.get_zooms()[:3]), "labels", {int(k): int(v) for k, v in zip(*np.unique(a, return_counts=True))})

mu = sorted(glob.glob(str(ROOT / "data/tcia/MU-Glioma-Post/*/Timepoint_*/*_tumorMask.nii.gz")))
print("MU masks", len(mu))
for f in mu[:1] + mu[len(mu)//2:len(mu)//2+1] + mu[-1:]: hdr(nib.load(f), "MU " + os.path.basename(f))
zs = set()
for f in mu[::25]: zs.add(tuple(round(float(x), 2) for x in nib.load(f).header.get_zooms()[:3]) + nib.load(f).shape)
print("MU distinct zoom+shape in sample", zs)
uc = sorted(glob.glob(str(ROOT / "data/external/UCSD-PTGBM-v3/**/*_BraTS_tumor_seg.nii.gz"), recursive=True))
print("UCSD segs", len(uc))
for f in uc[:1] + uc[-1:]: hdr(nib.load(f), "UCSD " + os.path.basename(f))
allseg = sorted(glob.glob(str(ROOT / "data/external/UCSD-PTGBM-v3/**/*seg*.nii.gz"), recursive=True))
print("UCSD seg name kinds", sorted({os.path.basename(f).split("_", 3)[-1] if f.count("_") else f for f in allseg})[:12])
zs = set()
for f in uc[::10]: zs.add(tuple(round(float(x), 2) for x in nib.load(f).header.get_zooms()[:3]) + nib.load(f).shape)
print("UCSD distinct zoom+shape in sample", zs)

for name in ["track_b_v8_ucsd_external_20261114", "track_b_v8_ucsf_alptdg_20261114"]:
    d = pd.read_csv(ROOT / "output" / name / "predictions.csv")
    if "used_before" in d: d = d[~d.used_before]
    e = (d.persistence - d.y_out).abs(); z = d.y_out == 0
    print(name, "pairs", len(d), "pers pair-mean", round(e.mean(), 3), "median", round(e.median(), 3), "y_out==0", int(z.sum()), "share err", round(e[z].sum() / e.sum(), 3),
          "v8 err on zero-out", round((d.v8 - d.y_out).abs()[z].mean(), 3) if z.any() else None,
          "v8-pers mean signed", round((d.v8 - d.persistence).mean(), 3), "y_in median", round(d.persistence.median(), 2),
          "true change median", round((d.y_out - d.persistence).median(), 3), "phase", d.phase.value_counts().to_dict())
d = pd.read_csv(ROOT / "output/track_b_mu_v8_20261114/predictions.csv")
print("MU v8 preds cols", list(d.columns), "NaN v8", int(d.v8.isna().sum()))

Z = Path(os.environ.get("UCSF_ALPTDG_ZIP", r"C:\Users\vihan\Downloads\UCSF_POSTOP_GLIOMA_DATASET_FINAL_v1.0.zip"))
print("UCSF zip exists", Z.exists())
if Z.exists():
    z = zipfile.ZipFile(Z)
    names = z.namelist()
    print("zip n", len(names), "sample", [n for n in names if "/100" in n or n.count("/") == 1][:14])
    x = pd.read_excel(io.BytesIO(z.read("UCSF_PostopGlioma_Table S1 R1 V5.0_UNBLINDED_FINAL.xlsx")), sheet_name=None)
    for k, v in x.items(): print("SHEET", k, v.shape, [str(c).strip() for c in v.columns][:40])
    c = x["Clinical Info"]; sid = int(c["SubjectID"].iloc[0])
    for tp in (1, 2):
        img = nib.Nifti1Image.from_bytes(gzip.decompress(z.read(f"{sid}/{sid}_time{tp}_seg.nii.gz")))
        hdr(img, f"UCSF {sid} t{tp}")
cu = pd.read_excel(ROOT / "research_notes/trackB_confirmatory_audit_ucsd/UCSD_PTGBM-clinical-information_v3_2026-12-Mar.xlsx", sheet_name="Cinical info")
print("UCSD clin cols", [str(c).strip() for c in cu.columns])
print(cu.iloc[:3, :12].to_string())
