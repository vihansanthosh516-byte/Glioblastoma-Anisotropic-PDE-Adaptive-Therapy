"""Track B image-feature build. Pre-registration: paper/trackB_imgfeat_prereg_2026-11-19.md. Writes pairs_with_features.csv (six cohorts, 1,271 pairs)
with log_edema, nec_frac, enh_frac, log_cavity at the INPUT scan. Volumes in mm3. Reads masks and label tables only."""
import glob, gzip, io, os, re, sys, warnings, zipfile
from pathlib import Path
import numpy as np, pandas as pd, nibabel as nib

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent
AUDIT = ROOT / "output" / "track_b_code_audit_20261114"; TEST = ROOT / "output" / "track_b_vig_cfb_waw_20261115"; OUT = ROOT / "output" / "track_b_imgfeat_dev_20261119"
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(AUDIT))
import track_b_vig_lumiere as vig, track_b_v8_ucsd_external as ext, track_b_v8_ucsf_alptdg as uf


def counts(a, vox=1.0): return [float((a == k).sum()) * vox for k in (1, 2, 3, 4)]


def mu():
    c = pd.read_csv(AUDIT / "mu_mask_counts.csv"); c = c[c.has_mask]
    return pd.DataFrame({"patient_id": "MU_" + c.patient_id, "n1": c.n1, "n2": c.n2, "n3": c.n3, "n4": c.n4})


def ucsd():
    rows = []
    for p in sorted(glob.glob(str(ext.UCSD / "**" / "*_BraTS_tumor_seg.nii.gz"), recursive=True)):
        s = os.path.basename(p).split("_BraTS")[0]; im = nib.load(p); n = counts(np.asarray(im.dataobj), float(np.prod(im.header.get_zooms()[:3])))
        rows.append({"patient_id": "UCSD_" + s.rsplit("_", 1)[0], "n1": n[0], "n2": n[1], "n3": n[2], "n4": n[3]})
    return pd.DataFrame(rows)


def ucsf():
    rows = []
    with zipfile.ZipFile(uf.ZIP) as z:
        for n in z.namelist():
            m = re.match(r"^(\d+)/\1_time1_seg\.nii\.gz$", n)
            if m:
                im = nib.Nifti1Image.from_bytes(gzip.decompress(z.read(n))); k = counts(np.asarray(im.dataobj), float(np.prod(im.header.get_zooms()[:3])))
                rows.append({"patient_id": "UCSF_" + m.group(1), "n1": k[0], "n2": k[1], "n3": k[2], "n4": k[3]})
    return pd.DataFrame(rows)


def lumiere():
    s = vig.lumiere_scans()
    return pd.DataFrame({"patient_id": "LUMIERE_" + s["Patient"], "n1": s["Necrosis"], "n2": s["Edema"], "n3": s["Contrast-enhancing"], "n4": np.nan})


def cfb():
    rows = []
    for f in sorted(glob.glob(str(ROOT / "kaggle_run" / "out_g2" / "pred" / "CFB*_t?.nii*")) + glob.glob(str(ROOT / "kaggle_run" / "out_cfb" / "pred" / "CFB*_t?.nii*"))):
        m = re.match(r"CFB(\d+)_t\d\.nii", os.path.basename(f)); im = nib.load(f); k = counts(np.asarray(im.dataobj), float(np.prod(im.header.get_zooms()[:3])))
        rows.append({"patient_id": f"CFB_{int(m.group(1))}", "n1": k[0], "n2": k[1], "n3": k[2], "n4": k[3]})
    return pd.DataFrame(rows)


def waw():
    rows = []
    with zipfile.ZipFile(ROOT / "data" / "external" / "waw" / "3_2_postop_masks.zip") as z:
        for n in sorted(x for x in z.namelist() if x.endswith(".nii.gz")):
            m = re.match(r"^(\d+)_(\d+)_", Path(n).name); im = nib.Nifti1Image.from_bytes(gzip.decompress(z.read(n)))
            k = counts(np.asarray(im.dataobj), float(np.prod(im.header.get_zooms()[:3]))); rows.append({"patient_id": "WAW_" + m.group(1), "n1": k[0], "n2": k[1], "n3": k[2], "n4": k[3]})
    return pd.DataFrame(rows)


def main():
    a = pd.read_csv(AUDIT / "pooled4_pairs.csv"); b = pd.read_csv(TEST / "test_pairs.csv"); d = pd.concat([a, b[a.columns]], ignore_index=True)
    sc = pd.concat([f() for f in (mu, ucsd, ucsf, lumiere, cfb, waw)], ignore_index=True)
    sc["core"] = sc.n1 + sc.n3; sc["k"] = sc.core.round(3)
    sc["log_edema"] = np.log(sc.n2 + 1.0); sc["nec_frac"] = sc.n1 / sc.core.where(sc.core > 0); sc["enh_frac"] = sc.n3 / sc.core.where(sc.core > 0); sc["log_cavity"] = np.log(sc.n4 + 1.0)
    dup = sc.duplicated(["patient_id", "k"], keep=False); print("scans", len(sc), "with duplicate (patient, core) keys", int(dup.sum()))
    sc = sc.drop_duplicates(["patient_id", "k"], keep=False)          # ambiguous keys get NaN features
    d["k"] = d.core_vox_in.round(3); m = d.merge(sc[["patient_id", "k", "log_edema", "nec_frac", "enh_frac", "log_cavity"]], on=["patient_id", "k"], how="left").drop(columns="k")
    assert len(m) == len(d)
    for c in d.cohort.unique():
        g = m[m.cohort == c]; print(c, "pairs", len(g), "matched", int(g.log_edema.notna().sum()), "unmatched", int(g.log_edema.isna().sum()), "| cavity present", int(g.log_cavity.notna().sum()))
    m.to_csv(OUT / "pairs_with_features.csv", index=False)


if __name__ == "__main__":
    main()
