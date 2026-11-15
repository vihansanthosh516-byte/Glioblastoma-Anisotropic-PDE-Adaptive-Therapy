"""Code and data audit 1: MU-Glioma-Post masks, labels, days and the pair table. Read-only.

Rebuilds label counts and scan days from the raw files and compares them with
data/manifests/forecast_pairs_mu.csv. Then lists the worst pairs with all four labels.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, nibabel as nib

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import track_b_mu_tcks as t

HERE = Path(__file__).resolve().parent
MU = ROOT / "data" / "tcia" / "MU-Glioma-Post"
CACHE = HERE / "mu_mask_counts.csv"
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)


def counts():
    if CACHE.exists():
        return pd.read_csv(CACHE)
    rows = []
    for pdir in sorted(MU.iterdir()):
        for tdir in sorted(pdir.glob("Timepoint_*")):
            tp = int(tdir.name.split("_")[1])
            files = sorted(f.name.replace(f"{pdir.name}_Timepoint_{tp}_", "") for f in tdir.iterdir())
            m = tdir / f"{pdir.name}_Timepoint_{tp}_tumorMask.nii.gz"
            r = {"patient_id": pdir.name, "tp": tp, "files": "|".join(files), "has_mask": m.exists()}
            if m.exists():
                img = nib.load(m); a = np.asarray(img.dataobj)
                u, c = np.unique(a, return_counts=True)
                d = {int(k): int(v) for k, v in zip(u, c)}
                r.update({"n1": d.get(1, 0), "n2": d.get(2, 0), "n3": d.get(3, 0), "n4": d.get(4, 0), "other_labels": "|".join(str(k) for k in d if k not in (0, 1, 2, 3, 4)),
                          "zoom": str(tuple(round(float(z), 3) for z in img.header.get_zooms()[:3])), "shape": str(img.shape)})
            rows.append(r)
    c = pd.DataFrame(rows); c.to_csv(CACHE, index=False)
    return c


def main():
    c = counts()
    print("A. masks:", len(c), "patients", c.patient_id.nunique(), "| missing mask", int((~c.has_mask).sum()), "| zoom", c.zoom.value_counts().to_dict(), "| shape", c["shape"].value_counts().to_dict())
    print("   file sets:", c.files.value_counts().to_dict())
    print("   unexpected labels:", c.other_labels.fillna("").replace("", np.nan).dropna().value_counts().to_dict())
    c["core"] = c.n1 + c.n3
    print("   label present in share of masks:", {k: round(float((c[k] > 0).mean()), 3) for k in ["n1", "n2", "n3", "n4"]})
    print("   median voxels when present:", {k: float(c.loc[c[k] > 0, k].median()) for k in ["n1", "n2", "n3", "n4"]})
    print("   masks with core == 0:", int((c.core == 0).sum()), "| of those with label 4 > 0:", int(((c.core == 0) & (c.n4 > 0)).sum()), "| with label 2 > 0:", int(((c.core == 0) & (c.n2 > 0)).sum()), "| fully empty:", int(((c.n1 + c.n2 + c.n3 + c.n4) == 0).sum()))
    print("   label 3 > 0 but label 1 == 0:", int(((c.n3 > 0) & (c.n1 == 0)).sum()), "| label 1 > 0 but label 3 == 0:", int(((c.n1 > 0) & (c.n3 == 0)).sum()), "| label 4 > 0 and label 3 == 0:", int(((c.n4 > 0) & (c.n3 == 0)).sum()))

    cl = pd.read_excel(t.CLIN, sheet_name="MU Glioma Post"); cl.columns = [str(k).strip() for k in cl.columns]
    cl["patient_id"] = cl["Patient_ID"].astype(str).str.strip()
    dcol = {k: next(x for x in cl.columns if f"(Timepoint_{k})" in x) for k in range(1, 7)}
    days = cl.melt(id_vars="patient_id", value_vars=list(dcol.values()), var_name="col", value_name="day")
    days["tp"] = days.col.map({v: k for k, v in dcol.items()}); days["day"] = pd.to_numeric(days.day, errors="coerce")
    c = c.merge(days[["patient_id", "tp", "day"]], on=["patient_id", "tp"], how="left")
    print("\nB. scan days: masks with no day in the clinical table:", int(c.day.isna().sum()), "| days in the table with no mask folder:", int((days.day.notna()).sum() - c.day.notna().sum()))
    bad = []
    for pid, g in c.sort_values(["patient_id", "tp"]).groupby("patient_id"):
        d = g.day.values
        if np.any(np.diff(d[~np.isnan(d)]) <= 0): bad.append((pid, g.tp.tolist(), d.tolist()))
    print("   patients whose days do not increase with timepoint number:", len(bad), bad[:8])

    p = pd.read_csv(t.PAIRS)
    k = p.merge(c.add_suffix("_a").rename(columns={"patient_id_a": "patient_id", "tp_a": "tp_in"}), on=["patient_id", "tp_in"], how="left")
    k = k.merge(c.add_suffix("_b").rename(columns={"patient_id_b": "patient_id", "tp_b": "tp_out"}), on=["patient_id", "tp_out"], how="left")
    print("\nC. pair table vs raw files: pairs", len(k), "| core_in mismatch", int((k.core_vox_in != k.core_a).sum()), "| core_out mismatch", int((k.core_vox_out != k.core_b).sum()),
          "| day_in mismatch", int((np.abs(k.day_in - k.day_a) > 0.5).sum()), "| day_out mismatch", int((np.abs(k.day_out - k.day_b) > 0.5).sum()), "| dt mismatch", int((np.abs(k.dt_days - (k.day_b - k.day_a)) > 0.5).sum()))
    nxt = {}
    for pid, g in c.dropna(subset=["day"]).sort_values(["patient_id", "day"]).groupby("patient_id"):
        tps = g.tp.tolist(); nxt.update({(pid, a): b for a, b in zip(tps[:-1], tps[1:])})
    k["next_tp"] = [nxt.get((a, b), -1) for a, b in zip(k.patient_id, k.tp_in)]
    print("   pairs whose tp_out is not the next scan in time:", int((k.next_tp != k.tp_out).sum()), "| possible consecutive pairs from raw files:", len(nxt), "| pairs in table:", len(k))
    print("   in_primary false:", int((~k.in_primary).sum()), "| reasons: dt<14", int((k.dt_days < 14).sum()), "dt>365", int((k.dt_days > 365).sum()), "empty input", int((k.core_vox_in == 0).sum()))

    prim = k[k.in_primary & k.input_core_nonempty].copy()
    prim["ch"] = np.log(prim.core_vox_out + 1) - np.log(prim.core_vox_in + 1)
    prim["whole_a"], prim["whole_b"] = prim.n1_a + prim.n2_a + prim.n3_a, prim.n1_b + prim.n2_b + prim.n3_b
    prim["d_whole"] = np.log(prim.whole_b + 1) - np.log(prim.whole_a + 1)
    prim["d_cav"] = prim.n4_b - prim.n4_a
    cols = ["patient_id", "tp_in", "tp_out", "dt_days", "n1_a", "n3_a", "n2_a", "n4_a", "n1_b", "n3_b", "n2_b", "n4_b", "ch", "d_whole"]
    big = prim[prim.ch.abs() > 3].sort_values("ch")
    print("\nD. worst pairs (|log change| > 3):", len(big), "| share of total abs change", round(float(big.ch.abs().sum() / prim.ch.abs().sum()), 3))
    print(big[cols].round(2).to_string(index=False))
    dn, up = big[big.ch < 0], big[big.ch > 0]
    print("   falls:", len(dn), "| whole tumor also fell > 1 log:", int((dn.d_whole < -1).sum()), "| cavity grew > 1000 vox:", int((dn.d_cav > 1000).sum()), "| cavity appeared (0 -> >0):", int(((dn.n4_a == 0) & (dn.n4_b > 0)).sum()), "| edema grew by more than the lost core:", int(((dn.n2_b - dn.n2_a) > 0.5 * (dn.n1_a + dn.n3_a)).sum()))
    print("   rises:", len(up), "| whole tumor also rose > 1 log:", int((up.d_whole > 1).sum()), "| cavity shrank > 1000 vox:", int((up.d_cav < -1000).sum()), "| edema fell by more than half the gained core:", int(((up.n2_a - up.n2_b) > 0.5 * (up.n1_b + up.n3_b)).sum()))
    # how much error is label re-assignment between core and edema with a stable whole tumor?
    stable = (prim.d_whole.abs() < 0.3)
    print("\nE. pairs where the whole tumor (1+2+3) changed < 0.3 log:", int(stable.sum()), "| their mean |core change|", round(float(prim.ch.abs()[stable].mean()), 3), "| others", round(float(prim.ch.abs()[~stable].mean()), 3))
    print("   corr(core change, whole change)", round(float(np.corrcoef(prim.ch, prim.d_whole)[0, 1]), 3), "| corr(core change, edema log change)", round(float(np.corrcoef(prim.ch, np.log(prim.n2_b + 1) - np.log(prim.n2_a + 1))[0, 1]), 3))
    prim[["patient_id", "tp_in", "tp_out", "dt_days", "ch", "d_whole", "d_cav"] + [x for x in prim.columns if x[:2] in ("n1", "n2", "n3", "n4")]].to_csv(HERE / "mu_pairs_with_labels.csv", index=False)


if __name__ == "__main__":
    main()
