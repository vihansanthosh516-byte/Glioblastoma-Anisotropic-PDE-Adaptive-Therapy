"""Frozen v8 test on UCSF-ALPTDG time1 -> time2. Run once.

Protocol: paper/trackB_v8_ucsf_alptdg_protocol_2026-11-14.md (locked before scoring).
Reads segs straight from the downloaded zip; nothing is extracted to the repo.
"""
from __future__ import annotations

import gzip
import io
import json
import os
import zipfile
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

import track_b_mu_tcks as t
import track_b_mu_v3 as v3
import track_b_mu_v5 as v5
import track_b_v8_ucsd_external as ext

ROOT = Path(__file__).resolve().parent.parent
ZIP = Path(os.environ.get("UCSF_ALPTDG_ZIP", r"C:\Users\vihan\Downloads\UCSF_POSTOP_GLIOMA_DATASET_FINAL_v1.0.zip"))
TABLE = "UCSF_PostopGlioma_Table S1 R1 V5.0_UNBLINDED_FINAL.xlsx"
USED = ROOT / "output" / "predictgbm" / "h3_final" / "ucsf_test_pairs.csv"
OUT = ROOT / "output" / "track_b_v8_ucsf_alptdg_20261114"
V5C = ["log_edema", "nec_frac", "enh_frac", "log_whole", "compact"]


def num(x):
    return pd.to_numeric(x, errors="coerce")


def seg_from_zip(z, sid, tp):
    img = nib.Nifti1Image.from_bytes(gzip.decompress(z.read(f"{sid}/{sid}_time{tp}_seg.nii.gz")))
    return img


def build(z):
    x = pd.read_excel(io.BytesIO(z.read(TABLE)), sheet_name=["Clinical Info", "UCSF-LPTDG", "UCSF-PDGM Mapping"])
    c = x["Clinical Info"]
    c.columns = [str(k).strip() for k in c.columns]
    col = lambda key: next(k for k in c.columns if key.lower() in k.lower())
    lp = x["UCSF-LPTDG"]
    age = lp[lp.Timepoint == 1].set_index("SubjectID")["Patient Age"]
    used = set(pd.read_csv(USED).pid.astype(int)) | set(x["UCSF-PDGM Mapping"]["UCSF-LPTDG ID"].astype(int))
    rows = []
    for _, r in c.iterrows():
        sid = int(r["SubjectID"])
        day_in, dt = num(r[col("1st surgery/DX")]), num(r[col("1st scan to 2nd scan")])
        a1 = np.asarray(seg_from_zip(z, sid, 1).dataobj)
        a2 = np.asarray(seg_from_zip(z, sid, 2).dataobj)
        n1, n2, n3 = (int((a1 == k).sum()) for k in (1, 2, 3))
        core_in, core_out = n1 + n3, int(((a2 == 1) | (a2 == 3)).sum())
        box = 1
        if core_in:
            idx = np.argwhere((a1 == 1) | (a1 == 3))
            box = int(np.prod(idx.max(0) - idx.min(0) + 1))
        rt_s = day_in + num(r[col("1st RT start")])
        ad_s = day_in - num(r[col("1st chemo start to 1st scan")])
        rows.append({"patient_id": str(sid), "used_before": sid in used, "day_in": day_in, "day_out": day_in + dt,
                     "dt_days": dt, "core_vox_in": core_in, "core_vox_out": core_out,
                     "r_hist": np.nan, "dt_prior": np.nan, "rt_s": rt_s, "rt_e": rt_s + 42, "ad_s": ad_s, "ad_e": ad_s + 180,
                     "age": num(age.get(sid)), "mgmt": np.nan,
                     "gbm": float("glioblastoma" in str(r[col("WHO 2021 Diagnosis")]).lower()),
                     "log_edema": np.log(n2 + 1.0), "nec_frac": n1 / core_in if core_in else 0.0,
                     "enh_frac": n3 / core_in if core_in else 0.0, "log_whole": np.log(n1 + n2 + n3 + 1.0),
                     "compact": core_in / box if core_in else 0.0})
    d = pd.DataFrame(rows)
    n_all = len(d)
    d = d[(d.dt_days >= 14) & (d.dt_days <= 365) & (d.core_vox_in > 0)].reset_index(drop=True)
    d["y_in"], d["y_out"] = np.log(d.core_vox_in + 1.0), np.log(d.core_vox_out + 1.0)
    m = (d.day_in + d.day_out) / 2
    ph = np.full(len(d), "other", dtype=object)
    ph[((d.ad_s <= m) & (m <= d.ad_e)).values] = "adjuvant"
    ph[((d.rt_s <= m) & (m <= d.rt_e)).values] = "chemoradiation"
    d["phase"] = ph
    d["r_obs"] = (d.y_out - d.y_in) / d.dt_days
    return d, n_all


def score(d, preds, rng):
    y = d.y_out.values
    pat = pd.DataFrame({k: np.abs(v - y) for k, v in preds.items()}).assign(pid=d.patient_id.values).groupby("pid").mean()
    vs = {}
    for b in [k for k in pat.columns if k != "v8"]:
        diff = (pat.v8 - pat[b]).values
        vs[b] = {"mean": float(diff.mean()), "ci95": v3.boot_ci(diff, rng, 0.95)}
    return {"n_patients": int(len(pat)), "patient_mae_log_volume": pat.mean().to_dict(), "v8_minus": vs}


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "result.json").exists():
        raise SystemExit("locked test already run; refusing to overwrite")
    z = zipfile.ZipFile(ZIP)
    d, n_all = build(z)
    mu, _ = t.load()
    comp_mu = v5.features_table(mu).drop(columns=["patient_id", "tp_in"]).reset_index(drop=True)
    preds, p = ext.fit_predict(mu, d, comp_mu, d[V5C])
    rng = np.random.default_rng(t.SEED)
    prim = (~d.used_before).values
    res = {"script": "src/track_b_v8_ucsf_alptdg.py", "protocol": "paper/trackB_v8_ucsf_alptdg_protocol_2026-11-14.md",
           "n_patients_in_table": n_all, "n_eligible_pairs": int(len(d)),
           "phase_counts": d.phase.value_counts().to_dict(), "n_predicted_collapse": int((p > 0.5).sum())}
    pr = score(d[prim].reset_index(drop=True), {k: v[prim] for k, v in preds.items()}, rng)
    ci = pr["v8_minus"]["pop_rate"]["ci95"]
    res["primary_untouched"] = pr
    res["primary_success"] = bool(pr["n_patients"] >= 25 and ci[1] < 0)
    res["secondary_used_before"] = score(d[~prim].reset_index(drop=True), {k: v[~prim] for k, v in preds.items()}, rng)
    (OUT / "result.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    pd.DataFrame(preds).assign(patient_id=d.patient_id, used_before=d.used_before, y_out=d.y_out, phase=d.phase).to_csv(OUT / "predictions.csv", index=False)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    run()
