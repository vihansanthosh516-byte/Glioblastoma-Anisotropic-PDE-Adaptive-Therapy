#!/usr/bin/env python3
"""Script 94: canonical patient manifest, patient-level split, CONSORT counts (Phase 1.1, 1.2).

One source of truth for who is eligible and in which fold (analysis_plan_v1.md A1.1-A1.3).
Every later experiment must read these files, never re-derive eligibility.

Inputs : output/mu_glioma_cohort.json (+ masks under data/tcia/MU-Glioma-Post),
         data/external/lumiere/*.csv, data/external/rhuh/clinical_data_TCIA_RHUH-GBM.csv,
         data/tcga_gbm_clinical.csv, data/external/cgga/*clinical*.txt
Outputs: data/manifests/patient_manifest.csv   one row per patient (all datasets)
         data/manifests/scan_manifest.csv      one row per scan (MU, LUMIERE)
         data/manifests/forecast_pairs_mu.csv  one row per consecutive-scan forecast (MU)
         data/manifests/split_mu.csv           patient -> fold (5 folds, patient level)
         data/manifests/consort_mu.json        counts and exclusion reasons
Rules  : target = labels {1, 3} (core), as script 81. Interval window 14-365 d for the primary
         population; other pairs are kept with in_primary=False. Nothing is dropped silently.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
DATA = PROJECT_ROOT / "data"
MAN = DATA / "manifests"
COHORT = OUTPUT_DIR / "mu_glioma_cohort.json"
MU_CLINICAL = DATA / "tcia" / "MU-Glioma-Post_ClinicalData-July2025.xlsx"
N_FOLDS = 5
SPLIT_SEED = "20261004"
DT_MIN, DT_MAX = 14.0, 365.0
CORE_LABELS = (1, 3)


def core_voxels(mask_path: _Path):
    """Native-resolution core voxel count, or None if the file is missing/unreadable."""
    import nibabel as nib
    try:
        arr = np.asarray(nib.load(str(mask_path)).dataobj)
    except Exception:
        return None
    return int(np.isin(arr, CORE_LABELS).sum())


def fold_of(pid: str) -> int:
    """Deterministic patient-level fold from a hash, independent of row order."""
    h = hashlib.sha256(f"{SPLIT_SEED}|{pid}".encode()).hexdigest()
    return int(h, 16) % N_FOLDS


def mu_clinical() -> pd.DataFrame:
    """Per-patient clinical fields used for subgroup and shift tables (MGMT codes per the data dictionary:
    0 none, 1 methylated, 2 indeterminate, 3 unable, 4 unknown)."""
    c = pd.read_excel(MU_CLINICAL, sheet_name="MU Glioma Post")
    out = pd.DataFrame({
        "patient_id": c["Patient_ID"], "age": c["Age at diagnosis"], "sex": c["Sex at Birth"],
        "primary_diagnosis": c["Primary Diagnosis"], "who_grade": c["Grade of Primary Brain Tumor"].astype(str),
        "mgmt_code": c["MGMT methylation"], "progression": c["Progression"],
        "died": c["Overall Survival (Death)"]})
    out["is_gbm"] = out["primary_diagnosis"].isin(["GBM", "Glioma w/ GBM features"])
    return out


def build_mu():
    cohort = json.loads(COHORT.read_text())
    scans, pairs, pats = [], [], []
    for p in cohort:
        pid = p["patient_id"]
        tps = p["timepoints"]
        ts = p.get("treatment_schedule") or {}
        reasons = []
        days = [t["day_from_diagnosis"] for t in tps]
        if any(d is None for d in days):
            reasons.append("missing scan day")
        known = [d for d in days if d is not None]
        if any(b <= a for a, b in zip(known, known[1:])):
            reasons.append("non-increasing scan days")
        vox = {}
        for t in tps:
            mp = PROJECT_ROOT / t["mask_path"]
            v = core_voxels(mp) if mp.exists() else None
            vox[t["number"]] = v
            scans.append({"patient_id": pid, "dataset": "MU", "scan_number": t["number"],
                          "day": t["day_from_diagnosis"], "mask_exists": mp.exists(),
                          "core_voxels_native": v, "volume_mm3_all_labels": t["volume_mm3"]})
        n_pairs = n_primary = 0
        for a, b in zip(tps, tps[1:]):
            da, db = a["day_from_diagnosis"], b["day_from_diagnosis"]
            if da is None or db is None:
                continue
            dt = db - da
            if dt <= 0:
                continue
            va, vb = vox[a["number"]], vox[b["number"]]
            ok_input = va is not None and va > 0
            in_primary = ok_input and (DT_MIN <= dt <= DT_MAX) and vb is not None
            n_pairs += 1
            n_primary += int(in_primary)
            pairs.append({"patient_id": pid, "tp_in": a["number"], "tp_out": b["number"],
                          "day_in": da, "day_out": db, "dt_days": dt,
                          "core_vox_in": va, "core_vox_out": vb,
                          "has_prior_scan": a is not tps[0],
                          "input_core_nonempty": ok_input, "in_primary": in_primary})
        eligible = (not reasons) and n_primary > 0
        if not eligible and not reasons:
            reasons.append("no pair in primary window with non-empty input core")
        pats.append({"patient_id": pid, "dataset": "MU", "eligible": eligible,
                     "reason_excluded": "; ".join(reasons),
                     "n_scans": len(tps), "n_forecast_pairs": n_pairs, "n_primary_pairs": n_primary,
                     "scan_ids": ",".join(str(t["number"]) for t in tps),
                     "scan_days": ",".join("" if d is None else f"{d:g}" for d in days),
                     "radiation_end_day": ts.get("radiation_end_day"),
                     "treatment_known": bool(ts), "segment_available": all(vox[t["number"]] is not None for t in tps),
                     "molecular_available": False, "fold": fold_of(pid)})
    pats = pd.DataFrame(pats).merge(mu_clinical(), on="patient_id", how="left")
    return pats, pd.DataFrame(scans), pd.DataFrame(pairs)


def build_lumiere():
    d = DATA / "external" / "lumiere"
    demo = pd.read_csv(d / "LUMIERE-Demographics_Pathology.csv")
    rate = pd.read_csv(d / "LUMIERE-ExpertRating-v202211.csv")
    rcol = [c for c in rate.columns if c.startswith("Rating")][0]
    rate = rate.rename(columns={rcol: "rating"})
    rate["week"] = rate["Date"].str.extract(r"week-(\d+)").astype(float)
    scans, pats = [], []
    for pid, g in rate.groupby("Patient"):
        post = g[~g["rating"].isin(["Pre-Op"])]
        follow = g[g["rating"].isin(["PD", "SD", "PR", "CR"])]
        scans += [{"patient_id": pid, "dataset": "LUMIERE", "scan_number": i, "day": w * 7,
                   "mask_exists": False, "core_voxels_native": None, "volume_mm3_all_labels": None,
                   "rating": r} for i, (w, r) in enumerate(zip(g["week"], g["rating"]))]
        pats.append({"patient_id": pid, "dataset": "LUMIERE", "eligible": len(follow) >= 3,
                     "reason_excluded": "" if len(follow) >= 3 else "fewer than 3 rated follow-up scans",
                     "n_scans": len(g), "n_forecast_pairs": max(len(post) - 1, 0), "n_primary_pairs": None,
                     "scan_ids": ",".join(g["Date"]), "scan_days": ",".join(f"{w * 7:g}" for w in g["week"]),
                     "radiation_end_day": None, "treatment_known": False,
                     "segment_available": False, "molecular_available": bool(
                         demo.loc[demo["Patient"] == pid, "MGMT qualitative"].astype(str).isin(
                             ["methylated", "not methylated"]).any()),
                     "fold": -1})
    return pd.DataFrame(pats), pd.DataFrame(scans)


def build_other():
    rows = []
    r = pd.read_csv(DATA / "external" / "rhuh" / "clinical_data_TCIA_RHUH-GBM.csv")
    for pid in r.iloc[:, 0]:
        rows.append({"patient_id": pid, "dataset": "RHUH", "eligible": None,
                     "reason_excluded": "not yet screened", "n_scans": None,
                     "segment_available": None, "molecular_available": False, "fold": -1})
    t = pd.read_csv(DATA / "tcga_gbm_clinical.csv", encoding="utf-8-sig")
    expr_ids = None
    for pid in t["sample"]:
        rows.append({"patient_id": pid, "dataset": "TCGA-GBM", "eligible": None,
                     "reason_excluded": "survival cohort; see script 36-37", "n_scans": 0,
                     "segment_available": False, "molecular_available": None, "fold": -1})
    for f in ["CGGA.mRNAseq_693_clinical.20200506.txt", "CGGA.mRNAseq_325_clinical.20200506.txt"]:
        c = pd.read_csv(DATA / "external" / "cgga" / f, sep="\t")
        for pid in c["CGGA_ID"]:
            rows.append({"patient_id": f"{f.split('_')[1]}:{pid}", "dataset": "CGGA", "eligible": None,
                         "reason_excluded": "survival cohort; see script 90", "n_scans": 0,
                         "segment_available": False, "molecular_available": True, "fold": -1})
    return pd.DataFrame(rows)


def main():
    MAN.mkdir(parents=True, exist_ok=True)
    mu_p, mu_s, mu_pairs = build_mu()
    lu_p, lu_s = build_lumiere()
    other = build_other()
    pd.concat([mu_p, lu_p, other], ignore_index=True).to_csv(MAN / "patient_manifest.csv", index=False)
    pd.concat([mu_s, lu_s], ignore_index=True).to_csv(MAN / "scan_manifest.csv", index=False)
    mu_pairs.to_csv(MAN / "forecast_pairs_mu.csv", index=False)
    mu_p[["patient_id", "fold"]].to_csv(MAN / "split_mu.csv", index=False)

    prim = mu_pairs[mu_pairs["in_primary"]]
    per_pat = prim.groupby("patient_id").size()
    first_pairs = mu_pairs[~mu_pairs["has_prior_scan"]]
    consort = {
        "script": "94_patient_manifest",
        "split": {"seed": SPLIT_SEED, "n_folds": N_FOLDS, "unit": "patient"},
        "mu": {
            "patients_in_cohort_file": int(len(mu_p)),
            "excluded_by_reason": mu_p.loc[~mu_p["eligible"], "reason_excluded"].value_counts().to_dict(),
            "eligible_patients": int(mu_p["eligible"].sum()),
            "patients_with_ge2_scans": int((mu_p["n_scans"] >= 2).sum()),
            "patients_with_ge3_scans": int((mu_p["n_scans"] >= 3).sum()),
            "eligible_with_ge2_primary_pairs": int((per_pat >= 2).sum()),
            "eligible_with_ge3_primary_pairs": int((per_pat >= 3).sum()),
            "forecast_pairs_all": int(len(mu_pairs)),
            "forecast_pairs_primary": int(len(prim)),
            "pairs_dropped_dt_lt14": int((mu_pairs["dt_days"] < DT_MIN).sum()),
            "pairs_dropped_dt_gt365": int((mu_pairs["dt_days"] > DT_MAX).sum()),
            "pairs_empty_input_core": int((~mu_pairs["input_core_nonempty"]).sum()),
            "pairs_with_prior_scan_primary": int(prim["has_prior_scan"].sum()),
            "folds_patient_counts": mu_p["fold"].value_counts().sort_index().to_dict(),
            "gbm_by_primary_diagnosis_all": int(mu_p["is_gbm"].sum()),
            "gbm_among_eligible": int((mu_p["is_gbm"] & mu_p["eligible"]).sum()),
            "non_gbm_among_eligible": int((~mu_p["is_gbm"] & mu_p["eligible"]).sum()),
            "primary_diagnosis_counts_eligible": mu_p.loc[mu_p["eligible"], "primary_diagnosis"].value_counts().to_dict(),
        },
        "rolling_origin_rule_A1_1": ("primary" if int((per_pat >= 2).sum()) >= 30 else "secondary"),
        "lumiere": {"patients": int(len(lu_p)), "eligible_ge3_rated_followups": int(lu_p["eligible"].sum()),
                    "masks_available": False},
        "rhuh_patients": int((other["dataset"] == "RHUH").sum()),
        "tcga_rows": int((other["dataset"] == "TCGA-GBM").sum()),
        "cgga_rows": int((other["dataset"] == "CGGA").sum()),
    }
    (MAN / "consort_mu.json").write_text(json.dumps(consort, indent=1, default=int))
    print(json.dumps(consort, indent=1, default=int))


if __name__ == "__main__":
    main()
