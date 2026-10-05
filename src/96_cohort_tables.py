#!/usr/bin/env python3
"""Script 96: missingness, eligible-vs-excluded, and cohort-shift tables (Phase 1.2).

Reads data/manifests/ (script 94) and the LUMIERE demographics. Writes
output/cohort_tables.json and output/cohort_shift.csv. Shift = standardised mean difference
(SMD) for continuous variables, difference in share for binary variables. These are descriptive;
no p-values, no tuning. LUMIERE is the external cohort: nothing here is used to change a model.

Not compared: tumour volume. MU volumes in the manifest are all-label volumes in mm3; LUMIERE
volumes are contrast-enhancing pyradiomics volumes. They are different quantities.
"""
from __future__ import annotations

import json
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
MAN = PROJECT_ROOT / "data" / "manifests"
LUM = PROJECT_ROOT / "data" / "external" / "lumiere"
OUT_JSON = PROJECT_ROOT / "output" / "cohort_tables.json"
OUT_CSV = PROJECT_ROOT / "output" / "cohort_shift.csv"


def smd(a: pd.Series, b: pd.Series) -> float:
    a, b = a.dropna().astype(float), b.dropna().astype(float)
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    sp = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return float((a.mean() - b.mean()) / sp) if sp > 0 else float("nan")


def summary(s: pd.Series) -> dict:
    s = s.dropna().astype(float)
    return {"n": int(len(s)), "mean": float(s.mean()), "median": float(s.median()),
            "q25": float(s.quantile(0.25)), "q75": float(s.quantile(0.75))} if len(s) else {"n": 0}


def main():
    pm = pd.read_csv(MAN / "patient_manifest.csv", low_memory=False)
    pairs = pd.read_csv(MAN / "forecast_pairs_mu.csv")
    mu = pm[pm["dataset"] == "MU"].copy()
    lu = pm[pm["dataset"] == "LUMIERE"].copy()
    demo = pd.read_csv(LUM / "LUMIERE-Demographics_Pathology.csv").rename(columns={"Patient": "patient_id"})
    demo = demo.replace("na", np.nan)
    lu = lu.merge(demo, on="patient_id", how="left")
    lu["age"] = pd.to_numeric(lu["Age at surgery (years)"], errors="coerce")
    lu["female"] = (lu["Sex"] == "female").astype(float).where(lu["Sex"].notna())
    lu["mgmt_meth"] = lu["MGMT qualitative"].map({"methylated": 1.0, "not methylated": 0.0})
    lu["idh_wt"] = lu["IDH (WT: wild type)"].map({"WT": 1.0}).where(lu["IDH (WT: wild type)"].notna())
    mu["female"] = (mu["sex"] == "Female").astype(float).where(mu["sex"].notna())
    mu["mgmt_meth"] = mu["mgmt_code"].map({1: 1.0, 0: 0.0})
    mu["gbm"] = mu["is_gbm"].astype(float)

    # scan intervals, same definition for both: gaps between consecutive scans, days
    mu_ids = set(mu.loc[mu["eligible"] == True, "patient_id"])  # noqa: E712
    mu_dt = pairs[pairs["in_primary"] & pairs["patient_id"].isin(mu_ids)]["dt_days"]
    sc = pd.read_csv(MAN / "scan_manifest.csv", low_memory=False)
    ls = sc[(sc["dataset"] == "LUMIERE") & sc["rating"].isin(["PD", "SD", "PR", "CR"])].sort_values(["patient_id", "day"])
    lu_ids = set(lu.loc[lu["eligible"] == True, "patient_id"])  # noqa: E712
    lu_dt = ls[ls["patient_id"].isin(lu_ids)].groupby("patient_id")["day"].diff().dropna()

    mu_e = mu[mu["eligible"] == True]  # noqa: E712
    lu_e = lu[lu["eligible"] == True]  # noqa: E712
    rows = []

    def add(name, a, b, kind):
        rows.append({"variable": name, "kind": kind, "MU_eligible": summary(a)["mean"] if kind == "binary" and len(a.dropna()) else summary(a).get("median"),
                     "MU_n": int(a.notna().sum()), "LUMIERE_eligible": summary(b)["mean"] if kind == "binary" and len(b.dropna()) else summary(b).get("median"),
                     "LUMIERE_n": int(b.notna().sum()),
                     "shift": smd(a, b) if kind == "continuous" else (float(a.mean() - b.mean()) if a.notna().any() and b.notna().any() else float("nan"))})

    add("age (years; median)", mu_e["age"], lu_e["age"], "continuous")
    add("female (share)", mu_e["female"], lu_e["female"], "binary")
    add("MGMT methylated (share of known)", mu_e["mgmt_meth"], lu_e["mgmt_meth"], "binary")
    add("scans per patient (median)", mu_e["n_scans"], lu_e["n_scans"], "continuous")
    add("interval between consecutive scans (days; median)", mu_dt, lu_dt, "continuous")
    shift = pd.DataFrame(rows)
    shift.to_csv(OUT_CSV, index=False)

    def miss(df, cols):
        return {c: {"missing": int(df[c].isna().sum()), "n": int(len(df))} for c in cols}

    # selection check inside MU: eligible vs excluded
    sel = {}
    for name, col, kind in [("age", "age", "continuous"), ("female", "female", "binary"),
                            ("gbm", "gbm", "binary"), ("mgmt_methylated", "mgmt_meth", "binary"),
                            ("died", "died", "binary"), ("progression", "progression", "binary")]:
        a, b = mu.loc[mu["eligible"] == True, col], mu.loc[mu["eligible"] == False, col]  # noqa: E712
        sel[name] = {"eligible_mean": float(a.mean()), "excluded_mean": float(b.mean()),
                     "smd": smd(a, b), "n_eligible": int(a.notna().sum()), "n_excluded": int(b.notna().sum())}

    out = {
        "script": "96_cohort_tables",
        "mu_missingness": {
            **miss(mu, ["age", "sex", "primary_diagnosis", "mgmt_code", "radiation_end_day"]),
            "mgmt_indeterminate_or_unknown_codes_2_3_4": int(mu["mgmt_code"].isin([2, 3, 4]).sum()),
            "scan_day_missing_patients": int((mu["reason_excluded"].fillna("").str.contains("missing scan day")).sum()),
            "treatment_schedule_unknown": int((~mu["treatment_known"].astype(bool)).sum()),
            "segmentation_file_missing": int((~mu["segment_available"].astype(bool)).sum()),
        },
        "lumiere_missingness": {
            **miss(lu, ["age", "female", "mgmt_meth", "idh_wt"]),
            "masks_available": False,
        },
        "mu_eligible_vs_excluded": sel,
        "mu_diagnosis_mix_eligible": mu_e["primary_diagnosis"].value_counts().to_dict(),
        "shift_table_MU_eligible_vs_LUMIERE_eligible": shift.to_dict(orient="records"),
        "notes": ["Shift is descriptive. LUMIERE age/sex are for all patients with 3+ rated follow-ups.",
                  "MU mixes WHO grades; see mu_diagnosis_mix_eligible. A GBM-only sensitivity analysis is needed.",
                  "Tumour volume is not compared (different quantities)."],
    }
    OUT_JSON.write_text(json.dumps(out, indent=1, default=lambda o: None if o != o else str(o)))
    print(json.dumps(out, indent=1, default=str)[:6000])


if __name__ == "__main__":
    main()
