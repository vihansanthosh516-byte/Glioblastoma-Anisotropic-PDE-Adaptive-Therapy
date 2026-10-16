#!/usr/bin/env python3
"""Script 129: RANO-aware check X5, exploratory, no new model runs
(Amendment 10 A10.6).

(a) Burdenko open CSV: pseudo-progression and distant-progression shares per
    follow-up and per patient; timing vs RT end (<=84 d vs later).
(b) LUMIERE: frozen H-3v test errors (output/h3v_test_pairs_*.csv, which carry
    the target-scan RANO rating) stratified PD vs non-PD; rationale pseudo
    mentions counted from LUMIERE-ExpertRating-v202211.csv.
(c) PREDICT-GBM/MU have no RANO field: contextual only.
Output: output/rano_check.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def burdenko() -> dict:
    df = pd.read_csv(PROJECT_ROOT / "data" / "external" / "burdenko_clinical.csv")
    resp = [c for c in df.columns if c.startswith("Response_")]
    datec = [c for c in df.columns if c.startswith("Anonym") and "fup" in c.lower()
             and "esponse" not in c]
    rt_end = pd.to_datetime(df["AnonymEnd of RT"], errors="coerce")
    per_fup, per_pt_pseudo, per_pt_dist = [], set(), set()
    early = {"pseudo": 0, "n": 0}
    late = {"pseudo": 0, "n": 0}
    for i, col in enumerate(resp):
        # date columns keep the ordinal: Response_1st_fup -> Anonym1st_fup
        dcol = col.replace("Response_", "Anonym")
        dd = pd.to_datetime(df.get(dcol, pd.Series([pd.NaT] * len(df))),
                            errors="coerce")
        dt = (dd - rt_end).dt.days
        for pid, val, days in zip(df["AnonymPatientID"], df[col], dt):
            if pd.isna(val) or str(val).strip() == "" or str(val) == "nan":
                continue
            v = str(val).strip().lower()
            pseudo = "pseudo" in v
            dist = "distant" in v
            per_fup.append({"pseudo": pseudo, "distant": dist,
                            "early84": bool(pd.notna(days) and days <= 84)})
            if pseudo:
                per_pt_pseudo.add(pid)
            if dist:
                per_pt_dist.add(pid)
            if pd.notna(days):
                bucket = early if days <= 84 else late
                bucket["n"] += 1
                bucket["pseudo"] += int(pseudo)
    F = pd.DataFrame(per_fup)
    return {
        "n_rated_followups": int(len(F)),
        "pseudo_share_followups": float(F["pseudo"].mean()),
        "distant_share_followups": float(F["distant"].mean()),
        "n_patients": int(df["AnonymPatientID"].nunique()),
        "pseudo_ever_share_patients": float(len(per_pt_pseudo) / df["AnonymPatientID"].nunique()),
        "distant_ever_share_patients": float(len(per_pt_dist) / df["AnonymPatientID"].nunique()),
        "early_le84d": {"n": early["n"],
                        "pseudo_share": early["pseudo"] / early["n"] if early["n"] else None},
        "late_gt84d": {"n": late["n"],
                       "pseudo_share": late["pseudo"] / late["n"] if late["n"] else None},
    }


def lumiere() -> dict:
    out = {}
    rat = pd.read_csv(PROJECT_ROOT / "data" / "external" / "lumiere" /
                      "LUMIERE-ExpertRating-v202211.csv")
    rat.columns = ["Patient", "tp", "lt3", "nonmeas", "rating", "rationale"]
    r = rat["rationale"].astype(str).str.lower()
    out["rationale_pseudo_mentions"] = int(r.str.contains("pseudo").sum())
    out["rationale_rows"] = int(len(rat))
    out["rating_values"] = sorted(rat["rating"].astype(str).str.strip().unique().tolist())
    for name in ("hdglioauto", "deepbratumia"):
        p = pd.read_csv(PROJECT_ROOT / "output" / f"h3v_test_pairs_{name}.csv")
        p["pd"] = (p["rating_target"] == "PD").astype(int)
        strata = {}
        for lab, sub in (("PD", p[p["pd"] == 1]), ("nonPD", p[p["pd"] == 0])):
            if len(sub) == 0:
                strata[lab] = {"n": 0}
                continue
            strata[lab] = {"n": int(len(sub)),
                           "mean_ae_mix": float(sub["ae_mix"].mean()),
                           "mean_ae_persist": float(sub["ae_persist"].mean()),
                           "mean_gain_persist_minus_mix": float(
                               (sub["ae_persist"] - sub["ae_mix"]).mean())}
        out[name] = strata
    return out


def main() -> None:
    res = {"script": "129_rano_check", "burdenko": burdenko(), "lumiere": lumiere(),
           "note": "PREDICT-GBM and MU have no RANO field; their nulls are "
                   "contextualised with these rates, not re-run."}
    (PROJECT_ROOT / "output" / "rano_check.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
