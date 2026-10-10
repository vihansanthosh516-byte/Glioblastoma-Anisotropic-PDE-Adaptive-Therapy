#!/usr/bin/env python3
"""Script 118: GRAND_PLAN 14 Y5 (age adjusted for IDH and treatment) and Y2 (MGMT within TMZ-treated), CGGA GBM.

Declared in GRAND_PLAN section 14 before this script ran. Same cohort rules as script 90 (WHO IV, OS > 0, censor flag,
duplicate samples kept in batch 1 only); batches pooled with a batch indicator.
Y5: one Cox model on overall survival: age (per year) + IDH wild-type + radiotherapy + TMZ + recurrent (vs primary)
    + batch; complete cases. Reported even if age stays weak.
Y2: MGMT methylated vs un-methylated, (a) univariate and (b) adjusted for age, IDH, radiotherapy, recurrent, batch,
    in TMZ-treated patients only; and the interaction MGMT x TMZ in all patients (the literature's predictive claim,
    meta-analysis PFS HR 0.48 with TMZ). Endpoint here is OS (CGGA has no PFS), so this is a consistency check.
Output: output/cgga_adjusted_age_mgmt.json
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path as _Path

import numpy as np
from lifelines import CoxPHFitter

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT_JSON = PROJECT_ROOT / "output" / "cgga_adjusted_age_mgmt.json"
_s = importlib.util.spec_from_file_location("s90", PROJECT_ROOT / "src" / "90_cgga_validation.py")
s90 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(s90)


def cohort():
    import pandas as pd
    b1, b2 = s90.gbm(s90.load("693")), s90.gbm(s90.load("325"))
    b2 = b2[~b2["CGGA_ID"].isin(b1["CGGA_ID"])]
    b1["batch2"], b2["batch2"] = 0.0, 1.0
    d = pd.concat([b1, b2], ignore_index=True)
    d["radio"] = pd.to_numeric(d["Radio_status (treated=1;un-treated=0)"], errors="coerce")
    d["tmz"] = pd.to_numeric(d["Chemo_status (TMZ treated=1;un-treated=0)"], errors="coerce")
    d["recurrent"] = (d["PRS_type"] == "Recurrent").astype(float)
    return d


def cox(d, cols):
    x = d[["OS", "event"] + cols].dropna()
    f = CoxPHFitter().fit(x, "OS", "event")
    s = f.summary
    return {"n": int(len(x)), "events": int(x["event"].sum()),
            "terms": {c: {"hr": float(s.loc[c, "exp(coef)"]), "ci95": [float(s.loc[c, "exp(coef) lower 95%"]),
                                                                         float(s.loc[c, "exp(coef) upper 95%"])],
                          "p": float(s.loc[c, "p"])} for c in cols}}


def main():
    d = cohort()
    res = {"script": "118_cgga_adjusted_age_mgmt", "declared": "GRAND_PLAN 14 Y5, Y2", "n_gbm_pooled": int(len(d)),
           "tcga_reference_age_hr": 1.0307,
           "Y5_age_adjusted": cox(d, ["age", "idh_wt", "radio", "tmz", "recurrent", "batch2"]),
           "Y5_age_univariate_same_rows": None}
    rows = d[["OS", "event", "age", "idh_wt", "radio", "tmz", "recurrent", "batch2"]].dropna()
    res["Y5_age_univariate_same_rows"] = cox(rows, ["age"])
    t = d[d["tmz"] == 1]
    res["Y2_mgmt_tmz_univariate"] = cox(t, ["mgmt_meth"])
    res["Y2_mgmt_tmz_adjusted"] = cox(t, ["mgmt_meth", "age", "idh_wt", "radio", "recurrent", "batch2"])
    d["mgmt_x_tmz"] = d["mgmt_meth"] * d["tmz"]
    res["Y2_mgmt_by_tmz_interaction"] = cox(d, ["mgmt_meth", "tmz", "mgmt_x_tmz", "age", "idh_wt", "radio", "recurrent", "batch2"])
    res["weakest_points"] = ["OS, not PFS/TTP; treatment recorded as yes/no only (no dose, timing).",
                             "Complete-case analysis; CGGA treatment and MGMT have missing values.",
                             "Retrospective; treatment is not randomised, so the TMZ effects are confounded by who got TMZ."]
    OUT_JSON.write_text(json.dumps(res, indent=1))
    for k in ("Y5_age_adjusted", "Y5_age_univariate_same_rows", "Y2_mgmt_tmz_univariate", "Y2_mgmt_tmz_adjusted",
              "Y2_mgmt_by_tmz_interaction"):
        r = res[k]
        print(k, "n", r["n"], "ev", r["events"], {c: (round(v["hr"], 3), [round(x, 3) for x in v["ci95"]], round(v["p"], 4))
                                                   for c, v in r["terms"].items() if c in ("age", "mgmt_meth", "mgmt_x_tmz")})


if __name__ == "__main__":
    main()
