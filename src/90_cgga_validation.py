"""Script 90: external check of the TCGA-GBM survival findings on CGGA (limit L7).

Tests fixed BEFORE looking at any CGGA result (same covariates as script 36 on TCGA-GBM):
  T1  age (per year), univariate Cox on overall survival. TCGA: HR 1.031.  Expected: HR > 1.
  T2  sex (male vs female), univariate Cox. TCGA: HR 1.15, null.
  T3  MGMT promoter (methylated vs un-methylated), univariate Cox. Not tested on TCGA here.
  Holm correction over T1-T3 inside each cohort.
Also reported: Kaplan-Meier median OS (TCGA: 432 d), a multivariate Cox (age + sex + MGMT + IDH), and the
same tests on IDH-wild-type primary GBM only.

Cohort: WHO grade IV rows with OS > 0 and a censor flag. Batch 1 = mRNAseq_693, batch 2 = mRNAseq_325.
Samples that appear in both batches are kept in batch 1 only (checked by CGGA_ID).
Not done here: the gene-expression survival models (needs the TCGA gene weights).

Input:  data/external/cgga/CGGA.mRNAseq_{693,325}_clinical.20200506.txt
Output: output/cgga_validation.json
"""
import json
from pathlib import Path as _Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
DATA = PROJECT_ROOT / "data" / "external" / "cgga"
OUT_JSON = PROJECT_ROOT / "output" / "cgga_validation.json"
CENSOR = "Censor (alive=0; dead=1)"


def load(batch):
    d = pd.read_csv(DATA / f"CGGA.mRNAseq_{batch}_clinical.20200506.txt", sep="\t")
    for c in ("Gender", "Grade", "IDH_mutation_status", "MGMTp_methylation_status", "PRS_type", "Histology"):
        d[c] = d[c].astype(str).str.strip()
    d["OS"] = pd.to_numeric(d["OS"], errors="coerce")
    d["event"] = pd.to_numeric(d[CENSOR], errors="coerce")
    d["age"] = pd.to_numeric(d["Age"], errors="coerce")
    d["male"] = (d["Gender"] == "Male").astype(float)
    d["mgmt_meth"] = d["MGMTp_methylation_status"].map({"methylated": 1.0, "un-methylated": 0.0})
    d["idh_wt"] = d["IDH_mutation_status"].map({"Wildtype": 1.0, "Mutant": 0.0})
    return d


def gbm(d):
    g = d[(d["Grade"] == "WHO IV") & (d["OS"] > 0) & d["event"].notna() & d["age"].notna()].copy()
    return g


def holm(ps):
    order = np.argsort(ps)
    m = len(ps)
    adj = np.empty(m)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, (m - rank) * ps[i])
        adj[i] = min(1.0, run)
    return adj


def cox_uni(d, col):
    x = d[["OS", "event", col]].dropna()
    cph = CoxPHFitter().fit(x, "OS", "event")
    r = cph.summary.loc[col]
    return {"covariate": col, "n": int(len(x)), "events": int(x["event"].sum()), "hr": float(r["exp(coef)"]),
            "ci_lo": float(r["exp(coef) lower 95%"]), "ci_hi": float(r["exp(coef) upper 95%"]), "p": float(r["p"])}


def km_median(d):
    k = KaplanMeierFitter().fit(d["OS"], d["event"])
    from lifelines.utils import median_survival_times
    ci = median_survival_times(k.confidence_interval_)
    return {"median_days": float(k.median_survival_time_), "ci_lo": float(ci.iloc[0, 0]), "ci_hi": float(ci.iloc[0, 1])}


def analyse(d):
    out = {"n": int(len(d)), "events": int(d["event"].sum()), "km": km_median(d), "tests": {}}
    res = [cox_uni(d, c) for c in ("age", "male", "mgmt_meth")]
    adj = holm([r["p"] for r in res])
    for r, a in zip(res, adj):
        r["holm_p"] = float(a)
        out["tests"][r["covariate"]] = r
    cols = ["age", "male", "mgmt_meth", "idh_wt"]
    m = d[["OS", "event"] + cols].dropna()
    if m["idh_wt"].nunique() > 1 and len(m) > 30:
        cph = CoxPHFitter().fit(m, "OS", "event")
    else:
        cols = ["age", "male", "mgmt_meth"]
        cph = CoxPHFitter().fit(m[["OS", "event"] + cols], "OS", "event")
    s = cph.summary
    out["multivariate"] = {"n": int(len(m)), "features": cols, "c_index": float(cph.concordance_index_),
                           "hr": {c: float(s.loc[c, "exp(coef)"]) for c in cols},
                           "ci_lo": {c: float(s.loc[c, "exp(coef) lower 95%"]) for c in cols},
                           "ci_hi": {c: float(s.loc[c, "exp(coef) upper 95%"]) for c in cols},
                           "p": {c: float(s.loc[c, "p"]) for c in cols}}
    return out


def main():
    b1, b2 = load("693"), load("325")
    overlap = sorted(set(b1["CGGA_ID"]) & set(b2["CGGA_ID"]))
    b2u = b2[~b2["CGGA_ID"].isin(b1["CGGA_ID"])]
    res = {"script": "90_cgga_validation", "tcga_reference": {"age_hr": 1.0307, "gender_hr": 1.1457, "km_median_days": 432},
           "batch_overlap_ids": len(overlap), "cohorts": {}}
    sets = {"batch1_693_all_gbm": gbm(b1), "batch2_325_all_gbm": gbm(b2u)}
    sets["batch1_693_primary_idhwt"] = sets["batch1_693_all_gbm"].query("PRS_type == 'Primary' and idh_wt == 1")
    sets["batch2_325_primary_idhwt"] = sets["batch2_325_all_gbm"].query("PRS_type == 'Primary' and idh_wt == 1")
    both = pd.concat([sets["batch1_693_all_gbm"], sets["batch2_325_all_gbm"]])
    sets["pooled_all_gbm"] = both
    sets["pooled_primary_idhwt"] = both.query("PRS_type == 'Primary' and idh_wt == 1")
    for name, d in sets.items():
        res["cohorts"][name] = analyse(d)
        t = res["cohorts"][name]["tests"]
        print(name, "n", len(d), "ev", int(d["event"].sum()), "KM", round(res["cohorts"][name]["km"]["median_days"]),
              "age HR %.3f p %.2g" % (t["age"]["hr"], t["age"]["p"]),
              "sex HR %.2f p %.2g" % (t["male"]["hr"], t["male"]["p"]),
              "MGMT HR %.2f p %.2g" % (t["mgmt_meth"]["hr"], t["mgmt_meth"]["p"]))
    OUT_JSON.write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
