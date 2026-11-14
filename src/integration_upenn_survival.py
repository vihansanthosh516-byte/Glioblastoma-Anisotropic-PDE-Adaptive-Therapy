"""Supporting test: does the frozen TCGA MRI MES score predict survival in UPENN? Run once.

Protocol: paper/integration_protocol_backup_upenn_2026-11-12.md. Coding: paper/integration_protocol_bamd_note1_2026-11-12.md.
Score = locked model (neftel_MES, 'new' features, pls2) fit on all TCGA, frozen, applied to UPENN.
Primary: Cox time ~ score_z + age + logvol_WT; one-sided Wald on score_z (higher score = higher hazard).
Secondary (descriptive): Harrell C of age + volume with and without the score, 10 x 5-fold CV.
Volume = logvol_WT from output/integration_bamd/features.csv (same masks).

  py -3.14 -I src/integration_upenn_survival.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy.stats import norm

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import integration_methoddev_search as ms  # noqa: E402

OUT = PROJECT_ROOT / "output" / "integration_methoddev"
CLIN = PROJECT_ROOT / "data" / "external" / "upenn_meta" / "clin.csv"
SEED = 20261115


def cohort() -> pd.DataFrame:
    c = pd.read_csv(CLIN)
    c = c[c["IDH1"] != "Mutated"].copy()
    c["event"] = (c["Survival_Status"] == "Deceased").astype(int)
    c["time"] = np.where(c["event"] == 1, c["Survival_from_surgery_days_UPDATED"], c["Survival_Censor"])
    c["time"] = pd.to_numeric(c["time"], errors="coerce")
    return c.rename(columns={"ID": "patient_id", "Age_at_scan_years": "age"})[["patient_id", "age", "time", "event"]]


def main() -> int:
    d, T, _old, newcols = ms.load()
    Xtr = d[newcols].to_numpy(float)
    med = np.nanmedian(Xtr, axis=0)
    med = np.where(np.isnan(med), 0.0, med)
    model = ms.make_model("pls2").fit(np.where(np.isnan(Xtr), med, Xtr), T["neftel_MES"].to_numpy(float))

    f = pd.read_csv(OUT / "features_UPENN-GBM.csv")
    vol = pd.read_csv(PROJECT_ROOT / "output" / "integration_bamd" / "features.csv")[["brats_id", "logvol_WT"]]
    f = f.merge(vol, on="brats_id", validate="one_to_one")
    X = f[newcols].to_numpy(float)
    f["score"] = np.ravel(model.predict(np.where(np.isnan(X), med, X)))
    df = f[["patient_id", "score", "logvol_WT"]].merge(cohort(), on="patient_id", validate="one_to_one")
    df = df.dropna(subset=["score", "logvol_WT", "age", "time"])
    df = df[df["time"] > 0].reset_index(drop=True)
    df["score_z"] = (df["score"] - df["score"].mean()) / df["score"].std()

    cph = CoxPHFitter().fit(df[["time", "event", "score_z", "age", "logvol_WT"]], "time", "event")
    b, se = float(cph.params_["score_z"]), float(cph.standard_errors_["score_z"])
    p_one = float(norm.sf(b / se))  # H1: beta > 0

    rng = np.random.Generator(np.random.PCG64DXSM(SEED))
    cs = {"base": [], "with_score": []}
    for _ in range(10):
        idx = rng.permutation(len(df))
        for name, cols in (("base", ["age", "logvol_WT"]), ("with_score", ["age", "logvol_WT", "score_z"])):
            risk = np.empty(len(df))
            for k in range(5):
                te = idx[k::5]
                tr = np.setdiff1d(idx, te)
                m = CoxPHFitter(penalizer=1e-4).fit(df.loc[tr, ["time", "event"] + cols], "time", "event")
                risk[te] = np.ravel(m.predict_partial_hazard(df.loc[te, cols]))
            cs[name].append(concordance_index(df["time"], -risk, df["event"]))
    res = {"n": len(df), "events": int(df["event"].sum()),
           "cox_score_z": {"hr": float(np.exp(b)), "beta": b, "se": se, "p_one_sided": p_one},
           "cox_age_hr": float(np.exp(cph.params_["age"])), "cox_logvolWT_hr": float(np.exp(cph.params_["logvol_WT"])),
           "pass_rule": "one-sided p < 0.05, higher score = higher hazard", "passed": bool(p_one < 0.05),
           "cindex_cv_base": float(np.mean(cs["base"])), "cindex_cv_with_score": float(np.mean(cs["with_score"])),
           "label": "supporting analysis only; real patients; prognostic association, not treatment benefit",
           "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT / "upenn_survival.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
