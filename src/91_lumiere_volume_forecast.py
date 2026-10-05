"""Script 91: LUMIERE as new real patients: growth rates and a volume forecast (limits L1, L3).

LUMIERE (91 GBM patients, 638 study dates) gives tumour volumes but not masks here (masks are in the 30 GB zip,
not downloaded). So this is a VOLUME forecast, not the Dice forecast of scripts 79/81.

PRE-SPECIFIED DESIGN (fixed before any LUMIERE result was seen)
- Volume: contrast-enhancing label, HD-GLIO-AUTO, original_shape_VoxelVolume (mm3); one row per scan (CT1 row).
  Robustness: the same with DeepBraTumIA.
- Follow-up scans: study dates whose expert RANO rating is PD, SD, PR or CR (Pre-Op and Post-Op dates are dropped)
  and whose volume is not missing.
  Time: week-NNN is day 7*NNN; a "-k" suffix only orders studies inside the same week.
- Cohort: patients with >= 3 follow-up scans whose first three volumes are all > 0 (same gate as script 79).
- Forecast (primary): follow-up scan index 1 -> 2, from the volumes at index 0 and 1 (the script 79 design).
    no_change   V2_hat = V1
    exp_growth  V2_hat = V1 * exp(rho * (t2 - t1)),  rho = ln(V1/V0)/(t1 - t0)   (a growth-only model)
  Error: |ln(V_hat / V2)|. Test: paired two-sided Wilcoxon on the per-patient errors, exp_growth vs no_change.
- Secondary (reported in full): every rolling forecast k >= 2 with a patient-cluster bootstrap CI; strata by whether
  the tumour grew or shrank between scans k-1 and k (descriptive only, outcome-selected); DeepBraTumIA repeat;
  sensitivity with volumes >= 500 mm3 only.
- Growth rates: per-patient log-linear fit over all follow-up scans with volume > 0 (>= 2 scans), per day.
  Reported: median, range, share above 0.02 /day (the script 85 threshold), and the MU-Glioma values for comparison.

Input:  data/external/lumiere/LUMIERE-pyradiomics-{hdglioauto,deepbratumia}-features.csv, LUMIERE-ExpertRating-v202211.csv
Output: output/lumiere_volume_forecast.json
"""
import json
import re
from pathlib import Path as _Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
DATA = PROJECT_ROOT / "data" / "external" / "lumiere"
OUT_JSON = PROJECT_ROOT / "output" / "lumiere_volume_forecast.json"
KEEP = {"PD", "SD", "PR", "CR"}
RNG = np.random.default_rng(0)


def tp_key(tp):
    m = re.match(r"week-(\d+)(?:-(\d+))?$", tp)
    return int(m.group(1)), int(m.group(2) or 0)


def load_volumes(name, label="Contrast-enhancing"):
    p = pd.read_csv(DATA / f"LUMIERE-pyradiomics-{name}-features.csv",
                    usecols=["Patient", "Time point", "Label name", "Sequence", "original_shape_VoxelVolume"])
    p = p[(p["Label name"] == label) & (p["Sequence"] == "CT1")]
    return p.rename(columns={"Time point": "tp", "original_shape_VoxelVolume": "vol"})[["Patient", "tp", "vol"]]


def load_ratings():
    r = pd.read_csv(DATA / "LUMIERE-ExpertRating-v202211.csv")
    r.columns = ["Patient", "tp", "lt3", "nonmeas", "rating", "rationale"]
    r["rating"] = r["rating"].astype(str).str.strip()
    return r[["Patient", "tp", "rating"]]


def follow_up_series(vol, rat):
    d = vol.merge(rat, on=["Patient", "tp"], how="inner")
    d = d[d["rating"].isin(KEEP) & d["vol"].notna()].copy()
    d["wk"], d["sub"] = zip(*d["tp"].map(tp_key))
    d["day"] = 7.0 * d["wk"]
    d = d.sort_values(["Patient", "wk", "sub"])
    return {p: g[["day", "vol", "rating"]].reset_index(drop=True) for p, g in d.groupby("Patient")}


def forecast_pair(series, k):
    """Forecast scan k from scans k-1 and k-2. Returns (err_exp, err_nochange, grew)."""
    t0, t1, t2 = series["day"][k - 2], series["day"][k - 1], series["day"][k]
    v0, v1, v2 = series["vol"][k - 2], series["vol"][k - 1], series["vol"][k]
    if min(v0, v1, v2) <= 0 or t1 <= t0 or t2 <= t1:
        return None
    rho = np.log(v1 / v0) / (t1 - t0)
    vhat = v1 * np.exp(rho * (t2 - t1))
    return abs(np.log(vhat / v2)), abs(np.log(v1 / v2)), bool(v2 > v1)


def summarize_errors(rows):
    e_exp = np.array([r[0] for r in rows])
    e_nc = np.array([r[1] for r in rows])
    d = e_nc - e_exp  # positive: growth model better
    try:
        p = float(wilcoxon(e_exp, e_nc).pvalue)
    except ValueError:
        p = None
    return {"n": len(rows), "mean_abs_logerr_exp": float(e_exp.mean()), "mean_abs_logerr_nochange": float(e_nc.mean()),
            "median_abs_logerr_exp": float(np.median(e_exp)), "median_abs_logerr_nochange": float(np.median(e_nc)),
            "mean_gain_nochange_minus_exp": float(d.mean()), "share_exp_better": float((d > 0).mean()),
            "wilcoxon_p": p}


def primary(series_by_patient, min_vol=0.0):
    rows = []
    for s in series_by_patient.values():
        if len(s) < 3 or (s["vol"].iloc[:3] <= min_vol).any():
            continue
        r = forecast_pair(s, 2)
        if r:
            rows.append(r)
    return rows


def rolling(series_by_patient, min_vol=0.0, n_boot=2000):
    per_patient = []
    for s in series_by_patient.values():
        rs = []
        for k in range(2, len(s)):
            if min(s["vol"].iloc[k - 2:k + 1]) <= min_vol:
                continue
            r = forecast_pair(s, k)
            if r:
                rs.append(r)
        if rs:
            per_patient.append(rs)
    allr = [r for rs in per_patient for r in rs]
    out = summarize_errors(allr)
    out["n_patients"] = len(per_patient)
    diffs = np.array([r[1] - r[0] for r in allr])
    owner = np.concatenate([[i] * len(rs) for i, rs in enumerate(per_patient)])
    means = []
    for _ in range(n_boot):
        pick = RNG.integers(0, len(per_patient), len(per_patient))
        sel = np.concatenate([diffs[owner == i] for i in pick])
        means.append(sel.mean())
    out["gain_ci95_cluster_bootstrap"] = [float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))]
    for name, flag in (("grew", True), ("shrank_or_same", False)):
        sub = [r for r in allr if r[2] == flag]
        if len(sub) > 1:
            out[f"stratum_{name}"] = summarize_errors(sub)
    return out


def growth_rates(series_by_patient):
    rhos = []
    for s in series_by_patient.values():
        s = s[s["vol"] > 0]
        if len(s) >= 2 and s["day"].nunique() >= 2:
            rhos.append(float(np.polyfit(s["day"], np.log(s["vol"]), 1)[0]))
    a = np.array(rhos)
    return {"n": len(a), "median_rho_per_day": float(np.median(a)), "min": float(a.min()), "max": float(a.max()),
            "q25": float(np.percentile(a, 25)), "q75": float(np.percentile(a, 75)),
            "n_above_0.02": int((a > 0.02).sum()), "share_positive": float((a > 0).mean())}


def main():
    rat = load_ratings()
    res = {"script": "91_lumiere_volume_forecast", "mu_glioma_reference": {"median_rho_per_day_real_test": 0.001674},
           "segmenters": {}}
    for name in ("hdglioauto", "deepbratumia"):
        ser = follow_up_series(load_volumes(name), rat)
        lens = [len(s) for s in ser.values()]
        out = {"n_patients_with_followup": len(ser), "n_ge3_followup": int(sum(l >= 3 for l in lens)),
               "n_followup_scans": int(sum(lens)), "growth_rates": growth_rates(ser)}
        prim = primary(ser)
        out["primary_idx1_to_2"] = summarize_errors(prim)
        out["primary_ge500mm3"] = summarize_errors(primary(ser, 500.0)) if len(primary(ser, 500.0)) > 5 else None
        out["rolling_all_k"] = rolling(ser)
        res["segmenters"][name] = out
        print(name, "ge3:", out["n_ge3_followup"], "primary n", out["primary_idx1_to_2"]["n"],
              "err exp %.3f vs nochange %.3f p=%s" % (out["primary_idx1_to_2"]["mean_abs_logerr_exp"],
                                                      out["primary_idx1_to_2"]["mean_abs_logerr_nochange"],
                                                      out["primary_idx1_to_2"]["wilcoxon_p"]),
              "| median rho %.4f, >0.02: %d" % (out["growth_rates"]["median_rho_per_day"], out["growth_rates"]["n_above_0.02"]))
    OUT_JSON.write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
