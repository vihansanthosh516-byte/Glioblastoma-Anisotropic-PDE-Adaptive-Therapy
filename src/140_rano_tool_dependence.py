#!/usr/bin/env python3
"""Script 140: Amendment 15 tests 1 and 2. EXPLORATORY. Declared in analysis_plan_v1.md (commit 51fd053) before this file was written.

Operational details fixed here, before the first run (they only make A15.2-A15.4 executable):
  * Volumes: contrast-enhancing label, CT1, original_shape_VoxelVolume (mm3) from the two LUMIERE pyradiomics CSVs, exactly as script 104.
  * Post-operative scan = scan whose time point is after the patient's scan rated "Pre-Op" in LUMIERE-ExpertRating.
    Patients without a "Pre-Op" rating are excluded and counted.
  * Pair = two consecutive post-operative scans of a patient in the merged table of both tools (a scan with a missing tool
    row stays in the sequence). Eligible if all four volumes (2 tools x 2 scans) are >= 500 mm3.
  * Test 1 (A15.3): class by r = V1/V0 inside each tool: progression if r >= 1.40, partial response if r <= 0.35, else stable.
    flip = classes differ. Patient-cluster bootstrap, 10,000 draws, seed 20261018. Reading: lower CI >= 0.20 material;
    upper CI <= 0.10 robust; else inconclusive.
  * Test 2 (A15.4): d = dlnV(HD) - dlnV(DBT). Size band = tertile (edges from eligible pairs) of the smaller of the two scans'
    tool-geometric-mean volumes. Upper noise bound s_up = 1.4826 MAD(d) / sqrt(2). Lower noise bound from UPenn expert correction
    (region enhancing, rows without skip): median |ln(V_ai/V_rad)| per band (band by geometric mean of V_ai and V_rad, same
    edges) / 0.6745 = per-scan SD, times sqrt(2) = change SD. MDC95 = 1.96 s. Line = ln(1.40) = 0.3365.
    Reading: bracket = [min, max] of the two MDC95; max < line -> above the noise; min > line -> inside the noise; else undetermined.
Run once. No re-run with other thresholds, floors or definitions.
Output: output/rano_noise/rano_tool_dependence.json, rano_pairs.csv
"""
from __future__ import annotations

import json
import re
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
LUM = PROJECT_ROOT / "data" / "external" / "lumiere"
OUT = PROJECT_ROOT / "output" / "rano_noise"
UPENN = PROJECT_ROOT / "output" / "upenn_rater_noise_scans.csv"
LABEL = "Contrast-enhancing"
FLOOR = 500.0
PD_LINE, PR_LINE = 1.40, 0.35
LINE_LN = float(np.log(PD_LINE))
SEED, NBOOT = 20261018, 10000


def tp_key(tp):
    m = re.match(r"week-(\d+)(?:-(\d+))?$", tp)
    return int(m.group(1)), int(m.group(2) or 0)


def load(name):
    p = pd.read_csv(LUM / f"LUMIERE-pyradiomics-{name}-features.csv",
                    usecols=["Patient", "Time point", "Label name", "Sequence", "original_shape_VoxelVolume"])
    p = p[(p["Label name"] == LABEL) & (p["Sequence"] == "CT1")]
    return p.rename(columns={"Time point": "tp", "original_shape_VoxelVolume": "vol"})[["Patient", "tp", "vol"]]


def cls(r):
    return np.where(r >= PD_LINE, "PD", np.where(r <= PR_LINE, "PR", "SD"))


def boot_ci(df, fn, rng):
    pats = df["Patient"].unique()
    groups = {p: g for p, g in df.groupby("Patient")}
    vals = []
    for _ in range(NBOOT):
        pick = rng.choice(pats, len(pats), replace=True)
        vals.append(fn(pd.concat([groups[p] for p in pick], ignore_index=True)))
    v = np.array(vals, float)
    v = v[~np.isnan(v)]
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def main():
    hd, dbt = load("hdglioauto"), load("deepbratumia")
    m = hd.merge(dbt, on=["Patient", "tp"], how="outer", suffixes=("_hd", "_dbt"))
    rat = pd.read_csv(LUM / "LUMIERE-ExpertRating-v202211.csv")
    rat.columns = ["Patient", "tp", "lt3", "nonmeas", "rating", "rationale"]
    rat["rating"] = rat["rating"].astype(str).str.strip()
    pre = rat[rat["rating"] == "Pre-Op"].set_index("Patient")["tp"].map(tp_key).to_dict()
    m["k"] = m["tp"].map(tp_key)
    n_patients_all = m["Patient"].nunique()
    m = m[m["Patient"].isin(pre)]
    n_no_preop = n_patients_all - m["Patient"].nunique()
    m = m[[k > pre[p] for p, k in zip(m["Patient"], m["k"])]]
    rating_of = {(p, t): r for p, t, r in zip(rat["Patient"], rat["tp"], rat["rating"])}
    rows, excl = [], {"tool_missing_or_zero_or_below_500mm3": 0}
    n_consec = 0
    for pid, g in m.sort_values(["Patient", "k"]).groupby("Patient"):
        g = g.reset_index(drop=True)
        for i in range(1, len(g)):
            n_consec += 1
            v = [g.loc[i - 1, "vol_hd"], g.loc[i - 1, "vol_dbt"], g.loc[i, "vol_hd"], g.loc[i, "vol_dbt"]]
            if any((pd.isna(x) or x < FLOOR) for x in v):
                excl["tool_missing_or_zero_or_below_500mm3"] += 1
                continue
            rows.append({"Patient": pid, "tp0": g.loc[i - 1, "tp"], "tp1": g.loc[i, "tp"],
                         "hd0": v[0], "dbt0": v[1], "hd1": v[2], "dbt1": v[3],
                         "expert_later": rating_of.get((pid, g.loc[i, "tp"]), "")})
    df = pd.DataFrame(rows)
    res = {"script": "140_rano_tool_dependence", "amendment": "A15", "n_patients_in_lumiere_tables": int(n_patients_all),
           "n_patients_excluded_no_preop_rating": int(n_no_preop), "n_consecutive_postop_pairs": int(n_consec),
           "n_pairs_excluded": excl, "n_pairs_eligible": int(len(df)), "n_patients_eligible": int(df["Patient"].nunique()),
           "share_excluded_by_floor_or_missing": float(excl["tool_missing_or_zero_or_below_500mm3"] / max(n_consec, 1))}
    if len(df) == 0:
        (OUT).mkdir(parents=True, exist_ok=True)
        (OUT / "rano_tool_dependence.json").write_text(json.dumps(res, indent=1))
        print(json.dumps(res, indent=1))
        return
    df["r_hd"], df["r_dbt"] = df["hd1"] / df["hd0"], df["dbt1"] / df["dbt0"]
    df["cls_hd"], df["cls_dbt"] = cls(df["r_hd"]), cls(df["r_dbt"])
    df["flip"] = (df["cls_hd"] != df["cls_dbt"]).astype(int)
    df["one_pd_only"] = ((df["cls_hd"] == "PD") ^ (df["cls_dbt"] == "PD")).astype(int)
    rng = np.random.default_rng(SEED)

    res["TEST1_flip_rate"] = {"rate": float(df["flip"].mean()), "ci95": boot_ci(df, lambda d: d["flip"].mean(), rng)}
    lo, hi = res["TEST1_flip_rate"]["ci95"]
    res["TEST1_reading"] = ("material tool dependence" if lo >= 0.20 else "robust" if hi <= 0.10 else "inconclusive")
    res["TEST1a_one_tool_PD_other_not"] = {"rate": float(df["one_pd_only"].mean()),
                                           "ci95": boot_ci(df, lambda d: d["one_pd_only"].mean(), rng)}
    res["class_counts"] = {"hd": df["cls_hd"].value_counts().to_dict(), "dbt": df["cls_dbt"].value_counts().to_dict()}
    res["class_crosstab_hd_by_dbt"] = pd.crosstab(df["cls_hd"], df["cls_dbt"]).to_dict()

    # 1b: expert PD vs not PD at the later scan (only pairs whose later scan has a PD/SD/PR/CR rating)
    e = df[df["expert_later"].isin(["PD", "SD", "PR", "CR"])].copy()
    e["exp_pd"] = (e["expert_later"] == "PD").astype(int)
    b1 = {"n_pairs": int(len(e)), "n_patients": int(e["Patient"].nunique()), "n_expert_pd": int(e["exp_pd"].sum())}
    for tool in ("hd", "dbt"):
        e["call"] = (e[f"cls_{tool}"] == "PD").astype(int)
        sens = lambda d: float(d.loc[d["exp_pd"] == 1, "call"].mean()) if (d["exp_pd"] == 1).any() else np.nan
        spec = lambda d: float(1 - d.loc[d["exp_pd"] == 0, "call"].mean()) if (d["exp_pd"] == 0).any() else np.nan
        b1[tool] = {"sensitivity": sens(e), "sens_ci95": boot_ci(e, sens, rng), "specificity": spec(e), "spec_ci95": boot_ci(e, spec, rng)}
    res["TEST1b_vs_expert_rating"] = b1

    # TEST 2
    df["d"] = np.log(df["hd1"] / df["hd0"]) - np.log(df["dbt1"] / df["dbt0"])
    g0, g1 = np.sqrt(df["hd0"] * df["dbt0"]), np.sqrt(df["hd1"] * df["dbt1"])
    df["vsize"] = np.minimum(g0, g1)
    edges = list(df["vsize"].quantile([1 / 3, 2 / 3]))
    df["band"] = np.digitize(df["vsize"], edges)
    up = pd.read_csv(UPENN)
    up = up[(up["region"] == "enhancing") & (up["skip"].isna() | (up["skip"].astype(str) == ""))].copy()
    up = up[(up["v_ai_mm3"] > 0) & (up["v_rad_mm3"] > 0)]
    up["vg"] = np.sqrt(up["v_ai_mm3"] * up["v_rad_mm3"])
    up["band"] = np.digitize(up["vg"], edges)
    t2 = {"band_edges_mm3": edges, "line_ln_1.40": LINE_LN, "bands": {}}
    for b in (0, 1, 2):
        d = df.loc[df["band"] == b, "d"]
        u = up.loc[up["band"] == b, "abs_ln_ratio"]
        mad = float(np.median(np.abs(d - d.median())))
        s_up = 1.4826 * mad / np.sqrt(2)
        s_low = float(np.sqrt(2) * (u.median() / 0.6745)) if len(u) else float("nan")
        mdc_up, mdc_low = 1.96 * s_up, 1.96 * s_low
        lo_b, hi_b = min(mdc_up, mdc_low), max(mdc_up, mdc_low)
        reading = ("above the noise" if hi_b < LINE_LN else "inside the noise" if lo_b > LINE_LN else "undetermined")
        t2["bands"][str(b)] = {"n_pairs": int(len(d)), "n_patients": int(df.loc[df["band"] == b, "Patient"].nunique()),
                               "median_vsize_mm3": float(df.loc[df["band"] == b, "vsize"].median()),
                               "s_change_upper_bound": float(s_up), "mdc95_upper_bound_ln": float(mdc_up),
                               "mdc95_upper_bound_pct": float(np.expm1(mdc_up) * 100),
                               "upenn_n_scans": int(len(u)), "s_change_lower_bound": s_low, "mdc95_lower_bound_ln": float(mdc_low),
                               "mdc95_lower_bound_pct": float(np.expm1(mdc_low) * 100),
                               "bracket_ln": [float(lo_b), float(hi_b)], "reading": reading}
    res["TEST2_mdc_vs_progression_line"] = t2
    res["weakest_points"] = [
        "Two automated tools with different label definitions are not a scan-rescan pair; the upper bound includes definition differences that change over time.",
        "UPenn expert correction anchors on the AI mask, so the lower bound understates independent noise.",
        "LUMIERE is one centre; expert ratings are 2D RANO at selected follow-ups and are relative to baseline or nadir, not to the previous scan.",
        "Exploratory on public data; aggregate noise summaries (script 104) were seen before the declaration."]
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "rano_pairs.csv", index=False)
    (OUT / "rano_tool_dependence.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
