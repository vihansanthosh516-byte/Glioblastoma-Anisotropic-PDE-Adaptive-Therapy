#!/usr/bin/env python3
"""Script 130: H-3v frozen LUMIERE test, run ONCE after the freeze
(configs/h3v_frozen.yaml + git tag h3v-frozen). Amendment 10 A10.4.

Reads the frozen MixedLM params, applies them to LUMIERE enhancing volumes
(HD-GLIO-AUTO primary, DeepBraTumIA sensitivity; script-91 loaders and
week->day mapping reused unchanged), same cohort gate as script 91
(>=3 follow-up scans, first three volumes > 0; forecast scan idx 2 from 0,1).
Endpoints: mean |ln| error MixedLM vs persistence, 80/90% coverage with the
frozen quantiles, RANO-PD AUC using predicted ln ratio, abstain subgroup.
Output: output/h3v_test.json, output/h3v_test_pairs.csv
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path as _Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import wilcoxon

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
FROZEN = PROJECT_ROOT / "configs" / "h3v_frozen.yaml"
OUT_JSON = PROJECT_ROOT / "output" / "h3v_test.json"
OUT_CSV = PROJECT_ROOT / "output" / "h3v_test_pairs.csv"  # only a name stem: the CSVs written are h3v_test_pairs_<segmenter>.csv
# NOTE 2026-10-19 audit: is_frozen() checks that the config and tag exist; it does not verify the config equals the tagged
# version or enforce runs_allowed. Missing daysRT/MGMT are imputed to the frozen median / unmethylated + indicator.
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from src import patient_stats as ps  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "s91", PROJECT_ROOT / "src" / "91_lumiere_volume_forecast.py")
s91 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s91)


def is_frozen() -> bool:
    if not FROZEN.exists():
        return False
    r = subprocess.run(["git", "tag", "--list", "h3v-frozen"], cwd=PROJECT_ROOT,
                       capture_output=True, text=True)
    return r.stdout.strip() == "h3v-frozen"


def floor_for_vin(v_mm3: float, edges, floors) -> float:
    if v_mm3 < edges[0]:
        return floors[0]
    if v_mm3 < edges[1]:
        return floors[1]
    return floors[2]


def load_demo() -> pd.DataFrame:
    d = pd.read_csv(PROJECT_ROOT / "data" / "external" / "lumiere" /
                    "LUMIERE-Demographics_Pathology.csv")
    d = d.rename(columns={"Patient": "Patient", "Sex": "sex",
                          "Age at surgery (years)": "age",
                          "MGMT qualitative": "mgmt"})
    d["sex"] = d["sex"].astype(str).str.strip().str.lower()
    d["mgmt"] = d["mgmt"].astype(str).str.strip().str.lower()
    return d[["Patient", "sex", "age", "mgmt"]]


def run_segmenter(name: str, cfg: dict, demo: pd.DataFrame) -> dict:
    ser = s91.follow_up_series(s91.load_volumes(name), s91.load_ratings())
    fe = cfg["fixed_effects"]
    beta = np.array([cfg["params"][k] for k in fe])
    rows, skipped = [], 0
    for pid, s in ser.items():
        if len(s) < 3 or (s["vol"].iloc[:3] <= 0).any():
            continue
        t0, t1, t2 = s["day"].iloc[0], s["day"].iloc[1], s["day"].iloc[2]
        v0, v1, v2 = s["vol"].iloc[0], s["vol"].iloc[1], s["vol"].iloc[2]
        r2 = s["rating"].iloc[2]
        if min(v0, v1, v2) <= 0 or t1 <= t0 or t2 <= t1:
            skipped += 1
            continue
        m = demo[demo["Patient"] == pid]
        age = float(m["age"].iloc[0]) if len(m) and pd.notna(m["age"].iloc[0]) else np.nan
        sex = m["sex"].iloc[0] if len(m) else ""
        mg = m["mgmt"].iloc[0] if len(m) else ""
        if pd.isna(age):
            skipped += 1
            continue
        x = {"const": 1.0, "logVin": float(np.log(v1)), "dt_days": float(t2 - t1),
             "daysRT": float(cfg["rt_median_days"]), "rt_known": 0.0,
             "age": age, "sexM": 1.0 if sex == "male" else 0.0,
             "mgmt_meth": 1.0 if mg == "methylated" else 0.0,
             "mgmt_miss": 0.0 if mg in ("methylated", "not methylated") else 1.0}
        if x["mgmt_miss"]:
            x["mgmt_meth"] = 0.0
        yhat = float(np.dot(beta, [x[k] for k in fe]))
        y = float(np.log(v2 / v1))
        fl = floor_for_vin(float(v1), cfg["floor_edges_mm3"], cfg["floors"])
        rows.append({"pid": pid, "v0": float(v0), "v1": float(v1), "v2": float(v2),
                     "dt": float(t2 - t1), "y": y, "yhat": yhat,
                     "ae_mix": abs(y - yhat), "ae_persist": abs(y),
                     "rating_target": str(r2), "pd": int(str(r2) == "PD"),
                     "abstain": bool(abs(yhat) < fl)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV.with_name(f"h3v_test_pairs_{name}.csv"), index=False)
    d = (df["ae_mix"] - df["ae_persist"]).to_numpy()
    s = ps.summarize_delta(d)
    s["wilcoxon_p_two_sided"] = float(wilcoxon(
        df["ae_mix"], df["ae_persist"]).pvalue) if np.any(d != 0) else 1.0
    from sklearn.metrics import roc_auc_score
    try:
        auc = float(roc_auc_score(df["pd"], df["yhat"])) if df["pd"].nunique() > 1 else None
    except Exception:  # noqa: BLE001
        auc = None
    return {"n": int(len(df)), "skipped": int(skipped),
            "mean_ae_mix": float(df["ae_mix"].mean()),
            "mean_ae_persist": float(df["ae_persist"].mean()),
            "mix_minus_persist": s,
            "coverage_80": float((((df["y"] - df["yhat"]).abs() <= cfg["q80"])).mean()),
            "coverage_90": float((((df["y"] - df["yhat"]).abs() <= cfg["q90"])).mean()),
            "rano_pd_auc": auc, "pd_share": float(df["pd"].mean()),
            "abstain_share": float(df["abstain"].mean())}


def main() -> None:
    if not is_frozen():
        raise SystemExit("not frozen: configs/h3v_frozen.yaml + tag h3v-frozen required")
    if OUT_JSON.exists() and not os.environ.get("ALLOW_SPENT_RESCORE"):
        raise SystemExit("h3v_test.json already exists: the LUMIERE set is spent; set ALLOW_SPENT_RESCORE=1 only to overwrite deliberately")
    cfg = yaml.safe_load(FROZEN.read_text())
    demo = load_demo()
    res = {"script": "130_h3v_frozen_lumiere",
           "frozen_sha256": hashlib.sha256(FROZEN.read_bytes()).hexdigest(),
           "segmenters": {}}
    for name in ("hdglioauto", "deepbratumia"):
        res["segmenters"][name] = run_segmenter(name, cfg, demo)
        print(name, res["segmenters"][name]["n"],
              "mix %.3f vs persist %.3f" % (
                  res["segmenters"][name]["mean_ae_mix"],
                  res["segmenters"][name]["mean_ae_persist"]),
              "cov80 %.2f cov90 %.2f auc %s" % (
                  res["segmenters"][name]["coverage_80"],
                  res["segmenters"][name]["coverage_90"],
                  res["segmenters"][name]["rano_pd_auc"]))
    OUT_JSON.write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
