#!/usr/bin/env python3
"""Script 128: corpus-callosum subgroup, TUM development only, exploratory
(Amendment 7 A7.3 as implemented by Amendment 10 A10.5).

Proxy rule (no CC atlas label locally; labelled PROXY everywhere): a patient
is "callosum-adjacent" if the pre-op core has >= 1 voxel with
|x - midline| <= 5 (1 mm voxels) AND wm_pbmap > 0.5, where midline is the
middle x-index of the brain mask.
Contrast: per-patient M2(r1|lam2) minus M1(iso|lam2) enhancing coverage from
output/predictgbm/cache_112/*.json (selected cells, ledger 2al), by subgroup.
M3a-M3i subgroup is PENDING (no per-patient M3 file saved; A10.5).
Output: output/predictgbm/callosum_subgroup.json, callosum_subgroup.csv
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import predictgbm_io as pg  # noqa: E402
from src import patient_stats as ps  # noqa: E402

C1 = "M1_iso|lam=2"
C2 = "M2_aniso_r1|lam=2"


def flag_callosum(pid: str) -> dict:
    c = pg.load_case(pid, with_recurrence=False)
    core = pg.core(c["seg"])
    nx = c["brain"].shape[0]
    mid = nx // 2
    xs = np.arange(nx)[:, None, None]
    near_mid = np.abs(xs - mid) <= 5
    hit = core & near_mid & (c["wm"].astype(np.float32) > 0.5)
    return {"pid": pid, "callosum_adjacent": bool(hit.any()),
            "n_hit_vox": int(hit.sum()), "n_core_vox": int(core.sum())}


def main() -> None:
    rows = []
    for f in sorted(glob.glob(str(PROJECT_ROOT / "output" / "predictgbm"
                                          / "cache_112" / "*.json"))):
        d = json.load(open(f))
        if d.get("skip") or C1 not in d.get("runs", {}) or C2 not in d.get("runs", {}):
            continue
        pid = d["pid"]
        fl = flag_callosum(pid)
        fl["m1_enh"] = float(d["runs"][C1]["enh"])
        fl["m2_enh"] = float(d["runs"][C2]["enh"])
        fl["m2_minus_m1"] = fl["m2_enh"] - fl["m1_enh"]
        fl["n_rec_enh"] = int(d["n_rec_enh"])
        rows.append(fl)
    df = pd.DataFrame(rows)
    df.to_csv(PROJECT_ROOT / "output" / "predictgbm" / "callosum_subgroup.csv",
              index=False)
    out = {"script": "128_callosum_subgroup", "proxy_rule":
           "core voxel with |x-midline|<=5 AND wm_pbmap>0.5 (PROXY, no atlas label)",
           "cells": {"M1": C1, "M2": C2}, "subgroups": {}}
    for lab, sub in (("adjacent", df[df["callosum_adjacent"]]),
                     ("rest", df[~df["callosum_adjacent"]])):
        x = sub[sub["n_rec_enh"] > 0]["m2_minus_m1"].to_numpy()
        s = ps.summarize_delta(x)
        out["subgroups"][lab] = {"n_cache": int(len(sub)),
                                 "n_with_recurrence": int(len(x)), **s}
    out["n_adjacent"] = int(df["callosum_adjacent"].sum())
    out["n_total"] = int(len(df))
    (PROJECT_ROOT / "output" / "predictgbm" / "callosum_subgroup.json"
     ).write_text(json.dumps(out, indent=1))
    print(json.dumps({k: {"n": v["n_with_recurrence"],
                          "mean": v["mean"], "ci": v["mean_ci95"]}
                      for k, v in out["subgroups"].items()}, indent=1))


if __name__ == "__main__":
    main()
