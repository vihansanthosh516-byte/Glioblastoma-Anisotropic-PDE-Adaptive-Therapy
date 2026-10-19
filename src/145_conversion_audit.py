#!/usr/bin/env python3
"""Script 145: Amendment 21, audit of the A19 conversion from published repeatability to change-SD. SIMULATION.
Declared in analysis_plan_v1.md (commit 2d3e6a4) before this file was written. A19 stands as run.

Published values (Abu Khalaf 2021, PMC8191681), enhancing tumour:
  BraTumIA 1.2: RC% 46 (33-67), RC 6.9 mL (4.9-10), mean volumes 24.2 and 26.3 mL
  BraTumIA 2.0: RC% 39 (28-57), RC 5.2 mL (3.7-7.5), mean volumes 22.9 and 24.9 mL
Readings: R1 s = RC%/100/1.96 ; R2 s = ln(1+RC%/100)/1.96 ; R3 s = (RC_mL/mean_volume_mL)/1.96.
s is the SD of the CHANGE of ln volume between two scans (per-scan SD s/sqrt(2) in the simulator, as in A16).
Same model, seed, draw order, estimands and materiality rules as A16/A18/A19 (module 142 supplies shared functions).
Run once. No re-run with other values.
Output: output/conversion_audit/conversion_audit_cells.csv, conversion_audit_summary.json
"""
from __future__ import annotations

import importlib.util
import json
import math
import time
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "conversion_audit"


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, PROJECT_ROOT / "src" / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


M = _load("sens142", "142_rano_sensitivity.py")
T = M.T

PUB = {
    "BraTumIA_1.2": dict(rcpct=(46, 33, 67), rcml=(6.9, 4.9, 10.0), mean_ml=(24.2 + 26.3) / 2),
    "BraTumIA_2.0": dict(rcpct=(39, 28, 57), rcml=(5.2, 3.7, 7.5), mean_ml=(22.9 + 24.9) / 2),
}
GRID = [0.08, 0.111, 0.126, 0.139, 0.168, 0.193, 0.262]
NS, DS = (50, 100), (8, 12)


def conversions():
    rows = []
    for ver, p in PUB.items():
        for k, lab in enumerate(("point", "ci_low", "ci_high")):
            rc = p["rcpct"][k] / 100.0
            rows.append(dict(version=ver, which=lab, R1_linear=rc / 1.96, R2_log=math.log1p(rc) / 1.96,
                             R3_absolute=(p["rcml"][k] / p["mean_ml"]) / 1.96))
    return pd.DataFrame(rows)


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    conv = conversions()
    conv.to_csv(OUT / "conversion_table.csv", index=False)
    rows = []
    for D in DS:
        for N in NS:
            M.N, M.D = N, D
            times = np.arange(0, M.HORIZON + 1, D, dtype=float)
            J = len(times) - 1
            W = M.make_world(times, J)
            for rule in T.CONFIRMS:
                base = M.evaluate(W, J, times, 0, None, rule, False)
                for s in GRID:
                    r = M.summarise(M.evaluate(W, J, times, s, None, rule, False), base)
                    rows.append(dict(N=N, D=D, rule=rule, s=s, **r, E1_material=M.mat_e1(r), E2_material=M.mat_e2(r)))
    cells = pd.DataFrame(rows)
    cells.to_csv(OUT / "conversion_audit_cells.csv", index=False)

    c = cells[(cells.N == 100) & (cells.D == 8) & (cells.rule == "none")].set_index("s")
    big = [s for s in GRID if s >= 0.111]
    pu1 = bool(all(c.loc[s, "E1_material"] for s in big))
    pu2 = bool(all(c.loc[s, "E2_material"] for s in big))
    pu3 = bool(not c.loc[0.08, "E1_material"])
    a19 = pd.read_csv(PROJECT_ROOT / "output" / "anchored_sim" / "anchored_sim_cells.csv")
    a19 = a19[(a19.N == 100) & (a19.D == 8) & (a19.rule == "none")].set_index("s")
    # PU4: R1 central vs R3 central per version, single-scan E1
    pu4_detail = {}
    for ver, r1_lab, r3_s in (("BraTumIA_1.2", 0.235, 0.139), ("BraTumIA_2.0", 0.199, 0.111)):
        e_r1 = float(a19.loc[r1_lab, "E1"]); e_r3 = float(c.loc[r3_s, "E1"])
        pu4_detail[ver] = dict(E1_R1=e_r1, E1_R3=e_r3, abs_diff=abs(e_r1 - e_r3), more_than_4=bool(abs(e_r1 - e_r3) > 4.0))
    pu4 = bool(all(v["more_than_4"] for v in pu4_detail.values()))

    def rng_for(col):
        out = {}
        for ver in PUB:
            v = conv[conv.version == ver][col]
            lo, hi = float(v.min()), float(v.max())
            out[ver] = dict(s_point=float(conv[(conv.version == ver) & (conv.which == "point")][col].iloc[0]), s_range=[lo, hi])
        return out

    summary = dict(
        label="simulation; audit of an assumed conversion from one published repeatability study; independent normal error; volumetric component only; not evidence about any real trial",
        amendment="A21 in analysis_plan_v1.md", seed=M.SEED, R=M.R, runtime_s=round(time.time() - t0, 1),
        conversion_table=conv.round(4).to_dict("records"),
        readings={k: rng_for(k) for k in ("R1_linear", "R2_log", "R3_absolute")},
        union_range=[float(conv[["R1_linear", "R2_log", "R3_absolute"]].min().min()), float(conv[["R1_linear", "R2_log", "R3_absolute"]].max().max())],
        N100_D8_single_scan={str(s): dict(E1=float(c.loc[s, "E1"]), E1_lo=float(c.loc[s, "E1_lo"]), E1_hi=float(c.loc[s, "E1_hi"]),
                                           E2=float(c.loc[s, "E2"]), E2_hi=float(c.loc[s, "E2_hi"]),
                                           E1_material=bool(c.loc[s, "E1_material"]), E2_material=bool(c.loc[s, "E2_material"])) for s in GRID},
        N100_D8_two_scan={str(r.s): dict(E1=float(r.E1), E2=float(r.E2), E1_material=bool(r.E1_material), E2_material=bool(r.E2_material))
                          for r in cells[(cells.N == 100) & (cells.D == 8) & (cells.rule == "two-scan")].itertuples()},
        PU1_supported=pu1, PU2_supported=pu2, PU3_supported=pu3, PU4_supported=pu4, PU4_detail=pu4_detail)
    (OUT / "conversion_audit_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(conv.round(4).to_string())
    print(json.dumps({k: summary[k] for k in ("union_range", "N100_D8_single_scan", "PU1_supported", "PU2_supported", "PU3_supported", "PU4_supported", "PU4_detail", "runtime_s")}, indent=2, default=float))


if __name__ == "__main__":
    main()
