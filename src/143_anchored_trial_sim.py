#!/usr/bin/env python3
"""Script 143: Amendment 19, trial simulation anchored to a published scan-rescan repeatability.
SIMULATION ONLY. Declared in analysis_plan_v1.md (commit d9b924b) before this file was written.

Anchor: Abu Khalaf et al., AJNR 2021 (PMC8191681), enhancing-tumour RC% 46% (33-67) and 39% (28-57), large recurrent
tumours, BraTumIA, 20 patients, 2 days apart. Assumed conversion s = RC% / 1.96. Anchors s in {0.143, 0.199, 0.235, 0.342}.
Same model, seed, draw order, estimands and materiality rules as A16 and A18 (module 142 supplies the shared functions).
Run once. No re-run with other anchors or conversions.
Output: output/anchored_sim/anchored_sim_cells.csv, anchored_sim_summary.json
"""
from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "anchored_sim"


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, PROJECT_ROOT / "src" / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


M = _load("sens142", "142_rano_sensitivity.py")
T = M.T

ANCHORS = {"0.143": 0.143, "0.199": 0.199, "0.235": 0.235, "0.342": 0.342}
NS, DS = (50, 100), (8, 12)


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for D in DS:
        for N in NS:
            M.N, M.D = N, D
            times = np.arange(0, M.HORIZON + 1, D, dtype=float)
            J = len(times) - 1
            W = M.make_world(times, J)
            for rule in T.CONFIRMS:
                base = M.evaluate(W, J, times, 0, None, rule, False)
                for lab, s in ANCHORS.items():
                    r = M.summarise(M.evaluate(W, J, times, s, None, rule, False), base)
                    rows.append(dict(N=N, D=D, rule=rule, s_label=lab, s=s,
                                     ctrl_median_nonoise=float(np.nanmean(base["med"])), **r,
                                     E1_material=M.mat_e1(r), E2_material=M.mat_e2(r)))
    cells = pd.DataFrame(rows)
    cells.to_csv(OUT / "anchored_sim_cells.csv", index=False)

    c = cells[(cells.N == 100) & (cells.D == 8)]
    single = c[c.rule == "none"]; two = c[c.rule == "two-scan"]
    px1 = bool(single.E1_material.all() and (single.E1 < 0).all())
    px2 = bool(single.E2_material.all())
    two_c = two[two.s_label.isin(["0.199", "0.235"])]
    px3 = bool(two_c.E1_material.all())
    px4 = bool(((c.type1 >= 0.04) & (c.type1 <= 0.06)).all())
    keep = ["rule", "s_label", "E1", "E1_lo", "E1_hi", "E2", "E2_lo", "E2_hi", "E4", "type1", "ctrl_median", "E1_material", "E2_material"]
    summary = dict(
        label="simulation anchored to one published repeatability study; assumed conversion s = RC%/1.96; large lesions only; not evidence about any real trial",
        amendment="A19 in analysis_plan_v1.md", seed=M.SEED, R=M.R, runtime_s=round(time.time() - t0, 1),
        PX1_supported=px1, PX2_supported=px2, PX3_supported=px3, PX4_supported=px4,
        type1_outside_band_N100_D8=c[(c.type1 < 0.04) | (c.type1 > 0.06)][["rule", "s_label", "type1"]].to_dict("records"),
        N100_D8=c[keep].to_dict("records"),
        n_cells=int(len(cells)),
        material_E1_cells=int(cells.E1_material.sum()), material_E2_cells=int(cells.E2_material.sum()),
    )
    (OUT / "anchored_sim_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    pd.set_option("display.width", 220)
    print(c[keep].round(3).to_string())
    print({k: summary[k] for k in ("PX1_supported", "PX2_supported", "PX3_supported", "PX4_supported", "n_cells", "material_E1_cells", "material_E2_cells", "runtime_s")})


if __name__ == "__main__":
    main()
