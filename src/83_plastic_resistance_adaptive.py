#!/usr/bin/env python3
"""
Script 83: Drug-induced, reversible (plastic) resistance and adaptive therapy
(Track C negative 3, model-class follow-up).
=============================================================================
Scripts 76 and 82 modelled resistance as FIXED pre-existing clones. In GBM, temozolomide
resistance is partly drug-induced and reversible (MGMT up-regulation, a transient slow-growth
tolerant state: PMC6944699; internal-variable adaptation models that favour intermittent TMZ:
J R Soc Interface 17:20190722; Comput Biol Med 2024 S001048252400951X). Theory (Kuosmanen et
al. 2021, PLoS Comput Biol) says drug-induced resistance favours less aggressive treatment.
Question: when resistance is plastic, does an adaptive rule beat the same-window paced
non-adaptive arm, by more than in the matching no-plasticity control?

Model = script 76 TwoPopSim plus switching between sensitive (S) and tolerant (T) cells:
   dS/dt = rho S (1 - N/K) - k S - sig_i x S + sig_r (1 - x) T
   dT/dt = rho (1 - c) T (1 - N/K) - k eps T + sig_i x S - sig_r (1 - x) T
with x = 1 on a day with drug (k > 0), 0 otherwise; the switching terms conserve cells.

PRE-SPECIFIED DESIGN (fixed before the first run; not changed after seeing data)
-------------------------------------------------------------------------------
Fixed     f_r0 = 0.001 (script 76's smallest: plasticity, not pre-existing clones, must do
          the work), cost c = 0.25, eps = 0, script-59 seed, window 365, Stupp budget,
          script-76 arms, evaluation sets and RANO TTP to day 365.
Grid      sig_i in {0.005, 0.02, 0.1} /day  (induction while on drug)
          sig_r in {0.005, 0.05}  /day      (reversion while off drug)
          occupancy: seed59 (PRIMARY) and near_capacity (secondary)  -> 12 cells.
          No literature value for sig_i / sig_r is available; the grid spans slow to fast
          switching and every cell is reported, so the finding is a regime map.
Control   script 76's cell f_r0 = 0.001, cost 0.25, eps 0, same occupancy, window 365
          (sig_i = sig_r = 0). Validity check: with sig = 0 this script reproduces it exactly.
Endpoint  TTP (RANO, >= 40% above nadir) restricted to day 365; per set, AT50 / AT80 minus the
          same-window paced arm (95% bootstrap CI) and the plasticity-attributable part
          (difference-in-differences vs the control). Also AT vs best non-adaptive.
POSITIVE (regime) iff, on real_test (the only real-patient set) and the primary occupancy,
          for the same arm (AT50 or AT80), criteria (a) AT - paced CI lower bound > 0 and
          (b) plasticity-attributable days > 0 hold in >= 3 of the 6 cells, and those cells
          form a contiguous region of the (sig_i, sig_r) grid. An isolated cell is reported
          as a corner, not a finding (24 cell x arm tests, no multiplicity protection).
          Otherwise C3 stays negative. All sets, cells and arms are reported either way.
Output: output/plastic_resistance/{results.json, cells/}   (resumable per cell)
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import torch

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
OUT_DIR = PROJECT_ROOT / "output" / "plastic_resistance"

_spec = spec_from_file_location("s76", PROJECT_ROOT / "src" / "76_resistance_adaptive_equal_budget.py")
s76 = module_from_spec(_spec)
_spec.loader.exec_module(s76)
K, N_SUB, DT_SUB = s76.K, s76.N_SUB, s76.DT_SUB

SIG_I = (0.005, 0.02, 0.1)
SIG_R = (0.005, 0.05)
OCC = ("seed59", "near_capacity")
BASE = {"f_r0": 0.001, "cost": 0.25, "eps": 0.0, "window": 365}


class PlasticSim(s76.TwoPopSim):
    sig_i = 0.0
    sig_r = 0.0

    def step(self, actions: torch.Tensor) -> torch.Tensor:
        k = self.kill.gather(1, actions.view(-1, 1))
        x = (k > 0).to(k.dtype)
        for _ in range(N_SUB):
            tot = self.us + self.ur
            to_t = self.sig_i * x * self.us
            to_s = self.sig_r * (1.0 - x) * self.ur
            gs = self.rho * self.us * (1.0 - tot / K) - k * self.us - to_t + to_s
            gr = self.rho * (1.0 - self.cost) * self.ur * (1.0 - tot / K) - k * self.eps * self.ur + to_t - to_s
            self.us = torch.clamp(self.us + DT_SUB * gs, 0.0, K)
            self.ur = torch.clamp(self.ur + DT_SUB * gr, 0.0, K)
        self.u = self.us + self.ur
        return self.volume()


def cfg_for(occ: str) -> dict:
    return {**BASE, "occupancy": occ}


def run_plastic(idx: int, sig_i: float, sig_r: float, occ: str, sets: dict) -> dict:
    PlasticSim.sig_i, PlasticSim.sig_r = sig_i, sig_r
    s76.TwoPopSim = PlasticSim                      # run_cell / best_fixed look this name up in s76
    return s76.run_cell(idx, cfg_for(occ), sets, {})


def main() -> None:
    (OUT_DIR / "cells").mkdir(parents=True, exist_ok=True)
    sets = s76.s68.reaction_sets()
    t0 = time.time()
    # validity: sig = 0 reproduces script 76's control cell exactly
    ctrl = {occ: json.loads(s76.cell_path(cfg_for(occ)).read_text()) for occ in OCC}
    chk = run_plastic(0, 0.0, 0.0, "seed59", {"real_test": sets["real_test"]})
    diff = max(abs(chk["sets"]["real_test"]["arms"][a]["rmst_days"] - ctrl["seed59"]["sets"]["real_test"]["arms"][a]["rmst_days"])
               for a in ("stupp", "paced", "AT50", "AT80", "best_fixed"))
    assert diff < 1e-9, f"PlasticSim(sigma=0) differs from script-76 control by {diff} days"
    print(f"validation OK: sigma = 0 reproduces script 76 (max RMST diff {diff})", flush=True)

    rows, cells = [], {}
    for n, (si, sr, occ) in enumerate(itertools.product(SIG_I, SIG_R, OCC), start=1):
        tag = f"sigi{si}_sigr{sr}_{occ}"
        path = OUT_DIR / "cells" / f"{tag}.json"
        if path.exists():
            cell = json.loads(path.read_text())
        else:
            cell = run_plastic(100 + n, si, sr, occ, sets)
            path.write_text(json.dumps(cell))
        cells[tag] = cell
        rs = cell["sets"]["real_test"]["arms"]
        print(f"[{time.time() - t0:6.0f}s] {tag}  real_test RMST stupp {rs['stupp']['rmst_days']:.1f} paced "
              f"{rs['paced']['rmst_days']:.1f} AT50 {rs['AT50']['rmst_days']:.1f} AT80 {rs['AT80']['rmst_days']:.1f}", flush=True)
        for set_name, st in cell["sets"].items():
            for at in ("AT50", "AT80"):
                a, a0 = st["arms"][at], ctrl[occ]["sets"][set_name]["arms"][at]
                rows.append({"sig_i": si, "sig_r": sr, "occupancy": occ, "set": set_name, "arm": at,
                             "vs_paced": a["vs_paced"]["mean_diff_days"], "vs_paced_ci95": a["vs_paced"]["ci95"],
                             "control_vs_paced": a0["vs_paced"]["mean_diff_days"],
                             "plasticity_attributable_days": a["vs_paced"]["mean_diff_days"] - a0["vs_paced"]["mean_diff_days"],
                             "vs_best_nonadaptive": a["vs_best_nonadaptive"]["mean_diff_days"],
                             "vs_best_nonadaptive_ci95": a["vs_best_nonadaptive"]["ci95"],
                             "tolerant_fraction_day365_median": a["resistant_fraction_day365_median"]})
    # decision rule on real_test, primary occupancy
    verdict = {}
    for at in ("AT50", "AT80"):
        hit = {(r["sig_i"], r["sig_r"]) for r in rows if r["set"] == "real_test" and r["occupancy"] == "seed59"
               and r["arm"] == at and r["vs_paced_ci95"][0] > 0 and r["plasticity_attributable_days"] > 0}
        contiguous = False
        if len(hit) >= 3:
            ii, jj = {s: i for i, s in enumerate(SIG_I)}, {s: j for j, s in enumerate(SIG_R)}
            idxs = {(ii[a], jj[b]) for a, b in hit}
            seen, stack = set(), [next(iter(idxs))]
            while stack:
                c = stack.pop()
                if c in seen:
                    continue
                seen.add(c)
                stack += [(c[0] + d0, c[1] + d1) for d0, d1 in ((1, 0), (-1, 0), (0, 1), (0, -1))
                          if (c[0] + d0, c[1] + d1) in idxs and (c[0] + d0, c[1] + d1) not in seen]
            contiguous = len(seen) == len(idxs)
        verdict[at] = {"n_cells_meeting_criteria_of_6": len(hit), "cells": sorted(hit), "contiguous": contiguous,
                       "positive_regime": len(hit) >= 3 and contiguous}
    res = {"script": "83_plastic_resistance_adaptive", "design": {"base": BASE, "sig_i": SIG_I, "sig_r": SIG_R,
                                                                    "occupancy": OCC},
           "validation_sigma0_max_rmst_diff_days": diff, "rows": rows, "verdict_real_test_seed59": verdict,
           "positive_overall": any(v["positive_regime"] for v in verdict.values())}
    (OUT_DIR / "results.json").write_text(json.dumps(res, indent=2))
    print(json.dumps({k: v for k, v in res.items() if k != "rows"}, indent=2))
    print(f"[saved] {OUT_DIR / 'results.json'}")


if __name__ == "__main__":
    main()
