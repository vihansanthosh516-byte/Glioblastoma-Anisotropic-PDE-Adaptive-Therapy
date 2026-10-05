"""Script 93: more real patients, and kill rates that are not one fixed table (limits L4, L14).

Problem: `real_test` (script 67/75) is 21 MU-Glioma patients with REAL growth rates but the SAME assumed kill table
[0, 0.05, 0.08, 0.13] for every patient. Gains there are small (+2.1 d pacing).
This script (a) uses more real growth rates and (b) draws each patient's kill rates from the script-67 training
ranges instead of one fixed table, to see whether the conclusions depend on the fixed table.

PRE-SPECIFIED (before any result)
Rho sets:   mu91     all 91 growing MU-Glioma patients (script 66 loader; sorted, positive rho)
            lumiere_hd   positive per-patient growth rates from LUMIERE, HD-GLIO-AUTO contrast-enhancing volume (script 91 rule)
            lumiere_deep same with DeepBraTumIA
Kill sets:  fixed      script 66 table, same for all patients
            random_s{0,1,2}  chemo ~ U(0.01, 0.10), rad ~ U(0.01, 0.13), combo = chemo + rad (script 67 ranges),
                       drawn per patient with rng seed 93 + s
Arms (no training involved, so they are valid on the 70 MU patients that script 67's networks saw):
            stupp (reference), TRUEKILL_efficiency_rule (one-line rule given the true kills),
            paced_heuristic59, probe_rule L=1 at noise sigma 0, 0.05 (3 noise seeds)
Metrics:    exactly script 75: day-90 win rate vs Stupp, TTP minus Stupp (days, bootstrap CI). Same equal drug budget.
Network arms (kill-conditioned PPO/DAgger) are NOT run here: they were trained on a pool that resampled MU training rho.
Kill rates stay ASSUMED. This tests sensitivity to the assumption, not the truth of any kill rate.

Output: output/real_rho_kill_sensitivity.json
"""
import json
import sys
import time
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path as _Path

import numpy as np
import torch

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
OUT_JSON = PROJECT_ROOT / "output" / "real_rho_kill_sensitivity.json"


def load(name, file):
    spec = spec_from_file_location(name, PROJECT_ROOT / "src" / file)
    mod = module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


s75 = load("s75", "75_probe_and_paced_policies.py")
s91 = load("s91", "91_lumiere_volume_forecast.py")
s68, s67, s66 = s75.s68, s75.s67, s75.s66
ARMS = ["TRUEKILL_efficiency_rule", "paced_heuristic59", "probe_rule_L1_sig0.00_s0",
        "probe_rule_L1_sig0.05_s0", "probe_rule_L1_sig0.05_s1", "probe_rule_L1_sig0.05_s2"]


def lumiere_rho(name):
    ser = s91.follow_up_series(s91.load_volumes(name), s91.load_ratings())
    out = []
    for s in ser.values():
        s = s[s["vol"] > 0]
        if len(s) >= 2 and s["day"].nunique() >= 2:
            r = float(np.polyfit(s["day"], np.log(s["vol"]), 1)[0])
            if r > 0:
                out.append(r)
    return out


def mu_rho():
    rows = s66.real_cohort()
    return [p["rho"] for p in rows["train"] + rows["test"]]


def kill_table(n, mode, seed):
    if mode == "fixed":
        return s66.KILL.view(1, 4).repeat(n, 1)
    rng = np.random.default_rng(93 + seed)
    kc = rng.uniform(0.01, 0.10, n)
    kr = rng.uniform(0.01, 0.13, n)
    return s67.kill_table(kc, kr)


def main():
    rho_sets = {"mu91": mu_rho(), "lumiere_hd": lumiere_rho("hdglioauto"), "lumiere_deep": lumiere_rho("deepbratumia")}
    kill_modes = [("fixed", 0), ("random_s0", 0), ("random_s1", 1), ("random_s2", 2)]
    rng = np.random.default_rng(93)
    res = {"script": "93_real_rho_kill_sensitivity", "budget": s75.BUDGET, "arms": ARMS, "sets": {}}
    for rname, rho in rho_sets.items():
        for kname, seed in kill_modes:
            t0 = time.time()
            kill = kill_table(len(rho), "fixed" if kname == "fixed" else "random", seed)
            facs = s75.arms(kill)
            sim = lambda: s67.KillReactionSim(rho, kill)  # noqa: E731
            runs = {"stupp": s68.rollout_ttp(sim(), facs["stupp"]())}
            for a in ARMS:
                runs[a] = s68.rollout_ttp(sim(), facs[a]())
            assert all(float(r["auc"].max()) <= s75.BUDGET + 1e-9 for r in runs.values())
            block = {a: s75.stats_vs(runs[a], runs["stupp"], rng) for a in ARMS}
            kr = kill[:, 3] / torch.tensor(rho, dtype=kill.dtype)
            res["sets"][f"{rname}__{kname}"] = {
                "n": len(rho), "rho_median": float(np.median(rho)), "kill_over_rho_median": float(kr.median()),
                "arms": block}
            b = block
            print(f"{rname:12s} {kname:10s} n={len(rho):3d} kill/rho med {float(kr.median()):6.1f} | "
                  f"rule d90 win {b['TRUEKILL_efficiency_rule']['day90_win_rate_vs_stupp_pct']:5.1f}% "
                  f"TTP {b['TRUEKILL_efficiency_rule']['ttp_minus_stupp_days_mean']:+6.1f} | "
                  f"paced d90 {b['paced_heuristic59']['day90_win_rate_vs_stupp_pct']:5.1f}% "
                  f"TTP {b['paced_heuristic59']['ttp_minus_stupp_days_mean']:+6.1f} "
                  f"CI {b['paced_heuristic59']['ttp_minus_stupp_days_ci95'][0]:+.1f}..{b['paced_heuristic59']['ttp_minus_stupp_days_ci95'][1]:+.1f} "
                  f"({time.time() - t0:.0f}s)", flush=True)
    OUT_JSON.write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
