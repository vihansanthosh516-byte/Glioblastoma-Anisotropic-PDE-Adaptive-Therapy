"""Script 89: fair baseline for the script-47 dual-agent arm (limit L10).

Script 47 compares dual-agent MPC (TMZ plus a secondary drug on every TMZ-holiday day) against single-drug
MTD (TMZ 5-on/23-off). The dual arm gets a second drug that the MTD arm never receives. This script adds
  mtd_plus_secondary : TMZ 5-on/23-off plus the secondary drug on every off day (a2 = 1 when TMZ is off),
using the same reduced ODE, patients and initial states as script 47. No MPC is needed.
Drug exposure is reported per drug (TMZ AUC and secondary-drug AUC) as well as combined.

Input:  src/47_optimal_control.py (imported for the model), output/dual_drug_comparison.json
Output: output/dual_mtd_baseline.json
"""
import importlib.util
import json
from pathlib import Path as _Path

import numpy as np

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUT_JSON = OUTPUT_DIR / "dual_mtd_baseline.json"

spec = importlib.util.spec_from_file_location("s47", PROJECT_ROOT / "src" / "47_optimal_control.py")
s47 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s47)


def run_mtd_plus_secondary(state0, rho_s, rho_r):
    state = state0.copy()
    n = s47.N_SIM_STEPS
    Ms, Mr, C, C2 = (np.zeros(n) for _ in range(4))
    for step in range(n):
        a = 1.0 if int(step * s47.DT) % 28 < 5 else 0.0
        a2 = 1.0 if a < 0.01 else 0.0
        state = s47.step_reduced_ode(state, a, a2, rho_s=rho_s, rho_r=rho_r)
        Ms[step], Mr[step], C[step], C2[step] = state
    return {"M_s": Ms, "M_r": Mr, "C": C, "C2": C2}


def main():
    prev = {p["patient_id"]: p for p in json.loads((OUTPUT_DIR / "dual_drug_comparison.json").read_text())["patients"]}
    rows = []
    for pid in s47.COHORT_PATIENTS:
        s0 = s47.get_patient_state(pid)
        rs = s47._params[pid]["rho_per_day"]
        rr = rs * 0.75
        r = run_mtd_plus_secondary(s0, rs, rr)
        traj = s47._trajectory_from_result(r)
        tot = r["M_s"] + r["M_r"]
        row = {
            "patient_id": pid, "rho_per_day": rs,
            "mtd_plus_secondary": {
                "ttp_days": s47.compute_ttp(tot),
                "resistant_fraction": s47.resistant_fraction(r["M_s"], r["M_r"]),
                "auc_total": s47.compute_auc(traj),
                "auc_tmz": float(np.trapezoid(r["C"], dx=s47.DT)),
                "auc_secondary": float(np.trapezoid(r["C2"], dx=s47.DT)),
            },
            "from_script47": prev[pid],
        }
        rows.append(row)
        d, m = prev[pid]["dual_adaptive"], row["mtd_plus_secondary"]
        print(f"{pid} rho={rs:.4f}  dual-adaptive TTP {d['ttp_days']:.0f} Rf {d['resistant_fraction']:.4f} AUC {d['auc']:.0f}"
              f"  | MTD+secondary TTP {m['ttp_days']:.0f} Rf {m['resistant_fraction']:.4f} AUC {m['auc_total']:.0f}")
    OUT_JSON.write_text(json.dumps({"script": "89_dual_mtd_baseline", "patients": rows}, indent=2))


if __name__ == "__main__":
    main()
