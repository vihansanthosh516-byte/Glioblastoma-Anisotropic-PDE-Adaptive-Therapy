"""Script 92: does the adaptive-therapy result depend on the assumed kill scale E_MAX_RATIO? (limit L8)

Script 44 sets the TMZ kill rate per patient as E_MAX_RATIO x that patient's growth rate (default 1000, chosen by a
sweep that was only recorded in a code comment). Here script 44 is re-run with ratio 500, 1000 (reproduction check)
and 2000 (`ADAPTIVE_E_MAX_RATIO`, `ADAPTIVE_OUT_DIR`; runs in output/emax_sweep/r{ratio}/), and the same cohort
metrics are tabulated.

PRE-SPECIFIED (before the 500 and 2000 results were read)
Reproduction check: ratio 1000 must match the committed output/adaptive_cohort_summary.json
  (10 earlier, 22 any progression, 39 both censored, mean drug 32.36% of MTD, mean VR 2.391).
Reported per ratio, same 61 patients: n where MTD ends with more resistant cells than adaptive (and the reverse);
  median resistant fraction in each arm; mean drug % of MTD; mean final-mass ratio (adaptive / MTD);
  n where adaptive progresses earlier / later / same; TTP-censored counts; and the rho > 0.02 split
  (how many of the patients above 0.02 /day progress earlier, and how many below).
No threshold or conclusion is changed after seeing the numbers.

Input:  output/emax_sweep/r{500,1000,2000}/adaptive_geometry_metrics.json, adaptive_cohort_summary.json
Output: output/emax_sweep_summary.json
"""
import json
from pathlib import Path as _Path

import numpy as np

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
RATIOS = (500, 1000, 2000)


def summarise(ratio):
    d = OUT / "emax_sweep" / f"r{ratio}"
    m = json.loads((d / "adaptive_geometry_metrics.json").read_text())
    s = json.loads((d / "adaptive_cohort_summary.json").read_text())
    rm = np.array([p["resistant_fraction_mtd"] for p in m])
    ra = np.array([p["resistant_fraction_adaptive"] for p in m])
    tm = np.array([p["ttp_mtd"] for p in m], float)
    ta = np.array([p["ttp_adaptive"] for p in m], float)
    rho = np.array([p["rho_per_day"] for p in m])
    hi = rho > 0.02
    earlier, later = ta < tm, ta > tm
    return {
        "e_max_ratio_in_summary": s["e_max_ratio"], "n": len(m),
        "n_mtd_more_resistant": int((rm > ra).sum()), "n_adaptive_more_resistant": int((ra > rm).sum()),
        "n_resistance_ties": int((rm == ra).sum()),
        "median_resistant_fraction_mtd": float(np.median(rm)), "median_resistant_fraction_adaptive": float(np.median(ra)),
        "mean_drug_percent_of_mtd": float(np.mean([p["drug_percent_of_mtd"] for p in m])),
        "mean_vr_ratio": float(np.mean([p["vr_ratio"] for p in m])),
        "mean_final_mass_mtd": s["mean_final_mass_mtd"], "mean_final_mass_adaptive": s["mean_final_mass_adaptive"],
        "ttp_adaptive_earlier": int(earlier.sum()), "ttp_adaptive_later": int(later.sum()),
        "ttp_same": int((ta == tm).sum()), "n_ttp_both_censored": s["n_ttp_both_censored"],
        "n_ttp_any_progression": s["n_ttp_any_progression"],
        "rho_gt_0.02": {"n": int(hi.sum()), "adaptive_earlier": int(earlier[hi].sum())},
        "rho_le_0.02": {"n": int((~hi).sum()), "adaptive_earlier": int(earlier[~hi].sum())},
    }


def main():
    res = {"script": "92_emax_sweep_summary", "ratios": {}}
    for r in RATIOS:
        res["ratios"][str(r)] = summarise(r)
    ref = json.loads((OUT / "adaptive_cohort_summary.json").read_text())
    r1000 = res["ratios"]["1000"]
    res["reproduction_check_ratio_1000"] = {
        "committed": {"earlier": ref["n_ttp_adaptive_earlier"], "any_progression": ref["n_ttp_any_progression"],
                      "both_censored": ref["n_ttp_both_censored"], "mean_drug": ref["mean_drug_percent_of_mtd"],
                      "mean_vr": ref["mean_vr_ratio"]},
        "rerun": {"earlier": r1000["ttp_adaptive_earlier"], "any_progression": r1000["n_ttp_any_progression"],
                  "both_censored": r1000["n_ttp_both_censored"], "mean_drug": r1000["mean_drug_percent_of_mtd"],
                  "mean_vr": r1000["mean_vr_ratio"]}}
    (OUT / "emax_sweep_summary.json").write_text(json.dumps(res, indent=2))
    for r, v in res["ratios"].items():
        print(r, v["n"], "MTD>adapt resist", v["n_mtd_more_resistant"], "rev", v["n_adaptive_more_resistant"],
              "ties", v["n_resistance_ties"], "| med res", round(v["median_resistant_fraction_mtd"], 3),
              round(v["median_resistant_fraction_adaptive"], 3), "| drug%", round(v["mean_drug_percent_of_mtd"], 1),
              "VR", round(v["mean_vr_ratio"], 2), "| TTP adapt earlier/later/same", v["ttp_adaptive_earlier"],
              v["ttp_adaptive_later"], v["ttp_same"], "| rho>.02", v["rho_gt_0.02"], "rho<=.02", v["rho_le_0.02"])
    print("repro", res["reproduction_check_ratio_1000"])


if __name__ == "__main__":
    main()
