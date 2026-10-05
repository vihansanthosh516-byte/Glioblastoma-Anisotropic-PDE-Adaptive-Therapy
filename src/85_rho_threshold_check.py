"""Script 85: does adaptive therapy lose time-to-progression above a growth-rate threshold?

Reads the per-patient adaptive-vs-MTD results of script 44 and reports, by rho bin, how many
patients progress earlier under adaptive, plus the dose sparing. No simulation is rerun.

Input:  output/adaptive_geometry_metrics.json (61 records)
Output: output/rho_threshold_check.json
"""
import json
from pathlib import Path as _Path

import numpy as np
from scipy.stats import spearmanr

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
IN_JSON = OUTPUT_DIR / "adaptive_geometry_metrics.json"
OUT_JSON = OUTPUT_DIR / "rho_threshold_check.json"
BINS = [(0.0, 0.002), (0.002, 0.01), (0.01, 0.02), (0.02, 0.03), (0.03, 1.0)]


def main():
    recs = json.loads(IN_JSON.read_text())
    rho = np.array([r["rho_per_day"] for r in recs])
    spare = np.array([r["drug_reduction"] for r in recs])
    t_mtd = np.array([r["ttp_mtd"] for r in recs], float)
    t_ad = np.array([r["ttp_adaptive"] for r in recs], float)
    earlier = t_ad < t_mtd

    bins = []
    for lo, hi in BINS:
        m = (rho >= lo) & (rho < hi)
        bins.append({
            "rho_range": [lo, hi],
            "n": int(m.sum()),
            "n_adaptive_progressed_earlier": int(earlier[m].sum()),
            "dose_sparing_min": float(spare[m].min()) if m.any() else None,
            "dose_sparing_mean": float(spare[m].mean()) if m.any() else None,
            "dose_sparing_max": float(spare[m].max()) if m.any() else None,
        })
    above = rho > 0.02
    rs = spearmanr(rho, spare)
    res = {
        "script": "85_rho_threshold_check",
        "input": str(IN_JSON.relative_to(PROJECT_ROOT)),
        "n_patients": len(recs),
        "n_adaptive_earlier_total": int(earlier.sum()),
        "rho_gt_0.02": {
            "n": int(above.sum()),
            "n_adaptive_earlier": int(earlier[above].sum()),
        },
        "rho_le_0.02": {
            "n": int((~above).sum()),
            "n_adaptive_earlier": int(earlier[~above].sum()),
        },
        "min_rho_among_earlier": float(rho[earlier].min()) if earlier.any() else None,
        "max_rho_among_not_earlier": float(rho[~earlier].max()),
        "spearman_rho_vs_dose_sparing": {"rho": float(rs.statistic), "p": float(rs.pvalue)},
        "rho_of_23rd_highest_patient": float(np.sort(rho)[::-1][22]),
        "bins": bins,
    }
    OUT_JSON.write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
