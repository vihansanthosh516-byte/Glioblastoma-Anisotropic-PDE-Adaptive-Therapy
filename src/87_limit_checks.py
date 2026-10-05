"""Script 87: three analyses that tighten limits without new data.

L16  The rho > 0.02 threshold (script 85): exact CIs, a permutation test for a perfect split by
     chance, and leave-one-out cutoff selection (does a cutoff chosen on 60 patients classify the 61st?).
L19  MGMT null (script 82): events needed and minimum detectable hazard ratio (Schoenfeld).
L14  Real vs synthetic parameters in the Track C evaluation sets (why real_test gains are small).

Inputs:  output/adaptive_geometry_metrics.json, mgmt_resistance/results.json, rl_kill_conditioned/evaluate.json
Output:  output/limit_checks.json
"""
import json
from pathlib import Path as _Path

import numpy as np
from scipy.stats import beta, norm

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUT_JSON = OUTPUT_DIR / "limit_checks.json"


def clopper_pearson(k, n, a=0.05):
    lo = 0.0 if k == 0 else beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - a / 2, k + 1, n - k)
    return float(lo), float(hi)


def threshold_checks():
    recs = json.loads((OUTPUT_DIR / "adaptive_geometry_metrics.json").read_text())
    rho = np.array([r["rho_per_day"] for r in recs])
    earlier = np.array([r["ttp_adaptive"] < r["ttp_mtd"] for r in recs])
    n, k = len(rho), int(earlier.sum())

    def best_cut(r, e):
        """Cutoff (midpoint between sorted values) with the fewest misclassifications of 'earlier = rho > cut'."""
        s = np.sort(r)
        cands = (s[:-1] + s[1:]) / 2
        errs = [int(((r > c) != e).sum()) for c in cands]
        i = int(np.argmin(errs))
        return float(cands[i]), errs[i]

    cut, err = best_cut(rho, earlier)
    rng = np.random.default_rng(0)
    perm_perfect = 0
    n_perm = 20000
    for _ in range(n_perm):
        if best_cut(rho, rng.permutation(earlier))[1] == 0:
            perm_perfect += 1
    # leave-one-out: choose the cutoff without patient i, classify patient i
    loo_wrong = 0
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        c, _ = best_cut(rho[m], earlier[m])
        loo_wrong += int((rho[i] > c) != earlier[i])
    above, below = rho > 0.02, rho <= 0.02
    return {
        "n": n, "n_earlier": k,
        "rho_gt_0.02": {"n": int(above.sum()), "k_earlier": int(earlier[above].sum()),
                        "exact_ci95_prop_earlier": clopper_pearson(int(earlier[above].sum()), int(above.sum()))},
        "rho_le_0.02": {"n": int(below.sum()), "k_earlier": int(earlier[below].sum()),
                        "exact_ci95_prop_earlier": clopper_pearson(int(earlier[below].sum()), int(below.sum()))},
        "best_cutoff_on_all": cut, "misclassified_at_best_cutoff": err,
        "permutation_p_perfect_split_by_chance": (perm_perfect + 1) / (n_perm + 1), "n_permutations": n_perm,
        "leave_one_out_misclassified": loo_wrong, "leave_one_out_accuracy": 1 - loo_wrong / n,
    }


def mgmt_power():
    m = json.loads((OUTPUT_DIR / "mgmt_resistance" / "results.json").read_text())["clinical"]
    d = m["events"]["mgmt0"] + m["events"]["mgmt1"]
    za, zb = norm.ppf(0.975), norm.ppf(0.80)
    p0 = m["n"]["mgmt0"] / (m["n"]["mgmt0"] + m["n"]["mgmt1"])
    # Schoenfeld: events = (za+zb)^2 / (p0 (1-p0) ln(HR)^2)
    def events_needed(hr):
        return (za + zb) ** 2 / (p0 * (1 - p0) * np.log(hr) ** 2)
    mdhr = float(np.exp((za + zb) / np.sqrt(d * p0 * (1 - p0))))
    return {
        "events_observed": d, "alloc_fraction_mgmt0": float(p0), "alpha": 0.05, "power": 0.80,
        "observed_hr_mgmt0_vs_mgmt1": m["cox_hr_mgmt0_vs_mgmt1"]["hr"],
        "events_needed_for_observed_hr": float(events_needed(m["cox_hr_mgmt0_vs_mgmt1"]["hr"])),
        "events_needed_for_hr_1.5": float(events_needed(1.5)),
        "events_needed_for_hr_2.0": float(events_needed(2.0)),
        "min_detectable_hr_at_observed_events": mdhr,
    }


def real_vs_synthetic():
    sets = json.loads((OUTPUT_DIR / "rl_kill_conditioned" / "evaluate.json").read_text())["sets"]
    out = {}
    for s, v in sets.items():
        rho = np.array(v["rho"])
        kill = np.array(v["kill"])  # columns: none, chemo, rad, combo
        out[s] = {"n": int(v["n"]),
                  "rho_median": float(np.median(rho)), "rho_min": float(rho.min()), "rho_max": float(rho.max()),
                  "kill_chemo_median": float(np.median(kill[:, 1])), "kill_rad_median": float(np.median(kill[:, 2])),
                  "kill_combo_median": float(np.median(kill[:, 3])),
                  "kill_over_rho_median": float(np.median(kill[:, 3] / rho))}
    return out


def main():
    res = {"script": "87_limit_checks", "L16_threshold": threshold_checks(), "L19_mgmt_power": mgmt_power(),
           "L14_real_vs_synthetic": real_vs_synthetic()}
    OUT_JSON.write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
