#!/usr/bin/env python3
"""
Script 202: Track C Stage 2, feedback rules against the fixed-schedule frontier at equal drug use
================================================================================================
SIMULATION ONLY. Exploratory and POST HOC: reads the saved output of script 194 (no simulation, no seed). The question
was chosen after the first regime-map result (config single) had been seen.

Question. Script 194 reports each arm's gain in restricted mean time to progression against a reference arm, but the arms
use different amounts of drug. A rule that gains less time than "give 12 courses at once" and uses fewer courses is not
necessarily worse. The fair comparison is at equal drug use: does a rule that reads scans reach more restricted time than a
fixed schedule that uses the same mean number of courses?

Method. For each world, endpoint and horizon, the fixed-schedule points are the arms of script 194 that ignore readings
(soc6 or std_capmix, front3, front9, front12, paced2, paced3, best_fixed and the best_fixed_lam arms; the best_fixed arms
were chosen on separate selection patients). Each is a point (mean courses used, restricted mean time). The frontier is
the upper concave envelope of these points: a coin-flip assignment of patients between two fixed schedules reaches every
point on the line between them, so the envelope is what fixed dosing can do. For every arm that reads scans (raw_*,
nadir_*, resp_*) the gap is its restricted mean time minus the frontier at its own mean number of courses; a positive gap
means it beats fixed dosing at equal drug use. The informed arms (informed_*; they know each patient's parameters and are
not realisable) are reported next to them as the headroom. Over worlds the script reports the median and quartiles of the
gap, the share of worlds above and below the frontier (a difference of 0.05 days or less counts as on the frontier), and
an exact two-sided sign test over the worlds that are not ties (descriptive; not corrected for the number of arms).

Input:  output/trackc_stage2_E_regime_map<suffix>/arms_by_world.csv and results.json   (script 194)
Output: output/trackc_stage2_E_frontier<suffix>/frontier.json and gaps_by_world.csv  (refuses to overwrite)

Run:  python src/202_trackc_frontier_analysis.py --config mixture
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_calibration as cal  # noqa: E402
import trackc_manifest as man  # noqa: E402

SCRIPT = "src/202_trackc_frontier_analysis.py"
FEEDBACK_PREFIXES = ("raw_", "nadir_", "resp_")
INFORMED_PREFIX = "informed"
TIE_DAYS = 0.05
WEAKEST = [
    "post hoc: the question was chosen after the first regime-map result was seen; the sign tests are descriptive, uncorrected for the number of arms and endpoints, and the worlds are a sample of an assumed family, not of a population",
    "the fixed-schedule frontier is built from about a dozen fixed arms (front3, front9, front12, paced2, paced3, the six-course or mixed-cap standard arm and the best_fixed arms), not from all 56 schedules of script 194; a richer fixed family would lift the frontier, so a gap above the frontier is an upper bound on what the rule gains over fixed dosing",
    "each gap uses 1,000 evaluation patients per world, so a single world carries a sampling error of a few days; the paired errors of a gap are not in the saved file (only paired differences to the reference arm are)",
    "the informed arms use each patient's true parameters and are not realisable; their gap is a headroom, not a result for any rule",
    "the best feedback arm of each world is a maximum over several arms and is biased upward by selection",
    "restricted mean time is a mean; it says nothing about the worst-case patients (see the paired-regret column of script 194)",
]


def _num(s):
    try:
        return float(s)
    except ValueError:
        return s


def read_csv(path):
    with open(path, encoding="utf-8") as fh:
        return [{k: _num(v) for k, v in row.items()} for row in csv.DictReader(fh)]


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def family(arm):
    """'feedback' (reads scans), 'informed' (knows the patient's parameters) or 'fixed' (ignores readings)."""
    if arm.startswith(FEEDBACK_PREFIXES):
        return "feedback"
    if arm.startswith(INFORMED_PREFIX):
        return "informed"
    return "fixed"


def upper_hull(points):
    """Upper concave envelope of (x, y) points: vertices sorted by x (for equal x the highest y is kept)."""
    best = {}
    for x, y in points:
        x = round(float(x), 9)
        if x not in best or y > best[x]:
            best[x] = float(y)
    hull = []
    for p in sorted(best.items()):
        while len(hull) >= 2 and ((hull[-1][1] - hull[-2][1]) * (p[0] - hull[-2][0])
                                  <= (p[1] - hull[-2][1]) * (hull[-1][0] - hull[-2][0])):
            hull.pop()
        hull.append(p)
    return hull


def frontier_at(hull, c):
    """Height of the envelope at c; NaN outside the range of the fixed points or with fewer than two vertices."""
    if len(hull) < 2 or c < hull[0][0] - 1e-9 or c > hull[-1][0] + 1e-9:
        return float("nan")
    return float(np.interp(c, [h[0] for h in hull], [h[1] for h in hull]))


def gaps_for(rows):
    """Long-format gap rows for every world, endpoint, horizon and non-fixed arm. Pure function."""
    out = []
    cells = sorted({(int(r["world"]), r["endpoint"], int(r["tau"])) for r in rows})
    for world, ep, tau in cells:
        sel = [r for r in rows if int(r["world"]) == world and r["endpoint"] == ep and int(r["tau"]) == tau]
        hull = upper_hull([(r["mean_courses"], r["rmst"]) for r in sel if family(r["arm"]) == "fixed"])
        for r in sel:
            fam = family(r["arm"])
            if fam == "fixed":
                continue
            fr = frontier_at(hull, r["mean_courses"])
            out.append({"world": world, "endpoint": ep, "tau": tau, "arm": r["arm"], "family": fam,
                        "mean_courses": r["mean_courses"], "rmst": r["rmst"], "frontier": fr,
                        "gap": r["rmst"] - fr if np.isfinite(fr) else float("nan")})
    return out


def summarise(gaps, courses):
    """Summary of the gaps of one arm over worlds (NaN gaps are outside the range of the fixed arms and are counted)."""
    g = np.asarray(gaps, float)
    ok = np.isfinite(g)
    n_out = int((~ok).sum())
    g = g[ok]
    if g.size == 0:
        return {"n_worlds": 0, "n_outside_fixed_range": n_out}
    above, below = int((g > TIE_DAYS).sum()), int((g < -TIE_DAYS).sum())
    p = float(binomtest(above, above + below, 0.5).pvalue) if above + below > 0 else None
    return {"n_worlds": int(g.size), "n_outside_fixed_range": n_out, "median_gap_days": float(np.median(g)),
            "q25_q75": [float(q) for q in np.percentile(g, [25, 75])], "min_max": [float(g.min()), float(g.max())],
            "n_above_frontier": above, "n_below_frontier": below, "n_on_frontier": int(g.size - above - below),
            "frac_above_frontier": above / g.size, "frac_below_frontier": below / g.size,
            "sign_test_p_two_sided": p,
            "median_mean_courses": float(np.median(np.asarray(courses, float)[ok]))}


def analyse(rows):
    """Summary over worlds for every endpoint, horizon and arm. Pure function, no I/O."""
    gap_rows = gaps_for(rows)
    res = {"tie_days": TIE_DAYS, "n_worlds": len({int(r["world"]) for r in rows}), "endpoints": {}}
    for ep in sorted({g["endpoint"] for g in gap_rows}):
        res["endpoints"][ep] = {}
        for tau in sorted({g["tau"] for g in gap_rows}):
            block = {"arms": {}}
            arms = sorted({g["arm"] for g in gap_rows if g["endpoint"] == ep and g["tau"] == tau})
            for arm in arms:
                sel = sorted([g for g in gap_rows if g["endpoint"] == ep and g["tau"] == tau and g["arm"] == arm],
                             key=lambda g: g["world"])
                entry = summarise([g["gap"] for g in sel], [g["mean_courses"] for g in sel])
                entry["family"] = family(arm)
                block["arms"][arm] = entry
            # best realisable-looking feedback arm in each world (a maximum over arms: biased upward by selection)
            per_world = {}
            for g in gap_rows:
                if g["endpoint"] == ep and g["tau"] == tau and g["family"] == "feedback" and np.isfinite(g["gap"]):
                    if g["world"] not in per_world or g["gap"] > per_world[g["world"]][0]:
                        per_world[g["world"]] = (g["gap"], g["arm"])
            if per_world:
                best = summarise([v[0] for v in per_world.values()], [0.0] * len(per_world))
                best.pop("median_mean_courses", None)
                best["most_often_best_arm"] = Counter(v[1] for v in per_world.values()).most_common(1)[0][0]
                block["best_feedback_arm_per_world"] = best
            res["endpoints"][ep][str(tau)] = block
    return res


def write_gaps(rows, out_dir):
    gap_rows = gaps_for(rows)
    cols = ["world", "endpoint", "tau", "arm", "family", "mean_courses", "rmst", "frontier", "gap"]
    with open(out_dir / "gaps_by_world.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for g in gap_rows:
            w.writerow([f"{g[c]:.6g}" if isinstance(g[c], float) else g[c] for c in cols])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="single", help="run configuration of src/trackc_calibration.py")
    args = ap.parse_args()
    cfg = cal.get_config(args.config)
    e_dir, out_dir = cal.stage_dir("E", cfg), cal.stage_dir("F", cfg)
    for need in (e_dir / "arms_by_world.csv", e_dir / "results.json"):
        if not need.exists():
            sys.exit(f"missing input {need}")
    if (out_dir / "frontier.json").exists():
        sys.exit(f"refusing to overwrite {out_dir / 'frontier.json'}")
    rows = read_csv(e_dir / "arms_by_world.csv")
    e_res = json.loads((e_dir / "results.json").read_text(encoding="utf-8"))
    res = analyse(rows)
    if res["n_worlds"] != e_res["n_controlled_worlds"]:
        sys.exit(f"world count differs from results.json: {res['n_worlds']} against {e_res['n_controlled_worlds']}")
    meta = {"script": "202_trackc_frontier_analysis", "config": cfg["name"], "family": cfg["family"],
            "reference_arm_of_script_194": e_res.get("reference"),
            "label": "SIMULATION; exploratory; post hoc analysis of saved outputs; no simulation and no seed used",
            "inputs_sha256": {str(p.relative_to(PROJECT_ROOT)): _sha256(p) for p in (e_dir / "arms_by_world.csv", e_dir / "results.json")},
            "git": man.git_state((SCRIPT,)), "code_sha256": man.code_hashes((SCRIPT,)), "weakest_points": WEAKEST}
    out_dir.mkdir(parents=True, exist_ok=True)
    write_gaps(rows, out_dir)
    (out_dir / "frontier.json").write_text(json.dumps({**meta, **res}, indent=2, default=float) + "\n", encoding="utf-8")
    print(f"wrote {out_dir / 'frontier.json'}")
    for ep, blocks in res["endpoints"].items():
        for tau, block in blocks.items():
            print(f"== {ep}, restricted mean {tau} d: gap to the fixed frontier (days; positive = above fixed dosing at equal drug use)")
            for arm, e in sorted(block["arms"].items(), key=lambda kv: -kv[1].get("median_gap_days", -1e9)):
                if e["n_worlds"]:
                    print(f"   {arm:22s} {e['family']:9s} median {e['median_gap_days']:+7.2f}  above in {e['n_above_frontier']:3d} / below in "
                          f"{e['n_below_frontier']:3d} of {e['n_worlds']} worlds  courses {e['median_mean_courses']:.1f}")


if __name__ == "__main__":
    main()
