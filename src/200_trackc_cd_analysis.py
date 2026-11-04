#!/usr/bin/env python3
"""
Script 200: Track C Stage 2, analysis of tests C and D (reads saved files only; no simulation, no seed).
=====================================================================================================
SIMULATION ONLY. Exploratory and POST HOC: the questions below were chosen after the first results of scripts 192 and
193 had been seen, so p-values are descriptive and are not corrected for the number of comparisons.

Why this script exists. The interim report of 2026-11-04 quotes numbers that came from ad hoc code (the pass count under
the plain step median, the share of high-kill worlds that pass, rank tests between worlds that pass the landmark control
and worlds that do not, the number of worlds inside a second published landmark interval). Every number in a report must
cite a file, so this script recomputes them from the saved outputs of scripts 192 and 193 and writes one JSON file.

Inputs (per configuration, see src/trackc_calibration.py):
    output/trackc_stage2_C_calibration<suffix>/worlds.csv, results.json         (script 192)
    output/trackc_stage2_D_negative_controls<suffix>/controls.csv, results.json (script 193)
Output: output/trackc_stage2_CD_analysis<suffix>/analysis.json (refuses to overwrite).

Run:  python src/200_trackc_cd_analysis.py
      python src/200_trackc_cd_analysis.py --config mixture
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
from scipy.stats import mannwhitneyu

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_calibration as cal  # noqa: E402
import trackc_manifest as man  # noqa: E402

SCRIPT = "src/200_trackc_cd_analysis.py"
SCALE_LINE = 0.5                                          # kill scale at which a world is "near the earlier design range"
BLUMENTHAL_HR = (0.80, 0.65, 0.98)                        # post hoc landmark analysis, targets file row P (Neuro-Oncol 2017)
WEAKEST = [
    "post hoc: the comparisons were chosen after the first C and D results were seen; p-values are not corrected for the number of axes and splits",
    "the worlds are a Latin hypercube sample of an assumed family, so shares of worlds are shares of that sample, not of any population",
    "41 against 40 worlds (single family) is a small sample; a rank test between them says little about axes with small effects",
    "the kill scale in the mixture family is the scale of the non-responders; the responders keep the earlier range, so scale-based statements need the responder fraction q beside them",
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


def kill_scale_key(keys):
    return "log_m_scale" if "log_m_scale" in keys else "log_m_scale_non"


def _q(x, qs=(5, 25, 50, 75, 95)):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return [float(v) for v in np.percentile(x, qs)] if x.size else None


def _med(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.median(x)) if x.size else None


def analyse_c(rows, keys, target=cal.TARGET_PFS_DAYS):
    """Counts and axis summaries from the rows of worlds.csv. Pure function."""
    lo, hi = target
    kill = kill_scale_key(keys)
    n = len(rows)
    col = lambda name: np.array([r[name] for r in rows], float)  # noqa: E731
    scale = np.exp(col(kill))
    med_i, med_s = col("base_median_interp"), col("base_median_step")
    pass_i = (med_i >= lo) & (med_i <= hi)
    pass_s = (med_s >= lo) & (med_s <= hi)
    both, either = int((pass_i & pass_s).sum()), int((pass_i | pass_s).sum())
    big = scale >= SCALE_LINE
    out = {"n_worlds": n, "kill_scale_axis": kill, "kill_scale_definition": f"exp({kill})",
           "target_days": [lo, hi],
           "interpolated_rule": {"n_pass": int(pass_i.sum())},
           "step_rule": {"n_pass": int(pass_s.sum()), "n_pass_in_both_rules": both,
                         "jaccard_with_interpolated": float(both / either) if either else None,
                         "median_kill_scale_of_step_passers": _med(scale[pass_s]),
                         "n_step_passers_with_scale_ge_0.5": int((pass_s & big).sum())},
           "step_median_of_interpolated_passers": {f"{k:g}": int(v) for k, v in sorted(Counter(med_s[pass_i]).items())},
           "kill_scale_line": SCALE_LINE,
           "worlds_with_scale_ge_line": {"n": int(big.sum()), "n_pass_interpolated": int((pass_i & big).sum()),
                                         "median_interp_of_passers": [float(v) for v in np.sort(med_i[pass_i & big])],
                                         "median_interp_q05_q25_q50_q75_q95": _q(med_i[big])},
           "axis_medians": {}, "variants": {}}
    for k in keys:
        x = col(k)
        out["axis_medians"][k] = {"all_worlds": _med(x), "interpolated_passers": _med(x[pass_i]), "step_passers": _med(x[pass_s])}
    out["axis_medians"]["kill_scale"] = {"all_worlds": _med(scale), "interpolated_passers": _med(scale[pass_i]),
                                         "step_passers": _med(scale[pass_s])}
    names = sorted({k[:-5] for k in rows[0] if isinstance(k, str) and k.endswith("_pass")}) if rows else []
    base = col("base_pass").astype(bool) if "base_pass" in rows[0] else pass_i
    for v in names:
        pv = col(f"{v}_pass").astype(bool)
        j_both, j_either = int((pv & base).sum()), int((pv | base).sum())
        out["variants"][v] = {"n_pass": int(pv.sum()), "jaccard_with_base": float(j_both / j_either) if j_either else None,
                              "median_kill_scale_of_passers": _med(scale[pv]),
                              "n_passers_with_scale_ge_0.5": int((pv & big).sum()),
                              "median_q_of_passers": _med(col("q")[pv]) if "q" in keys else None}
    return out


def _flag(rows, name):
    return np.array([bool(r[name]) for r in rows])


def mann_whitney_by_axis(rows, flag, keys):
    """Rank test of each world axis between rows with the flag set and rows without it. Uncorrected and descriptive."""
    res = {}
    f = _flag(rows, flag)
    for k in keys:
        x = np.array([r[k] for r in rows], float)
        a, b = x[f], x[~f]
        if a.size < 2 or b.size < 2 or np.ptp(x) == 0:
            res[k] = None
            continue
        u = mannwhitneyu(a, b, alternative="two-sided")
        res[k] = {"n_flag": int(a.size), "n_not_flag": int(b.size), "mean_flag": float(a.mean()), "mean_not_flag": float(b.mean()),
                  "median_flag": float(np.median(a)), "median_not_flag": float(np.median(b)),
                  "auc_flag_higher": float(u.statistic / (a.size * b.size)), "p_two_sided": float(u.pvalue)}
    n_tests = sum(v is not None for v in res.values())
    for v in res.values():
        if v is not None:
            v["p_bonferroni_over_axes"] = float(min(1.0, v["p_two_sided"] * n_tests))
    return res


def analyse_d(rows, keys, target=cal.TARGET_PFS_DAYS, mult="dd2"):
    """Counts, flips, hazard-ratio summaries and rank tests from the rows of controls.csv. Pure function."""
    lo, hi = target
    n = len(rows)
    out = {"n_kept_worlds": n, "primary_dose_column": mult}
    if n == 0:
        return out
    col = lambda name: np.array([r[name] for r in rows], float)  # noqa: E731
    re_pass = _flag(rows, "re_pass")
    out["n_still_in_target"] = int(re_pass.sum())
    out["flips"] = [{"world": int(r["world"]), "first_draw_median": float(r["c_median_interp"]),
                     "redraw_median": float(r["re_median_interp"]),
                     "first_draw_distance_to_nearest_edge_days": float(min(abs(r["c_median_interp"] - lo), abs(r["c_median_interp"] - hi)))}
                    for r in rows if not r["re_pass"]]
    for name in ("dd1.5_ok", "dd2_ok", "dd3_ok", "lm_ok", "controlled"):
        if name in rows[0]:
            out[f"n_{name}"] = int(_flag(rows, name).sum())
    dd, lm = col(f"{mult}_hr"), col("lm_hr")
    out["dose_dense_hr_q05_q25_q50_q75_q95"], out["landmark_hr_q05_q25_q50_q75_q95"] = _q(dd), _q(lm)
    out["landmark_hr_min_max"] = [float(np.nanmin(lm)), float(np.nanmax(lm))]
    out["n_landmark_hr_below_1"] = int((lm < 1.0).sum())
    out["frac_landmark_hr_below_published_lower_bound"] = float(np.mean(lm < cal.TWELVE_VS_SIX_HR[1]))
    h, hlo, hhi = BLUMENTHAL_HR
    in_b = (lm >= hlo) & (lm <= hhi)
    out["second_landmark_interval"] = {"source": "Blumenthal 2017 post hoc landmark analysis, PFS adjusted HR 0.80 (0.65-0.98), non-randomised",
                                       "interval": [hlo, hhi], "n_worlds_inside": int(in_b.sum()),
                                       "n_inside_and_dose_dense_ok": int((in_b & _flag(rows, f"{mult}_ok")).sum())}
    out["rank_tests"] = {flag: mann_whitney_by_axis(rows, flag, keys) for flag in ("lm_ok", f"{mult}_ok", "controlled") if flag in rows[0]}
    ctl = _flag(rows, "controlled") if "controlled" in rows[0] else np.zeros(n, bool)
    out["controlled_axis_ranges"] = ({k: [float(col(k)[ctl].min()), float(col(k)[ctl].max())] for k in keys} if ctl.any() else None)
    if ctl.any() and "mu" in keys and "nu" in keys:
        cells = Counter((float(r["mu"]), float(r["nu"])) for r, c in zip(rows, ctl) if c)
        out["controlled_mu_nu_cells"] = {f"mu={a:g},nu={b:g}": int(v) for (a, b), v in sorted(cells.items())}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="single", help="run configuration of src/trackc_calibration.py")
    args = ap.parse_args()
    cfg = cal.get_config(args.config)
    keys = cal.world_keys(cfg["family"])
    c_dir, d_dir, out_dir = cal.stage_dir("C", cfg), cal.stage_dir("D", cfg), cal.stage_dir("CD", cfg)
    for need in (c_dir / "worlds.csv", c_dir / "results.json", d_dir / "controls.csv", d_dir / "results.json"):
        if not need.exists():
            sys.exit(f"missing input {need}")
    if (out_dir / "analysis.json").exists():
        sys.exit(f"refusing to overwrite {out_dir / 'analysis.json'}")
    c_rows, d_rows = read_csv(c_dir / "worlds.csv"), read_csv(d_dir / "controls.csv")
    c_res = json.loads((c_dir / "results.json").read_text(encoding="utf-8"))
    d_res = json.loads((d_dir / "results.json").read_text(encoding="utf-8"))
    c, d = analyse_c(c_rows, keys), analyse_d(d_rows, keys)
    checks = {"c_n_pass_matches_results_json": c["interpolated_rule"]["n_pass"] == c_res["base"]["n_pass"],
              "d_n_controlled_matches_results_json": d.get("n_controlled") == d_res.get("n_controlled"),
              "d_n_still_in_target_matches_results_json": d.get("n_still_in_target") == d_res.get("n_still_in_target")}
    meta = {"script": "200_trackc_cd_analysis", "config": cfg["name"], "family": cfg["family"],
            "label": "SIMULATION; exploratory; post hoc analysis of saved outputs; no simulation and no seed used",
            "inputs_sha256": {str(p.relative_to(PROJECT_ROOT)): _sha256(p) for p in
                              (c_dir / "worlds.csv", c_dir / "results.json", d_dir / "controls.csv", d_dir / "results.json")},
            "git": man.git_state((SCRIPT, "src/trackc_calibration.py")),
            "code_sha256": man.code_hashes((SCRIPT, "src/trackc_calibration.py")),
            "consistency_with_results_json": checks, "weakest_points": WEAKEST}
    if not all(checks.values()):
        sys.exit(f"inconsistent with the saved results.json files: {checks}")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "analysis.json").write_text(json.dumps({**meta, "c": c, "d": d}, indent=2, default=float) + "\n", encoding="utf-8")
    print(f"wrote {out_dir / 'analysis.json'}")
    print(f"C: interpolated rule {c['interpolated_rule']['n_pass']} pass, step rule {c['step_rule']['n_pass']} pass "
          f"(both {c['step_rule']['n_pass_in_both_rules']}); worlds with scale >= {SCALE_LINE}: {c['worlds_with_scale_ge_line']['n']}, "
          f"passing {c['worlds_with_scale_ge_line']['n_pass_interpolated']}")
    print(f"D: kept {d['n_kept_worlds']}, still in target {d.get('n_still_in_target')}, flips {len(d.get('flips', []))}, "
          f"controlled {d.get('n_controlled')}, inside the second landmark interval {d['second_landmark_interval']['n_worlds_inside']}")


if __name__ == "__main__":
    main()
