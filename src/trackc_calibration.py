"""Track C Stage 2 helpers: the assumed world family, the open-loop standard-of-care arm, PFS summaries.

SIMULATION ONLY. Every range below is an ASSUMPTION, not a measurement.

A "world" is one setting of the parameter ranges of the old simulator (src/trackc_sim.py). Patients of a world are
drawn from those ranges. The published anchor is the RTOG 0525 standard arm (Gilbert 2013, J Clin Oncol
31:4085-4091): median PFS 5.5 months (95% CI 4.7 to 6.1), scans about every 3 months. Verification levels are in
paper/trackc_calibration_targets_2026-11-03.md. The anchor fixes one number per world; it does not say which
world is true.

World axes (all ASSUMED):
    log_r_shift   shifts the growth-rate range r in [0.001, 0.05] per day by exp(shift)
    log_m_scale   scales the kill-to-growth ratio range m in [6, 20] by exp(scale) (kill strength is NOT anchored)
    log_f0_shift  shifts the initial resistant fraction range f0 in [0.001, 0.1] by exp(shift), capped at 0.9
    B0_center     centre of the starting-burden range (half width 0.10); B0 is the crowding axis
    c             cost of resistance (constant inside a world)
    mu, nu        drug-induced switching and reversion rates per day, from the levels used in script 177

Mixture family (config "mixture", 2026-11-04; all ASSUMED). The family above gives every patient of a world the same
kill-strength range, so a short median PFS can only come from a weak drug. The mixture family tests that reading: a
fraction q of the patients ("responders") keep the kill range of the earlier design, m in [6, 20]; the others
("non-responders") have the weak range, scaled by exp(log_m_scale_non). The anchor and the controls are then asked of
this family with the same code. It is a robustness probe, not a fit, and it is not an MGMT model.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.stats import qmc

import trackc_sim as tc
import trackc_survival as sv

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DAYS_PER_MONTH = 30.44
TARGET_PFS_DAYS = (4.7 * DAYS_PER_MONTH, 6.1 * DAYS_PER_MONTH)       # 143.1 to 185.7 d (RTOG 0525 standard arm)
STD_COURSES = 6                                                       # label: 6 maintenance cycles
SCAN_EVERY = 3                                                        # cycles, i.e. every 84 d
MU_NU_LEVELS = (0.0, 0.01, 0.04)                                      # per day; script 177 stress values, ASSUMED
F0_CAP = 0.9
B0_HALF_WIDTH = 0.10
AXES = (("log_r_shift", -0.7, 0.9), ("log_m_scale", float(np.log(0.05)), 0.0),
        ("log_f0_shift", -0.7, 2.3), ("B0_center", 0.12, 0.55), ("c", 0.0, 0.25))
WORLD_KEYS = [a[0] for a in AXES] + ["mu", "nu"]
# Mixture family: the kill-scale axis is the scale of the NON-responders; q is the responder fraction.
MIX_AXES = (("log_r_shift", -0.7, 0.9), ("log_m_scale_non", float(np.log(0.02)), float(np.log(0.3))),
            ("log_f0_shift", -0.7, 2.3), ("B0_center", 0.12, 0.55), ("c", 0.0, 0.25), ("q", 0.1, 0.6))
FAMILY_AXES = {"single": AXES, "mixture": MIX_AXES}
SCAN_TIMES = (84, 168, 252, 336, 504)


def world_keys(family="single"):
    """Names of the world columns of one family, in a fixed order (the continuous axes, then mu and nu)."""
    return [a[0] for a in FAMILY_AXES[family]] + ["mu", "nu"]


# ---------------------------------------------------------------- worlds
def draw_worlds(n, seed, family="single"):
    """n worlds by Latin hypercube over the continuous axes; mu and nu take the three levels with equal weight."""
    axes = FAMILY_AXES[family]
    rng = np.random.default_rng(seed)
    u = qmc.LatinHypercube(d=len(axes) + 2, seed=rng).random(n)
    worlds = []
    for i in range(n):
        w = {"world": i}
        for k, (name, lo, hi) in enumerate(axes):
            w[name] = float(lo + u[i, k] * (hi - lo))
        for k, name in enumerate(("mu", "nu")):
            w[name] = MU_NU_LEVELS[min(int(u[i, len(axes) + k] * len(MU_NU_LEVELS)), len(MU_NU_LEVELS) - 1)]
        worlds.append(w)
    return worlds


def world_ranges(w):
    """The parameter-range dictionary of one world, in the format of tc.draw_scenarios.
    For a mixture world this is the range of the non-responders (the responders are handled in draw_patients)."""
    base = tc.RANGES
    log_m = w["log_m_scale"] if "log_m_scale" in w else w["log_m_scale_non"]
    s_r, s_m, s_f = np.exp(w["log_r_shift"]), np.exp(log_m), np.exp(w["log_f0_shift"])
    ranges = dict(base)
    ranges["r"] = (base["r"][0] * s_r, base["r"][1] * s_r, "log")
    ranges["m"] = (base["m"][0] * s_m, base["m"][1] * s_m, "log")
    ranges["f0"] = (base["f0"][0] * s_f, min(base["f0"][1] * s_f, F0_CAP), "log")
    ranges["B0"] = (w["B0_center"] - B0_HALF_WIDTH, w["B0_center"] + B0_HALF_WIDTH, "lin")
    ranges["c"] = (w["c"], w["c"], "lin")
    return ranges


def draw_patients(w, n, seed):
    """n simulated patients of one world (arrays of length n, ready for tc.simulate).
    A world with a responder fraction q (mixture family) gives each patient the earlier kill range m in [6, 20] with
    probability q; the draw for the responders and the responder flags use their own random streams, so the other
    parameters of the patients do not depend on q."""
    scen = tc.draw_scenarios(n, seed, ranges=world_ranges(w))
    scen["mu"] = np.full(n, float(w["mu"]))
    scen["nu"] = np.full(n, float(w["nu"]))
    if "q" in w:
        seq = [int(s) for s in seed] if np.iterable(seed) else [int(seed)]
        responder = np.random.default_rng(seq + [2]).random(n) < float(w["q"])
        m_resp = tc.draw_scenarios(n, seq + [1], ranges=dict(world_ranges(w), m=tc.RANGES["m"]))["m"]
        scen["m"] = np.where(responder, m_resp, scen["m"])
    return scen


def concat(scens):
    """Stack scenario dictionaries (one per world) into one batch; patients stay independent."""
    return {k: np.concatenate([s[k] for s in scens]) for k in scens[0]}


# ---------------------------------------------------------------- run configurations
# A configuration names a world family, its seed blocks (registry: src/trackc_seed_registry.py, purpose
# stage2_dev) and a suffix for the output folders. "single" is the configuration of the runs of 2026-11-04
# (scripts 192 to 194); "single_rep" repeats it with fresh seeds; "mixture" uses the mixture family. The scripts take
# --config NAME and hand the dictionary to their worker processes (module globals set in main() do not reach
# spawned workers). Only "single_rep" and "mixture" add the 3- and 4-course variants of the standard arm.
CONFIGS = {
    "single": dict(family="single", suffix="", extra_variants=False, world_seed=20261801, patient_seed=20261802,
                   noise_c=1_810_000, d_seed=20261803, noise_d=1_820_000, e_eval_seed=20261804, e_sel_seed=20261805,
                   noise_e=1_830_000, i_seed=20261806, noise_i=1_840_000, j_seed=20261807),
    "single_rep": dict(family="single", suffix="_rep", extra_variants=True, world_seed=20261810, patient_seed=20261811,
                       noise_c=1_850_000, d_seed=20261812, noise_d=1_860_000, e_eval_seed=20261813,
                       e_sel_seed=20261814, noise_e=1_870_000, i_seed=20261815, noise_i=1_880_000, j_seed=20261816),
    "mixture": dict(family="mixture", suffix="_mixture", extra_variants=True, world_seed=20261820,
                    patient_seed=20261821, noise_c=1_890_000, d_seed=20261822, noise_d=1_890_100,
                    e_eval_seed=20261823, e_sel_seed=20261824, noise_e=1_890_200, i_seed=20261825,
                    noise_i=1_890_300, j_seed=20261826, k_seed=20261827, noise_k=1_890_400),
}
STAGE_DIRS = {"C": "trackc_stage2_C_calibration", "D": "trackc_stage2_D_negative_controls",
              "E": "trackc_stage2_E_regime_map", "I": "trackc_stage2_I_numerics",
              "J": "trackc_stage2_J_response_depth", "CD": "trackc_stage2_CD_analysis",
              "K": "trackc_stage2_K_stratified_anchor"}


def get_config(name):
    """A copy of the named configuration, with its name added."""
    if name not in CONFIGS:
        raise KeyError(f"unknown config {name!r}; choose one of {sorted(CONFIGS)}")
    return dict(CONFIGS[name], name=name)


def stage_dir(stage, cfg, smoke=False):
    """Output folder of one stage ('C', 'D', 'E', 'I', 'J', 'CD') for one configuration."""
    return PROJECT_ROOT / "output" / (STAGE_DIRS[stage] + cfg["suffix"] + ("_smoke" if smoke else ""))


# ---------------------------------------------------------------- stratified anchor (mixture family only)
# RTOG 0525, both arms pooled, PFS from random assignment (Gilbert 2013, PMC3816958; page opened 2026-11-04):
# tumours with a methylated MGMT promoter 8.7 months (95% CI 6.6 to 11.2), unmethylated 5.7 months (5.1 to 6.1);
# 245 methylated and 517 unmethylated tumours among the 762 patients with a result. Script 201 asks whether a mixture
# world gives these stratified medians when its responders are mapped to the methylated tumours. That mapping is an
# ASSUMPTION of the check; no published data say that MGMT status is the response class of this model.
MGMT_UNMETHYLATED_PFS_DAYS = (5.1 * DAYS_PER_MONTH, 6.1 * DAYS_PER_MONTH)
MGMT_METHYLATED_PFS_DAYS = (6.6 * DAYS_PER_MONTH, 11.2 * DAYS_PER_MONTH)
MGMT_COUNTS = (245, 762)                                              # methylated, with a result


def wilson_interval(k, n, z=1.959964):
    """Wilson score interval for a proportion k / n."""
    p = k / n
    den = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return float(centre - half), float(centre + half)


MGMT_METHYLATED_FRACTION_CI = wilson_interval(*MGMT_COUNTS)


def responder_mask(scen, w):
    """True for the patients of a mixture world that were drawn from the responder range (m in 6 to 20).
    The non-responder range tops out at 20 * exp(log_m_scale_non), which is below 6 on the axis (scale below 0.3)."""
    if "q" not in w:
        raise ValueError("responder_mask needs a mixture world (a world with a responder fraction q)")
    return scen["m"] > tc.RANGES["m"][1] * np.exp(w["log_m_scale_non"]) * (1.0 + 1e-12)


# ---------------------------------------------------------------- arms
def run_arm(scen, courses=STD_COURSES, scan_every=SCAN_EVERY, dose=1.0, endpoint="rano", seed=(0, 0)):
    """Open-loop arm: a course every cycle from cycle 0 until `courses` are given or progression is called.
    Returns (time to progression in days, event flag, courses used). A patient without an event is censored at
    the horizon. dose multiplies the drug term of the true dynamics (see tc.simulate)."""
    s = dict(scen)
    n = len(s["r"])
    if scan_every is not None:
        s["scan_every"] = np.full(n, int(scan_every))
    if dose != 1.0:
        s["dose"] = np.full(n, float(dose))
    out = tc.simulate(s, tc.Periodic(1), list(seed), endpoint=endpoint, budget=courses)
    return out["ttp"], ~out["capped"], out["courses"]


def pfs_stats(ttp, event, times=SCAN_TIMES):
    """Median PFS (step and interpolated) and P(progression-free after t) at the given days."""
    d = {"median_step": sv.km_median(ttp, event), "median_interp": sv.km_median_interp(ttp, event),
         "event_fraction": float(np.mean(event))}
    for t in times:
        d[f"S{t}"] = sv.km_survival(ttp, event, t)
    return d


def in_target(median_days):
    return TARGET_PFS_DAYS[0] <= median_days <= TARGET_PFS_DAYS[1]


# ---------------------------------------------------------------- negative controls (test D)
# Published comparisons that found little or no benefit from more drug. A calibrated world should not show a much
# larger benefit. Intervals are 95% CIs of the published hazard ratios for PFS (targets file, rows P).
LANDMARK_DAYS = 6 * tc.CYCLE_DAYS                     # scan at the end of cycle 6, before cycle 7 would start
DOSE_DENSE_HR = (0.87, 0.75, 1.00)                    # RTOG 0525 dose-dense vs standard, PFS
TWELVE_VS_SIX_HR = (0.99, 0.71, 1.38)                 # GEINO 14-01, randomised after cycle 6, PFS
DOSE_DENSE_MULT = 2.0                                 # ASSUMED ratio of drug exposure per cycle; label doses not verified


def hazard_ratio(t_ref, e_ref, t_new, e_new):
    """Cox hazard ratio (new arm relative to reference) for two arms; NaN if an arm has no event."""
    t = np.r_[t_ref, t_new]
    e = np.r_[e_ref, e_new]
    g = np.r_[np.zeros(len(t_ref), int), np.ones(len(t_new), int)]
    try:
        hr, ci = sv.cox_hr(t, e, g)
    except ValueError:
        return float("nan"), (float("nan"), float("nan"))
    return hr, ci


def dose_dense_control(scen, mult=DOSE_DENSE_MULT, seed=(0, 0), ref=None):
    """Standard arm against the same arm with the drug effect multiplied by `mult`, same patients.
    ref: the standard-arm result of run_arm(scen, seed=seed) if it is already at hand."""
    ttp0, ev0, _ = run_arm(scen, seed=seed) if ref is None else ref
    ttp1, ev1, _ = run_arm(scen, dose=mult, seed=seed)
    hr, ci = hazard_ratio(ttp0, ev0, ttp1, ev1)
    return {"hr": hr, "ci": ci, "median_standard": sv.km_median_interp(ttp0, ev0),
            "median_dose_dense": sv.km_median_interp(ttp1, ev1)}


def landmark_control(scen, seed=(0, 0), ref=None):
    """12 versus 6 courses among patients without a progression call at the day-168 scan; time from the landmark.
    ref: the standard-arm result of run_arm(scen, seed=seed) if it is already at hand."""
    ttp0, ev0, _ = run_arm(scen, courses=STD_COURSES, seed=seed) if ref is None else ref
    ttp1, ev1, _ = run_arm(scen, courses=2 * STD_COURSES, seed=seed)
    keep = ttp0 > LANDMARK_DAYS
    if not np.array_equal(keep, ttp1 > LANDMARK_DAYS):
        raise RuntimeError("the two arms must agree before the landmark")
    hr, ci = hazard_ratio(ttp0[keep] - LANDMARK_DAYS, ev0[keep], ttp1[keep] - LANDMARK_DAYS, ev1[keep])
    return {"hr": hr, "ci": ci, "n_landmark": int(keep.sum()), "n_total": int(keep.size)}


def control_flags(dd_hr, lm_hr):
    """Pass flags: the model hazard ratio lies inside the published 95% interval."""
    dd_ok = DOSE_DENSE_HR[1] <= dd_hr <= DOSE_DENSE_HR[2]
    lm_ok = TWELVE_VS_SIX_HR[1] <= lm_hr <= TWELVE_VS_SIX_HR[2]
    return bool(dd_ok), bool(lm_ok)
