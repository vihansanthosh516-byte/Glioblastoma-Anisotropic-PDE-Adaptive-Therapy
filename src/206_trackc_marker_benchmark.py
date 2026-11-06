#!/usr/bin/env python3
"""Track C candidate 1: fixed schedules conditioned on an imperfect binary marker.

SIMULATION ONLY. The binary marker is hypothetical and labels the simulator's responder class;
it is not an MGMT test. Policies receive marker results or scan readings, never true tumour
parameters. Selection and evaluation patients are disjoint. Every arm uses the same evaluation
patients within a world.

The marker policy selects one of the existing 56 fixed schedules for each observed marker group.
It selects schedules on the selection patients by averaging each patient's worse endpoint-specific
restricted utility: min(min(RMST_RANO_730, 730) - lambda * courses_RANO,
min(RMST_baseline_730, 730) - lambda * courses_baseline). Lambda is 0, 5, or 10 restricted-mean
days per delivered course. Exact ties prefer fewer mean courses, then earlier catalog order. A
marker group with fewer than 25 selection patients uses the global selection.

For each endpoint and tau (365, 730 days), output reports RMST, delivered courses, paired
worst-decile difference from the configuration's standard reference, and gap from the upper
concave frontier of all 56 fixed schedules at the same mean courses. Standard and dose-dense
proxies, Stage 2 scan rules, and three current particle-filter rules are comparators.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import csv
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_calibration as cal  # noqa: E402
import trackc_manifest as man  # noqa: E402
import trackc_policies as pol  # noqa: E402
import trackc_regime as rg  # noqa: E402
import trackc_sim as tc  # noqa: E402
from trackc_integrated_policy_utils import patient_robust_utility  # noqa: E402

SCRIPT = "src/206_trackc_marker_benchmark.py"
TEST_SCRIPT = "tests/test_trackc_marker_benchmark.py"
ENDPOINTS = {"3b_rano": "rano", "3a_baseline": "baseline"}
TAUS = (365.0, 730.0)
LAMBDAS = (0.0, 5.0, 10.0)  # assumed RMST days per delivered course
MARKER_LEVELS = (0.50, 0.70, 0.80, 0.90, 0.95, 1.00)
MIN_GROUP_SELECTION = 25
N_EVAL, N_SEL = 1000, 400
PF_PARTICLES = 500

SEED_CONFIGS = {
    "mixture": {
        "calibration_config": "mixture",
        "eval_seed": 20261880,
        "selection_seed": 20261881,
        "noise_eval": 1_892_500,
        "noise_selection": 1_892_501,
    },
    "mixture_mix": {
        "calibration_config": "mixture_mix",
        "eval_seed": 20261882,
        "selection_seed": 20261883,
        "noise_eval": 1_892_600,
        "noise_selection": 1_892_601,
    },
}

LIMITATIONS = [
    "The marker is hypothetical. It is a noisy label for a simulator-defined response class, not MGMT status or an MGMT assay.",
    "The responder class is assigned by the assumed mixture generator. No clinical marker accuracy is estimated here.",
    "The policy selector receives marker results and simulated schedule outcomes only. It does not receive r, m, f0, c, or other true tumour parameters.",
    "The existing particle-filter implementation receives scan readings and simulator-supplied measurement-noise scale; it does not receive true tumour dynamics.",
    "Scan noise follows the simulator's per-patient assumed measurement model. This benchmark does not sweep scan-noise levels.",
    "Lambda is an assumed RMST-day penalty per delivered course, not a dollar price or toxicity estimate.",
    "Selection uses the same fixed schedule library for all marker accuracy settings. The 56-schedule upper frontier uses evaluation outcomes, so it is a descriptive held-out upper envelope, not a deployable policy or proven optimum.",
    "The controlled worlds are Stage 2 development worlds selected by earlier assumed calibration checks. They do not represent a patient population.",
    "Worst-decile paired differences describe the lowest 10 percent of patient-level RMST differences within each simulated world; they are not clinical risk bounds.",
]

MANIFEST_PATHS = tuple(dict.fromkeys(man.CORE_CODE + (
    "src/trackc_integrated_policy_utils.py",
    "src/trackc_calibration.py",
    "src/trackc_regime.py",
    "src/trackc_policies.py",
    "src/trackc_survival.py",
    SCRIPT,
    TEST_SCRIPT,
)))


def marker_from_uniforms(responder, uniforms, sensitivity, specificity):
    """Return positive marker calls with the requested sensitivity and specificity.

    Shared uniforms across the sweep couple the marker draws and reduce Monte Carlo noise when
    comparing accuracy settings. The selector never receives the responder array.
    """
    responder = np.asarray(responder, bool)
    uniforms = np.asarray(uniforms, float)
    if responder.shape != uniforms.shape:
        raise ValueError("responder and uniforms must have the same shape")
    if not (0.0 <= sensitivity <= 1.0 and 0.0 <= specificity <= 1.0):
        raise ValueError("sensitivity and specificity must lie in [0, 1]")
    return np.where(responder, uniforms < sensitivity, uniforms < 1.0 - specificity)


def marker_uniforms(n, noise_key, world, partition):
    """Common random numbers for marker calls, separate for selection and evaluation patients."""
    rng = np.random.default_rng(np.random.SeedSequence([int(noise_key), int(world), int(partition)]))
    return rng.random(int(n))


def robust_training_scores(ttp_rano, courses_rano, ttp_baseline, courses_baseline,
                           patient_indices, lam, tau=730.0):
    """Score schedules by mean patient-wise worst-endpoint utility minus course cost.

    Inputs have shape (n_schedules, n_patients). For each patient, cap each endpoint's RMST,
    subtract its endpoint-specific course cost, then take the lower endpoint utility. Average
    those patient utilities. Returned course use is the endpoint-averaged mean, used for ties.
    """
    ix = np.asarray(patient_indices, int)
    if ix.size == 0:
        raise ValueError("cannot score an empty selection group")
    tr = np.asarray(ttp_rano, float)[:, ix]
    tb = np.asarray(ttp_baseline, float)[:, ix]
    cr = np.asarray(courses_rano, float)[:, ix]
    cb = np.asarray(courses_baseline, float)[:, ix]
    utility = patient_robust_utility(
        tb, tr, cb, cr, tau=tau, course_penalty=lam
    )
    course_cost = 0.5 * (cr + cb)
    return utility.mean(axis=1), course_cost.mean(axis=1)


def _best_schedule_index(scores, mean_courses):
    """Maximize utility; break exact ties by fewer courses, then catalog order."""
    scores = np.asarray(scores, float)
    mean_courses = np.asarray(mean_courses, float)
    if scores.ndim != 1 or scores.shape != mean_courses.shape or scores.size == 0:
        raise ValueError("scores and mean_courses must be nonempty vectors of equal shape")
    if not np.isfinite(scores).all() or not np.isfinite(mean_courses).all():
        raise ValueError("scores and mean_courses must be finite")
    return int(np.lexsort((np.arange(scores.size), mean_courses, -scores))[0])


def select_global_fixed_schedule(selection_outcomes, schedule_names, lam):
    """Select one fixed schedule on all selection patients for a global comparator."""
    tr, cr = selection_outcomes["3b_rano"]
    tb, cb = selection_outcomes["3a_baseline"]
    n = np.asarray(tr).shape[1]
    scores, courses = robust_training_scores(tr, cr, tb, cb, np.arange(n), lam)
    idx = _best_schedule_index(scores, courses)
    return {
        "schedule_index": idx,
        "schedule": str(schedule_names[idx]),
        "training_robust_utility": float(scores[idx]),
        "training_mean_courses": float(courses[idx]),
    }


def choose_marker_schedules(marker_selection, selection_outcomes, schedule_names, lam,
                            min_group_n=MIN_GROUP_SELECTION):
    """Select marker-negative and marker-positive schedules using selection outcomes only.

    selection_outcomes maps endpoint labels to (schedule-major TTP, schedule-major courses).
    Groups smaller than min_group_n use the schedule selected from the full selection set.
    """
    marker_selection = np.asarray(marker_selection, bool)
    n = marker_selection.size
    if len(schedule_names) == 0:
        raise ValueError("schedule catalog must not be empty")
    for endpoint in ("3b_rano", "3a_baseline"):
        if endpoint not in selection_outcomes:
            raise ValueError(f"missing selection endpoint {endpoint}")
        ttp, courses = selection_outcomes[endpoint]
        if np.asarray(ttp).shape != (len(schedule_names), n):
            raise ValueError(f"unexpected TTP shape for {endpoint}")
        if np.asarray(courses).shape != (len(schedule_names), n):
            raise ValueError(f"unexpected course shape for {endpoint}")

    tr, cr = selection_outcomes["3b_rano"]
    tb, cb = selection_outcomes["3a_baseline"]
    all_ix = np.arange(n)
    global_scores, global_courses = robust_training_scores(tr, cr, tb, cb, all_ix, lam)
    global_idx = _best_schedule_index(global_scores, global_courses)

    chosen = {}
    for value, label in ((False, "negative"), (True, "positive")):
        ix = np.flatnonzero(marker_selection == value)
        fallback = ix.size < int(min_group_n)
        if fallback:
            idx, scores, costs = global_idx, global_scores, global_courses
        else:
            scores, costs = robust_training_scores(tr, cr, tb, cb, ix, lam)
            idx = _best_schedule_index(scores, costs)
        chosen[label] = {
            "schedule_index": idx,
            "schedule": str(schedule_names[idx]),
            "selection_n": int(ix.size),
            "used_global_fallback": bool(fallback),
            "training_robust_utility": float(scores[idx]),
            "training_mean_courses": float(costs[idx]),
            "global_schedule_index": global_idx,
            "global_schedule": str(schedule_names[global_idx]),
            "global_training_robust_utility": float(global_scores[global_idx]),
        }
    return chosen


def fit_marker_schedule_policy(marker_selection, selection_outcomes, schedule_names, lam,
                               min_group_n=MIN_GROUP_SELECTION):
    """Public pure fitting API for a shared selection stream and precomputed schedule outcomes."""
    return choose_marker_schedules(marker_selection, selection_outcomes, schedule_names, lam, min_group_n)


def marker_policy_outcomes(marker_evaluation, positive_schedule_index, negative_schedule_index,
                           schedule_ttp, schedule_courses):
    """Apply a two-schedule marker rule to evaluation patients; no latent class is used."""
    marker_evaluation = np.asarray(marker_evaluation, bool)
    ttp = np.asarray(schedule_ttp, float)
    courses = np.asarray(schedule_courses, float)
    if ttp.shape != courses.shape or ttp.ndim != 2 or ttp.shape[1] != marker_evaluation.size:
        raise ValueError("schedule outcomes must have shape (n_schedules, n_evaluation_patients)")
    schedule_index = np.where(marker_evaluation, int(positive_schedule_index), int(negative_schedule_index))
    patient_index = np.arange(marker_evaluation.size)
    return ttp[schedule_index, patient_index], courses[schedule_index, patient_index], schedule_index


def apply_marker_schedule_policy(marker_evaluation, fitted_policy, schedule_ttp, schedule_courses):
    """Public pure application API; caller supplies marker calls and common schedule outcomes."""
    return marker_policy_outcomes(
        marker_evaluation,
        fitted_policy["positive"]["schedule_index"],
        fitted_policy["negative"]["schedule_index"],
        schedule_ttp,
        schedule_courses,
    )


def upper_hull(points):
    """Upper concave envelope of (mean courses, RMST) points."""
    best = {}
    for x, y in points:
        x, y = float(x), float(y)
        x = round(x, 9)
        if x not in best or y > best[x]:
            best[x] = y
    hull = []
    for point in sorted(best.items()):
        while len(hull) >= 2 and (
            (hull[-1][1] - hull[-2][1]) * (point[0] - hull[-2][0])
            <= (point[1] - hull[-2][1]) * (hull[-1][0] - hull[-2][0])
        ):
            hull.pop()
        hull.append(point)
    return hull


def frontier_at(hull, mean_courses):
    """Return fixed-frontier RMST at mean_courses, or NaN outside observed fixed range."""
    c = float(mean_courses)
    if len(hull) < 2 or c < hull[0][0] - 1e-9 or c > hull[-1][0] + 1e-9:
        return float("nan")
    return float(np.interp(c, [point[0] for point in hull], [point[1] for point in hull]))


def comparison_policy_specs():
    """Return established scan policies plus selected current particle-filter rules."""
    arms = {
        name: (factory, budget, "scan_policy", 1.0)
        for name, (factory, budget) in rg.policy_arms().items()
    }
    arms.update({
        # Existing wide-prior risk containment and later-trigger candidates from scripts 167/173.
        "pf_wide_containment": (
            lambda: tc.RiskContainment(0.10, 2, 0.0, True, n_particles=PF_PARTICLES),
            tc.BUDGET, "particle_filter", 1.0,
        ),
        "pf_tail_trigger_candidate": (
            lambda: tc.RiskContainment(0.05, 4, 0.0, True, n_particles=PF_PARTICLES, front=4),
            tc.BUDGET, "particle_filter", 1.0,
        ),
        "pf_holiday_containment": (
            lambda: tc.HolidayContainment(0.10, 2, front=0, max_run=2, n_particles=PF_PARTICLES),
            tc.BUDGET, "particle_filter", 1.0,
        ),
    })
    return arms


def read_controlled(config):
    """Read worlds retained by the matching Stage 2 negative-control run."""
    cfg = cal.get_config(SEED_CONFIGS[config]["calibration_config"])
    source = cal.stage_dir("D", cfg) / "controlled_worlds.csv"
    with source.open(newline="", encoding="utf-8") as handle:
        return [
            {key: (int(value) if key == "world" else float(value)) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]


def _reference_name(cfg):
    return "std_capmix" if cfg.get("std_courses", cal.STD_COURSES) == "mix36" else "soc6"


def _world_policy_results(evaluation, cfg, world):
    """Run comparator policies on one common evaluation cohort for each endpoint."""
    policy_specs = comparison_policy_specs()
    result = {}
    reference = _reference_name(cfg)
    caps = cal.std_caps(cfg["std_courses"], len(evaluation["r"]))
    for endpoint_label, endpoint in ENDPOINTS.items():
        seed = [SEED_CONFIGS[cfg["benchmark_config"]]["noise_eval"], int(world["world"])]
        endpoint_arms = {}
        for name, (factory, budget, family, dose_mult) in policy_specs.items():
            ttp, courses = rg.run_policy(evaluation, factory, budget, endpoint, seed)
            endpoint_arms[name] = (ttp, courses, family, dose_mult, "")
        if reference == "std_capmix":
            ttp, _, courses = cal.run_arm(
                evaluation, courses=caps, scan_every=rg.SCAN_EVERY,
                endpoint=endpoint, seed=seed,
            )
            endpoint_arms["std_capmix"] = (ttp, courses, "standard", 1.0, "")
        dense_ttp, _, dense_courses = cal.run_arm(
            evaluation, courses=caps, scan_every=rg.SCAN_EVERY,
            dose=cal.DOSE_DENSE_MULT, endpoint=endpoint, seed=seed,
        )
        endpoint_arms["dose_dense_proxy"] = (
            dense_ttp, dense_courses, "dose_dense", cal.DOSE_DENSE_MULT, ""
        )
        if reference == "std_capmix":
            soc6_ttp, _, soc6_courses = cal.run_arm(
                evaluation, courses=cal.STD_COURSES, scan_every=rg.SCAN_EVERY,
                endpoint=endpoint, seed=seed,
            )
            endpoint_arms["dose_dense_soc6"] = (
                soc6_ttp, soc6_courses, "dose_dense", cal.DOSE_DENSE_MULT, ""
            )
        result[endpoint_label] = endpoint_arms
    return result


def _metric_row(world_id, endpoint_label, tau, arm, family, ttp, courses, dose_multiplier,
                reference_ttp, reference_courses, reference_name, frontier, metadata=None):
    ttp = np.asarray(ttp, float)
    courses = np.asarray(courses, float)
    stats = rg.paired_stats(ttp, reference_ttp, tau)
    rmst_value = rg.rmst(ttp, tau)
    mean_courses = float(courses.mean())
    frontier_basis = "excluded_dose_intensity_differs" if float(dose_multiplier) != 1.0 else "delivered_courses"
    frontier_rmst = (
        float("nan") if float(dose_multiplier) != 1.0
        else frontier_at(frontier, mean_courses)
    )
    metadata = metadata or {}
    row = {
        "world": int(world_id),
        "endpoint": endpoint_label,
        "tau": int(tau),
        "arm": arm,
        "family": family,
        "reference_standard": reference_name,
        "selection_lambda": metadata.get("selection_lambda", ""),
        "marker_sensitivity": metadata.get("marker_sensitivity", ""),
        "marker_specificity": metadata.get("marker_specificity", ""),
        "marker_positive_eval_prevalence": metadata.get("marker_positive_eval_prevalence", ""),
        "selected_marker_positive_schedule": metadata.get("selected_marker_positive_schedule", ""),
        "selected_marker_negative_schedule": metadata.get("selected_marker_negative_schedule", ""),
        "rmst_days": rmst_value,
        "diff_vs_standard_days": stats["diff"],
        "ci95_lo_days": stats["ci_lo"],
        "ci95_hi_days": stats["ci_hi"],
        "paired_worst10_diff_vs_standard_days": stats["regret_cvar"],
        "fraction_worse_than_standard": stats["frac_worse"],
        "mean_courses": mean_courses,
        "mean_dose_equivalent_courses": float(mean_courses * float(dose_multiplier)),
        "mean_courses_diff_vs_standard": float(mean_courses - np.asarray(reference_courses, float).mean()),
        "fixed_frontier_rmst_at_same_mean_courses_days": frontier_rmst,
        "gap_vs_fixed_frontier_days": float(rmst_value - frontier_rmst) if np.isfinite(frontier_rmst) else float("nan"),
        "frontier_in_range": bool(np.isfinite(frontier_rmst)),
        "frontier_comparison_basis": frontier_basis,
        "course_cost_lam0_days": 0.0 * mean_courses,
        "course_cost_lam5_days": 5.0 * mean_courses,
        "course_cost_lam10_days": 10.0 * mean_courses,
        "dose_equivalent_cost_lam0_days": 0.0 * mean_courses * float(dose_multiplier),
        "dose_equivalent_cost_lam5_days": 5.0 * mean_courses * float(dose_multiplier),
        "dose_equivalent_cost_lam10_days": 10.0 * mean_courses * float(dose_multiplier),
        "net_utility_lam0_days": rmst_value,
        "net_utility_lam5_days": rmst_value - 5.0 * mean_courses,
        "net_utility_lam10_days": rmst_value - 10.0 * mean_courses,
        "net_dose_equivalent_utility_lam0_days": rmst_value,
        "net_dose_equivalent_utility_lam5_days": rmst_value - 5.0 * mean_courses * float(dose_multiplier),
        "net_dose_equivalent_utility_lam10_days": rmst_value - 10.0 * mean_courses * float(dose_multiplier),
        "schedule": metadata.get("schedule", ""),
        "dose_multiplier": float(dose_multiplier),
    }
    return row


def evaluate_world(world, config, n_eval=N_EVAL, n_sel=N_SEL):
    """Evaluate all fixed, marker, scan-policy, and particle-filter arms in one controlled world."""
    seed_cfg = SEED_CONFIGS[config]
    cfg = cal.get_config(seed_cfg["calibration_config"])
    cfg["benchmark_config"] = config
    world_id = int(world["world"])
    evaluation = cal.draw_patients(world, int(n_eval), [seed_cfg["eval_seed"], world_id])
    selection = cal.draw_patients(world, int(n_sel), [seed_cfg["selection_seed"], world_id])
    responder_eval = cal.responder_mask(evaluation, world)
    responder_selection = cal.responder_mask(selection, world)
    uniforms_eval = marker_uniforms(n_eval, seed_cfg["noise_eval"], world_id, 1)
    uniforms_selection = marker_uniforms(n_sel, seed_cfg["noise_selection"], world_id, 2)

    schedules, schedule_names = rg.candidate_schedules()
    selection_outcomes = {}
    evaluation_schedule_outcomes = {}
    for endpoint_label, endpoint in ENDPOINTS.items():
        selection_outcomes[endpoint_label] = rg.evaluate_schedules(selection, schedules, endpoint)
        evaluation_schedule_outcomes[endpoint_label] = rg.evaluate_schedules(evaluation, schedules, endpoint)
    global_fixed = {
        lam: select_global_fixed_schedule(selection_outcomes, schedule_names, lam)
        for lam in LAMBDAS
    }

    comparator_outcomes = _world_policy_results(evaluation, cfg, world)
    reference_name = _reference_name(cfg)
    rows = []
    marker_rows = []
    for lam, selected in global_fixed.items():
        marker_rows.append({
            "world": world_id,
            "arm": f"best_fixed_lam{lam:g}",
            "selection_lambda": float(lam),
            "selected_global_schedule": selected["schedule"],
            "global_training_robust_utility": selected["training_robust_utility"],
            "global_training_mean_courses": selected["training_mean_courses"],
            "n_selection": int(n_sel),
        })
    marker_outcomes_by_arm = {}
    marker_meta_by_arm = {}

    # Open-loop 56-schedule family and established comparators share the same evaluation patients.
    for endpoint_label in ENDPOINTS:
        ttp_sched, courses_sched = evaluation_schedule_outcomes[endpoint_label]
        endpoint_arms = comparator_outcomes[endpoint_label]
        if reference_name not in endpoint_arms:
            raise RuntimeError(f"reference arm {reference_name} was not evaluated")
        ref_ttp, ref_courses = endpoint_arms[reference_name][:2]
        for tau in TAUS:
            fixed_points = [
                (float(np.mean(courses_sched[i])), rg.rmst(ttp_sched[i], tau))
                for i in range(len(schedule_names))
            ]
            if reference_name == "std_capmix":
                t_ref_capmix, c_ref_capmix = endpoint_arms[reference_name][:2]
                fixed_points.append((float(np.mean(c_ref_capmix)), rg.rmst(t_ref_capmix, tau)))
            frontier = upper_hull(fixed_points)
            for i, schedule_name in enumerate(schedule_names):
                rows.append(_metric_row(
                    world_id, endpoint_label, tau, f"schedule_{schedule_name}", "fixed_schedule",
                    ttp_sched[i], courses_sched[i], 1.0, ref_ttp, ref_courses, reference_name,
                    frontier, {"schedule": schedule_name},
                ))
            for lam, selected in global_fixed.items():
                i = selected["schedule_index"]
                rows.append(_metric_row(
                    world_id, endpoint_label, tau, f"best_fixed_lam{lam:g}", "selected_fixed_schedule",
                    ttp_sched[i], courses_sched[i], 1.0, ref_ttp, ref_courses, reference_name,
                    frontier, {"selection_lambda": lam, "schedule": selected["schedule"]},
                ))
            for arm_name, (ttp, courses, family, dose_mult, _) in endpoint_arms.items():
                rows.append(_metric_row(
                    world_id, endpoint_label, tau, arm_name, family, ttp, courses, dose_mult,
                    ref_ttp, ref_courses, reference_name, frontier,
                ))

            # Build marker-conditioned outcomes from one pre-generated marker draw per patient.
            for sensitivity in MARKER_LEVELS:
                for specificity in MARKER_LEVELS:
                    marker_eval = marker_from_uniforms(responder_eval, uniforms_eval, sensitivity, specificity)
                    marker_sel = marker_from_uniforms(responder_selection, uniforms_selection, sensitivity, specificity)
                    for lam in LAMBDAS:
                        tag = f"marker_se{sensitivity:.2f}_sp{specificity:.2f}_lam{lam:g}"
                        if tag not in marker_outcomes_by_arm:
                            chosen = fit_marker_schedule_policy(marker_sel, selection_outcomes, schedule_names, lam)
                            marker_meta_by_arm[tag] = {
                                "sensitivity": float(sensitivity),
                                "specificity": float(specificity),
                                "lambda": float(lam),
                                "chosen": chosen,
                                "n_positive_selection": int(marker_sel.sum()),
                                "n_negative_selection": int((~marker_sel).sum()),
                                "n_positive_evaluation": int(marker_eval.sum()),
                                "n_negative_evaluation": int((~marker_eval).sum()),
                                "marker_positive_eval_prevalence": float(marker_eval.mean()),
                                "selection_sensitivity_observed": float(marker_sel[responder_selection].mean()),
                                "selection_specificity_observed": float((~marker_sel[~responder_selection]).mean()),
                                "evaluation_sensitivity_observed": float(marker_eval[responder_eval].mean()),
                                "evaluation_specificity_observed": float((~marker_eval[~responder_eval]).mean()),
                                "responder_prevalence_selection": float(responder_selection.mean()),
                                "responder_prevalence_evaluation": float(responder_eval.mean()),
                            }
                            marker_rows.append({
                                "world": world_id,
                                "marker_sensitivity": float(sensitivity),
                                "marker_specificity": float(specificity),
                                "selection_lambda": float(lam),
                                "n_positive_selection": int(marker_sel.sum()),
                                "n_negative_selection": int((~marker_sel).sum()),
                                "n_positive_evaluation": int(marker_eval.sum()),
                                "n_negative_evaluation": int((~marker_eval).sum()),
                                "marker_positive_eval_prevalence": float(marker_eval.mean()),
                                "responder_prevalence_selection": float(responder_selection.mean()),
                                "responder_prevalence_evaluation": float(responder_eval.mean()),
                                "selection_sensitivity_observed": float(marker_sel[responder_selection].mean()),
                                "selection_specificity_observed": float((~marker_sel[~responder_selection]).mean()),
                                "evaluation_sensitivity_observed": float(marker_eval[responder_eval].mean()),
                                "evaluation_specificity_observed": float((~marker_eval[~responder_eval]).mean()),
                                "selected_marker_positive_schedule": chosen["positive"]["schedule"],
                                "selected_marker_negative_schedule": chosen["negative"]["schedule"],
                                "positive_group_used_global_fallback": chosen["positive"]["used_global_fallback"],
                                "negative_group_used_global_fallback": chosen["negative"]["used_global_fallback"],
                                "positive_group_training_robust_utility": chosen["positive"]["training_robust_utility"],
                                "negative_group_training_robust_utility": chosen["negative"]["training_robust_utility"],
                                "positive_group_training_n": chosen["positive"]["selection_n"],
                                "negative_group_training_n": chosen["negative"]["selection_n"],
                            })
                            # Store both endpoint outcomes now; all later tau rows reuse same evaluation patients.
                            marker_outcomes_by_arm[tag] = {}
                            for ep in ENDPOINTS:
                                ep_ttp, ep_courses = evaluation_schedule_outcomes[ep]
                                marker_outcomes_by_arm[tag][ep] = apply_marker_schedule_policy(
                                    marker_eval, chosen, ep_ttp, ep_courses
                                )
                        metadata = marker_meta_by_arm[tag]
                        ttp_marker, courses_marker, _ = marker_outcomes_by_arm[tag][endpoint_label]
                        rows.append(_metric_row(
                            world_id, endpoint_label, tau, tag, "marker_conditioned_fixed",
                            ttp_marker, courses_marker, 1.0, ref_ttp, ref_courses,
                            reference_name, frontier,
                            {
                                "selection_lambda": lam,
                                "marker_sensitivity": sensitivity,
                                "marker_specificity": specificity,
                                "marker_positive_eval_prevalence": metadata["marker_positive_eval_prevalence"],
                                "selected_marker_positive_schedule": metadata["chosen"]["positive"]["schedule"],
                                "selected_marker_negative_schedule": metadata["chosen"]["negative"]["schedule"],
                            },
                        ))
    return rows, marker_rows


def _worker(args):
    return evaluate_world(*args)


def _protocol(config, n_eval=N_EVAL, n_sel=N_SEL):
    seed = SEED_CONFIGS[config]
    cfg = cal.get_config(seed["calibration_config"])
    return {
        "script": SCRIPT,
        "label": "SIMULATION ONLY; hypothetical marker benchmark; no clinical meaning",
        "config": config,
        "world_source": str(cal.stage_dir("D", cfg) / "controlled_worlds.csv"),
        "selection_seed": seed["selection_seed"],
        "evaluation_seed": seed["eval_seed"],
        "selection_marker_noise_key": seed["noise_selection"],
        "evaluation_marker_and_scan_noise_key": seed["noise_eval"],
        "n_selection_per_world": int(n_sel),
        "n_evaluation_per_world": int(n_eval),
        "se_grid": list(MARKER_LEVELS),
        "sp_grid": list(MARKER_LEVELS),
        "marker_positive_definition": "simulator responder-class membership (m-range); hypothetical marker, not MGMT",
        "selection_rule": "for each observed marker group, choose one of the existing 56 schedules by min(mean RMST_RANO_730, mean RMST_baseline_730) - lambda * mean delivered courses; course count is averaged over endpoint-specific simulations",
        "selection_lambda_days_per_course": list(LAMBDAS),
        "selection_tie_break": "prefer fewer mean delivered courses, then candidate catalog order",
        "minimum_marker_group_selection_n": MIN_GROUP_SELECTION,
        "evaluation_endpoints": dict(ENDPOINTS),
        "evaluation_tau_days": list(TAUS),
        "standard_reference": _reference_name(cfg),
        "fixed_frontier": "upper concave envelope of all 56 fixed schedules; mixed-cap standard point is also included when present; dose-dense proxy excluded because dose intensity differs",
        "paired_tail": "mean of lowest 10 percent patient-level RMST differences from configuration standard, using same evaluation patients",
        "dose_dense_proxy": {"course_cap": cfg["std_courses"], "dose_multiplier": cal.DOSE_DENSE_MULT},
        "course_cost": "lambda penalises delivered courses; dose-equivalent cost utilities also reported for dose-dense arms; no dollar price or toxicity model",
        "threshold_comparators": {
            "raw_zhang_1.0_0.5": "Zhang/Gatenby-style baseline on/off threshold in the existing Track C policy catalog",
            "raw_0.6_0.6": "additional raw-observation threshold comparator",
        },
        "particle_filter_comparators": ["pf_wide_containment", "pf_tail_trigger_candidate", "pf_holiday_containment"],
        "particle_filter_particles": PF_PARTICLES,
        "limits": LIMITATIONS,
    }


def _clean_json(value):
    if isinstance(value, dict):
        return {str(key): _clean_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean_json(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def _summary(rows):
    groups = {}
    for row in rows:
        key = (row["arm"], row["endpoint"], int(row["tau"]), row["selection_lambda"],
               row["marker_sensitivity"], row["marker_specificity"])
        groups.setdefault(key, []).append(row)
    summaries = []
    fields = (
        "rmst_days", "diff_vs_standard_days", "paired_worst10_diff_vs_standard_days",
        "mean_courses", "mean_dose_equivalent_courses", "gap_vs_fixed_frontier_days",
        "net_utility_lam5_days", "net_dose_equivalent_utility_lam5_days",
    )
    for key, block in sorted(groups.items(), key=lambda item: tuple(str(x) for x in item[0])):
        entry = {
            "arm": key[0], "endpoint": key[1], "tau": key[2],
            "selection_lambda": key[3], "marker_sensitivity": key[4], "marker_specificity": key[5],
            "n_worlds": len(block),
        }
        for field in fields:
            vals = np.asarray([row[field] for row in block], float)
            vals = vals[np.isfinite(vals)]
            entry[field] = {
                "n_finite": int(vals.size),
                "median": float(np.median(vals)) if vals.size else float("nan"),
                "mean": float(vals.mean()) if vals.size else float("nan"),
                "q25_q75": [float(x) for x in np.percentile(vals, [25, 75])] if vals.size else [float("nan"), float("nan")],
            }
        summaries.append(entry)
    return summaries


def _write_csv(path, rows):
    if not rows:
        raise ValueError(f"refusing to write empty result file {path}")
    fields = list(dict.fromkeys(field for row in rows for field in row))
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: ("" if isinstance(value, float) and not np.isfinite(value) else value)
                             for key, value in row.items()})


def run_config(config, workers):
    """Write a manifest first, then per-world CSVs and read-back-friendly JSON summary."""
    seed = SEED_CONFIGS[config]
    cfg = cal.get_config(seed["calibration_config"])
    out_dir = PROJECT_ROOT / "output" / f"trackc_stage2_F_marker_benchmark_{config}"
    if out_dir.exists() and any(out_dir.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty output directory {out_dir}")
    worlds = read_controlled(config)
    protocol = _protocol(config)
    manifest = man.build_manifest(
        f"stage2_F_marker_benchmark_{config}", protocol, purpose="stage2_dev",
        scenario=(seed["eval_seed"], seed["selection_seed"]),
        noise=(seed["noise_eval"], seed["noise_selection"]),
        script=SCRIPT, require_clean=True, paths=MANIFEST_PATHS, root=PROJECT_ROOT,
    )
    man.write_manifest(out_dir, manifest)
    start = time.time()
    jobs = [(world, config, N_EVAL, N_SEL) for world in worlds]
    if int(workers) > 1:
        with ProcessPoolExecutor(max_workers=int(workers)) as pool:
            outcomes = list(pool.map(_worker, jobs))
    else:
        outcomes = [_worker(job) for job in jobs]
    arm_rows = [row for result, _ in outcomes for row in result]
    selection_rows = [row for _, result in outcomes for row in result]
    _write_csv(out_dir / "arms_by_world.csv", arm_rows)
    _write_csv(out_dir / "marker_selection_by_world.csv", selection_rows)
    result_doc = {
        "manifest": "manifest.json",
        "protocol": protocol,
        "n_controlled_worlds": len(worlds),
        "arm_rows": len(arm_rows),
        "marker_selection_rows": len(selection_rows),
        "elapsed_seconds": time.time() - start,
        "results_csv": "arms_by_world.csv",
        "marker_selection_csv": "marker_selection_by_world.csv",
        "summaries": _summary(arm_rows),
        "limits": LIMITATIONS,
    }
    (out_dir / "results.json").write_text(
        json.dumps(_clean_json(result_doc), indent=2, sort_keys=False) + "\n", encoding="utf-8"
    )
    return out_dir


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", choices=("mixture", "mixture_mix", "both"), required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    configs = tuple(SEED_CONFIGS) if args.config == "both" else (args.config,)
    outputs = [str(run_config(config, args.workers)) for config in configs]
    print(json.dumps({"outputs": outputs}, indent=2))


if __name__ == "__main__":
    main()
