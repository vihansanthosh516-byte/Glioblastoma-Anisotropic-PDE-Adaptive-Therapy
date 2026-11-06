"""Paired development harness for Track C signal-to-treatment policies.

SIMULATION ONLY. The signal-to-resistance and drug-response links are assumptions.
This module scores no patient data and writes no files. Callers must register fresh
development seeds and freeze the run protocol before evaluating policy outcomes.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import importlib.util
from pathlib import Path
import sys
from typing import Mapping, Sequence

import numpy as np

import trackc_regime as rg
import trackc_sim as tc
from trackc_integrated_signal_model import (
    F0_REFERENCE_BOUNDS,
    SignalModelParameters,
    SignalObservations,
)


def _load_source_module(module_name: str, filename: str):
    if module_name in sys.modules:
        return sys.modules[module_name]
    source = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(module_name, source)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {filename} from {source}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


@dataclass(frozen=True)
class PolicyBenchmarkSettings:
    scan_every: int = 3
    course_budget: int = tc.BUDGET
    tau_days: float = 730.0
    course_penalty: float = 5.0
    signal_f0_bounds: tuple[float, float] = tuple(map(float, tc.RANGES["f0"][:2]))
    n_f0_nodes: int = 256
    mpc_particles: int = 500
    prediction_particles: int = 16
    forecast_cycles: int = 6
    q_episodes: int = 50
    q_learning_rate: float = 0.1
    q_epsilon: float = 0.25
    q_epsilon_decay: float = 0.98
    q_min_epsilon: float = 0.02

    def __post_init__(self):
        _validate_int(self.scan_every, "scan_every", 1, 4)
        _validate_int(self.course_budget, "course_budget", 1, tc.BUDGET)
        _validate_int(self.n_f0_nodes, "n_f0_nodes", 8, 100_000)
        _validate_int(self.mpc_particles, "mpc_particles", 2, 100_000)
        _validate_int(self.prediction_particles, "prediction_particles", 1, self.mpc_particles)
        _validate_int(self.forecast_cycles, "forecast_cycles", self.scan_every, 12)
        _validate_int(self.q_episodes, "q_episodes", 1, 1_000_000)
        if self.forecast_cycles % self.scan_every:
            raise ValueError("forecast_cycles must be a multiple of scan_every")
        if not np.isfinite(self.tau_days) or not 0 < float(self.tau_days) <= tc.HORIZON:
            raise ValueError("tau_days must be finite and in (0, simulator horizon]")
        if not np.isfinite(self.course_penalty) or self.course_penalty < 0:
            raise ValueError("course_penalty must be finite and non-negative")
        try:
            f0_bounds = tuple(map(float, self.signal_f0_bounds))
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("signal_f0_bounds must contain two finite bounds") from exc
        if (
            len(f0_bounds) != 2
            or not np.isfinite(f0_bounds).all()
            or f0_bounds[0] <= 0.0
            or f0_bounds[1] <= f0_bounds[0]
            or f0_bounds[0] < F0_REFERENCE_BOUNDS[0]
            or f0_bounds[1] > F0_REFERENCE_BOUNDS[1]
        ):
            raise ValueError("signal_f0_bounds must be ordered inside the signal reference bounds")
        object.__setattr__(self, "signal_f0_bounds", f0_bounds)
        if not np.isfinite(self.q_learning_rate) or not 0 < self.q_learning_rate <= 1:
            raise ValueError("q_learning_rate must be in (0, 1]")
        if not 0 <= self.q_min_epsilon <= self.q_epsilon <= 1:
            raise ValueError("Q epsilon values must satisfy 0 <= min <= initial <= 1")
        if not 0 <= self.q_epsilon_decay <= 1:
            raise ValueError("q_epsilon_decay must be in [0, 1]")


def _validate_int(value, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer")
    try:
        numeric = float(value)
        integer = int(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if not np.isfinite(numeric) or numeric != integer or not minimum <= integer <= maximum:
        raise ValueError(f"{name} must be an integer in [{minimum}, {maximum}]")
    return integer


def _validate_seed(seed, name: str) -> tuple[int, ...]:
    raw = np.atleast_1d(seed).tolist()
    values = []
    for value in raw:
        if isinstance(value, (bool, np.bool_)):
            raise ValueError(f"{name} must contain non-negative integers")
        try:
            integer = int(value)
            numeric = float(value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"{name} must contain non-negative integers") from exc
        if not np.isfinite(numeric) or integer != numeric or integer < 0:
            raise ValueError(f"{name} must contain non-negative integers")
        values.append(integer)
    if not values:
        raise ValueError(f"{name} must not be empty")
    return tuple(values)


def _validate_scenarios(scenarios: Mapping, name: str) -> int:
    if not isinstance(scenarios, Mapping) or "r" not in scenarios or "sigma" not in scenarios:
        raise ValueError(f"{name} must be a scenario mapping with r and sigma")
    r = np.asarray(scenarios["r"])
    sigma = np.asarray(scenarios["sigma"], dtype=float)
    if r.ndim != 1 or r.size == 0 or sigma.shape != r.shape:
        raise ValueError(f"{name} must contain matching nonempty one-dimensional r and sigma arrays")
    if not np.isfinite(sigma).all() or (sigma <= 0).any():
        raise ValueError(f"{name} sigma values must be finite and positive")
    for key, value in scenarios.items():
        array = np.asarray(value)
        if array.ndim > 0 and array.shape[0] != r.size:
            raise ValueError(f"{name} scenario field {key!r} has a different patient count")
    return int(r.size)


def _copy_scenarios(scenarios: Mapping) -> dict:
    return {key: np.array(value, copy=True) if isinstance(value, np.ndarray) else deepcopy(value)
            for key, value in scenarios.items()}


def _repeat_scenarios(scenarios: Mapping, n_repeats: int, n_patients: int) -> dict:
    repeated = {}
    for key, value in scenarios.items():
        array = np.asarray(value)
        if array.ndim > 0 and array.shape[0] == n_patients:
            repeated[key] = np.tile(array, (n_repeats,) + (1,) * (array.ndim - 1))
        else:
            repeated[key] = deepcopy(value)
    return repeated


class FixedSchedulePolicy:
    """Expose fixed boolean cycle schedules through the simulator policy interface."""

    uses_filter = False

    def __init__(self, schedules: np.ndarray):
        value = np.asarray(schedules, dtype=bool)
        if value.ndim == 1:
            value = value[None, :]
        if value.ndim != 2 or value.shape[1] != tc.N_CYCLES or value.shape[0] == 0:
            raise ValueError(f"schedules must have shape (patients, {tc.N_CYCLES}) or ({tc.N_CYCLES},)")
        if np.any(value.sum(axis=1) > tc.BUDGET):
            raise ValueError(f"no fixed schedule may exceed the {tc.BUDGET}-course simulator budget")
        self.schedules = value.copy()

    def reset(self, scen, rng):
        del rng
        n_patients = len(np.asarray(scen["r"]))
        if self.schedules.shape[0] == 1:
            self._active = np.broadcast_to(self.schedules, (n_patients, tc.N_CYCLES))
        elif self.schedules.shape[0] == n_patients:
            self._active = self.schedules
        else:
            raise ValueError("fixed schedule rows must be one or match the scenario patient count")

    def act(self, cycle: int, readings) -> np.ndarray:
        del readings
        if isinstance(cycle, bool) or int(cycle) != cycle or not 0 <= int(cycle) < tc.N_CYCLES:
            raise ValueError("cycle must be a simulator cycle index")
        return np.asarray(self._active[:, int(cycle)], dtype=bool)

    def advance(self, cycle: int, treatment):
        del cycle, treatment


def conditional_independence_grid(model_parameters: Sequence[SignalModelParameters]):
    """Make a matched d=0 grid, preserving the candidate's q/a/b prior support."""
    models = tuple(model_parameters)
    if not models or any(not isinstance(row, SignalModelParameters) for row in models):
        raise ValueError("model_parameters must be a nonempty sequence of SignalModelParameters")
    if len(set(models)) != len(models):
        raise ValueError("model_parameters must not contain duplicate settings")
    grouped: dict[tuple[float, float, float], set[float]] = {}
    for model in models:
        key = (float(model.q), float(model.a), float(model.b))
        grouped.setdefault(key, set()).add(float(model.d))
    if not any(d != 0.0 for ds in grouped.values() for d in ds):
        raise ValueError("candidate model grid must include at least one nonzero d value")
    counts = {len(ds) for ds in grouped.values()}
    if len(counts) != 1:
        raise ValueError("candidate grid must have the same number of d settings per q/a/b tuple")
    grid = tuple(
        SignalModelParameters(q=q, a=a, b=b, d=0.0)
        for q, a, b in sorted(grouped)
    )
    if len(set(grid)) != len(grid):
        raise ValueError("conditional-independence grid contains duplicate q/a/b settings")
    return grid


def robust_utility(outcome: Mapping, *, tau_days: float, course_penalty: float) -> np.ndarray:
    """Patient utility: worse restricted time minus courses delivered before tau."""
    try:
        baseline = np.asarray(outcome["ttp_baseline"], dtype=float)
        rano = np.asarray(outcome["ttp_rano"], dtype=float)
        courses = np.asarray(outcome["courses_through_cutoff"], dtype=float)
    except (KeyError, TypeError) as exc:
        raise ValueError("outcome must include both endpoints and courses_through_cutoff") from exc
    if baseline.ndim != 1 or rano.shape != baseline.shape or courses.shape != baseline.shape:
        raise ValueError("endpoint times and course counts must be matching one-dimensional arrays")
    if not all(np.isfinite(a).all() for a in (baseline, rano, courses)) or (courses < 0).any():
        raise ValueError("outcome arrays must be finite with non-negative course counts")
    tau = float(tau_days)
    cost = float(course_penalty)
    if not np.isfinite(tau) or tau <= 0 or not np.isfinite(cost) or cost < 0:
        raise ValueError("tau_days must be positive and course_penalty non-negative")
    return np.minimum(np.minimum(baseline, tau), np.minimum(rano, tau)) - cost * courses


def _checked_outcome(outcome: Mapping, n_patients: int, course_budget: int) -> dict:
    needed = ("ttp_baseline", "ttp_rano", "courses", "courses_through_cutoff", "stop_endpoint")
    if not isinstance(outcome, Mapping) or not all(key in outcome for key in needed):
        raise ValueError("simulator result must include both endpoints, courses, and stop_endpoint")
    clean = {}
    for key in ("ttp_baseline", "ttp_rano", "courses", "courses_through_cutoff"):
        value = np.asarray(outcome[key], dtype=float)
        if value.shape != (n_patients,) or not np.isfinite(value).all() or (value < 0).any():
            raise ValueError(f"simulator result {key!r} must be finite, non-negative, and match the cohort")
        clean[key] = value.copy()
    if (clean["courses"] > course_budget).any() or not np.equal(clean["courses"], np.floor(clean["courses"])).all():
        raise ValueError("simulator course counts violate the shared budget")
    if ((clean["courses_through_cutoff"] > clean["courses"]).any()
            or not np.equal(clean["courses_through_cutoff"], np.floor(clean["courses_through_cutoff"])).all()):
        raise ValueError("restricted course counts must be integers no greater than total courses")
    if outcome["stop_endpoint"] != "either":
        raise ValueError("all arms must use stop_endpoint='either'")
    clean["stop_endpoint"] = "either"
    return clean


def select_fixed_schedule(
    selection_scenarios: Mapping,
    *,
    seed,
    settings: PolicyBenchmarkSettings,
    schedules: np.ndarray | None = None,
    schedule_names: Sequence[str] | None = None,
) -> dict:
    """Select one fixed schedule on selection patients under the shared-path utility."""
    n_patients = _validate_scenarios(selection_scenarios, "selection_scenarios")
    seed_key = _validate_seed(seed, "selection seed")
    using_default_catalog = schedules is None and schedule_names is None
    if schedules is None or schedule_names is None:
        catalog_schedules, catalog_names = rg.candidate_schedules()
        schedules = catalog_schedules if schedules is None else schedules
        schedule_names = catalog_names if schedule_names is None else schedule_names
    schedule_array = np.asarray(schedules, dtype=bool)
    names = tuple(schedule_names)
    if schedule_array.ndim != 2 or schedule_array.shape[1] != tc.N_CYCLES or schedule_array.shape[0] == 0:
        raise ValueError("schedules must be a nonempty schedule catalog with one column per cycle")
    if len(names) != schedule_array.shape[0] or any(not isinstance(name, str) or not name for name in names):
        raise ValueError("schedule_names must name every schedule")
    if using_default_catalog:
        eligible = schedule_array.sum(axis=1) <= settings.course_budget
        schedule_array = schedule_array[eligible]
        names = tuple(name for name, keep in zip(names, eligible) if keep)
        if schedule_array.shape[0] == 0:
            raise ValueError("default schedule catalog has no schedule within the configured budget")
    if np.any(schedule_array.sum(axis=1) > settings.course_budget):
        raise ValueError("schedule catalog exceeds the configured course budget")

    flat_scenarios = _repeat_scenarios(selection_scenarios, len(schedule_array), n_patients)
    flat_scenarios["scan_every"] = np.full(len(schedule_array) * n_patients, settings.scan_every)
    schedule_rows = np.repeat(schedule_array, n_patients, axis=0)
    outcome = tc.simulate_two_endpoints(
        flat_scenarios,
        FixedSchedulePolicy(schedule_rows),
        seed=list(seed_key),
        stop_endpoint="either",
        budget=settings.course_budget,
        course_cutoff_days=settings.tau_days,
    )
    checked = _checked_outcome(outcome, len(schedule_rows), settings.course_budget)
    utilities = robust_utility(
        checked, tau_days=settings.tau_days, course_penalty=settings.course_penalty
    ).reshape(len(schedule_array), n_patients)
    course_means = checked["courses_through_cutoff"].reshape(len(schedule_array), n_patients).mean(axis=1)
    total_course_means = checked["courses"].reshape(len(schedule_array), n_patients).mean(axis=1)
    utility_means = utilities.mean(axis=1)
    best_index = max(
        range(len(schedule_array)),
        key=lambda index: (
            float(utility_means[index]), -float(course_means[index]),
            -float(total_course_means[index]), -index,
        ),
    )
    return {
        "schedule_index": int(best_index),
        "schedule_name": names[best_index],
        "schedule": schedule_array[best_index].copy(),
        "selection_mean_robust_utility": float(utility_means[best_index]),
        "selection_mean_courses": float(course_means[best_index]),
        "selection_mean_total_courses": float(total_course_means[best_index]),
        "selection_seed": list(seed_key),
        "n_selection_patients": int(n_patients),
    }


def _evaluate_policy(policy, scenarios: Mapping, seed: tuple[int, ...], settings: PolicyBenchmarkSettings) -> dict:
    n_patients = _validate_scenarios(scenarios, "evaluation_scenarios")
    outcome = tc.simulate_two_endpoints(
        _copy_scenarios(scenarios), policy, seed=list(seed), stop_endpoint="either",
        budget=settings.course_budget, course_cutoff_days=settings.tau_days,
    )
    clean = _checked_outcome(outcome, n_patients, settings.course_budget)
    clean["rmst_baseline"] = np.minimum(clean["ttp_baseline"], settings.tau_days)
    clean["rmst_rano"] = np.minimum(clean["ttp_rano"], settings.tau_days)
    clean["robust_utility"] = robust_utility(
        clean, tau_days=settings.tau_days, course_penalty=settings.course_penalty
    )
    return clean


def run_shared_policy_world(
    *,
    world_id: str,
    selection_scenarios: Mapping,
    evaluation_scenarios: Mapping,
    evaluation_observations: SignalObservations,
    model_parameters: Sequence[SignalModelParameters],
    selection_seed,
    q_training_seed,
    evaluation_seed,
    settings: PolicyBenchmarkSettings = PolicyBenchmarkSettings(),
) -> dict:
    """Run a paired, shared-path development comparison in one virtual world.

    Training and fixed-schedule selection use selection patients only. The LOO signal
    controllers use only unlabeled baseline signals from the evaluation cohort and leave
    each patient's own signal out of that patient's model-weight fit. All arms share the
    same evaluation cohort, evaluation seed, scan interval, endpoint-stop rule, and budget.
    This function returns in-memory synthetic outcomes and writes no files.
    """
    if not isinstance(world_id, str) or not world_id.strip():
        raise ValueError("world_id must be a nonempty string")
    n_selection = _validate_scenarios(selection_scenarios, "selection_scenarios")
    n_evaluation = _validate_scenarios(evaluation_scenarios, "evaluation_scenarios")
    if settings.course_budget < 6:
        raise ValueError("course_budget must be at least 6 to include the standard six-course comparator")
    f0_low, f0_high = settings.signal_f0_bounds
    for name, scenarios in (
        ("selection_scenarios", selection_scenarios),
        ("evaluation_scenarios", evaluation_scenarios),
    ):
        if "f0" not in scenarios:
            raise ValueError(f"{name} must include f0 values covered by signal_f0_bounds")
        f0 = np.asarray(scenarios["f0"], dtype=float)
        if (
            f0.shape != np.asarray(scenarios["r"]).shape
            or not np.isfinite(f0).all()
            or (f0 < f0_low).any()
            or (f0 > f0_high).any()
        ):
            raise ValueError(f"{name} contains f0 values outside the declared signal_f0_bounds")
    selection_seed_key = _validate_seed(selection_seed, "selection_seed")
    q_seed_key = _validate_seed(q_training_seed, "q_training_seed")
    eval_seed_key = _validate_seed(evaluation_seed, "evaluation_seed")
    if len({selection_seed_key, q_seed_key, eval_seed_key}) != 3:
        raise ValueError("selection, Q-training, and evaluation seeds must be distinct")
    observations = evaluation_observations
    if not all(hasattr(observations, key) for key in ("molecular_report", "mri_score")):
        raise TypeError("evaluation_observations must provide molecular_report and mri_score")
    report = np.asarray(observations.molecular_report)
    mri_score = np.asarray(observations.mri_score, dtype=float)
    if report.shape != (n_evaluation,) or mri_score.shape != (n_evaluation,):
        raise ValueError("evaluation signals must have one value per evaluation patient")
    if not np.isin(report, (False, True, 0, 1)).all() or not np.isfinite(mri_score).all():
        raise ValueError("evaluation signals must be binary and finite")
    observations = SignalObservations(report.astype(bool, copy=True), mri_score.copy())
    models = tuple(model_parameters)
    ci_models = conditional_independence_grid(models)

    chosen_fixed = select_fixed_schedule(
        selection_scenarios, seed=selection_seed_key, settings=settings
    )
    q_module = _load_source_module("trackc_qlearning_baseline", "207_trackc_qlearning_baseline.py")
    loo_module = _load_source_module("trackc_loo_model_average", "214_trackc_loo_model_average.py")
    mpc = loo_module._MPC
    q_policy = q_module.train_q_learning_shared_endpoints(
        _copy_scenarios(selection_scenarios),
        episodes=settings.q_episodes,
        course_cap=settings.course_budget,
        course_penalty=settings.course_penalty,
        scan_every=settings.scan_every,
        tau_days=settings.tau_days,
        learning_rate=settings.q_learning_rate,
        discount=1.0,
        epsilon=settings.q_epsilon,
        epsilon_decay=settings.q_epsilon_decay,
        min_epsilon=settings.q_min_epsilon,
        seed=q_seed_key,
    )

    common_mpc = {
        "scan_every": settings.scan_every,
        "budget": settings.course_budget,
        "course_penalty": settings.course_penalty,
        "n_particles": settings.mpc_particles,
        "prediction_particles": settings.prediction_particles,
        "forecast_cycles": settings.forecast_cycles,
        "stop_endpoint": "either",
        "ranges": mpc.INTEGRATED_BELIEF_RANGES,
    }
    full_signal_mpc = loo_module.LOOModelAveragedBeliefMPC(
        observations, model_parameters=models, f0_bounds=settings.signal_f0_bounds,
        n_f0_nodes=settings.n_f0_nodes, **deepcopy(common_mpc)
    )
    ci_signal_mpc = loo_module.LOOModelAveragedBeliefMPC(
        observations, model_parameters=ci_models, f0_bounds=settings.signal_f0_bounds,
        n_f0_nodes=settings.n_f0_nodes, **deepcopy(common_mpc)
    )
    no_signal_mpc = mpc.RobustBeliefMPC(**deepcopy(common_mpc))
    standard_schedule = np.zeros(tc.N_CYCLES, dtype=bool)
    standard_schedule[:6] = True
    fixed_schedule = chosen_fixed["schedule"]
    policies = {
        "loo_dependence_aware_mpc": full_signal_mpc,
        "loo_conditional_independence_mpc": ci_signal_mpc,
        "no_signal_mpc": no_signal_mpc,
        "q_learning": q_policy,
        "selection_best_fixed": FixedSchedulePolicy(fixed_schedule),
        "standard_six_course": FixedSchedulePolicy(standard_schedule),
    }
    arms = {
        name: _evaluate_policy(policy, evaluation_scenarios, eval_seed_key, settings)
        for name, policy in policies.items()
    }
    candidate = arms["loo_dependence_aware_mpc"]
    control = arms["loo_conditional_independence_mpc"]
    primary = {
        "contrast": "loo_dependence_aware_mpc minus loo_conditional_independence_mpc",
        "robust_utility_difference": candidate["robust_utility"] - control["robust_utility"],
        "rmst_baseline_difference": candidate["rmst_baseline"] - control["rmst_baseline"],
        "rmst_rano_difference": candidate["rmst_rano"] - control["rmst_rano"],
    }
    if primary["robust_utility_difference"].shape != (n_evaluation,):
        raise ArithmeticError("primary paired contrast does not match the evaluation cohort")
    return {
        "analysis_type": "synthetic_virtual_world_development",
        "world_id": world_id,
        "settings": settings.__dict__.copy(),
        "seeds": {
            "selection": list(selection_seed_key),
            "q_training": list(q_seed_key),
            "evaluation": list(eval_seed_key),
        },
        "n_selection_patients": n_selection,
        "n_evaluation_patients": n_evaluation,
        "selected_fixed_schedule": chosen_fixed,
        "arm_outcomes": arms,
        "primary_paired_contrast": primary,
        "path_rule": "all arms use one shared treatment path; stop later treatment at the first baseline or RANO-like event; both endpoint times share one course count",
        "signal_rule": "LOO arms use only unlabeled evaluation-cohort baseline signals and leave each patient out of their own peer-weight fit",
        "limits": [
            "Development simulation only; no patient treatment effect is estimated.",
            "Track A and Duke data are not joined to these virtual patients.",
            "The signal-to-resistance and treatment-response relationships are assumed by the simulator.",
            "This function does not establish novelty or global method superiority.",
        ],
    }


__all__ = [
    "FixedSchedulePolicy",
    "PolicyBenchmarkSettings",
    "conditional_independence_grid",
    "robust_utility",
    "run_shared_policy_world",
    "select_fixed_schedule",
]
