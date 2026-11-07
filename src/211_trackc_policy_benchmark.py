"""Paired development harness for Track C signal-to-treatment policies.

SIMULATION ONLY. The signal-to-resistance and drug-response links are assumptions.
This module scores no patient data and writes no files. Callers must register fresh
development seeds and freeze the run protocol before evaluating policy outcomes.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import hashlib
import importlib.util
import itertools
import json
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
from trackc_equal_course_frontier import randomize_schedule_mix, select_schedule_mix


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


def _validate_patient_ids(values, n_patients: int, name: str) -> tuple[str, ...]:
    """Require ordered, unique synthetic patient IDs for cohort-split auditing."""
    if isinstance(values, (str, bytes)):
        raise ValueError(f"{name} must be a sequence of nonempty strings")
    try:
        ids = tuple(values)
    except TypeError as exc:
        raise ValueError(f"{name} must be a sequence of nonempty strings") from exc
    if len(ids) != n_patients:
        raise ValueError(f"{name} must contain one ID per patient")
    if any(not isinstance(value, str) or not value.strip() for value in ids):
        raise ValueError(f"{name} must contain nonempty strings")
    if any(value != value.strip() for value in ids):
        raise ValueError(f"{name} values must not have leading or trailing whitespace")
    if len(set(ids)) != len(ids):
        raise ValueError(f"{name} must contain unique IDs")
    return ids


def _validate_disjoint_patient_ids(selection_ids, evaluation_ids, n_selection, n_evaluation):
    selection = _validate_patient_ids(selection_ids, n_selection, "selection_patient_ids")
    evaluation = _validate_patient_ids(evaluation_ids, n_evaluation, "evaluation_patient_ids")
    overlap = set(selection).intersection(evaluation)
    if overlap:
        raise ValueError("selection and evaluation patient IDs must be disjoint")
    return selection, evaluation


def _validate_three_disjoint_patient_ids(
    selection_ids, calibration_ids, evaluation_ids,
    n_selection: int, n_calibration: int, n_evaluation: int,
):
    selection = _validate_patient_ids(selection_ids, n_selection, "selection_patient_ids")
    calibration = _validate_patient_ids(calibration_ids, n_calibration, "calibration_patient_ids")
    evaluation = _validate_patient_ids(evaluation_ids, n_evaluation, "evaluation_patient_ids")
    if set(selection).intersection(calibration):
        raise ValueError("selection and calibration patient IDs must be disjoint")
    if set(selection).intersection(evaluation) or set(calibration).intersection(evaluation):
        raise ValueError("selection, calibration, and evaluation patient IDs must be pairwise disjoint")
    return selection, calibration, evaluation


def _copy_signal_observations(observations, n_patients: int, name: str) -> SignalObservations:
    if not all(hasattr(observations, key) for key in ("molecular_report", "mri_score")):
        raise TypeError(f"{name} must provide molecular_report and mri_score")
    report = np.asarray(observations.molecular_report)
    mri_score = np.asarray(observations.mri_score, dtype=float)
    if report.shape != (n_patients,) or mri_score.shape != (n_patients,):
        raise ValueError(f"{name} must contain one molecular and MRI signal per patient")
    if not np.isin(report, (False, True, 0, 1)).all() or not np.isfinite(mri_score).all():
        raise ValueError(f"{name} signals must be binary and finite")
    return SignalObservations(report.astype(bool, copy=True), mri_score.copy())


def _named_assignment_seed(seed: int, label: str) -> int:
    """Derive stable independent schedule-assignment streams from one master seed."""
    payload = f"trackc-schedule-mixture-v1:{seed}:{label}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def _fingerprint_array_mapping(values: Mapping) -> str:
    """Hash named array inputs with stable key order, dtype, shape, and bytes."""
    if not isinstance(values, Mapping):
        raise TypeError("fingerprint input must be a mapping")
    digest = hashlib.sha256()
    for key in sorted(values):
        if not isinstance(key, str):
            raise TypeError("fingerprint mapping keys must be strings")
        array = np.asarray(values[key])
        if array.dtype.hasobject:
            raise TypeError(f"fingerprint value {key!r} has unsupported object dtype")
        if array.dtype.kind in "biufc":
            array = np.asarray(array, dtype=array.dtype.newbyteorder("<"))
            payload = np.ascontiguousarray(array).tobytes(order="C")
        else:
            payload = json.dumps(
                array.tolist(), sort_keys=True, separators=(",", ":"), ensure_ascii=False
            ).encode("utf-8")
        header = json.dumps(
            {"key": key, "dtype": array.dtype.str, "shape": list(array.shape)},
            sort_keys=True, separators=(",", ":"),
        ).encode("utf-8")
        digest.update(len(header).to_bytes(8, "little"))
        digest.update(header)
        digest.update(len(payload).to_bytes(8, "little"))
        digest.update(payload)
    return digest.hexdigest()


def _validate_signal_generation_metadata(
    parameters: SignalModelParameters,
    seed,
    mri_template_sha256: str | None,
) -> tuple[int, str | None]:
    if not isinstance(parameters, SignalModelParameters):
        raise TypeError("signal_generation_parameters must be SignalModelParameters")
    seed_key = _validate_seed(seed, "signal_generation_seed")
    if len(seed_key) != 1:
        raise ValueError("signal_generation_seed must be one integer")
    if mri_template_sha256 is not None:
        if (
            not isinstance(mri_template_sha256, str)
            or len(mri_template_sha256) != 64
            or any(char not in "0123456789abcdefABCDEF" for char in mri_template_sha256)
        ):
            raise ValueError("signal_generation_mri_template_sha256 must be a 64-character hex digest")
        mri_template_sha256 = mri_template_sha256.lower()
    return seed_key[0], mri_template_sha256


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


def _policy_patient_count(context: Mapping, name: str) -> int:
    """Validate the observable reset context and return its patient count."""
    if not isinstance(context, Mapping) or "sigma" not in context:
        raise ValueError(f"{name} must contain observable sigma values")
    sigma = np.asarray(context["sigma"], dtype=float)
    if sigma.ndim != 1 or sigma.size == 0:
        raise ValueError(f"{name} sigma must be a nonempty one-dimensional array")
    if not np.isfinite(sigma).all() or (sigma <= 0).any():
        raise ValueError(f"{name} sigma values must be finite and positive")
    return int(sigma.size)


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
        n_patients = _policy_patient_count(scen, "policy reset context")
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


def _fixed_schedule_control_catalog(settings: PolicyBenchmarkSettings) -> dict:
    """Build explicit fixed-course references allowed by the configured global budget.

    These are simulator schedule controls, not clinical regimens. The simulator's
    treatment cycle is 28 days; the every-other-cycle arm is a 56-day pacing control,
    not a dose-dense clinical schedule.
    """
    controls = {}

    six_course = np.zeros(tc.N_CYCLES, dtype=bool)
    six_course[:6] = True
    controls["standard_six_course"] = {
        "schedule": six_course,
        "planned_courses": 6,
        "cycle_indices": list(range(6)),
        "nominal_interval_days": tc.CYCLE_DAYS,
        "description": "six courses at the first six 28-day cycle starts",
    }

    if settings.course_budget >= 12:
        front_loaded = np.zeros(tc.N_CYCLES, dtype=bool)
        front_loaded[:12] = True
        controls["standard_twelve_course_front_loaded"] = {
            "schedule": front_loaded,
            "planned_courses": 12,
            "cycle_indices": list(range(12)),
            "nominal_interval_days": tc.CYCLE_DAYS,
            "description": "twelve courses at the first twelve 28-day cycle starts",
        }

        paced = np.zeros(tc.N_CYCLES, dtype=bool)
        paced[0:24:2] = True
        controls["paced_twelve_course_every_other_cycle"] = {
            "schedule": paced,
            "planned_courses": 12,
            "cycle_indices": list(range(0, 24, 2)),
            "nominal_interval_days": 2 * tc.CYCLE_DAYS,
            "description": "twelve courses at every other 28-day cycle start (56-day spacing)",
        }
    return controls


class CourseCappedPolicy:
    """Apply an arm-specific course cap while preserving its policy interface."""

    def __init__(self, policy, max_courses: int):
        self.policy = policy
        self.max_courses = _validate_int(max_courses, "max_courses", 0, tc.BUDGET)
        self.uses_filter = bool(getattr(policy, "uses_filter", False))
        self._used = None

    def reset(self, scenarios: Mapping, rng: np.random.Generator):
        n_patients = _policy_patient_count(scenarios, "policy reset context")
        self._used = np.zeros(n_patients, dtype=np.int16)
        self.policy.reset(scenarios, rng)

    def act(self, cycle: int, observations: np.ndarray) -> np.ndarray:
        if self._used is None:
            raise RuntimeError("policy must be reset before act")
        proposed = np.asarray(self.policy.act(cycle, observations), dtype=bool)
        if proposed.shape != self._used.shape:
            raise ValueError("wrapped policy must return one treatment flag per patient")
        return proposed & (self._used < self.max_courses)

    def advance(self, cycle: int, treatment: np.ndarray):
        if self._used is None:
            raise RuntimeError("policy must be reset before advance")
        actual = np.asarray(treatment, dtype=bool)
        if actual.shape != self._used.shape:
            raise ValueError("advance must receive one delivered-treatment flag per patient")
        self.policy.advance(cycle, actual)
        self._used += actual.astype(np.int16)


def baseline_policy_catalog(settings: PolicyBenchmarkSettings) -> dict:
    """Return fixed scan-rule and state-estimation comparators with their provenance."""
    scan_configurations = {
        "raw_zhang_1.0_0.5": "RawThreshold(a_on=1.0, a_off=0.5)",
        "raw_0.6_0.6": "RawThreshold(a_on=0.6, a_off=0.6)",
        "nadir_k3_on1.25": "NadirContainment(induction=3, a_on=1.25, a_off=1.0)",
        "nadir_k3_on1.15": "NadirContainment(induction=3, a_on=1.15, a_off=1.0)",
        "nadir_k6_on1.25": "NadirContainment(induction=6, a_on=1.25, a_off=1.0)",
        "resp_k3_th0.5": "ResponseGuided(induction=3, theta=0.5)",
        "resp_k3_th0.75": "ResponseGuided(induction=3, theta=0.75)",
        "resp_k3_th1.0": "ResponseGuided(induction=3, theta=1.0)",
        "resp_k6_th0.75": "ResponseGuided(induction=6, theta=0.75)",
    }
    established_scan_arms = rg.policy_arms()
    missing = set(scan_configurations) - set(established_scan_arms)
    if missing:
        raise RuntimeError(f"trackc_regime is missing baseline arms: {sorted(missing)}")
    catalog = {}
    for name, configuration in scan_configurations.items():
        factory, max_courses = established_scan_arms[name]
        catalog[name] = {
            "factory": factory,
            "max_courses": max_courses,
            "family": "scan_rule",
            "configuration": configuration,
            "provenance": "src/trackc_regime.py:policy_arms",
            "evaluate_all": True,
        }

    catalog.update({
        "filtered_threshold_pilot2": {
            "factory": lambda: tc.FilteredThreshold(
                0.6, 0.3, 0.6, 0.05, n_particles=settings.mpc_particles
            ),
            "max_courses": tc.BUDGET,
            "family": "state_estimation",
            "configuration": (
                f"FilteredThreshold(a_on=0.6, a_off=0.3, a_R=0.6, b=0.05, "
                f"n_particles={settings.mpc_particles})"
            ),
            "provenance": (
                "output/trackc_pilot2/pilot_results.json; fixed historical configuration "
                "selected on its training draw"
            ),
            "evaluate_all": True,
        },
        "pf_wide_containment": {
            "factory": lambda: tc.RiskContainment(
                0.10, 2, 0.0, True, n_particles=settings.mpc_particles
            ),
            "max_courses": tc.BUDGET,
            "family": "state_estimation",
            "configuration": (
                f"RiskContainment(p_star=0.10, horizon=2, wide=True, "
                f"n_particles={settings.mpc_particles})"
            ),
            "provenance": "src/206_trackc_marker_benchmark.py:comparison_policy_specs",
            "evaluate_all": True,
        },
        "pf_tail_trigger_candidate": {
            "factory": lambda: tc.RiskContainment(
                0.05, 4, 0.0, True, n_particles=settings.mpc_particles, front=4
            ),
            "max_courses": tc.BUDGET,
            "family": "state_estimation",
            "configuration": (
                f"RiskContainment(p_star=0.05, horizon=4, wide=True, front=4, "
                f"n_particles={settings.mpc_particles})"
            ),
            "provenance": "src/206_trackc_marker_benchmark.py:comparison_policy_specs",
            "evaluate_all": True,
        },
        "pf_holiday_containment": {
            "factory": lambda: tc.HolidayContainment(
                0.10, 2, front=0, max_run=2, n_particles=settings.mpc_particles
            ),
            "max_courses": tc.BUDGET,
            "family": "state_estimation",
            "configuration": (
                f"HolidayContainment(p_star=0.10, horizon=2, max_run=2, "
                f"n_particles={settings.mpc_particles})"
            ),
            "provenance": "src/206_trackc_marker_benchmark.py:comparison_policy_specs",
            "evaluate_all": True,
        },
    })
    on_off_grid = (
        (0.6, 0.3), (0.6, 0.5), (0.8, 0.3), (0.8, 0.5),
        (0.8, 0.7), (1.0, 0.3), (1.0, 0.5), (1.0, 0.7),
    )
    resistance_thresholds = (0.6, 0.8, 1.0, 1.1)
    trigger_probabilities = (0.05, 0.10, 0.20, 0.35)
    for (a_on, a_off), a_resistance, p_trigger in itertools.product(
        on_off_grid, resistance_thresholds, trigger_probabilities
    ):
        name = (
            f"filtered_threshold_grid_on{a_on:g}_off{a_off:g}_"
            f"ar{a_resistance:g}_ptrig{p_trigger:g}"
        )
        catalog[name] = {
            "factory": lambda a_on=a_on, a_off=a_off, a_resistance=a_resistance, p_trigger=p_trigger: (
                tc.FilteredThreshold(
                    a_on, a_off, a_resistance, p_trigger,
                    n_particles=settings.mpc_particles,
                )
            ),
            "max_courses": tc.BUDGET,
            "family": "filtered_threshold_grid",
            "configuration": (
                f"FilteredThreshold(a_on={a_on:g}, a_off={a_off:g}, "
                f"a_R={a_resistance:g}, b={p_trigger:g}, "
                f"n_particles={settings.mpc_particles})"
            ),
            "provenance": "src/150_trackc_pilot.py:hyb_specs; selection rerun on current selection patients",
            "evaluate_all": False,
        }
    raw_on_off_grid = (
        (0.6, 0.3), (0.6, 0.4), (0.6, 0.5), (0.8, 0.3),
        (0.8, 0.5), (0.8, 0.7), (0.9, 0.5), (1.0, 0.3),
        (1.0, 0.5), (1.0, 0.7), (1.1, 0.7), (1.2, 0.8),
    )
    raw_single_grid = (0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.18)
    raw_configurations = [
        (f"raw_threshold_grid_on{a_on:g}_off{a_off:g}", a_on, a_off)
        for a_on, a_off in raw_on_off_grid
    ] + [
        (f"raw_threshold_grid_single{threshold:g}", threshold, threshold)
        for threshold in raw_single_grid
    ]
    for name, a_on, a_off in raw_configurations:
        catalog[name] = {
            "factory": lambda a_on=a_on, a_off=a_off: tc.RawThreshold(a_on, a_off),
            "max_courses": tc.BUDGET,
            "family": "raw_threshold_grid",
            "configuration": f"RawThreshold(a_on={a_on:g}, a_off={a_off:g})",
            "provenance": "src/167_trackc_sealed_containment.py:raw_onoff and raw_single families",
            "evaluate_all": False,
        }
    return catalog


def select_baseline_policy_families(
    selection_scenarios: Mapping,
    *,
    seed,
    settings: PolicyBenchmarkSettings,
    catalog: Mapping | None = None,
) -> dict:
    """Select one comparator per predeclared family using selection patients only."""
    n_patients = _validate_scenarios(selection_scenarios, "selection_scenarios")
    seed_key = _validate_seed(seed, "baseline selection seed")
    candidates = baseline_policy_catalog(settings) if catalog is None else dict(catalog)
    if not candidates:
        raise ValueError("baseline policy catalog must not be empty")
    scores = {}
    for name, spec in candidates.items():
        family = spec.get("family")
        if not isinstance(family, str) or not family:
            raise ValueError(f"baseline {name!r} must have a named family")
        cap = min(
            _validate_int(spec["max_courses"], f"{name}.max_courses", 0, tc.BUDGET),
            settings.course_budget,
        )
        policy = CourseCappedPolicy(spec["factory"](), cap)
        outcome = _evaluate_policy(policy, selection_scenarios, seed_key, settings)
        scores[name] = float(np.mean(outcome["robust_utility"]))

    winners = {}
    for name, spec in candidates.items():
        family = spec["family"]
        if family not in winners or scores[name] > scores[winners[family]["name"]]:
            winners[family] = {"name": name, "selection_mean_robust_utility": scores[name]}
    return {
        "objective": "mean robust utility on selection patients",
        "tie_rule": "first catalog entry wins exact selection-score ties",
        "selection_seed": list(seed_key),
        "n_selection_patients": n_patients,
        "scores": scores,
        "winners": winners,
    }


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
    needed = (
        "ttp_baseline", "capped_baseline", "ttp_rano", "capped_rano",
        "courses", "courses_through_cutoff", "stop_endpoint",
    )
    if not isinstance(outcome, Mapping) or not all(key in outcome for key in needed):
        raise ValueError("simulator result must include both endpoints, censor flags, courses, and stop_endpoint")
    clean = {}
    for key in ("ttp_baseline", "ttp_rano", "courses", "courses_through_cutoff"):
        value = np.asarray(outcome[key], dtype=float)
        if value.shape != (n_patients,) or not np.isfinite(value).all() or (value < 0).any():
            raise ValueError(f"simulator result {key!r} must be finite, non-negative, and match the cohort")
        clean[key] = value.copy()
    for key in ("capped_baseline", "capped_rano"):
        value = np.asarray(outcome[key])
        if value.shape != (n_patients,) or value.dtype.kind != "b":
            raise ValueError(f"simulator result {key!r} must be a boolean vector matching the cohort")
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
    raw_schedules = np.asarray(schedules)
    if not np.isin(raw_schedules, (False, True, 0, 1)).all():
        raise ValueError("schedules must contain only binary actions")
    schedule_array = raw_schedules.astype(bool, copy=True)
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
        policy_context="observable",
    )
    checked = _checked_outcome(outcome, len(schedule_rows), settings.course_budget)
    utilities = robust_utility(
        checked, tau_days=settings.tau_days, course_penalty=settings.course_penalty
    ).reshape(len(schedule_array), n_patients)
    course_means = checked["courses_through_cutoff"].reshape(len(schedule_array), n_patients).mean(axis=1)
    total_course_means = checked["courses"].reshape(len(schedule_array), n_patients).mean(axis=1)
    utility_means = utilities.mean(axis=1)
    worst_endpoint_days = np.minimum(
        np.minimum(checked["ttp_baseline"], settings.tau_days),
        np.minimum(checked["ttp_rano"], settings.tau_days),
    ).reshape(len(schedule_array), n_patients)
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
        "schedule_metrics": [
            {
                "schedule_index": int(index),
                "schedule_name": names[index],
                "mean_courses": float(course_means[index]),
                "mean_total_courses": float(total_course_means[index]),
                "mean_worst_endpoint_days": float(worst_endpoint_days[index].mean()),
                "mean_robust_utility": float(utility_means[index]),
                "schedule": schedule_array[index].copy(),
            }
            for index in range(len(schedule_array))
        ],
        "selection_seed": list(seed_key),
        "n_selection_patients": int(n_patients),
    }


def _evaluate_policy(policy, scenarios: Mapping, seed: tuple[int, ...], settings: PolicyBenchmarkSettings) -> dict:
    n_patients = _validate_scenarios(scenarios, "evaluation_scenarios")
    evaluation_scenarios = _copy_scenarios(scenarios)
    # The simulator uses all cycles for progression scans when this field is absent.
    # Set it explicitly so each policy and endpoint uses the configured observation schedule.
    evaluation_scenarios["scan_every"] = np.full(n_patients, settings.scan_every, dtype=int)
    outcome = tc.simulate_two_endpoints(
        evaluation_scenarios, policy, seed=list(seed), stop_endpoint="either",
        budget=settings.course_budget, course_cutoff_days=settings.tau_days,
        policy_context="observable",
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
    calibration_scenarios: Mapping,
    evaluation_scenarios: Mapping,
    selection_patient_ids: Sequence[str],
    calibration_patient_ids: Sequence[str],
    evaluation_patient_ids: Sequence[str],
    calibration_observations: SignalObservations,
    evaluation_observations: SignalObservations,
    model_parameters: Sequence[SignalModelParameters],
    signal_generation_parameters: SignalModelParameters,
    signal_generation_seed,
    calibration_signal_generation_seed,
    signal_generation_mri_template_sha256: str | None,
    selection_scenario_seed,
    calibration_scenario_seed,
    evaluation_scenario_seed,
    selection_seed,
    q_training_seed,
    calibration_seed,
    evaluation_seed,
    schedule_mixture_seed,
    settings: PolicyBenchmarkSettings = PolicyBenchmarkSettings(),
    include_ambiguity_controller: bool = False,
) -> dict:
    """Run a three-cohort development comparison in one virtual world.

    Selection patients choose schedules and baseline families. Calibration patients
    estimate feedback-policy course use. Their target doses choose randomized mixtures
    on the selection-only fixed-schedule frontier. Evaluation patients score all arms.
    Cohort IDs and scenario, signal, policy, and assignment seeds must be disjoint.
    Evaluation outcomes never choose a schedule or mixture. This returns synthetic
    outcomes in memory and writes no files.
    """
    if not isinstance(world_id, str) or not world_id.strip():
        raise ValueError("world_id must be a nonempty string")
    n_selection = _validate_scenarios(selection_scenarios, "selection_scenarios")
    n_calibration = _validate_scenarios(calibration_scenarios, "calibration_scenarios")
    n_evaluation = _validate_scenarios(evaluation_scenarios, "evaluation_scenarios")
    selection_ids, calibration_ids, evaluation_ids = _validate_three_disjoint_patient_ids(
        selection_patient_ids, calibration_patient_ids, evaluation_patient_ids,
        n_selection, n_calibration, n_evaluation,
    )
    generation_seed, mri_template_sha256 = _validate_signal_generation_metadata(
        signal_generation_parameters,
        signal_generation_seed,
        signal_generation_mri_template_sha256,
    )
    if settings.course_budget < 6:
        raise ValueError("course_budget must be at least 6 to include the standard six-course comparator")
    f0_low, f0_high = settings.signal_f0_bounds
    for name, scenarios in (
        ("selection_scenarios", selection_scenarios),
        ("calibration_scenarios", calibration_scenarios),
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
    selection_scenario_seed_key = _validate_seed(
        selection_scenario_seed, "selection_scenario_seed"
    )
    calibration_scenario_seed_key = _validate_seed(
        calibration_scenario_seed, "calibration_scenario_seed"
    )
    evaluation_scenario_seed_key = _validate_seed(
        evaluation_scenario_seed, "evaluation_scenario_seed"
    )
    calibration_signal_seed_key = _validate_seed(
        calibration_signal_generation_seed, "calibration_signal_generation_seed"
    )
    if len(calibration_signal_seed_key) != 1:
        raise ValueError("calibration_signal_generation_seed must be one integer")
    selection_seed_key = _validate_seed(selection_seed, "selection_seed")
    q_seed_key = _validate_seed(q_training_seed, "q_training_seed")
    calibration_seed_key = _validate_seed(calibration_seed, "calibration_seed")
    eval_seed_key = _validate_seed(evaluation_seed, "evaluation_seed")
    mixture_seed_key = _validate_seed(schedule_mixture_seed, "schedule_mixture_seed")
    if len(mixture_seed_key) != 1:
        raise ValueError("schedule_mixture_seed must be one integer")
    all_seeds = {
        selection_scenario_seed_key,
        calibration_scenario_seed_key,
        evaluation_scenario_seed_key,
        selection_seed_key,
        q_seed_key,
        calibration_seed_key,
        eval_seed_key,
        mixture_seed_key,
        (generation_seed,),
        calibration_signal_seed_key,
    }
    if len(all_seeds) != 10:
        raise ValueError(
            "selection, calibration, and evaluation scenario, signal, policy, and mixture seeds must be distinct"
        )
    calibration_observations = _copy_signal_observations(
        calibration_observations, n_calibration, "calibration_observations"
    )
    observations = _copy_signal_observations(
        evaluation_observations, n_evaluation, "evaluation_observations"
    )
    models = tuple(model_parameters)
    input_fingerprints = {
        "selection_scenarios_sha256": _fingerprint_array_mapping(selection_scenarios),
        "calibration_scenarios_sha256": _fingerprint_array_mapping(calibration_scenarios),
        "evaluation_scenarios_sha256": _fingerprint_array_mapping(evaluation_scenarios),
        "selection_patient_ids_sha256": _fingerprint_array_mapping({
            "patient_ids": np.asarray(selection_ids, dtype=str),
        }),
        "evaluation_patient_ids_sha256": _fingerprint_array_mapping({
            "patient_ids": np.asarray(evaluation_ids, dtype=str),
        }),
        "calibration_patient_ids_sha256": _fingerprint_array_mapping({
            "patient_ids": np.asarray(calibration_ids, dtype=str),
        }),
        "calibration_observations_sha256": _fingerprint_array_mapping({
            "molecular_report": calibration_observations.molecular_report,
            "mri_score": calibration_observations.mri_score,
        }),
        "evaluation_observations_sha256": _fingerprint_array_mapping({
            "molecular_report": observations.molecular_report,
            "mri_score": observations.mri_score,
        }),
    }
    ci_models = conditional_independence_grid(models)

    chosen_fixed = select_fixed_schedule(
        selection_scenarios, seed=selection_seed_key, settings=settings
    )
    baseline_catalog = baseline_policy_catalog(settings)
    baseline_selection = select_baseline_policy_families(
        selection_scenarios, seed=selection_seed_key, settings=settings,
        catalog=baseline_catalog,
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
    def make_signal_policies(signal_observations):
        policies_for_signals = {
            "loo_dependence_aware_mpc": loo_module.LOOModelAveragedBeliefMPC(
                signal_observations, model_parameters=models,
                f0_bounds=settings.signal_f0_bounds, n_f0_nodes=settings.n_f0_nodes,
                **deepcopy(common_mpc),
            ),
            "loo_conditional_independence_mpc": loo_module.LOOModelAveragedBeliefMPC(
                signal_observations, model_parameters=ci_models,
                f0_bounds=settings.signal_f0_bounds, n_f0_nodes=settings.n_f0_nodes,
                **deepcopy(common_mpc),
            ),
            "no_signal_mpc": mpc.RobustBeliefMPC(**deepcopy(common_mpc)),
        }
        if include_ambiguity_controller:
            ambiguity_module = _load_source_module(
                "trackc_robust_signal_ambiguity",
                "216_trackc_robust_signal_ambiguity.py",
            )
            policies_for_signals["finite_set_minimax_mpc"] = (
                ambiguity_module.SignalAmbiguityBeliefMPC(
                    signal_observations, model_parameters=models,
                    **deepcopy(common_mpc),
                )
            )
            policies_for_signals["finite_set_d0_minimax_mpc"] = (
                ambiguity_module.SignalAmbiguityBeliefMPC(
                    signal_observations, model_parameters=ci_models,
                    **deepcopy(common_mpc),
                )
            )
        return policies_for_signals

    fixed_control_catalog = _fixed_schedule_control_catalog(settings)
    fixed_schedule = chosen_fixed["schedule"]
    policies = make_signal_policies(observations)
    policies["q_learning"] = q_policy
    policies["selection_best_fixed"] = FixedSchedulePolicy(fixed_schedule)
    policies.update({
        name: FixedSchedulePolicy(spec["schedule"])
        for name, spec in fixed_control_catalog.items()
    })
    selected_baseline_names = {
        row["name"] for row in baseline_selection["winners"].values()
    }
    for name, spec in baseline_catalog.items():
        if not spec["evaluate_all"] and name not in selected_baseline_names:
            continue
        cap = min(int(spec["max_courses"]), settings.course_budget)
        policies[name] = CourseCappedPolicy(spec["factory"](), cap)

    calibration_policies = make_signal_policies(calibration_observations)
    calibration_policies["q_learning"] = deepcopy(q_policy)
    for name, spec in baseline_catalog.items():
        if not spec["evaluate_all"] and name not in selected_baseline_names:
            continue
        calibration_policies[name] = deepcopy(policies[name])
    calibration_arms = {
        name: _evaluate_policy(policy, calibration_scenarios, calibration_seed_key, settings)
        for name, policy in calibration_policies.items()
    }
    arms = {
        name: _evaluate_policy(policy, evaluation_scenarios, eval_seed_key, settings)
        for name, policy in policies.items()
    }

    schedule_by_index = {
        int(row["schedule_index"]): row for row in chosen_fixed["schedule_metrics"]
    }
    matched_schedule_controls = {}
    for policy_name, calibration_outcome in calibration_arms.items():
        target_courses = float(np.mean(calibration_outcome["courses_through_cutoff"]))
        match = select_schedule_mix(chosen_fixed["schedule_metrics"], target_courses)
        match_record = {key: value for key, value in match.items() if key != "frontier_points"}
        match_record["frontier_points"] = [
            {key: value for key, value in point.items() if key != "schedule"}
            for point in match["frontier_points"]
        ]
        if match["status"] != "matched":
            matched_schedule_controls[policy_name] = {
                "status": match["status"],
                "calibration_mean_courses": target_courses,
                "selection_mix": match_record,
                "evaluation_arm": None,
            }
            continue
        lower_schedule = schedule_by_index[match["lower_schedule_index"]]["schedule"]
        upper_schedule = schedule_by_index[match["upper_schedule_index"]]["schedule"]
        assignment_seed = _named_assignment_seed(mixture_seed_key[0], policy_name)
        mixed_schedules, assignment = randomize_schedule_mix(
            lower_schedule, upper_schedule,
            match["upper_mix_probability"], n_evaluation, assignment_seed,
        )
        arm_name = f"matched_fixed_mix_for_{policy_name}"
        mixture_outcome = _evaluate_policy(
            FixedSchedulePolicy(mixed_schedules), evaluation_scenarios,
            eval_seed_key, settings,
        )
        arms[arm_name] = mixture_outcome
        adaptive_worst_days = np.minimum(arms[policy_name]["rmst_baseline"], arms[policy_name]["rmst_rano"])
        mixture_worst_days = np.minimum(
            mixture_outcome["rmst_baseline"], mixture_outcome["rmst_rano"]
        )
        matched_schedule_controls[policy_name] = {
            "status": "evaluated",
            "calibration_mean_courses": target_courses,
            "selection_mix": match_record,
            "assignment": assignment,
            "evaluation_arm": arm_name,
            "evaluation_mean_courses": float(np.mean(mixture_outcome["courses_through_cutoff"])),
            "evaluation_course_use_difference": float(
                np.mean(mixture_outcome["courses_through_cutoff"]) - target_courses
            ),
            "paired_worst_endpoint_days_difference": adaptive_worst_days - mixture_worst_days,
            "paired_robust_utility_difference": (
                arms[policy_name]["robust_utility"] - mixture_outcome["robust_utility"]
            ),
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
        "signal_generation": {
            "generator": "trackc_integrated_signal_model.generate_integrated_signals",
            "parameters": {
                "q": float(signal_generation_parameters.q),
                "a": float(signal_generation_parameters.a),
                "b": float(signal_generation_parameters.b),
                "d": float(signal_generation_parameters.d),
            },
            "seed": generation_seed,
            "f0_source": "evaluation_scenarios['f0']",
            "mri_template_sha256": mri_template_sha256,
        },
        "calibration_signal_generation": {
            "generator": "trackc_integrated_signal_model.generate_integrated_signals",
            "parameters": {
                "q": float(signal_generation_parameters.q),
                "a": float(signal_generation_parameters.a),
                "b": float(signal_generation_parameters.b),
                "d": float(signal_generation_parameters.d),
            },
            "seed": calibration_signal_seed_key[0],
            "f0_source": "calibration_scenarios['f0']",
        },
        "inference_model_parameters": [
            {"q": float(model.q), "a": float(model.a), "b": float(model.b), "d": float(model.d)}
            for model in models
        ],
        "signal_generation_in_inference_grid": signal_generation_parameters in models,
        "input_fingerprints": input_fingerprints,
        "seeds": {
            "selection_scenario_generation": list(selection_scenario_seed_key),
            "calibration_scenario_generation": list(calibration_scenario_seed_key),
            "evaluation_scenario_generation": list(evaluation_scenario_seed_key),
            "calibration_signal_generation": list(calibration_signal_seed_key),
            "selection": list(selection_seed_key),
            "q_training": list(q_seed_key),
            "calibration_policy_simulation": list(calibration_seed_key),
            "evaluation": list(eval_seed_key),
            "schedule_mixture_master": list(mixture_seed_key),
        },
        "scenario_seed_provenance": (
            "caller-declared; realized scenario arrays are independently fingerprinted"
        ),
        "cohorts": {
            "selection": {"n_patients": n_selection, "patient_ids_sha256": input_fingerprints["selection_patient_ids_sha256"]},
            "calibration": {"n_patients": n_calibration, "patient_ids_sha256": input_fingerprints["calibration_patient_ids_sha256"]},
            "evaluation": {"n_patients": n_evaluation, "patient_ids_sha256": input_fingerprints["evaluation_patient_ids_sha256"]},
        },
        "n_selection_patients": n_selection,
        "n_calibration_patients": n_calibration,
        "n_evaluation_patients": n_evaluation,
        "selected_fixed_schedule": chosen_fixed,
        "fixed_schedule_control_catalog": {
            name: {
                "planned_courses": int(spec["planned_courses"]),
                "cycle_indices": list(spec["cycle_indices"]),
                "nominal_interval_days": int(spec["nominal_interval_days"]),
                "description": spec["description"],
                "role": "explicit synthetic schedule reference; not a clinical regimen",
            }
            for name, spec in fixed_control_catalog.items()
        },
        "baseline_selection": baseline_selection,
        "baseline_arm_metadata": {
            name: {
                "family": spec["family"],
                "max_courses": min(int(spec["max_courses"]), settings.course_budget),
                "configuration": spec["configuration"],
                "provenance": spec["provenance"],
                "evaluated_on_holdout": bool(
                    spec["evaluate_all"] or name in selected_baseline_names
                ),
            }
            for name, spec in baseline_catalog.items()
        },
        "selected_baseline_by_family": baseline_selection["winners"],
        "arm_outcomes": arms,
        "calibration_course_use": {
            name: float(np.mean(outcome["courses_through_cutoff"]))
            for name, outcome in calibration_arms.items()
        },
        "matched_fixed_schedule_controls": matched_schedule_controls,
        "primary_paired_contrast": primary,
        "exploratory_controller": {
            "name": "finite_set_minimax_mpc",
            "included": bool(include_ambiguity_controller),
            "model_set_size": len(models),
            "role": "additional development arm; excluded from frozen primary contrast",
            "interpretation": (
                "not evaluated" if not include_ambiguity_controller else
                "synthetic development outcome only; no clinical or superiority claim"
            ),
        },
        "exploratory_d0_controller": {
            "name": "finite_set_d0_minimax_mpc",
            "included": bool(include_ambiguity_controller),
            "model_set_size": len(ci_models),
            "dependence_rule": "same q/a/b support as the full set, with d fixed to 0",
            "role": "matched d=0 minimax control; excluded from frozen primary contrast",
            "interpretation": (
                "not evaluated" if not include_ambiguity_controller else
                "synthetic development outcome only; no clinical or superiority claim"
            ),
        },
        "path_rule": "all arms share initial virtual patients, precomputed observation-noise stream, scan schedule, and budget; each arm has its own action-dependent path; both endpoint times share that arm's path and course count; treatment stops at the first baseline or RANO-like event",
        "policy_reset_context": {
            "mode": "observable",
            "fields_exposed": ["sigma"],
            "latent_scenario_fields_exposed": False,
        },
        "signal_rule": "LOO controllers use only unlabeled baseline signals from the cohort being simulated and leave each patient out of their own peer-weight fit; calibration signals set course-use targets, evaluation signals do not set mixtures",
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
    "_fixed_schedule_control_catalog",
    "conditional_independence_grid",
    "robust_utility",
    "run_shared_policy_world",
    "select_fixed_schedule",
]
