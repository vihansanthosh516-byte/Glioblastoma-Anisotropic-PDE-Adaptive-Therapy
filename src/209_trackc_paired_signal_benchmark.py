"""Paired synthetic benchmark for the integrated Track C signal policies.

The primary comparison pairs a dependence-aware joint likelihood with the same
likelihood under conditional independence. All arms receive the same baseline
signals, virtual world, scan schedule, random seed, treatment budget, and utility.
Only the signal channel set or the declared residual-dependence coefficient changes.
This module runs no benchmark at import time and reads no patient data.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass
import importlib.util
from pathlib import Path
import sys
from typing import Mapping

import numpy as np

import trackc_sim as tc
from trackc_integrated_signal_model import SignalModelParameters, SignalObservations


@dataclass(frozen=True)
class SignalArm:
    """One fixed inference likelihood and its visible signal channels."""

    parameters: SignalModelParameters
    channels: tuple[str, ...]


def matched_signal_arms(parameters: SignalModelParameters) -> dict[str, SignalArm]:
    """Build the primary matched pair and channel-removal controls.

    ``joint`` and ``conditional_independence`` share q, a, b, and both observed
    channels. They differ only in d. The other arms remove signal channels while
    keeping the same q, a, and b settings.
    """
    if not isinstance(parameters, SignalModelParameters):
        raise TypeError("parameters must be SignalModelParameters")
    if parameters.d == 0.0:
        raise ValueError("joint primary model requires nonzero residual dependence")
    independent = SignalModelParameters(
        q=parameters.q, a=parameters.a, b=parameters.b, d=0.0
    )
    both = ("molecular", "mri")
    return {
        "joint": SignalArm(parameters, both),
        "conditional_independence": SignalArm(independent, both),
        "molecular_only": SignalArm(parameters, ("molecular",)),
        "mri_only": SignalArm(parameters, ("mri",)),
        "no_signal": SignalArm(parameters, ()),
    }


class ObservationOnlyPolicy:
    """Hide latent simulator fields from a policy's ``reset`` method."""

    def __init__(self, arm_name: str, policy):
        self.arm_name = str(arm_name)
        self.policy = policy
        self.uses_filter = bool(getattr(policy, "uses_filter", False))
        self.reset_context_keys: tuple[str, ...] = ()

    def reset(self, scenarios, rng):
        if "sigma" not in scenarios:
            raise ValueError("policy context must include scan precision sigma")
        sigma = np.asarray(scenarios["sigma"], dtype=float)
        if sigma.ndim != 1 or sigma.size == 0 or not np.isfinite(sigma).all() or (sigma <= 0).any():
            raise ValueError("scan precision sigma must be a finite positive vector")
        self.reset_context_keys = ("sigma",)
        self.policy.reset({"sigma": sigma.copy()}, rng)

    def act(self, cycle, readings):
        return self.policy.act(cycle, readings)

    def advance(self, cycle, treatment):
        return self.policy.advance(cycle, treatment)

    def __getattr__(self, name):
        return getattr(self.policy, name)


def _load_mpc_module():
    module_name = "trackc_robust_belief_mpc"
    if module_name in sys.modules:
        return sys.modules[module_name]
    path = Path(__file__).with_name("208_trackc_robust_belief_mpc.py")
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load Track C MPC module at {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _as_int(value, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer")
    try:
        numeric = float(value)
        result = int(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if not np.isfinite(numeric) or numeric != result or not minimum <= result <= maximum:
        raise ValueError(f"{name} must be an integer in [{minimum}, {maximum}]")
    return result


def _validate_observations(observations, n_patients: int):
    if not all(hasattr(observations, name) for name in ("molecular_report", "mri_score")):
        raise TypeError("observations must provide molecular_report and mri_score")
    report = np.asarray(observations.molecular_report)
    mri_score = np.asarray(observations.mri_score, dtype=float)
    if report.shape != (n_patients,) or mri_score.shape != (n_patients,):
        raise ValueError("each baseline signal must have one value per virtual patient")
    if not np.isin(report, (False, True, 0, 1)).all():
        raise ValueError("molecular reports must be binary")
    if not np.isfinite(mri_score).all():
        raise ValueError("MRI scores must be finite")
    return SignalObservations(report.astype(bool, copy=True), mri_score.copy())


def _validate_arms(arms: Mapping[str, SignalArm]) -> dict[str, SignalArm]:
    if not isinstance(arms, Mapping) or not arms:
        raise ValueError("arms must be a nonempty mapping")
    normalized = dict(arms)
    required = {"joint", "conditional_independence"}
    if not required.issubset(normalized):
        raise ValueError("arms must include joint and conditional_independence")
    for name, arm in normalized.items():
        if not isinstance(name, str) or not name:
            raise ValueError("arm names must be nonempty strings")
        if not isinstance(arm, SignalArm) or not isinstance(arm.parameters, SignalModelParameters):
            raise TypeError("each arm must contain SignalArm and SignalModelParameters")
        if not isinstance(arm.channels, tuple):
            raise TypeError(f"arm {name} channels must be a tuple")
        if len(set(arm.channels)) != len(arm.channels) or set(arm.channels) - {"molecular", "mri"}:
            raise ValueError(f"arm {name} has invalid signal channels")

    joint = normalized["joint"]
    independent = normalized["conditional_independence"]
    if joint.channels != ("molecular", "mri") or independent.channels != joint.channels:
        raise ValueError("primary pair must use the same joint molecular and MRI channels")
    if joint.parameters.d == 0.0 or independent.parameters.d != 0.0:
        raise ValueError("primary pair requires nonzero joint d and conditional-independence d=0")
    if (joint.parameters.q, joint.parameters.a, joint.parameters.b) != (
        independent.parameters.q, independent.parameters.a, independent.parameters.b
    ):
        raise ValueError("q, a, and b must match in the primary pair")
    for arm in normalized.values():
        if (arm.parameters.q, arm.parameters.a, arm.parameters.b) != (
            joint.parameters.q, joint.parameters.a, joint.parameters.b
        ):
            raise ValueError("all arms must use matched q, a, and b")
    return normalized


def _observations_for_arm(observations: SignalObservations, channels: tuple[str, ...]):
    """Give each ablation only its declared channels; remove unused values."""
    if not channels:
        return None
    n = np.asarray(observations.molecular_report).size
    report = (
        np.asarray(observations.molecular_report, dtype=bool).copy()
        if "molecular" in channels
        else np.zeros(n, dtype=bool)
    )
    mri_score = (
        np.asarray(observations.mri_score, dtype=float).copy()
        if "mri" in channels
        else np.zeros(n, dtype=float)
    )
    return SignalObservations(report, mri_score)


def _validate_outcome(outcome, n_patients: int, budget: int):
    required = {
        "ttp_baseline", "ttp_rano", "capped_baseline", "capped_rano",
        "courses", "stop_endpoint",
    }
    if not isinstance(outcome, Mapping) or not required.issubset(outcome):
        raise ValueError("simulator result must include both endpoints, courses, and stop_endpoint")
    ttp_baseline = np.asarray(outcome["ttp_baseline"], dtype=float)
    ttp_rano = np.asarray(outcome["ttp_rano"], dtype=float)
    capped_baseline = np.asarray(outcome["capped_baseline"])
    capped_rano = np.asarray(outcome["capped_rano"])
    courses = np.asarray(outcome["courses"], dtype=float)
    if (
        ttp_baseline.shape != (n_patients,)
        or ttp_rano.shape != (n_patients,)
        or capped_baseline.shape != (n_patients,)
        or capped_rano.shape != (n_patients,)
        or courses.shape != (n_patients,)
    ):
        raise ValueError("simulator results must contain one value per virtual patient")
    if capped_baseline.dtype.kind != "b" or capped_rano.dtype.kind != "b":
        raise ValueError("simulator censoring flags must be boolean")
    if not all(np.isfinite(values).all() for values in (ttp_baseline, ttp_rano, courses)):
        raise ValueError("simulator results must be finite")
    if (
        (ttp_baseline < 0).any()
        or (ttp_rano < 0).any()
        or (courses < 0).any()
        or (courses > budget).any()
    ):
        raise ValueError("simulator results violate nonnegative time or course budget")
    if (ttp_baseline > tc.HORIZON).any() or (ttp_rano > tc.HORIZON).any():
        raise ValueError("simulator result exceeds the configured horizon")
    if not np.equal(courses, np.floor(courses)).all():
        raise ValueError("course counts must be integers")
    stop_endpoint = outcome["stop_endpoint"]
    if stop_endpoint not in ("baseline", "rano", "either"):
        raise ValueError("stop_endpoint must be baseline, rano, or either")
    return {
        "ttp_baseline": ttp_baseline.copy(),
        "ttp_rano": ttp_rano.copy(),
        "capped_baseline": capped_baseline.copy(),
        "capped_rano": capped_rano.copy(),
        "courses": courses.copy(),
        "stop_endpoint": stop_endpoint,
    }


def run_paired_signal_world(
    *,
    world_id: str,
    scenarios: Mapping,
    observations: SignalObservations,
    arms: Mapping[str, SignalArm],
    generating_parameters: SignalModelParameters,
    seed: int,
    scan_every: int,
    budget: int,
    tau: float = 730.0,
    course_penalty: float = 5.0,
    policy_kwargs: Mapping | None = None,
    stop_endpoint: str = "either",
) -> dict:
    """Evaluate fixed signal arms on one synthetic world under both endpoints.

    This is a paired evaluation helper, not a training or model-selection routine.
    It uses the same world, baseline reports, scan interval, seed, and course budget
    for every arm. Each arm runs once and scores both endpoint times on one treatment
    trajectory. By default, treatment stops at the first of baseline or RANO-like
    progression. ``trackc_sim.simulate_two_endpoints`` pre-draws the scan-noise stream
    before policy initialization; repeated arms with the same seed therefore use common
    scan noise. Each policy reset receives only the observable ``sigma`` field.
    """
    if not isinstance(world_id, str) or not world_id.strip():
        raise ValueError("world_id must be a nonempty string")
    if isinstance(seed, (bool, np.bool_)) or int(seed) != seed or int(seed) < 0:
        raise ValueError("seed must be a nonnegative integer")
    seed = int(seed)
    scan_every = _as_int(scan_every, "scan_every", 1, 4)
    budget = _as_int(budget, "budget", 1, tc.BUDGET)
    if stop_endpoint not in ("baseline", "rano", "either"):
        raise ValueError("stop_endpoint must be baseline, rano, or either")
    tau = float(tau)
    course_penalty = float(course_penalty)
    if not np.isfinite(tau) or tau <= 0:
        raise ValueError("tau must be finite and positive")
    if not np.isfinite(course_penalty) or course_penalty < 0:
        raise ValueError("course_penalty must be finite and nonnegative")
    if not isinstance(scenarios, Mapping) or "sigma" not in scenarios or "r" not in scenarios:
        raise ValueError("scenarios must include sigma and r")
    scenario_sigma = np.asarray(scenarios["sigma"], dtype=float)
    n_patients = scenario_sigma.size
    if (
        scenario_sigma.ndim != 1
        or n_patients == 0
        or not np.isfinite(scenario_sigma).all()
        or (scenario_sigma <= 0).any()
    ):
        raise ValueError("scan precision sigma must be a nonempty finite positive vector")
    required_vectors = ("r", "m", "f0", "c", "eps", "B0", "a_sr", "a_rs")
    for key in required_vectors:
        if key not in scenarios:
            raise ValueError(f"scenarios must include {key}")
        value = np.asarray(scenarios[key], dtype=float)
        if value.shape != (n_patients,):
            raise ValueError(f"{key} must have shape ({n_patients},)")
        if not np.isfinite(value).all():
            raise ValueError(f"{key} must be finite")
    if not isinstance(generating_parameters, SignalModelParameters):
        raise TypeError("generating_parameters must be SignalModelParameters")
    if "scan_every" in scenarios:
        existing_scan_every = np.asarray(scenarios["scan_every"])
        if existing_scan_every.ndim != 0 or existing_scan_every.item() != scan_every:
            raise ValueError("scenario scan schedule must match the shared scan_every setting")
    observations = _validate_observations(observations, n_patients)
    arms = _validate_arms(arms)
    policy_kwargs = {} if policy_kwargs is None else dict(policy_kwargs)
    if {"scan_every", "budget", "stop_endpoint"} & set(policy_kwargs):
        raise ValueError("scan_every, budget, and stop_endpoint are shared runner settings")

    scenario_template = deepcopy(dict(scenarios))
    scenario_template["scan_every"] = scan_every
    fixed_policy_kwargs = dict(policy_kwargs)
    fixed_policy_kwargs.update(
        scan_every=scan_every,
        budget=budget,
        course_penalty=course_penalty,
        stop_endpoint=stop_endpoint,
    )
    mpc = _load_mpc_module()
    endpoint_results = {"baseline": {}, "rano": {}}
    common_courses = {}

    for arm_name, arm in arms.items():
        policy = mpc.SignalConditionedBeliefMPC(
            _observations_for_arm(observations, arm.channels),
            arm.parameters,
            channels=arm.channels,
            **fixed_policy_kwargs,
        )
        guarded_policy = ObservationOnlyPolicy(arm_name, policy)
        outcome = tc.simulate_two_endpoints(
            deepcopy(scenario_template),
            guarded_policy,
            seed=seed,
            stop_endpoint=stop_endpoint,
            budget=budget,
        )
        checked = _validate_outcome(outcome, n_patients, budget)
        if checked["stop_endpoint"] != stop_endpoint:
            raise ValueError("simulator used a different treatment stop rule")
        common_courses[arm_name] = checked["courses"]
        endpoint_results["baseline"][arm_name] = {
            "ttp": checked["ttp_baseline"],
            "capped": checked["capped_baseline"],
        }
        endpoint_results["rano"][arm_name] = {
            "ttp": checked["ttp_rano"],
            "capped": checked["capped_rano"],
        }

    robust_utility = {}
    for arm_name in arms:
        baseline = endpoint_results["baseline"][arm_name]
        rano = endpoint_results["rano"][arm_name]
        robust_utility[arm_name] = mpc.robust_endpoint_utility(
            baseline["ttp"],
            rano["ttp"],
            common_courses[arm_name],
            common_courses[arm_name],
            tau=tau,
            course_penalty=course_penalty,
        )

    arm_specs = {
        name: {
            "parameters": asdict(arm.parameters),
            "channels": list(arm.channels),
        }
        for name, arm in arms.items()
    }
    return {
        "analysis_type": "synthetic_virtual_world",
        "world_id": world_id,
        "seed": seed,
        "scan_every": scan_every,
        "budget": budget,
        "tau": tau,
        "course_penalty": course_penalty,
        "stop_endpoint": stop_endpoint,
        "generating_signal_parameters": asdict(generating_parameters),
        "arm_specs": arm_specs,
        "joint_likelihood_matches_generator": (
            arms["joint"].parameters == generating_parameters
        ),
        "policy_settings": {
            **deepcopy(policy_kwargs),
            "scan_every": scan_every,
            "budget": budget,
            "course_penalty": course_penalty,
            "stop_endpoint": stop_endpoint,
        },
        "common_randomness": (
            "Same simulator seed and virtual world per arm; scan-noise matrix is pre-drawn "
            "before policy initialization. One run per arm scores both endpoints."
        ),
        "endpoint_treatment_stop_rule": (
            f"One shared treatment path per arm; treatment stops at {stop_endpoint!r} "
            "progression. Both endpoint times use that path and one shared course count."
        ),
        "endpoints": endpoint_results,
        "courses": common_courses,
        "robust_utility": robust_utility,
    }
