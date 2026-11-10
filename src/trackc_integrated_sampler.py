"""Independent-world sampler for the integrated Track ABC development design.

This module is deliberately separate from ``trackc_calibration`` and
``trackc_sim.draw_scenarios`` so historical Stage 2 generators keep their
original outputs. It generates synthetic parameters only; it does not run a
controller or read patient data.
"""
from __future__ import annotations

from numbers import Integral
from typing import Mapping

import numpy as np
from scipy.stats import qmc

import trackc_sim as tc


WORLD_MU_NU_LEVELS = (0.0, 0.01, 0.04)
WORLD_Q_SIGNAL_LEVELS = (0.6, 0.7, 0.8, 0.9)
SIGMA_LEVELS = (0.05, 0.10, 0.20)
COHORT_ROLES = ("selection", "calibration", "evaluation")
PATIENT_DIMENSION_ORDER = ("r", "m", "f0", "B0", "c", "epsilon")
SIGMA_ASSIGNMENT_CHILD = 0x5349474D  # ASCII "SIGM", fixed substream of the world seed.
B0_CENTER_HALF_WIDTH = 0.10
WORLD_AXIS_BOUNDS = {
    "single": (
        ("log_r_shift", -0.7, 0.9),
        ("log_m_scale", float(np.log(0.05)), 0.0),
        ("log_f0_shift", -0.7, 2.3),
        ("B0_center", 0.12, 0.55),
    ),
    "mixture": (
        ("log_r_shift", -0.7, 0.9),
        ("log_m_scale_non", float(np.log(0.02)), float(np.log(0.30))),
        ("log_f0_shift", -0.7, 2.3),
        ("B0_center", 0.12, 0.55),
        ("pi_response", 0.1, 0.6),
    ),
}


def _seed(value: int, name: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def make_generator(seed: int) -> np.random.Generator:
    """Return the explicit bit generator used by the integrated sampler."""
    return np.random.Generator(np.random.PCG64DXSM(_seed(seed, "seed")))


def draw_world_parameters(family: str, seed: int) -> dict:
    """Draw one independent world from the proposed integrated parameter law."""
    if family not in ("single", "mixture"):
        raise ValueError("family must be 'single' or 'mixture'")
    rng = make_generator(seed)
    axes = WORLD_AXIS_BOUNDS[family]
    point = qmc.LatinHypercube(d=len(axes), scramble=True, seed=rng).random(1)[0]
    world = {"family": family}
    for coordinate, (name, low, high) in zip(point, axes):
        world[name] = float(low + coordinate * (high - low))
    world["mu"] = WORLD_MU_NU_LEVELS[int(rng.integers(0, len(WORLD_MU_NU_LEVELS)))]
    world["nu"] = WORLD_MU_NU_LEVELS[int(rng.integers(0, len(WORLD_MU_NU_LEVELS)))]
    world["q_signal"] = WORLD_Q_SIGNAL_LEVELS[int(rng.integers(0, len(WORLD_Q_SIGNAL_LEVELS)))]
    return world


def draw_patient_cohort(world: Mapping, n: int, seed: int) -> dict[str, np.ndarray]:
    """Draw one independent synthetic cohort using the declared patient-axis order.

    The simulator's ``eps`` field implements the protocol's ``epsilon`` axis.
    For mixture-response worlds, ``m`` is drawn from the world-scaled weak range
    for nonresponders and from the base range for responders. The responder flag
    is not included in the returned scenario, so it cannot be given to a policy.
    """
    seed = _seed(seed, "seed")
    if not isinstance(world, Mapping):
        raise TypeError("world must be a mapping returned by draw_world_parameters")
    if world.get("family") not in ("single", "mixture"):
        raise ValueError("world family must be single or mixture")
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, Integral) or n < 1:
        raise ValueError("n must be a positive integer")

    r_scale = float(np.exp(world["log_r_shift"]))
    f0_scale = float(np.exp(world["log_f0_shift"]))
    f0_low = tc.RANGES["f0"][0] * f0_scale
    f0_high = min(tc.RANGES["f0"][1] * f0_scale, 0.9)
    if not (0.0 < f0_low < f0_high < 1.0):
        raise ValueError("world produced invalid f0 bounds")
    b0_low = float(world["B0_center"]) - B0_CENTER_HALF_WIDTH
    b0_high = float(world["B0_center"]) + B0_CENTER_HALF_WIDTH
    if not (0.0 < b0_low < b0_high < 1.0):
        raise ValueError("world produced invalid B0 bounds")

    if world["family"] == "single":
        m_scale = float(np.exp(world["log_m_scale"]))
    else:
        m_scale = float(np.exp(world["log_m_scale_non"]))
    m_low, m_high = tc.RANGES["m"][0] * m_scale, tc.RANGES["m"][1] * m_scale
    rng = make_generator(seed)
    unit = qmc.LatinHypercube(d=len(PATIENT_DIMENSION_ORDER), scramble=True, seed=rng).random(int(n))

    ranges = {
        "r": (tc.RANGES["r"][0] * r_scale, tc.RANGES["r"][1] * r_scale),
        "m": (m_low, m_high),
        "f0": (f0_low, f0_high),
        "B0": (b0_low, b0_high),
        "c": tc.RANGES["c"][:2],
        "epsilon": tc.RANGES["eps"][:2],
    }
    scenarios: dict[str, np.ndarray] = {}
    for column, name in enumerate(PATIENT_DIMENSION_ORDER):
        low, high = map(float, ranges[name])
        if name in ("r", "m", "f0"):
            scenarios[name if name != "epsilon" else "eps"] = np.exp(
                np.log(low) + unit[:, column] * (np.log(high) - np.log(low))
            )
        else:
            output_name = "eps" if name == "epsilon" else name
            scenarios[output_name] = low + unit[:, column] * (high - low)

    if world["family"] == "mixture":
        pi_response = float(world["pi_response"])
        if not 0.0 <= pi_response <= 1.0:
            raise ValueError("pi_response must be in [0, 1]")
        responders = rng.random(int(n)) < pi_response
        m_response = np.exp(
            np.log(tc.RANGES["m"][0])
            + rng.random(int(n)) * (np.log(tc.RANGES["m"][1]) - np.log(tc.RANGES["m"][0]))
        )
        scenarios["m"] = np.where(responders, m_response, scenarios["m"])

    scenarios["mu"] = np.full(int(n), float(world["mu"]))
    scenarios["nu"] = np.full(int(n), float(world["nu"]))
    scenarios["a_sr"] = np.ones(int(n), dtype=float)
    scenarios["a_rs"] = np.ones(int(n), dtype=float)
    return scenarios


def assign_balanced_sigma(scenarios: Mapping, patient_ids: Mapping, *, world_seed: int) -> tuple[dict, dict]:
    """Assign sigma levels once across all sorted IDs from the three world cohorts."""
    world_seed = _seed(world_seed, "world_seed")
    if set(scenarios) != set(COHORT_ROLES) or set(patient_ids) != set(COHORT_ROLES):
        raise ValueError(f"scenarios and patient_ids must contain exactly {COHORT_ROLES}")
    rows = []
    seen = set()
    for role in COHORT_ROLES:
        ids = tuple(patient_ids[role])
        n = len(np.asarray(scenarios[role]["r"]))
        if len(ids) != n:
            raise ValueError(f"{role} IDs must match its scenario length")
        for index, patient_id in enumerate(ids):
            if not isinstance(patient_id, str) or not patient_id.strip() or patient_id in seen:
                raise ValueError("patient IDs must be unique, nonempty strings across the world")
            seen.add(patient_id)
            rows.append((patient_id, role, index))

    ordered = sorted(rows, key=lambda row: row[0])
    base = np.resize(np.asarray(SIGMA_LEVELS, dtype=float), len(ordered))
    child = np.random.SeedSequence([world_seed, SIGMA_ASSIGNMENT_CHILD])
    rng = np.random.Generator(np.random.PCG64DXSM(child))
    assigned = rng.permutation(base)
    by_id = {row[0]: float(value) for row, value in zip(ordered, assigned)}

    result = {}
    for role in COHORT_ROLES:
        result[role] = {key: np.asarray(value).copy() for key, value in scenarios[role].items()}
        result[role]["sigma"] = np.asarray([by_id[patient_id] for patient_id in patient_ids[role]], dtype=float)
    counts = {str(level): int(np.count_nonzero(assigned == level)) for level in SIGMA_LEVELS}
    record = {
        "levels": list(SIGMA_LEVELS),
        "counts": counts,
        "assignment_unit": "all sorted patient IDs across selection, calibration, and evaluation cohorts within one world",
        "generator": "numpy.random.Generator(PCG64DXSM)",
        "child_entropy": [world_seed, SIGMA_ASSIGNMENT_CHILD],
        "balanced_count_difference_max": max(counts.values()) - min(counts.values()),
    }
    return result, record
