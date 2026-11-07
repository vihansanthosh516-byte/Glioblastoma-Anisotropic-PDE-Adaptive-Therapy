"""Development-only orchestration prototype; no seeds are reserved here."""
from __future__ import annotations

import hashlib
import json
from numbers import Integral

import trackc_calibration as calibration
from trackc_integrated_signal_model import generate_integrated_signals
from trackc_world_schedule import validate_development_schedule

_SEED_MAP = {
    "selection_scenario": "selection_scenario_generation",
    "calibration_scenario": "calibration_scenario_generation",
    "evaluation_scenario": "evaluation_scenario_generation",
    "calibration_signal": "calibration_signal_generation",
    "selection_policy": "selection",
    "q_training": "q_training",
    "calibration_policy": "calibration_policy_simulation",
    "evaluation_policy": "evaluation",
    "schedule_mixture": "schedule_mixture_master",
}


class WorldRunError(RuntimeError):
    """A world failed; no partial result collection is returned."""


def prepare_synthetic_world(world, *, family, patients_per_cohort, signal_parameters,
                            model_parameters, settings, entries=None):
    """Build one development world's benchmark inputs without scoring policies.

    The returned world seed is caller-declared provenance, not proof of generation.
    Seed preflight is not an atomic reservation. No files are read or written.
    """
    validate_development_schedule([world], entries=entries)
    if family not in ("single", "mixture"):
        raise ValueError("family must be single or mixture")
    if isinstance(patients_per_cohort, bool) or not isinstance(patients_per_cohort, Integral) or patients_per_cohort < 2:
        raise ValueError("patients_per_cohort must be an integer of at least 2")
    n = int(patients_per_cohort)
    latent_world = calibration.draw_worlds(1, seed=world["world"], family=family)[0]
    scenarios = {
        role: calibration.draw_patients(latent_world, n, seed=[world[f"{role}_scenario"], 1])
        for role in ("selection", "calibration", "evaluation")
    }
    observations = {}
    for role in ("calibration", "evaluation"):
        observations[role], _ = generate_integrated_signals(
            n, q=signal_parameters.q, a=signal_parameters.a, b=signal_parameters.b,
            d=signal_parameters.d, seed=world[f"{role}_signal"], f0=scenarios[role]["f0"],
        )
    kwargs = {
        "world_id": world["world_id"],
        **{f"{role}_scenarios": scenarios[role] for role in scenarios},
        **{f"{role}_patient_ids": tuple(f"{world['world_id']}-{role}-{i:06d}" for i in range(n)) for role in scenarios},
        "calibration_observations": observations["calibration"],
        "evaluation_observations": observations["evaluation"],
        "model_parameters": tuple(model_parameters), "signal_generation_parameters": signal_parameters,
        "signal_generation_seed": world["evaluation_signal"],
        "calibration_signal_generation_seed": world["calibration_signal"],
        "signal_generation_mri_template_sha256": None,
        **{f"{role}_scenario_seed": world[f"{role}_scenario"] for role in scenarios},
        "selection_seed": [world["selection_policy"]], "q_training_seed": [world["q_training"]],
        "calibration_seed": [world["calibration_policy"]], "evaluation_seed": [world["evaluation_policy"]],
        "schedule_mixture_seed": [world["schedule_mixture"]], "settings": settings,
        "include_ambiguity_controller": True,
    }
    world_bytes = json.dumps(
        {"family": family, "parameters": latent_world}, sort_keys=True,
        separators=(",", ":"), allow_nan=False,
    ).encode("utf-8")
    return {"world_generation_seed": world["world"], "family": family,
            "world_parameters_sha256": hashlib.sha256(world_bytes).hexdigest(),
            "benchmark_kwargs": kwargs}


def execute_prepared_world(prepared, benchmark):
    """Run a caller-supplied benchmark and attach declared world provenance.

    This helper does not reserve seeds, seal inputs, or authorize outcome scoring.
    The attached hash identifies prepared parameters; it is not independent proof.
    """
    if not callable(benchmark):
        raise TypeError("benchmark must be callable")
    result = benchmark(**prepared["benchmark_kwargs"])
    if not isinstance(result, dict) or result.get("world_id") != prepared["benchmark_kwargs"]["world_id"]:
        raise ValueError("benchmark result does not match prepared world ID")
    if "world_provenance" in result:
        raise ValueError("benchmark result already contains world provenance")
    return {**result, "world_provenance": {
        "world_generation_seed": prepared["world_generation_seed"],
        "family": prepared["family"],
        "world_parameters_sha256": prepared["world_parameters_sha256"],
        "verification": "caller-declared builder provenance; not independently verified",
    }}


def run_development_worlds(worlds, executor, *, entries=None):
    """Preflight and execute supplied worlds sequentially in memory.

    The caller must freeze a protocol and manage seed reservation separately.
    Validation is not an atomic reservation and cannot prevent concurrent reuse.
    The executor is responsible for generating genuinely independent worlds.
    No files are written and no sealed execution is authorized by this helper.
    """
    worlds = tuple(dict(world) for world in worlds)
    validate_development_schedule(worlds, entries=entries)
    if not callable(executor):
        raise TypeError("executor must be callable")
    results = []
    for world in worlds:
        try:
            result = executor(dict(world))
            if result["world_id"] != world["world_id"]:
                raise ValueError("world ID does not match schedule")
            if result["analysis_type"] != "synthetic_virtual_world_development":
                raise ValueError("only synthetic development records are supported")
            for role, output_key in _SEED_MAP.items():
                if list(result["seeds"][output_key]) != [world[role]]:
                    raise ValueError(f"recorded seed for {role} does not match schedule")
            if result["signal_generation"]["seed"] != world["evaluation_signal"]:
                raise ValueError("recorded evaluation signal seed does not match schedule")
            results.append(result)
        except Exception as exc:
            raise WorldRunError(f"world {world['world_id']} failed after {len(results)} completed worlds") from exc
    return results
