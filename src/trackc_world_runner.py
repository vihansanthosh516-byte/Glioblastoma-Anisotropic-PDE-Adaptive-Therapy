"""Development-only orchestration prototype; no seeds are reserved here."""
from __future__ import annotations

import hashlib
import json
from numbers import Integral

import numpy as np
import trackc_integrated_sampler as integrated_sampler
from trackc_filter_priors import PF_PRIOR_INTEGRATED_V3
from trackc_integrated_signal_model import SignalModelParameters, generate_integrated_signals
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
                            model_parameters, settings, scan_loss=None, entries=None):
    """Build one development world's benchmark inputs without scoring policies.

    The returned world seed is caller-declared provenance, not proof of generation.
    Seed preflight is not an atomic reservation. No files are read or written.
    Integrated v2 requires the observed two-scan treatment-stop rule; reject any
    other setting before generating world or patient inputs.
    """
    validate_development_schedule([world], entries=entries)
    if family not in ("single", "mixture"):
        raise ValueError("family must be single or mixture")
    if getattr(settings, "stop_endpoint", None) != "observed_two_scan":
        raise ValueError("integrated v2 requires stop_endpoint='observed_two_scan'")
    if isinstance(patients_per_cohort, bool) or not isinstance(patients_per_cohort, Integral) or patients_per_cohort < 2:
        raise ValueError("patients_per_cohort must be an integer of at least 2")
    n = int(patients_per_cohort)
    latent_world = integrated_sampler.draw_world_parameters(family, seed=world["world"])
    scenarios = {
        role: integrated_sampler.draw_patient_cohort(
            latent_world, n, seed=world[f"{role}_scenario"]
        )
        for role in integrated_sampler.COHORT_ROLES
    }
    patient_ids = {
        role: tuple(f"{world['world_id']}-{role}-{i:06d}" for i in range(n))
        for role in integrated_sampler.COHORT_ROLES
    }
    scenarios, sigma_assignment = integrated_sampler.assign_balanced_sigma(
        scenarios, patient_ids, world_seed=world["world"]
    )
    loss = None
    if scan_loss is not None:
        if not isinstance(scan_loss, dict) or "miss_p" not in scan_loss or set(scan_loss) - {"miss_p", "miss_p_hi"}:
            raise ValueError("scan_loss must contain miss_p and optionally miss_p_hi")
        loss = {}
        for name, value in scan_loss.items():
            if isinstance(value, (bool, np.bool_)):
                raise ValueError(f"{name} must be a finite probability")
            try:
                probability = float(value)
            except (TypeError, ValueError, OverflowError) as exc:
                raise ValueError(f"{name} must be a finite probability") from exc
            if not np.isfinite(probability) or not 0.0 <= probability <= 1.0:
                raise ValueError(f"{name} must be a finite probability")
            loss[name] = probability
        for scenario in scenarios.values():
            for name, probability in loss.items():
                scenario[name] = np.full(n, probability, dtype=float)
            scenario["miss_visible"] = np.ones(n, dtype=bool)
    observations = {}
    effective_signal_parameters = SignalModelParameters(
        q=latent_world["q_signal"], a=signal_parameters.a,
        b=signal_parameters.b, d=signal_parameters.d,
    )
    for role in ("calibration", "evaluation"):
        observations[role], _ = generate_integrated_signals(
            n, q=effective_signal_parameters.q, a=effective_signal_parameters.a,
            b=effective_signal_parameters.b, d=effective_signal_parameters.d,
            seed=None, rng=integrated_sampler.make_generator(world[f"{role}_signal"]),
            f0=scenarios[role]["f0"],
        )
    kwargs = {
        "world_id": world["world_id"],
        **{f"{role}_scenarios": scenarios[role] for role in scenarios},
        **{f"{role}_patient_ids": patient_ids[role] for role in scenarios},
        "calibration_observations": observations["calibration"],
        "evaluation_observations": observations["evaluation"],
        "model_parameters": tuple(model_parameters), "signal_generation_parameters": effective_signal_parameters,
        "signal_generation_seed": world["evaluation_signal"],
        "calibration_signal_generation_seed": world["calibration_signal"],
        "signal_generation_mri_template_sha256": None,
        **{f"{role}_scenario_seed": world[f"{role}_scenario"] for role in scenarios},
        "selection_seed": [world["selection_policy"]], "q_training_seed": [world["q_training"]],
        "calibration_seed": [world["calibration_policy"]], "evaluation_seed": [world["evaluation_policy"]],
        "schedule_mixture_seed": [world["schedule_mixture"]], "settings": settings,
        "include_ambiguity_controller": True,
        "primary_contrast": "minimax",
        "filter_ranges": PF_PRIOR_INTEGRATED_V3,
        "include_single_signal_controls": True,
    }
    world_bytes = json.dumps(
        {"family": family, "parameters": latent_world}, sort_keys=True,
        separators=(",", ":"), allow_nan=False,
    ).encode("utf-8")
    return {"world_generation_seed": world["world"], "family": family, "scan_loss": loss,
            "world_parameters_sha256": hashlib.sha256(world_bytes).hexdigest(),
            "world_parameters": latent_world, "sigma_assignment": sigma_assignment,
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
        "scan_loss": prepared["scan_loss"],
        "world_parameters_sha256": prepared["world_parameters_sha256"],
        "world_parameters": prepared["world_parameters"],
        "sigma_assignment": prepared["sigma_assignment"],
        "verification": "caller-declared builder provenance; not independently verified",
    }}


def make_synthetic_executor(*, family, patients_per_cohort, signal_parameters,
                            model_parameters, settings, benchmark, scan_loss=None, entries=None):
    """Compose input preparation and a supplied development benchmark.

    No protocol is adopted, seeds reserved, or results written by this factory.
    """
    if not callable(benchmark):
        raise TypeError("benchmark must be callable")

    def execute(world):
        prepared = prepare_synthetic_world(
            world, family=family, patients_per_cohort=patients_per_cohort,
            signal_parameters=signal_parameters, model_parameters=model_parameters,
            settings=settings, scan_loss=scan_loss, entries=entries,
        )
        return execute_prepared_world(prepared, benchmark)

    return execute


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
            if "world_provenance" in result and result["world_provenance"]["world_generation_seed"] != world["world"]:
                raise ValueError("recorded world generation seed does not match schedule")
            results.append(result)
        except Exception as exc:
            raise WorldRunError(f"world {world['world_id']} failed after {len(results)} completed worlds") from exc
    return results
