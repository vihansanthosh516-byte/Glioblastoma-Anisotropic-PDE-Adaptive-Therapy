"""Development-only orchestration prototype; no seeds are reserved here."""
from __future__ import annotations

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
