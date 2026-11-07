"""Development-only test of whether paired synthetic signals change minimax plans.

This script records first-step action plans only. It does not simulate treatment
outcomes, estimate patient benefit, or validate a biological signal link.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import itertools
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy

import trackc_sim
from trackc_integrated_signal_model import (
    SignalModelParameters,
    SignalObservations,
    generate_integrated_signals,
)


ROOT = Path(__file__).resolve().parent.parent
SCRIPT_ID = "220"
PROJECT_TITLE = "Integrating Tumor Biology and Longitudinal MRI for Glioblastoma Progression Modeling and In-Silico Therapy Optimization"
INTERPRETATION_RULE = (
    "Report whether action patterns vary across generated paired-signal replicates and whether full_81 "
    "and matched_d0_27 agree. No pass rule, method selection, treatment outcome, or positive claim is allowed."
)
OUT_DIR = ROOT / "output" / "trackc_signal_action_sensitivity_20261106"
MANIFEST_PATH = OUT_DIR / "manifest.json"
RESULT_PATH = OUT_DIR / "diagnostic.json"
PROTOCOL_PATH = ROOT / "paper" / "trackc_signal_action_sensitivity_protocol_2026-11-06.md"

STATE_SEEDS = (20262000, 20262001, 20262002, 20262003)
SIGNAL_SEEDS = tuple(range(20262004, 20262068))
FILTER_INIT_SEED = 930001001
PREDICTION_SEED = 930001002
SIGNAL_REPLICATES_PER_STATE = 16
N_PARTICLES = 128
PREDICTION_PARTICLES = 64
SCAN_EVERY = 3
FORECAST_CYCLES = 3
TRUTH_PARAMETERS = SignalModelParameters(q=0.7, a=0.75, b=0.75, d=0.75)
MODEL_LEVELS = (-0.5, 0.0, 0.5)
MODEL_Q = (0.5, 0.7, 0.9)
SOURCE_FILES = (
    "src/220_trackc_signal_action_sensitivity.py",
    "src/216_trackc_robust_signal_ambiguity.py",
    "src/208_trackc_robust_belief_mpc.py",
    "src/trackc_integrated_signal_model.py",
    "src/trackc_sim.py",
)


def _load_sibling(name: str, filename: str):
    path = ROOT / "src" / filename
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {filename}: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


_AMBIGUITY = _load_sibling(
    "trackc_ambiguity_for_signal_action_sensitivity",
    "216_trackc_robust_signal_ambiguity.py",
)


def build_design() -> list[dict]:
    """Return the frozen 4-state by 16-signal-replicate design."""
    rows = []
    seed_index = 0
    for state_index in range(len(STATE_SEEDS)):
        for signal_replicate in range(SIGNAL_REPLICATES_PER_STATE):
            rows.append({
                "state_index": state_index,
                "signal_replicate": signal_replicate,
                "signal_seed": SIGNAL_SEEDS[seed_index],
            })
            seed_index += 1
    if seed_index != len(SIGNAL_SEEDS):
        raise ArithmeticError("signal seed count does not match the frozen design")
    return rows


def build_model_sets():
    """Return the frozen 81-model set and its matched d=0 subset."""
    full = tuple(
        SignalModelParameters(q, a, b, d)
        for q, a, b, d in itertools.product(MODEL_Q, MODEL_LEVELS, MODEL_LEVELS, MODEL_LEVELS)
    )
    matched_d0 = tuple(model for model in full if model.d == 0.0)
    if len(full) != 81 or len(matched_d0) != 27:
        raise ArithmeticError("frozen signal model grids have unexpected sizes")
    return full, matched_d0


def _source_hashes() -> dict[str, str]:
    paths = (*SOURCE_FILES, PROTOCOL_PATH.relative_to(ROOT).as_posix())
    hashes = {}
    for relative in paths:
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(f"required diagnostic source is missing: {relative}")
        hashes[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def _canonical_sha256(value) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _seed_registry_module():
    if str(ROOT / "src") not in sys.path:
        sys.path.insert(0, str(ROOT / "src"))
    import trackc_seed_registry
    return trackc_seed_registry


def build_manifest() -> dict:
    """Build a pre-run manifest and verify fresh Stage 2 development seeds."""
    registry = _seed_registry_module()
    scenario_range = (min((*STATE_SEEDS, *SIGNAL_SEEDS)), max((*STATE_SEEDS, *SIGNAL_SEEDS)))
    policy_range = (min(FILTER_INIT_SEED, PREDICTION_SEED), max(FILTER_INIT_SEED, PREDICTION_SEED))
    seed_problems = registry.check_free(scenario=scenario_range, purpose="stage2_dev")
    seed_problems += registry.check_free(policy=policy_range)
    if seed_problems:
        raise ValueError(f"seed allocation failed: {seed_problems}")
    if any(entry["script"] == SCRIPT_ID for entry in registry.ENTRIES):
        raise ValueError(f"script {SCRIPT_ID} already appears in the seed registry")
    full, matched_d0 = build_model_sets()
    return {
        "record_type": "pre_run_synthetic_action_sensitivity_manifest",
        "protocol_status": "prepared",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script_id": SCRIPT_ID,
        "project_title": PROJECT_TITLE,
        "source_hashes": _source_hashes(),
        "settings": {
            "state_seeds": list(STATE_SEEDS),
            "signal_seeds": list(SIGNAL_SEEDS),
            "filter_initialization_seed": FILTER_INIT_SEED,
            "prediction_seed": PREDICTION_SEED,
            "signal_replicates_per_state": SIGNAL_REPLICATES_PER_STATE,
            "synthetic_state_count": len(STATE_SEEDS),
            "posterior_particles": N_PARTICLES,
            "prediction_particles": PREDICTION_PARTICLES,
            "scan_every": SCAN_EVERY,
            "forecast_cycles": FORECAST_CYCLES,
            "truth_signal_parameters": {
                "q": TRUTH_PARAMETERS.q,
                "a": TRUTH_PARAMETERS.a,
                "b": TRUTH_PARAMETERS.b,
                "d": TRUTH_PARAMETERS.d,
            },
            "full_model_count": len(full),
            "matched_d0_model_count": len(matched_d0),
            "full_model_grid": {
                "q": list(MODEL_Q),
                "a": list(MODEL_LEVELS),
                "b": list(MODEL_LEVELS),
                "d": list(MODEL_LEVELS),
            },
            "matched_d0_grid": {"q": list(MODEL_Q), "a": list(MODEL_LEVELS), "b": list(MODEL_LEVELS), "d": [0.0]},
            "baseline_burden_reading": "exact log(B0); no burden measurement noise",
        },
        "seed_allocation_check": {
            "purpose": "stage2_dev",
            "scenario_range_free_at_prepare": True,
            "policy_range_free_at_prepare": True,
        },
        "interpretation_rule": INTERPRETATION_RULE,
        "outcome_metrics_included": False,
    }


def summarize_action_rows(rows: list[dict]) -> dict:
    """Summarize action diversity and paired model-set agreement without selecting a winner."""
    if not rows:
        raise ValueError("rows must be nonempty")
    required = {"state_id", "signal_replicate", "model_set", "pattern"}
    seen = set()
    grouped = {}
    for row in rows:
        if not required <= set(row):
            raise ValueError("row is missing required action-sensitivity fields")
        pattern = row["pattern"]
        if not isinstance(pattern, str) or not pattern or set(pattern) - {"0", "1"}:
            raise ValueError("pattern must be a nonempty binary string")
        if row["model_set"] not in ("full_81", "matched_d0_27"):
            raise ValueError("unknown model set")
        key = (str(row["state_id"]), int(row["signal_replicate"]), str(row["model_set"]))
        if key in seen:
            raise ValueError("duplicate state/signal/model-set row")
        seen.add(key)
        grouped.setdefault((key[0], key[1]), {})[key[2]] = pattern
    if any(set(pair) != {"full_81", "matched_d0_27"} for pair in grouped.values()):
        raise ValueError("every signal replicate must have paired model-set rows")
    if len({len(row["pattern"]) for row in rows}) != 1:
        raise ValueError("all patterns must have equal length")

    by_state = []
    for state_id in sorted({key[0] for key in grouped}):
        pairs = [(key, value) for key, value in grouped.items() if key[0] == state_id]
        pattern_counts = {}
        for model_set in ("full_81", "matched_d0_27"):
            counts = {}
            for _key, pair in pairs:
                pattern = pair[model_set]
                counts[pattern] = counts.get(pattern, 0) + 1
            pattern_counts[model_set] = dict(sorted(counts.items()))
        by_state.append({
            "state_id": state_id,
            "n_signal_replicates": len(pairs),
            "pattern_counts": pattern_counts,
            "unique_pattern_count": {
                model_set: len(pattern_counts[model_set])
                for model_set in ("full_81", "matched_d0_27")
            },
        })
    agreements = [
        pair["full_81"] == pair["matched_d0_27"]
        for pair in grouped.values()
    ]
    return {
        "n_rows": len(rows),
        "n_paired_signal_cases": len(grouped),
        "by_state": by_state,
        "full_vs_d0_agreement": {
            "agreements": int(sum(agreements)),
            "n_pairs": len(agreements),
            "fraction": float(np.mean(agreements)),
        },
        "outcome_metrics_included": False,
    }


def json_safe(value):
    """Convert NumPy values to strict JSON values; map non-finite floats to null."""
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError(f"cannot convert {type(value).__name__} to JSON")


def environment_record() -> dict:
    """Record runtime versions for numerical packages used by the diagnostic."""
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
    }


def _capture_plan(model_parameters, observations, scenario) -> dict:
    policy = _AMBIGUITY.SignalAmbiguityBeliefMPC(
        observations,
        model_parameters=model_parameters,
        scan_every=SCAN_EVERY,
        n_particles=N_PARTICLES,
        prediction_particles=PREDICTION_PARTICLES,
        forecast_cycles=FORECAST_CYCLES,
    )
    policy.reset(scenario, np.random.default_rng(FILTER_INIT_SEED))
    policy._prediction_seed = PREDICTION_SEED
    baseline_reading = np.log(np.maximum(np.asarray(scenario["B0"], dtype=float), 1e-9))
    captured = {}
    original = _AMBIGUITY.maximin_action_index

    def capture_action_values(model_action_values):
        values = np.asarray(model_action_values, dtype=float)
        captured["worst_case_values"] = np.min(values, axis=0).tolist()
        return original(values)

    _AMBIGUITY.maximin_action_index = capture_action_values
    started = time.perf_counter()
    try:
        policy.act(0, baseline_reading)
    finally:
        elapsed = time.perf_counter() - started
        _AMBIGUITY.maximin_action_index = original
    if "worst_case_values" not in captured:
        raise RuntimeError("controller did not produce a plan")
    ordered = sorted(captured["worst_case_values"], reverse=True)
    return {
        "pattern": "".join("1" if value else "0" for value in policy.plan[0, :FORECAST_CYCLES]),
        "worst_case_value_margin": ordered[0] - ordered[1] if len(ordered) > 1 else None,
        "elapsed_seconds": elapsed,
    }


def _registered_seed_check(registry) -> None:
    registered = [entry for entry in registry.ENTRIES if entry["script"] == SCRIPT_ID]
    if len(registered) != 1 or registered[0]["status"] != "dev_used":
        raise ValueError("script 220 development seeds must be registered before running")
    if registered[0]["scenario_seeds"] != [[STATE_SEEDS[0], SIGNAL_SEEDS[-1]]]:
        raise ValueError("registered scenario seed range does not match the frozen design")
    expected_policy = [[min(FILTER_INIT_SEED, PREDICTION_SEED), max(FILTER_INIT_SEED, PREDICTION_SEED)]]
    if registered[0]["policy_seeds"] != expected_policy:
        raise ValueError("registered policy seeds do not match the frozen design")


def run_diagnostic(manifest: dict) -> dict:
    """Run frozen action-sensitivity design; do not calculate treatment outcomes."""
    validate_prepared_manifest(manifest)
    registry = _seed_registry_module()
    _registered_seed_check(registry)

    full, matched_d0 = build_model_sets()
    model_sets = (("full_81", full), ("matched_d0_27", matched_d0))
    scenarios = [trackc_sim.draw_scenarios(1, seed=seed) for seed in STATE_SEEDS]
    rows = []
    for design_row in build_design():
        state_index = design_row["state_index"]
        replicate = design_row["signal_replicate"]
        scenario = scenarios[state_index]
        f0 = np.asarray(scenario["f0"], dtype=float)
        observations, _truth = generate_integrated_signals(
            1,
            q=TRUTH_PARAMETERS.q,
            a=TRUTH_PARAMETERS.a,
            b=TRUTH_PARAMETERS.b,
            d=TRUTH_PARAMETERS.d,
            seed=design_row["signal_seed"],
            f0=f0,
        )
        state_id = f"state_{state_index}"
        for model_name, models in model_sets:
            result = _capture_plan(models, observations, scenario)
            rows.append({
                "state_id": state_id,
                "state_seed": STATE_SEEDS[state_index],
                "f0": float(f0[0]),
                "signal_replicate": replicate,
                "signal_seed": design_row["signal_seed"],
                "molecular_report": bool(observations.molecular_report[0]),
                "mri_score": float(observations.mri_score[0]),
                "model_set": model_name,
                **result,
            })
    if len(rows) != len(STATE_SEEDS) * SIGNAL_REPLICATES_PER_STATE * 2:
        raise ArithmeticError("diagnostic row count does not match frozen design")
    return {
        "record_type": "synthetic_signal_action_sensitivity_result",
        "protocol_status": "completed",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "project_title": manifest["project_title"],
        "manifest_sha256": _canonical_sha256(manifest),
        "source_hashes": _source_hashes(),
        "seed_registry_sha256_at_run": hashlib.sha256(
            (ROOT / "src" / "trackc_seed_registry.py").read_bytes()
        ).hexdigest(),
        "environment": environment_record(),
        "settings": manifest["settings"],
        "rows": rows,
        "summary": summarize_action_rows(rows),
        "outcome_metrics_included": False,
        "limits": [
            "Four synthetic virtual-patient states and one assumed signal-generating model.",
            "Exact baseline log(B0) has no measurement noise; only paired signal draws vary.",
            "Action sensitivity is not treatment benefit, biological validation, novelty, or clinical evidence.",
        ],
    }


def build_manifest_unchecked() -> dict:
    """Return frozen manifest fields after seed registration has occurred."""
    full, matched_d0 = build_model_sets()
    return {
        "source_hashes": _source_hashes(),
        "settings": {
            "state_seeds": list(STATE_SEEDS),
            "signal_seeds": list(SIGNAL_SEEDS),
            "filter_initialization_seed": FILTER_INIT_SEED,
            "prediction_seed": PREDICTION_SEED,
            "signal_replicates_per_state": SIGNAL_REPLICATES_PER_STATE,
            "synthetic_state_count": len(STATE_SEEDS),
            "posterior_particles": N_PARTICLES,
            "prediction_particles": PREDICTION_PARTICLES,
            "scan_every": SCAN_EVERY,
            "forecast_cycles": FORECAST_CYCLES,
            "truth_signal_parameters": {
                "q": TRUTH_PARAMETERS.q,
                "a": TRUTH_PARAMETERS.a,
                "b": TRUTH_PARAMETERS.b,
                "d": TRUTH_PARAMETERS.d,
            },
            "full_model_count": len(full),
            "matched_d0_model_count": len(matched_d0),
            "full_model_grid": {
                "q": list(MODEL_Q),
                "a": list(MODEL_LEVELS),
                "b": list(MODEL_LEVELS),
                "d": list(MODEL_LEVELS),
            },
            "matched_d0_grid": {"q": list(MODEL_Q), "a": list(MODEL_LEVELS), "b": list(MODEL_LEVELS), "d": [0.0]},
            "baseline_burden_reading": "exact log(B0); no burden measurement noise",
        },
    }


def validate_prepared_manifest(manifest: dict) -> None:
    """Reject edits to frozen design metadata before any synthetic planning call."""
    if not isinstance(manifest, dict):
        raise ValueError("prepared manifest must be an object")
    current = build_manifest_unchecked()
    expected = {
        "record_type": "pre_run_synthetic_action_sensitivity_manifest",
        "protocol_status": "prepared",
        "script_id": SCRIPT_ID,
        "project_title": PROJECT_TITLE,
        "source_hashes": current["source_hashes"],
        "settings": current["settings"],
        "seed_allocation_check": {
            "purpose": "stage2_dev",
            "scenario_range_free_at_prepare": True,
            "policy_range_free_at_prepare": True,
        },
        "interpretation_rule": INTERPRETATION_RULE,
        "outcome_metrics_included": False,
    }
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise ValueError(f"prepared manifest field changed: {key}")


def _write_json_new(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    path.write_text(json.dumps(json_safe(value), indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def prepare() -> None:
    _write_json_new(MANIFEST_PATH, build_manifest())
    print(f"Prepared manifest: {MANIFEST_PATH}")


def run() -> None:
    if not MANIFEST_PATH.is_file():
        raise FileNotFoundError("pre-run manifest is required")
    if RESULT_PATH.exists():
        raise FileExistsError(f"refusing to overwrite {RESULT_PATH}")
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    result = run_diagnostic(manifest)
    _write_json_new(RESULT_PATH, result)
    print(f"Wrote diagnostic: {RESULT_PATH}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare", action="store_true", help="write pre-run manifest; run no diagnostic")
    mode.add_argument("--run", action="store_true", help="run diagnostic from frozen pre-run manifest")
    args = parser.parse_args(argv)
    prepare() if args.prepare else run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
