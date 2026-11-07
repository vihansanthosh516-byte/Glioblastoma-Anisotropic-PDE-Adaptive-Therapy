"""Development-only action-sampling stability diagnostic for finite-set MPC.

This diagnostic records chosen plans and Monte Carlo score margins only. It does
not simulate treatment outcomes or compare patient benefit.
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


ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "output" / "trackc_minimax_action_stability_20261106"
MANIFEST_PATH = OUT_DIR / "manifest.json"
RESULT_PATH = OUT_DIR / "diagnostic.json"
PROTOCOL_PATH = ROOT / "paper" / "trackc_minimax_action_stability_2026-11-06.md"
SCENARIO_SEED = 20261990
FILTER_INIT_SEED = 930000004
PREDICTION_SEEDS = (930000001, 930000002, 930000003)
PREDICTION_LEVELS = (16, 32, 64)
N_PARTICLES = 128
SCAN_EVERY = 3
FORECAST_CYCLES = 3
ANCHORS = (
    ("both_positive", True, 0.75),
    ("molecular_only_positive", True, -0.75),
    ("mri_only_positive", False, 0.75),
    ("both_negative", False, -0.75),
)
SOURCE_FILES = (
    "src/219_trackc_minimax_action_stability.py",
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


_SIM = _load_sibling("trackc_sim_for_action_stability", "trackc_sim.py")
_AMBIGUITY = _load_sibling("trackc_ambiguity_for_action_stability", "216_trackc_robust_signal_ambiguity.py")


def build_model_sets():
    """Return declared 81-model grid and its matched 27-model d=0 subset."""
    levels = (-0.5, 0.0, 0.5)
    full = tuple(
        _AMBIGUITY.SignalModelParameters(q, a, b, d)
        for q, a, b, d in itertools.product((0.5, 0.7, 0.9), levels, levels, levels)
    )
    d0 = tuple(model for model in full if model.d == 0.0)
    if len(full) != 81 or len(d0) != 27:
        raise ArithmeticError("declared finite model grids have unexpected sizes")
    return full, d0


def _source_hashes():
    paths = (*SOURCE_FILES, str(PROTOCOL_PATH.relative_to(ROOT)).replace("\\", "/"))
    hashes = {}
    for relative in paths:
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(f"required protocol source is missing: {relative}")
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


def build_manifest():
    """Build pre-run manifest and check that all development seeds are free."""
    registry = _seed_registry_module()
    scenario_problems = registry.check_free(
        scenario=(SCENARIO_SEED, SCENARIO_SEED), purpose="stage2_dev"
    )
    policy_problems = registry.check_free(
        policy=(min((*PREDICTION_SEEDS, FILTER_INIT_SEED)), max((*PREDICTION_SEEDS, FILTER_INIT_SEED)))
    )
    if scenario_problems or policy_problems:
        raise ValueError(f"seed allocation failed: {scenario_problems + policy_problems}")
    if any(entry["script"] == "219" for entry in registry.ENTRIES):
        raise ValueError("script 219 seeds already appear in the registry")
    full, d0 = build_model_sets()
    return {
        "protocol_status": "prepared",
        "diagnostic": "synthetic controller action-sampling stability; no outcome comparison",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_hashes": _source_hashes(),
        "settings": {
            "scenario_seed": SCENARIO_SEED,
            "filter_initialization_policy_seed": FILTER_INIT_SEED,
            "prediction_seeds": list(PREDICTION_SEEDS),
            "prediction_particles": list(PREDICTION_LEVELS),
            "posterior_particles": N_PARTICLES,
            "scan_every": SCAN_EVERY,
            "forecast_cycles": FORECAST_CYCLES,
            "anchor_count": len(ANCHORS),
            "anchors": [
                {"name": name, "molecular_report": report, "mri_score": score}
                for name, report, score in ANCHORS
            ],
            "full_model_count": len(full),
            "matched_d0_model_count": len(d0),
            "model_grid": {
                "q": [0.5, 0.7, 0.9],
                "a": [-0.5, 0.0, 0.5],
                "b": [-0.5, 0.0, 0.5],
                "d": [-0.5, 0.0, 0.5],
                "matched_d0": [0.0],
            },
            "baseline_burden_reading": "exact log(B0) from one synthetic scenario; no observation noise",
        },
        "seed_allocation_check": {
            "purpose": "stage2_dev",
            "scenario_seed_free_at_prepare": True,
            "policy_seed_range_free_at_prepare": True,
        },
        "interpretation_rule": (
            "Report exact plan agreement between adjacent particle counts and across prediction seeds. "
            "Any plan change from 32 to 64 marks this anchor/model-set diagnostic unstable at the tested setting. "
            "No positive or negative treatment-benefit claim is allowed."
        ),
        "outcome_metrics_included": False,
    }


def summarize_action_rows(rows):
    """Summarize paired plan agreement across adjacent draw counts and RNG roots."""
    if not rows:
        raise ValueError("rows must be nonempty")
    required = {"model_set", "anchor", "prediction_seed", "prediction_particles", "pattern"}
    seen = set()
    grouped = {}
    levels = set()
    model_sets = set()
    for row in rows:
        if not required <= set(row):
            raise ValueError("row is missing required action-stability fields")
        key = (row["model_set"], row["anchor"], row["prediction_seed"], row["prediction_particles"])
        if key in seen:
            raise ValueError("duplicate model/anchor/seed/particle row")
        seen.add(key)
        pattern = row["pattern"]
        if not isinstance(pattern, str) or not pattern or set(pattern) - {"0", "1"}:
            raise ValueError("pattern must be a nonempty binary string")
        levels.add(int(row["prediction_particles"]))
        model_sets.add(str(row["model_set"]))
        grouped.setdefault((row["model_set"], row["anchor"], row["prediction_seed"]), {})[
            int(row["prediction_particles"])
        ] = pattern
    ordered_levels = sorted(levels)
    if len(ordered_levels) < 2:
        raise ValueError("at least two particle levels are required")
    if len(model_sets) < 1:
        raise ValueError("at least one model set is required")
    pattern_lengths = {len(row["pattern"]) for row in rows}
    if len(pattern_lengths) != 1:
        raise ValueError("all patterns must have the same length")

    by_model = {}
    for model_set in sorted(model_sets):
        comparisons = []
        for low, high in zip(ordered_levels, ordered_levels[1:]):
            paired = [
                pair[low] == pair[high]
                for (candidate, _anchor, _seed), pair in grouped.items()
                if candidate == model_set and low in pair and high in pair
            ]
            comparisons.append((sum(paired), len(paired)))
        by_model[model_set] = comparisons

    adjacent = []
    for model_set, comparisons in by_model.items():
        adjacent.append({
            "model_set": model_set,
            "particle_pairs": [[a, b] for a, b in zip(ordered_levels, ordered_levels[1:])],
            "agreements": [same for same, _n in comparisons],
            "n_pairs": [n for _same, n in comparisons],
            "fraction_by_step": [same / n if n else None for same, n in comparisons],
        })

    shared = {}
    for row in rows:
        shared.setdefault((row["anchor"], row["prediction_seed"], row["prediction_particles"]), {})[
            row["model_set"]
        ] = row["pattern"]
    full_vs_d0 = [
        (values["full_81"] == values["matched_d0_27"])
        for values in shared.values()
        if "full_81" in values and "matched_d0_27" in values
    ]
    return {
        "particle_levels": ordered_levels,
        "adjacent_agreement": adjacent,
        "full_vs_d0_pattern_agreement": {
            "agreements": sum(full_vs_d0),
            "n_pairs": len(full_vs_d0),
            "fraction": sum(full_vs_d0) / len(full_vs_d0) if full_vs_d0 else None,
        },
        "outcome_metrics_included": False,
    }


def _capture_plan(model_parameters, observations, scenario, prediction_particles, prediction_seed):
    policy = _AMBIGUITY.SignalAmbiguityBeliefMPC(
        observations,
        model_parameters=model_parameters,
        scan_every=SCAN_EVERY,
        n_particles=N_PARTICLES,
        prediction_particles=prediction_particles,
        forecast_cycles=FORECAST_CYCLES,
    )
    policy.reset(scenario, np.random.default_rng(FILTER_INIT_SEED))
    policy._prediction_seed = prediction_seed
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
    values = captured["worst_case_values"]
    ordered = sorted(values, reverse=True)
    return {
        "pattern": "".join("1" if value else "0" for value in policy.plan[0, :FORECAST_CYCLES]),
        "worst_case_value_margin": ordered[0] - ordered[1] if len(ordered) > 1 else None,
        "worst_case_values_by_pattern": values,
        "elapsed_seconds": elapsed,
    }


def run_diagnostic(manifest):
    """Run fixed synthetic posterior cases across model sets, K levels and draw roots."""
    expected = build_manifest_unchecked()
    if manifest.get("protocol_status") != "prepared":
        raise ValueError("pre-run manifest is not in prepared state")
    if manifest.get("source_hashes") != expected["source_hashes"]:
        raise ValueError("source or protocol hash changed after pre-run manifest")
    if manifest.get("settings") != expected["settings"]:
        raise ValueError("diagnostic settings changed after pre-run manifest")
    registry = _seed_registry_module()
    registered = [entry for entry in registry.ENTRIES if entry["script"] == "219"]
    if len(registered) != 1 or registered[0]["status"] != "dev_used":
        raise ValueError("script 219 development seeds must be registered before running")
    if registered[0]["scenario_seeds"] != [[SCENARIO_SEED, SCENARIO_SEED]]:
        raise ValueError("registered scenario seed does not match the frozen manifest")
    if registered[0]["policy_seeds"] != [[min((*PREDICTION_SEEDS, FILTER_INIT_SEED)), max((*PREDICTION_SEEDS, FILTER_INIT_SEED))]]:
        raise ValueError("registered policy seeds do not match the frozen manifest")

    full, d0 = build_model_sets()
    model_sets = (("full_81", full), ("matched_d0_27", d0))
    scenario = _SIM.draw_scenarios(1, seed=SCENARIO_SEED)
    rows = []
    for anchor, report, mri_score in ANCHORS:
        observations = _AMBIGUITY.SignalObservations(
            molecular_report=np.asarray([report], dtype=bool),
            mri_score=np.asarray([mri_score], dtype=float),
        )
        for model_name, models in model_sets:
            for prediction_seed in PREDICTION_SEEDS:
                for particle_count in PREDICTION_LEVELS:
                    result = _capture_plan(
                        models,
                        observations,
                        scenario,
                        particle_count,
                        prediction_seed,
                    )
                    rows.append({
                        "model_set": model_name,
                        "anchor": anchor,
                        "prediction_seed": prediction_seed,
                        "prediction_particles": particle_count,
                        **result,
                    })
    return {
        "protocol_status": "completed",
        "run_type": "synthetic_action_stability_only",
        "scenario_seed": SCENARIO_SEED,
        "filter_initialization_policy_seed": FILTER_INIT_SEED,
        "prediction_seeds": list(PREDICTION_SEEDS),
        "settings": manifest["settings"],
        "manifest_sha256": _canonical_sha256(manifest),
        "source_hashes": _source_hashes(),
        "seed_registry_sha256_at_run": hashlib.sha256(
            (ROOT / "src" / "trackc_seed_registry.py").read_bytes()
        ).hexdigest(),
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
        "rows": rows,
        "summary": summarize_action_rows(rows),
        "outcome_metrics_included": False,
        "limits": [
            "One synthetic patient state; four analyst-chosen signal anchors.",
            "The signal-to-resistance relationship is assumed by the simulator.",
            "This diagnostic checks Monte Carlo action stability only, not treatment benefit or model validity.",
        ],
    }


def build_manifest_unchecked():
    """Build manifest fields without rechecking seeds after their reservation is spent."""
    full, d0 = build_model_sets()
    return {
        "source_hashes": _source_hashes(),
        "settings": {
            "scenario_seed": SCENARIO_SEED,
            "filter_initialization_policy_seed": FILTER_INIT_SEED,
            "prediction_seeds": list(PREDICTION_SEEDS),
            "prediction_particles": list(PREDICTION_LEVELS),
            "posterior_particles": N_PARTICLES,
            "scan_every": SCAN_EVERY,
            "forecast_cycles": FORECAST_CYCLES,
            "anchor_count": len(ANCHORS),
            "anchors": [
                {"name": name, "molecular_report": report, "mri_score": score}
                for name, report, score in ANCHORS
            ],
            "full_model_count": len(full),
            "matched_d0_model_count": len(d0),
            "model_grid": {
                "q": [0.5, 0.7, 0.9],
                "a": [-0.5, 0.0, 0.5],
                "b": [-0.5, 0.0, 0.5],
                "d": [-0.5, 0.0, 0.5],
                "matched_d0": [0.0],
            },
            "baseline_burden_reading": "exact log(B0) from one synthetic scenario; no observation noise",
        },
    }


def _write_json_new(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def prepare():
    manifest = build_manifest()
    _write_json_new(MANIFEST_PATH, manifest)
    print(f"Prepared manifest: {MANIFEST_PATH}")


def run():
    if not MANIFEST_PATH.is_file():
        raise FileNotFoundError("pre-run manifest is required")
    if RESULT_PATH.exists():
        raise FileExistsError(f"refusing to overwrite {RESULT_PATH}")
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    result = run_diagnostic(manifest)
    _write_json_new(RESULT_PATH, result)
    print(f"Wrote diagnostic: {RESULT_PATH}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare", action="store_true", help="write pre-run manifest; run no diagnostic")
    mode.add_argument("--run", action="store_true", help="run diagnostic from frozen pre-run manifest")
    args = parser.parse_args(argv)
    prepare() if args.prepare else run()


if __name__ == "__main__":
    main()
