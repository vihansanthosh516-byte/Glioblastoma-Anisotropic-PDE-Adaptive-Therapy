from __future__ import annotations

import importlib.util
import hashlib
import json
import math
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "src"
sys.path.insert(0, str(SRC))
import trackc_calibration as calibration  # noqa: E402
import trackc_sim as tc  # noqa: E402
from trackc_integrated_signal_model import F0_REFERENCE_BOUNDS  # noqa: E402


SCRIPT_PATH = SRC / "218_trackc_robust_policy_timing_pilot.py"
SPEC = importlib.util.spec_from_file_location("trackc_robust_policy_timing_pilot", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
pilot = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = pilot
SPEC.loader.exec_module(pilot)


def test_reference_bounds_cover_both_calibration_world_f0_ranges():
    base_low, base_high = map(float, tc.RANGES["f0"][:2])
    expected_low = base_low * math.exp(-0.7)
    expected_high = min(base_high * math.exp(2.3), calibration.F0_CAP)
    assert F0_REFERENCE_BOUNDS == pytest.approx((expected_low, expected_high))
    for family_axes in calibration.FAMILY_AXES.values():
        shift_low, shift_high = next((lo, hi) for name, lo, hi in family_axes if name == "log_f0_shift")
        assert F0_REFERENCE_BOUNDS[0] == pytest.approx(base_low * math.exp(shift_low))
        assert F0_REFERENCE_BOUNDS[1] == pytest.approx(min(base_high * math.exp(shift_high), calibration.F0_CAP))
    assert pilot.SETTINGS.signal_f0_bounds == pytest.approx(F0_REFERENCE_BOUNDS)


def test_timing_seed_assignments_are_unique_and_stage2_only():
    values = []
    for row in pilot.SEED_FAMILIES:
        values.append(row["world_seed"])
        values.extend(row["scenario_and_signal_seeds"])
        values.extend(row["policy_seeds"])
    assert len(values) == len(set(values))
    assert all(20261800 <= seed <= 20262099 for seed in values)
    assert min(values) == 20261952 and max(values) == 20261973


def test_pilot_manifest_contract_excludes_policy_outcomes():
    required = {
        "record_type", "status", "environment", "amendment_sha256", "runner_script_sha256",
        "source_sha256", "settings", "families", "timing_seconds", "outcome_scores_saved",
        "outcome_scores_inspected", "returned_outcome_objects_discarded",
    }
    manifest_source = SCRIPT_PATH.read_text(encoding="utf-8")
    assert all(f'"{field}"' in manifest_source for field in required)
    for forbidden in (""""policy_scores""", """"utility""", """"endpoint_times""", """"course_counts""", """"policy_ranks"""):
        assert forbidden not in manifest_source


def test_completed_manifest_has_timings_and_no_outcome_fields():
    manifest_path = REPO / "output" / "trackc_robust_policy_timing_pilot_20261106_followup" / "manifest.json"
    amendment_path = REPO / "paper" / "trackc_robust_policy_timing_amendment_2026-11-06.md"
    if not manifest_path.exists() or not amendment_path.exists():
        pytest.skip("completed timing record or amendment is not present")
    record = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert record["status"] == "completed"
    assert set(record["timing_seconds"]) == {"single", "mixture"}
    assert all(0 < value < 7200 for value in record["timing_seconds"].values())
    assert record["signal_f0_bounds"] == pytest.approx(F0_REFERENCE_BOUNDS)
    assert record["candidate_model_grid_size"] == 81
    assert record["matched_d0_model_grid_size"] == 27
    assert record["outcome_scores_saved"] is False
    assert record["outcome_scores_inspected"] is False
    assert record["returned_outcome_objects_discarded"] is True
    forbidden = {"policy_scores", "utilities", "endpoint_times", "course_counts", "policy_ranks"}
    assert not forbidden.intersection(record)

    def file_hash(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    assert record["runner_script_sha256"] == file_hash(SCRIPT_PATH)
    assert record["amendment_sha256"] == file_hash(amendment_path)
