from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "trackc_signal_fisher_grid", REPO / "src" / "212_trackc_signal_fisher_grid.py"
)
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)


def test_parameter_grid_contains_all_81_frozen_points_once():
    grid = _MODULE.parameter_grid()
    assert len(grid) == 81
    assert len(set(grid)) == 81
    assert all(len(point) == 4 for point in grid)


def test_grid_confirms_q_half_d_negative_control_and_normal_alias():
    result = _MODULE.build_diagnostic(n_f0_nodes=32, n_mri=201)
    assert len(result["all_points"]) == 324
    for rows in result["summary_by_f0_range"].values():
        assert rows["n_points"] == 81
        assert rows["q_0_5_all_d_unidentified"] is True


def test_validated_resolution_preserves_normal_control_rank_two():
    summary = _MODULE.expected_fisher_information(
        _MODULE.F0_RANGES["reference_normal_control"],
        parameters=(0.7, 0.5, 0.5, 0.5),
        n_f0_nodes=256,
        n_mri=1601,
    )
    assert summary["density_mass"] == pytest.approx(1.0, abs=1e-5)
    assert summary["effective_rank"] <= 2


def test_grid_rejects_no_parameter_selection_by_preserving_all_points():
    result = _MODULE.build_diagnostic(n_f0_nodes=16, n_mri=101)
    assert len(result["all_points"]) == 4 * 3**4
    assert all(row["density_mass"] == pytest.approx(1.0, abs=1e-5) for row in result["all_points"])
