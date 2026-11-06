from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

from trackc_integrated_signal_model import F0_REFERENCE_BOUNDS

_SPEC = importlib.util.spec_from_file_location(
    "trackc_signal_fisher_identifiability",
    REPO / "src" / "211_trackc_signal_fisher_identifiability.py",
)
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)
expected_fisher_information = _MODULE.expected_fisher_information


def test_standard_normal_reference_control_has_two_identifiable_directions():
    result = expected_fisher_information(
        F0_REFERENCE_BOUNDS, n_f0_nodes=64, n_mri=401
    )
    assert result["density_mass"] == pytest.approx(1.0, abs=1e-5)
    assert result["effective_rank"] == 2
    assert result["d_local_cramer_rao_sd_bound_per_sqrt_n"] is None
    assert np.all(np.asarray(result["eigenvalues"]) >= -1e-10)


def test_truncated_calibration_prior_has_full_local_rank_at_frozen_point():
    result = expected_fisher_information(
        (0.001, 0.1), n_f0_nodes=64, n_mri=401
    )
    matrix = np.asarray(result["fisher_information_per_pair"])
    assert result["density_mass"] == pytest.approx(1.0, abs=1e-5)
    assert result["effective_rank"] == 4
    assert result["condition_number"] > 1.0
    assert result["d_local_cramer_rao_sd_bound_per_sqrt_n"] > 0.0
    np.testing.assert_allclose(matrix, matrix.T, rtol=0.0, atol=1e-12)
    assert np.linalg.eigvalsh(matrix).min() >= -1e-10


def test_finite_difference_refuses_to_cross_parameter_boundary():
    with pytest.raises(ValueError, match="crosses q boundary"):
        expected_fisher_information(
            (0.001, 0.1), parameters=(0.00005, 0.5, 0.4, 0.6),
            finite_difference_step=1e-4, n_f0_nodes=16, n_mri=101,
        )
