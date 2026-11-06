from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

from trackc_integrated_signal_model import (
    F0_REFERENCE_BOUNDS,
    SignalModelParameters,
    SignalObservations,
    standard_normal_marginal_signature,
)
_SPEC = importlib.util.spec_from_file_location(
    "trackc_dependence_identifiability",
    REPO / "src" / "210_trackc_dependence_identifiability.py",
)
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)
log_uniform_f0_quadrature = _MODULE.log_uniform_f0_quadrature
marginal_signal_log_density = _MODULE.marginal_signal_log_density
compare_marginal_signal_models = _MODULE.compare_marginal_signal_models


def _same_signature_parameter_pair():
    first = SignalModelParameters(q=0.8, a=0.5, b=0.4, d=0.6)
    q, rho = standard_normal_marginal_signature(first)
    a, b = 0.8, 0.5
    d = (rho - a * b) / (np.sqrt(1.0 - a * a) * np.sqrt(1.0 - b * b))
    second = SignalModelParameters(q=q, a=a, b=b, d=d)
    assert standard_normal_marginal_signature(first) == pytest.approx(
        standard_normal_marginal_signature(second)
    )
    return first, second


def test_log_uniform_f0_quadrature_weights_are_normalized():
    f0_nodes, weights = log_uniform_f0_quadrature((0.001, 0.1), n_nodes=64)
    assert f0_nodes.shape == weights.shape == (64,)
    assert np.all((f0_nodes > 0.001) & (f0_nodes < 0.1))
    assert weights.sum() == pytest.approx(1.0)
    assert np.all(weights > 0.0)


def test_same_signature_tuples_have_same_marginal_density_for_normal_resistance():
    first, second = _same_signature_parameter_pair()
    reports = np.array([False, True, False, True, True, False])
    scores = np.array([-2.0, -1.0, -0.2, 0.0, 0.8, 2.0])
    observations = SignalObservations(reports, scores)
    log_first = marginal_signal_log_density(
        observations, first, f0_bounds=F0_REFERENCE_BOUNDS, n_nodes=256
    )
    log_second = marginal_signal_log_density(
        observations, second, f0_bounds=F0_REFERENCE_BOUNDS, n_nodes=256
    )
    np.testing.assert_allclose(log_first, log_second, rtol=0.0, atol=2e-5)


@pytest.mark.parametrize(
    "bounds",
    [
        (0.001 * np.exp(-0.7), 0.1 * np.exp(-0.7)),
        (0.001, 0.1),
        (0.001 * np.exp(2.3), 0.9),
    ],
)
def test_truncated_calibration_f0_ranges_can_distinguish_same_normal_signature(bounds):
    first, second = _same_signature_parameter_pair()
    reports = np.array([False, True, False, True, True, False])
    scores = np.array([-2.0, -1.0, -0.2, 0.0, 0.8, 2.0])
    observations = SignalObservations(reports, scores)
    log_first = marginal_signal_log_density(
        observations, first, f0_bounds=bounds, n_nodes=256
    )
    log_second = marginal_signal_log_density(
        observations, second, f0_bounds=bounds, n_nodes=256
    )
    assert np.isfinite(log_first).all()
    assert np.isfinite(log_second).all()
    assert np.max(np.abs(log_first - log_second)) > 1e-4


@pytest.mark.parametrize(
    "bounds",
    [
        F0_REFERENCE_BOUNDS,
        (0.001 * np.exp(-0.7), 0.1 * np.exp(-0.7)),
        (0.001, 0.1),
        (0.001 * np.exp(2.3), 0.9),
    ],
)
def test_pairwise_kl_is_nonnegative_and_marginal_densities_normalize(bounds):
    first, second = _same_signature_parameter_pair()
    summary = compare_marginal_signal_models(
        first, second, f0_bounds=bounds, n_nodes=256, n_mri=801
    )
    assert summary["p_mass"] == pytest.approx(1.0, abs=1e-5)
    assert summary["q_mass"] == pytest.approx(1.0, abs=1e-5)
    assert summary["kl_pq"] >= -1e-10
    assert summary["kl_qp"] >= -1e-10


def test_marginal_density_rejects_mismatched_or_empty_signal_arrays():
    parameters = SignalModelParameters(q=0.8, a=0.5, b=0.4, d=0.6)
    with pytest.raises(ValueError, match="nonempty one-dimensional"):
        marginal_signal_log_density(
            SignalObservations(np.array([], dtype=bool), np.array([], dtype=float)),
            parameters,
            f0_bounds=(0.001, 0.1),
        )
    with pytest.raises(ValueError, match="nonempty one-dimensional"):
        marginal_signal_log_density(
            SignalObservations(np.array([True]), np.array([0.0, 1.0])),
            parameters,
            f0_bounds=(0.001, 0.1),
        )
