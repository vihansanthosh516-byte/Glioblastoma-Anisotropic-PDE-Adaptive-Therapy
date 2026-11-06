from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "trackc_finite_sample_signal_recovery",
    REPO / "src" / "213_trackc_finite_sample_signal_recovery.py",
)
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)


def test_recovery_grid_is_the_fixed_81_point_grid_and_truths_are_on_it():
    candidates = _MODULE.parameter_grid()
    assert len(candidates) == 81
    assert len(set(candidates)) == 81
    assert all(truth in candidates for truth in _MODULE.TRUTH_ANCHORS)


def test_q_half_likelihood_is_invariant_to_a_and_d():
    observations = _MODULE.SignalObservations(
        molecular_report=np.asarray([True, False, True, True, False]),
        mri_score=np.asarray([-1.1, 0.2, 0.8, 1.4, -0.3]),
    )
    bounds = _MODULE.F0_RANGES["default"]
    likelihoods = _MODULE.marginal_log_likelihoods(observations, bounds, n_nodes=32)
    candidates = _MODULE.parameter_grid()
    same_q_and_b = [
        likelihoods[i]
        for i, point in enumerate(candidates)
        if point[0] == 0.5 and point[2] == 0.5
    ]
    assert len(same_q_and_b) == 9
    assert np.ptp(same_q_and_b) < 1e-12


def test_keyed_synthetic_observations_are_reproducible():
    truth = _MODULE.TRUTH_ANCHORS[0]
    bounds = _MODULE.F0_RANGES["default"]
    first = _MODULE._seeded_dataset(777001, 1, 0, 0, truth, bounds, 18)
    second = _MODULE._seeded_dataset(777001, 1, 0, 0, truth, bounds, 18)
    assert np.array_equal(first.molecular_report, second.molecular_report)
    assert np.array_equal(first.mri_score, second.mri_score)


def test_small_smoke_run_returns_every_fixed_condition_and_finite_scores():
    result = _MODULE.build_diagnostic(
        n_replicates=1,
        sample_sizes=(8,),
        n_nodes=16,
        seed_roots=(777002,),
        check_registered_seeds=False,
    )
    assert result["protocol"]["candidate_count"] == 81
    assert len(result["all_replicates"]) == 3 * len(_MODULE.TRUTH_ANCHORS)
    assert len(result["summary"]) == 3 * len(_MODULE.TRUTH_ANCHORS)
    assert all(np.isfinite(row["d_absolute_error"]) for row in result["all_replicates"])
    assert all(row["seed_root"] == 777002 for row in result["all_replicates"])


def test_registered_run_refuses_a_seed_outside_the_reserved_development_block():
    with pytest.raises(RuntimeError, match="registered development seed check failed"):
        _MODULE.build_diagnostic(
            n_replicates=1,
            sample_sizes=(8,),
            n_nodes=16,
            seed_roots=(777003,),
            check_registered_seeds=True,
        )
