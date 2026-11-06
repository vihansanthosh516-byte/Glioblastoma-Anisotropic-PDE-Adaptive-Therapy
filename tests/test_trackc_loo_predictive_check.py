import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest


_SOURCE = Path(__file__).parents[1] / "src" / "215_trackc_loo_predictive_check.py"
_SPEC = importlib.util.spec_from_file_location("trackc_loo_predictive_check", _SOURCE)
assert _SPEC is not None and _SPEC.loader is not None
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)


def test_protocol_constants_match_registered_synthetic_design():
    assert _MODULE.ROOT_SEEDS == tuple(range(20261905, 20261930))
    assert _MODULE.SAMPLE_SIZES == (18, 28, 100, 500)
    assert len(_MODULE.MODEL_GRID) == 81
    assert len(_MODULE.PRIMARY_TRUTHS) == 3


def test_smoke_run_reports_all_conditions_but_cannot_pass_protocol_gate():
    result = _MODULE.build_diagnostic(
        seed_roots=(12345, 12346),
        sample_sizes=(18, 28),
        n_f0_nodes=8,
        check_registered_seeds=False,
    )
    assert result["protocol_run"] is False
    assert result["primary_gate"]["status"] == "not_protocol_run"
    assert len(result["all_cohort_scores"]) == 3 * 7 * 2 * 2
    assert len(result["group_summaries"]) == 3 * 7 * 2
    assert all(
        np.isfinite(list(row["scores"].values())).all()
        for row in result["all_cohort_scores"]
    )
    assert set(result["primary_gate"]["comparisons"]) == {
        "n=18|vs_uniform_grid_mixture",
        "n=18|vs_matched_conditional_independence",
        "n=28|vs_uniform_grid_mixture",
        "n=28|vs_matched_conditional_independence",
    }


def test_spent_development_seeds_are_refused_before_diagnostic_work():
    with pytest.raises(RuntimeError, match="registered development seed check failed"):
        _MODULE.build_diagnostic(
            seed_roots=(20261880,),
            sample_sizes=(18,),
            n_f0_nodes=8,
            check_registered_seeds=True,
        )
