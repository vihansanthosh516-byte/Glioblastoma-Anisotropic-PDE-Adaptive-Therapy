import numpy as np
import pytest

from trackc_t_bonferroni_coverage_audit import bonferroni_t_intervals
from trackc_t_bonferroni_worldcount_extension import (
    DISTRIBUTION_PAIRS,
    WORLD_COUNTS,
    draw_standardized,
    run_coverage_screen,
)


def test_extension_grid_is_frozen_and_includes_bounded_skew():
    assert len(DISTRIBUTION_PAIRS) == 11
    assert WORLD_COUNTS == (30, 50, 100, 200)


def test_new_bounded_distributions_are_centered_and_standardized():
    rng = np.random.Generator(np.random.PCG64DXSM(7301))
    for name in ("beta_2_5", "beta_5_2", "bernoulli_0_1"):
        sample = draw_standardized(name, 200_000, rng)
        assert sample.shape == (200_000,)
        assert np.isfinite(sample).all()
        assert abs(sample.mean()) < 0.02
        assert abs(sample.var() - 1.0) < 0.04


def test_interval_candidate_is_same_fixed_bonferroni_t_method():
    values = np.array([-2.0, -1.0, 0.5, 2.5, 3.0])
    result = bonferroni_t_intervals({
        "single-response": values,
        "mixture-response": values[::-1],
    })
    assert result["method"] == "two-sided Student-t intervals with Bonferroni adjustment"
    assert result["family_count"] == 2
    assert result["intervals"]["single-response"]["n_worlds"] == 5


def test_full_extension_screen_is_complete_and_uses_no_project_outcomes():
    result = run_coverage_screen()
    assert result["cells_total"] == 44
    assert result["monte_carlo_repetitions_per_cell"] == 20_000
    assert result["world_counts_per_family"] == [30, 50, 100, 200]
    assert result["policy_outcomes_used"] is False
    assert result["patient_data_used"] is False
    assert result["inference_failures_total"] == 0
    assert len(result["cells"]) == 44
    assert all(cell["world_count_per_family"] in WORLD_COUNTS for cell in result["cells"])
