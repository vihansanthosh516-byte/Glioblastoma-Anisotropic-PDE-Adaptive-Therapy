"""Synthetic tests for the PROPOSED PDE-only cutoff option of the hybrid forecast. No patient data."""
import numpy as np
import pytest

from src.fk_hybrid_forecast import hybrid_posterior_predictive_forecast

SHAPE, H = (16, 16, 16), 2.0


def _run(**kw):
    g = np.ogrid[tuple(slice(0, s) for s in SHAPE)]
    r2 = sum(((g[i] - 8) * H) ** 2 for i in range(3))
    return hybrid_posterior_predictive_forecast(
        [r2 <= 6 ** 2, r2 <= 7 ** 2], [0.0, 42.0], 133.0, prior={"distribution": "log_uniform"}, h=H,
        boundary_covariance={"sigma_b_mm": 0.8, "tau_b_mm": 0.3, "ell_mm": 6.0},
        population_log_volume_slope_per_day=0.002, shrinkage_time_constant_days=90.0,
        volume_variance_fit={"intercept": 0.02, "slope_per_day": 0.0001}, domain=r2 <= 15 ** 2,
        n_coarse=5, n_refine=3, n_draws=6, n_anchors=16, seed=3, **kw)


def test_default_output_has_no_pde_only_keys_and_candidate_is_unchanged():
    base = _run()
    with_cut = _run(pde_only_concentration_cutoff=0.5)
    assert not any(k.startswith("pde_only") for k in base)
    assert np.array_equal(base["mask_draws"], with_cut["mask_draws"])
    assert np.array_equal(base["pre_boundary_probability"], with_cut["pre_boundary_probability"])


def test_pde_only_fields_are_consistent():
    res = _run(pde_only_concentration_cutoff=0.5)
    p = res["pde_only_inclusion_probability"]
    assert p.min() >= 0 and p.max() <= 1
    assert np.array_equal(res["pde_only_mask"], p >= 0.5)
    assert 0 <= res["pde_only_empty_draw_count"] <= 6


@pytest.mark.parametrize("bad", [0.0, 1.0, -0.1])
def test_cutoff_must_be_inside_unit_interval(bad):
    with pytest.raises(ValueError, match="cutoff"):
        _run(pde_only_concentration_cutoff=bad)


def test_cutoff_grid_leaves_candidate_unchanged_in_both_modes():
    grid = tuple(round(0.05 * i, 2) for i in range(1, 11))
    base, with_grid = _run(), _run(pde_only_cutoff_grid=grid)
    assert np.array_equal(base["mask_draws"], with_grid["mask_draws"])
    assert set(with_grid["pde_only_inclusion_probability_by_cutoff"]) == set(grid)
    cal = _run(calibration_only=True, pde_only_cutoff_grid=grid)
    cal0 = _run(calibration_only=True)
    assert np.array_equal(cal["pre_boundary_probability"], cal0["pre_boundary_probability"])
    p05 = with_grid["pde_only_inclusion_probability_by_cutoff"][0.5]
    single = _run(pde_only_concentration_cutoff=0.5)["pde_only_inclusion_probability"]
    assert np.array_equal(p05, single)
