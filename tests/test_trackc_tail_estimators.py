"""Tail estimators (src/trackc_tail_estimators.py) on synthetic data with a known answer.
No simulator is run; this checks the statistics, not any treatment policy."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_sim as tc  # noqa: E402
import trackc_tail_estimators as te  # noqa: E402

ALPHA = 0.10


def gaussian_world(n=4000, r_max=16, mu=500.0, sd_theta=100.0, sd_noise=60.0, seed=1):
    rng = np.random.default_rng(seed)
    theta = rng.normal(mu, sd_theta, n)
    X = theta[None, :] + rng.normal(0.0, sd_noise, (r_max, n))
    return theta, X


def test_weights_reproduce_tc_cvar_lower():
    rng = np.random.default_rng(3)
    for n in (7, 10, 33, 150, 400):
        x = rng.normal(size=n)
        w = te._weights(n, ALPHA)
        assert w.sum() == pytest.approx(ALPHA * n)
        assert (w * np.sort(x)).sum() / (ALPHA * n) == pytest.approx(tc.cvar_lower(x, ALPHA))


def test_plugin_is_biased_low_and_rises_with_repeats():
    theta, X = gaussian_world()
    truth = tc.cvar_lower(theta, ALPHA)
    vals = [te.plugin(X[:R]) for R in (1, 2, 4, 8, 16)]
    assert vals == sorted(vals)                         # rises with R
    assert vals[0] < truth - 5.0                        # clearly below at R = 1
    assert vals[-1] < truth                             # still below at R = 16
    assert truth - vals[-1] < 0.35 * (truth - vals[0])  # bias shrinks a lot by R = 16


def test_split_sample_is_biased_high_and_brackets_the_truth():
    theta, X = gaussian_world()
    truth = tc.cvar_lower(theta, ALPHA)
    rng = np.random.default_rng(5)
    s2 = te.split_sample(X[:2], rng=rng, n_splits=50)
    s16 = te.split_sample(X, rng=rng, n_splits=50)
    assert s2 > truth and s16 > truth                   # biased high
    assert s16 < s2                                     # shrinks with more repeats
    assert te.plugin(X[:16]) < truth < s16              # the two estimators bracket the truth
    assert np.isnan(te.split_sample(X[:1]))


def test_inverse_R_fit_extrapolates_close_to_the_truth():
    theta, X = gaussian_world()
    truth = tc.cvar_lower(theta, ALPHA)
    Rs = [2, 4, 8, 16]
    vals = [te.plugin(X[:R]) for R in Rs]
    c, a, rms = te.fit_inverse_R(Rs, vals)
    assert a > 0 and rms < 1.0
    plug1 = te.plugin(X[:1])
    assert abs(c - truth) < abs(plug1 - truth) / 3.0    # much closer than the R = 1 plug-in
    assert abs(c - truth) < 3.0


def test_inverse_R_fit_recovers_exact_curve():
    Rs = np.array([1, 2, 4, 8, 16])
    c, a, rms = te.fit_inverse_R(Rs, 177.9 - 6.8 / Rs)
    assert (c, a) == pytest.approx((177.9, 6.8), abs=1e-9) and rms < 1e-9


def test_joint_matches_gaussian_theory():
    theta, X = gaussian_world(n=6000, r_max=16)
    sd = float(np.hypot(100.0, 60.0))
    z = 1.2815515655446004                                   # normal quantile at 0.90
    phi = np.exp(-0.5 * z * z) / np.sqrt(2.0 * np.pi)
    theory = 500.0 - sd * phi / ALPHA
    assert te.joint(X) == pytest.approx(theory, abs=2.5)
    # joint is not selection-biased: flat in R up to sampling noise
    assert abs(te.joint(X[:2]) - te.joint(X)) < 4.0


def test_open_loop_arm_has_no_bias():
    theta = np.random.default_rng(9).normal(500, 100, 500)
    X = np.tile(theta, (8, 1))                          # repeats identical, like an open-loop rule
    base = tc.cvar_lower(theta, ALPHA)
    for R in (1, 2, 4, 8):
        assert te.plugin(X[:R]) == pytest.approx(base, abs=1e-9)
    assert te.joint(X) == pytest.approx(base, abs=1e-9)


def test_paired_regret_known_case():
    base = np.full((4, 100), 300.0)
    pol = base.copy()
    pol[:, :10] -= 10.0                                  # 10% of scenarios lose 10 days in every repeat
    out = te.paired_regret(pol, base)
    assert out["scenario_cvar"] == pytest.approx(-10.0) and out["joint_cvar"] == pytest.approx(-10.0)
    assert out["share_worse"] == pytest.approx(0.10) and out["mean_diff"] == pytest.approx(-1.0)
    with pytest.raises(ValueError):
        te.paired_regret(pol, base[:2])


def test_curve_uses_one_subset_at_full_size_and_is_reproducible():
    theta, X = gaussian_world(n=300, r_max=8)
    c1 = te.curve(X, [1, 2, 8], te.plugin, rng=np.random.default_rng(1), n_subsets=20)
    c2 = te.curve(X, [1, 2, 8], te.plugin, rng=np.random.default_rng(1), n_subsets=20)
    assert c1 == c2 and c1[8][1] == 0.0 and c1[1][1] > 0.0
    with pytest.raises(ValueError):
        te.curve(X, [9], te.plugin)


def test_scenario_bootstrap_interval_contains_point_and_zero_for_equal_arms():
    theta, X = gaussian_world(n=300, r_max=4)
    point, (lo, hi) = te.scenario_bootstrap(X, X, te.joint, n_boot=50, seed=1)
    assert point == 0.0 and lo == 0.0 and hi == 0.0
    Y = X - 20.0
    point, (lo, hi) = te.scenario_bootstrap(X, Y, te.joint, n_boot=100, seed=1)
    assert point == pytest.approx(20.0) and lo <= point <= hi
