"""Tests for src/forecast_metrics.py and src/patient_stats.py (Phase 2.6, 2.7)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src import forecast_metrics as fm  # noqa: E402
from src import patient_stats as ps  # noqa: E402


def cube(n=20, lo=5, hi=10):
    m = np.zeros((n, n, n), bool)
    m[lo:hi, lo:hi, lo:hi] = True
    return m


def test_dice_jaccard_known_values():
    a = np.zeros((10, 10), bool); a[:5] = True
    b = np.zeros((10, 10), bool); b[2:7] = True
    assert fm.dice(a, b) == pytest.approx(0.6)
    assert fm.jaccard(a, b) == pytest.approx(30 / 70)


def test_empty_mask_rules():
    z = np.zeros((5, 5, 5), bool)
    assert fm.dice(z, z) == 1.0 and fm.dice(z, cube(5, 1, 3)) == 0.0
    assert np.isnan(fm.hd95(z, cube(5, 1, 3))) and np.isnan(fm.centroid_displacement(z, cube(5, 1, 3)))


def test_hd95_and_asd_of_shifted_cube():
    a = cube()
    b = np.roll(a, 3, axis=0)
    assert fm.hd95(a, a) == 0.0 and fm.asd(a, a) == 0.0
    assert fm.hd95(a, b, spacing=1.0) == pytest.approx(3.0, abs=0.01)
    assert fm.hd95(a, b, spacing=2.0) == pytest.approx(6.0, abs=0.02)
    assert fm.centroid_displacement(a, b, spacing=2.0) == pytest.approx(6.0)


def test_volume_and_error_metrics():
    assert fm.volume_mm3(cube(), spacing=2.0) == pytest.approx(125 * 8)
    assert fm.symmetric_pct_error(100, 100) == 0 and fm.symmetric_pct_error(150, 50) == pytest.approx(100.0)
    assert fm.symmetric_pct_error(0, 0) == 0
    assert fm.log_volume_ratio(100, 100) == 0 and np.isfinite(fm.log_volume_ratio(0, 100))


def test_patient_means_one_value_per_patient():
    df = pd.DataFrame({"patient_id": ["a", "a", "b"], "d": [0.2, 0.4, 0.1]})
    m = ps.patient_means(df, "d")
    assert m["a"] == pytest.approx(0.3) and m["b"] == pytest.approx(0.1) and len(m) == 2


def test_bootstrap_ci_covers_mean_and_shrinks_with_n():
    rng = np.random.default_rng(0)
    x = rng.normal(0.1, 0.2, 400)
    lo, hi = ps.bootstrap_ci(x, n_boot=2000)
    assert lo < x.mean() < hi
    lo2, hi2 = ps.bootstrap_ci(x[:50], n_boot=2000)
    assert (hi2 - lo2) > (hi - lo)


def test_cluster_bootstrap_resamples_patients_not_rows():
    df = pd.DataFrame({"patient_id": np.repeat(list("abcdefghij"), 5), "v": np.repeat(np.arange(10.0), 5)})
    lo, hi = ps.cluster_bootstrap(df, lambda d: d["v"].mean(), n_boot=500)
    naive = ps.bootstrap_ci(df["v"].to_numpy(), n_boot=500)
    assert (hi - lo) > 0.9 * (naive[1] - naive[0])  # rows are not independent; width must not collapse


def test_sign_flip_p_null_and_effect():
    rng = np.random.default_rng(1)
    assert ps.paired_sign_flip_p(rng.normal(0, 1, 100), n_perm=3000) > 0.05
    assert ps.paired_sign_flip_p(rng.normal(1, 1, 100), n_perm=3000) < 0.01


def test_summarize_delta_flags():
    pos = np.full(60, 0.05) + np.random.default_rng(2).normal(0, 0.01, 60)
    s = ps.summarize_delta(pos, n_boot=2000)
    assert s["broadly_positive"] and not s["mixed"] and s["pct_improved"] > 95
    # a few large wins, median ~0: mean CI may exclude 0, median does not
    d = np.zeros(60)
    d[:6] = 0.4
    d[6:30] = -0.005
    s = ps.summarize_delta(d, n_boot=2000)
    assert not s["broadly_positive"]
    assert s["median"] <= 0


def test_holm_and_bh():
    p = {"a": 0.001, "b": 0.02, "c": 0.04}
    h = ps.holm(p)
    assert h["a"] == pytest.approx(0.003) and h["b"] == pytest.approx(0.04) and h["c"] == pytest.approx(0.04)
    q = ps.benjamini_hochberg(p)
    assert q["a"] == pytest.approx(0.003) and q["c"] == pytest.approx(0.04)
    assert all(q[k] <= h[k] + 1e-12 for k in p)


def test_power_increases_with_n():
    d = np.random.default_rng(3).normal(0.03, 0.1, 200)
    r = ps.precision_and_power(d, target_half_width=0.02, n_sim=300, n_grid=(25, 100, 400))
    pw = list(r["power"].values())[0]
    assert pw[25] <= pw[100] <= pw[400]
    assert r["n_for_target_half_width"] > 0 and r["half_width_at_n"][25] > r["half_width_at_n"][400]
