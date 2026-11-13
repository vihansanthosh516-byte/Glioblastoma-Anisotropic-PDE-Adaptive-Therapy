"""Synthetic tests for the BAMD survival tools (made-up data; no patient data)."""
import numpy as np
import pytest

import integration_bamd_survival as s


def sim(n=400, beta=(0.8, -0.5, 0.0), seed=0):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, len(beta)))
    t = rng.exponential(1 / np.exp(X @ np.array(beta)))
    c = rng.exponential(3.0, n)
    return X, np.minimum(t, c), (t <= c).astype(int)


def test_cox_matches_lifelines():
    lifelines = pytest.importorskip("lifelines")
    import pandas as pd
    X, t, e = sim()
    b = s.cox_fit(X, t, e, ridge=0.0)
    df = pd.DataFrame(X, columns=["a", "b", "c"]).assign(t=t, e=e)
    ref = lifelines.CoxPHFitter().fit(df, "t", "e").params_.to_numpy()
    assert np.allclose(b, ref, atol=1e-4)


def test_cox_recovers_truth():
    X, t, e = sim(n=3000)
    b = s.cox_fit(X, t, e)
    assert np.allclose(b, [0.8, -0.5, 0.0], atol=0.08)


def test_lasso_zeroes_null_feature():
    X, t, e = sim(n=2000)
    b = s.lasso_cox(X, t, e, lam=0.05 * 2000)
    assert abs(b[2]) < 0.05 and b[0] > 0.3


def test_cindex_extremes():
    t = np.array([1.0, 2, 3, 4])
    e = np.ones(4)
    assert s.cindex(-t, t, e) == 1.0
    assert s.cindex(t, t, e) == 0.0
    assert s.cindex(np.zeros(4), t, e) == 0.5
