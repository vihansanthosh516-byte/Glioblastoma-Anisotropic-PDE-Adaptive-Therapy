"""Synthetic tests for Gate 1b (made-up data)."""
import numpy as np

import integration_bamd_gate1b as g


def fake(n_pat=20, per=4, signal=True, seed=0):
    rng = np.random.default_rng(seed)
    groups = np.repeat(np.arange(n_pat).astype(str), per)
    X = rng.normal(size=(n_pat * per, 8))
    y1 = X[:, 0] + rng.normal(0, 0.5, len(X)) if signal else rng.normal(size=len(X))
    Y = np.column_stack([y1, rng.normal(size=len(X))])
    strata = (rng.random(len(X)) < 0.6).astype(float)
    return X, Y, groups, strata


def test_signal_passes_noise_does_not():
    X, Y, groups, strata = fake()
    obs, p, ph = g.gate(X, Y, groups, strata, 60, np.random.default_rng(1))
    assert obs[0] > 0.6 and ph[0] < 0.05
    assert ph[1] > 0.05


def test_lopo_never_uses_held_out_patient():
    X, Y, groups, _ = fake(n_pat=5)
    Y2 = Y.copy()
    Y2[groups == "0"] += 1000  # held-out patient's own labels must not affect its prediction
    p1, p2 = g.lopo_predict(X, Y, groups), g.lopo_predict(X, Y2, groups)
    assert np.allclose(p1[groups == "0"], p2[groups == "0"])
