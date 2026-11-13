"""Synthetic tests for BAMD Gate 1 (made-up numbers; no patient data)."""
import numpy as np
import pandas as pd

import integration_bamd_gate1 as g


def test_pick_samples_prefers_primary():
    cols = ["TCGA-AA-0001-02", "TCGA-AA-0001-01", "TCGA-AA-0002-11"]
    assert g.pick_samples(cols) == {"TCGA-AA-0001": "TCGA-AA-0001-01", "TCGA-AA-0002": "TCGA-AA-0002-11"}


def test_holm_matches_hand_values():
    assert np.allclose(g.holm([0.01, 0.04, 0.03, 0.5]), [0.04, 0.09, 0.09, 0.5])


def test_cv_rho_detects_signal_and_not_noise():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(90, 59))
    Y = np.column_stack([X[:, :3].sum(1) + rng.normal(0, 0.5, 90), rng.normal(size=90)])
    r = g.cv_rho(X, Y, np.random.default_rng(1), repeats=2)
    assert r[0] > 0.6 and abs(r[1]) < 0.3


def test_state_scores_rise_with_module_genes():
    genes = [f"G{i}" for i in range(200)]
    rng = np.random.default_rng(2)
    expr = pd.DataFrame(rng.normal(size=(200, 30)), index=genes)
    expr.iloc[:20, :15] += 3  # module A high in first 15 patients
    mods = {"A": genes[:20], "B": genes[20:40]}
    s = g.state_scores(expr, mods)
    assert s["A"].iloc[:15].mean() > s["A"].iloc[15:].mean() + 1


def test_freeze_residual_is_orthogonal_to_predicted_states():
    rng = np.random.default_rng(3)
    X = rng.normal(size=(80, 10))
    Y = X[:, :2] @ np.array([[1.0, 0.0], [0.5, 1.0]])
    fr = g.freeze(X, Y, [0, 1], np.random.default_rng(4))
    Xs = (X - np.array(fr["mu"])) / np.array(fr["sd"])
    P = Xs @ np.array(fr["W"])
    R = Xs - P @ np.array(fr["B"])
    assert np.allclose(P.T @ R, 0, atol=1e-8)
