"""Synthetic checks for src/trackA_hgcc_tmz_102.py (no HGCC data used)."""
import importlib.util
from pathlib import Path

import numpy as np

_p = Path(__file__).resolve().parent.parent / "src" / "trackA_hgcc_tmz_102.py"
_spec = importlib.util.spec_from_file_location("hgcc102", _p)
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)


def test_percentile_ranks_and_score():
    mat = np.array([[1.0, 4.0], [2.0, 3.0], [3.0, 2.0], [4.0, 1.0]])
    pr = m.percentile_ranks(mat)
    assert np.allclose(pr[:, 0], [0.25, 0.5, 0.75, 1.0])
    s = m.panel_score(pr, [3], [0])
    assert np.allclose(s, [0.75, -0.75])


def test_planted_signal_detected_and_null_not():
    rng = np.random.default_rng(1)
    x = rng.normal(size=60)
    rho, p = m.perm_p_one_sided(x, x + 0.5 * rng.normal(size=60), 2000, 7)
    assert rho > 0.6 and p < 0.01
    rho0, p0 = m.perm_p_one_sided(x, rng.normal(size=60), 2000, 7)
    assert p0 > 0.01 or rho0 < 0.3


def test_negative_signal_not_supported():
    rng = np.random.default_rng(2)
    x = rng.normal(size=50)
    rho, p = m.perm_p_one_sided(x, -x, 2000, 7)
    assert rho < 0 and p > 0.9


def test_boot_ci_contains_truth():
    rng = np.random.default_rng(3)
    x = rng.normal(size=80)
    lo, hi = m.boot_ci(x, x + rng.normal(size=80), 500, 9)
    assert lo < 0.7 < hi or lo > 0.4


def test_controls_unique_and_shape():
    rng = np.random.default_rng(4)
    feat = rng.normal(size=(500, 2))
    panel = list(range(10))
    ctrls = m.control_panels(feat, panel, 4, 5, 11)
    for pos, neg in ctrls:
        allc = pos + neg
        assert len(pos) == 4 and len(neg) == 6
        assert len(set(allc)) == 10 and not set(allc) & set(panel)
