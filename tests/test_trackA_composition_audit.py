"""Track A composition audit (src/trackA_composition_audit.py) on synthetic cells: a pure
cell-mix signal must be found by the cell-type-only model and NOT within a cell type; a true
within-type signal must be found within type; shuffled labels must give AUC near 0.5."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackA_composition_audit as ca  # noqa: E402


def make(within_signal, n_donors=5, n=600, n_genes=60, seed=0):
    rng = np.random.default_rng(seed)
    X, y, donor, ctype = [], [], [], []
    type_mean = {"tumor": rng.normal(0, 2, n_genes), "mg": rng.normal(0, 2, n_genes)}
    region_shift = np.zeros(n_genes)
    region_shift[:5] = within_signal
    for d in range(n_donors):
        dshift = rng.normal(0, 0.5, n_genes)
        for lab, p_tumor in ((0, 0.6), (1, 0.1)):         # periphery has fewer tumor cells
            t = np.where(rng.random(n) < p_tumor, "tumor", "mg")
            base = np.stack([type_mean[k] for k in t]) + dshift + rng.normal(0, 1, (n, n_genes))
            X.append(base + lab * region_shift)
            y.append(np.full(n, lab))
            donor.append(np.full(n, f"D{d}"))
            ctype.append(t)
    return (np.vstack(X).astype(np.float32), np.concatenate(y), np.concatenate(donor), np.concatenate(ctype))


def mean(res, key):
    return np.nanmean([r[key] for r in res.values()])


def test_pure_composition_signal():
    X, y, donor, ctype = make(within_signal=0.0)
    res = ca.lodo_audit(X, y, donor, ctype, n_genes=60, within_types=("mg", "tumor"))
    assert mean(res, "celltype_only") > 0.7
    assert abs(mean(res, "genes") - mean(res, "celltype_only")) < 0.05
    assert abs(mean(res, "within_mg") - 0.5) < 0.05
    assert abs(mean(res, "within_tumor") - 0.5) < 0.06


def test_within_type_signal_found():
    X, y, donor, ctype = make(within_signal=0.8, seed=1)
    res = ca.lodo_audit(X, y, donor, ctype, n_genes=60, within_types=("mg",))
    assert mean(res, "within_mg") > 0.8
    assert mean(res, "genes") > mean(res, "celltype_only") + 0.05


def test_shuffle_control_is_null():
    X, y, donor, ctype = make(within_signal=0.8, seed=2)
    res = ca.lodo_audit(X, y, donor, ctype, n_genes=60, within_types=("mg",), shuffle_seed=3)
    for key in ("celltype_only", "genes", "within_mg"):
        assert abs(mean(res, key) - 0.5) < 0.08


def test_sign_flip_p():
    assert ca.sign_flip_p([1, 1, 1, 1, 1, 1, 1]) == 1 / 128
    assert ca.sign_flip_p([1, -1]) > 0.4
