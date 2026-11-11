"""Synthetic mechanics checks for script 222 helpers. No expression data."""
import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("myeloid_origin_test", ROOT / "src" / "222_trackA_myeloid_origin.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_symbol_index_handles_common_id_formats():
    idx = mod.symbol_index(["ENSG1|P2RY12", "CX3CR1", "ENSG3_LYZ"])
    assert idx["P2RY12"] == 0 and idx["CX3CR1"] == 1 and idx["LYZ"] == 2
    found, cols = mod.set_columns(("P2RY12", "TMEM119"), idx)
    assert found == ["P2RY12"] and cols == [0]


def test_zscore_zero_variance_columns_are_zero():
    x = np.array([[1.0, 5.0], [3.0, 5.0]])
    z = mod.zscore(x)
    assert np.allclose(z[:, 0], [-1, 1]) and np.allclose(z[:, 1], 0)


def test_program_direction_on_planted_signal():
    rng = np.random.default_rng(0)
    y = np.r_[np.ones(30), np.zeros(30)].astype(int)
    mg = rng.normal(size=(60, 4)) + (y == 0)[:, None] * 1.5   # higher in periphery
    bm = rng.normal(size=(60, 4)) + (y == 1)[:, None] * 1.5   # higher in core
    score = mod.zscore(mg).mean(1) - mod.zscore(bm).mean(1)
    assert mod.auc(y, -score) > 0.9
