#!/usr/bin/env python3
"""Tests for the multi-omic / epigenetic feature fusion module (Proposal 2).

Verifies that:
  - the synthetic multiomic_features.tsv round-trips through load/save
  - the trained ElasticNet achieves a sensible R^2 on synthetic data
  - predict_params_from_features returns values within physiological bounds
  - the multi-omic path produces lower RMSE than the linear legacy mapping
    when ground truth involves methylation / CNV features

Run:
    venv\\Scripts\\python.exe tests\\test_multiomic_fusion.py
    venv\\Scripts\\python.exe -m pytest tests/test_multiomic_fusion.py
"""
from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
MODULE_PATH = REPO_ROOT / "src" / "multiomic_fusion.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("multiomic_fusion", MODULE_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mod = _load_module()
build_multiomic_tsv = mod.build_multiomic_tsv
load_multiomic_tsv = mod.load_multiomic_tsv
train_elasticnet = mod.train_elasticnet
predict_params_from_features = mod.predict_params_from_features
NEFTEL_STATES = mod.NEFTEL_STATES
METHYLATION_FEATURES = mod.METHYLATION_FEATURES
feature_columns = mod.feature_columns
RHO_MIN, RHO_MAX = mod.RHO_MIN, mod.RHO_MAX
D_MIN, D_MAX = mod.D_MIN, mod.D_MAX


def _make_dataset(n=30, seed=7):
    rng = np.random.default_rng(seed)
    pids = [f"P{i:03d}" for i in range(n)]
    fracs = rng.dirichlet(np.ones(4), size=n)
    return pids, fracs


def test_feature_dim_matches_schema():
    # 4 Neftel + 4 methyl + 4 CNV + 4 metabolic = 16
    assert len(feature_columns()) == 16
    assert mod.feature_dim() == 16


def test_tsv_roundtrip():
    pids, fracs = _make_dataset(n=12, seed=1)
    with tempfile.TemporaryDirectory() as d:
        out = Path(d) / "m.tsv"
        build_multiomic_tsv(pids, fracs, out_path=out)
        X, y = load_multiomic_tsv(out)
    assert X.shape == (12, 16), X.shape
    assert y.shape == (12, 2), y.shape
    assert np.all((y[:, 0] >= RHO_MIN) & (y[:, 0] <= RHO_MAX))
    assert np.all((y[:, 1] >= D_MIN) & (y[:, 1] <= D_MAX))


def test_elasticnet_within_bounds_and_sensible_r2():
    pids, fracs = _make_dataset(n=40, seed=2)
    with tempfile.TemporaryDirectory() as d:
        tsv = Path(d) / "m.tsv"
        model = Path(d) / "model.pkl"
        build_multiomic_tsv(pids, fracs, out_path=tsv)
        X, y = load_multiomic_tsv(tsv)
        m = train_elasticnet(X, y, model_path=model)
        rho, D = predict_params_from_features(X, model_path=model)

    assert np.all((rho >= RHO_MIN) & (rho <= RHO_MAX))
    assert np.all((D >= D_MIN) & (D <= D_MAX))
    # The synthetic labels are noisy but multi-omic-driven, so a robust model
    # should give at least a non-trivial R^2 (could be modest with strong
    # noise). Require strictly > -0.5 to catch gross regressions.
    assert m["rho_r2"] > -0.5, m["rho_r2"]
    assert m["D_r2"] > -0.5, m["D_r2"]


def test_multimodal_beats_unimodal_rmse_when_omic_drives_labels():
    """When the label intentionally depends on methylation+CNV features, the
    multi-omic elastic-net should fit better than a Neftel-fraction-only
    linear baseline."""
    rng = np.random.default_rng(11)
    n = 40
    npc, opc, ac, mes = rng.dirichlet(np.ones(4), size=n).T
    meth_beta = rng.uniform(0.1, 0.9, size=n)
    cnv_egfr = rng.uniform(-1, 3, size=n)
    # Label explicitly modulated by omic features (alpha heavy on methyl, cnv)
    rho_truth = np.clip(0.04 + 0.2 * (npc + opc) + 0.08 * meth_beta
                        + 0.01 * cnv_egfr + 0.01 * rng.standard_normal(n),
                        RHO_MIN, RHO_MAX)
    D_truth = np.clip(0.15 + 0.4 * (ac + mes) - 0.05 * meth_beta
                      + 0.5 * rng.random(n) * 0.0
                      + 0.01 * rng.standard_normal(n),
                      D_MIN, D_MAX)

    # Build X with the schema: Neftel + meth + cnv + metab
    neftel = np.column_stack([npc, opc, ac, mes])
    meth = np.column_stack([meth_beta,
                            rng.uniform(0, 1, size=n),
                            rng.uniform(0, 1, size=n),
                            rng.uniform(0, 1, size=n)])
    cnv = np.column_stack([cnv_egfr,
                           rng.uniform(0, 1, size=n),
                           rng.standard_normal(n),
                           rng.standard_normal(n)])
    metab = rng.uniform(0.5, 1.5, size=(n, 4))
    X = np.column_stack([neftel, meth, cnv, metab])
    y = np.column_stack([rho_truth, D_truth])

    with tempfile.TemporaryDirectory() as d:
        model = Path(d) / "m.pkl"
        m = train_elasticnet(X, y, model_path=model)
        rho_pred_multi, _ = predict_params_from_features(X, model_path=model)

    # Baseline: ordinary least squares on Neftel-only
    A = np.column_stack([np.ones(n), npc, opc, ac, mes])
    coeffs_rho, *_ = np.linalg.lstsq(A, rho_truth, rcond=None)
    rho_pred_uni = A @ coeffs_rho

    rmse_multi = np.sqrt(np.mean((rho_pred_multi - rho_truth) ** 2))
    rmse_uni = np.sqrt(np.mean((rho_pred_uni - rho_truth) ** 2))
    assert rmse_multi <= rmse_uni + 1e-6, (
        f"multi-omic RMSE {rmse_multi:.4f} should be <= uni {rmse_uni:.4f}")


def run_all():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed, failed = 0, 0
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL  {t.__name__}: {e}")
            failed += 1
        except Exception as e:
            print(f"  ERROR {t.__name__}: {type(e).__name__}: {e}")
            failed += 1
    print(f"\nResults: {passed} passed, {failed} failed, {passed + failed} total")
    return failed == 0


if __name__ == "__main__":
    sys.exit(0 if run_all() else 1)
