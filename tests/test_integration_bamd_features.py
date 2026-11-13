"""Synthetic tests for BAMD feature extraction (made-up volumes; no patient data)."""
import numpy as np
import pytest

import integration_bamd_features as f

nib = pytest.importorskip("nibabel")


def make_case(tmp_path, case="BraTS2021_99999"):
    d = tmp_path / case
    d.mkdir()
    seg = np.zeros((40, 40, 40), np.int16)
    seg[10:30, 10:30, 10:30] = 2          # edema cube 20^3
    seg[15:25, 15:25, 15:25] = 4          # enhancing 10^3
    seg[18:22, 18:22, 18:22] = 1          # necrosis 4^3
    rng = np.random.default_rng(0)
    for im in f.IMAGES:
        x = np.zeros((40, 40, 40), np.float32)
        x[5:35, 5:35, 5:35] = rng.normal(100, 10, (30, 30, 30))
        x[seg == 4] += 50 if im == "t1ce" else 0
        nib.save(nib.Nifti1Image(x, np.eye(4)), str(d / f"{case}_{im}.nii.gz"))
    nib.save(nib.Nifti1Image(seg, np.eye(4)), str(d / f"{case}_seg.nii.gz"))
    return case


def test_feature_list_has_59_names():
    assert len(f.feature_names()) == 59 == len(set(f.feature_names()))


def test_case_features_values(tmp_path):
    case = make_case(tmp_path)
    out = f.case_features(case, root=tmp_path)
    assert set(out) == {"brats_id", *f.feature_names()}
    et = 10**3 - 4**3
    assert out["logvol_ET"] == pytest.approx(np.log1p(et * 0.001))
    assert out["logvol_WT"] == pytest.approx(np.log1p(20**3 * 0.001))
    assert out["ratio_ncr_core"] == pytest.approx(64 / 1000)
    assert out["mean_ET_t1ce"] > out["mean_ED_t1ce"] + 2  # enhancement shows up after z-scoring
    assert 0 < out["sph_wt"] <= 1.0


def test_cube_sphericity_matches_formula():
    m = np.zeros((12, 12, 12), bool)
    m[2:10, 2:10, 2:10] = True
    sv, sph = f.shape(m)
    assert sv == pytest.approx(6 * 64 / 512)
    assert sph == pytest.approx((np.pi / 6) ** (1 / 3), rel=1e-6)  # cube sphericity
