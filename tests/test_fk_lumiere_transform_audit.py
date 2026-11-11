"""Synthetic tests: header affine parsing, ITK matrix convention, direction and composition. No patient data."""
import gzip

import nibabel as nib
import numpy as np
import pytest

from src import fk_lumiere_transform_audit as au
from src.fk_longitudinal_prep import apply_rigid
from src.track_b_transition_estimand import geometric_coverage


def _header_bytes(affine, shape, *, sform=True):
    img = nib.Nifti1Image(np.zeros(shape, np.int16), affine)
    img.header.set_qform(affine, code=1)
    if not sform:
        img.header.set_sform(None, code=0)
    return img.header.binaryblock[:348]


def test_header_affine_matches_nibabel_for_qform_and_sform():
    rot = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1.0]])
    aff = np.eye(4)
    aff[:3, :3] = rot @ np.diag([0.5, 0.5, 6.0])
    aff[:3, 3] = (95.3, 176.8, -2.0)
    for sform in (True, False):
        g = au.parse_header_geometry(_header_bytes(aff, (8, 9, 4), sform=sform))
        assert g["affine_used"] == ("sform" if sform else "qform")
        assert np.allclose(g["affine"], aff, atol=1e-4)
        assert g["shape"] == [8, 9, 4]


def test_itk_matrix_offset_centre_convention():
    A = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1.0]])
    c, t = np.array([1.0, 2.0, 3.0]), np.array([10.0, 0.0, -5.0])
    leaf = {"type": "MatrixOffsetTransformBase_double_3_3", "parameters": list(A.ravel()) + list(t),
            "fixed_parameters": list(c)}
    m = au.itk_leaf_to_4x4(leaf)
    x = np.array([4.0, -1.0, 7.0])
    assert np.allclose((m @ np.r_[x, 1])[:3], A @ (x - c) + c + t)


def test_non_matrix_itk_transform_is_refused_not_approximated():
    for kind in ("BSplineTransform_double_3_3", "DisplacementFieldTransform_double_3_3"):
        with pytest.raises(ValueError, match="unsupported"):
            au.itk_leaf_to_4x4({"type": kind, "parameters": [0.0] * 12, "fixed_parameters": [0.0] * 3})


def test_parse_tfm_composite_with_one_leaf():
    text = ("#Insight Transform File V1.0\n#Transform 0\nTransform: CompositeTransform_double_3_3\n"
            "#Transform 1\nTransform: MatrixOffsetTransformBase_double_3_3\n"
            "Parameters: 1 0 0 0 1 0 0 0 1 1 2 3\nFixedParameters: 0 0 0\n")
    tf = au.parse_itk_tfm(text)
    assert [b["type"] for b in tf["leaves"]] == ["MatrixOffsetTransformBase_double_3_3"]
    ras = au.tfm_to_ras_4x4(tf["leaves"])
    assert np.allclose(ras[:3, 3], [-1, -2, 3])          # LPS shift (1,2,3) is RAS (-1,-2,3)


def test_rigid_voxel_map_equals_apply_rigid():
    shape, sp = (12, 14, 10), (1.0, 1.0, 1.0)
    vol = np.zeros(shape)
    vol[3:6, 4:9, 2:5] = 1
    params = [4.0, -3.0, 6.0, 1.5, -2.0, 0.5]
    ref = apply_rigid(vol, params, spacing_mm=sp, order=0)
    m = au.rigid_voxel_map(params, shape=shape, spacing_mm=sp)
    from scipy.ndimage import affine_transform
    got = affine_transform(vol, m[:3, :3], offset=m[:3, 3], order=0, mode="constant", cval=0.0)
    assert np.array_equal(ref, got)


def test_chain_direction_and_composition_on_known_geometry():
    """Native grid covers only atlas slices z in [0, 10); every factor known; chain gives exactly that."""
    atlas_aff = np.eye(4)
    native_aff = np.diag([1.0, 1.0, 2.0, 1.0])            # 5 native slices of 2 mm -> 10 atlas slices
    shift = np.eye(4)
    shift[:3, 3] = (0, 0, 0)
    out = au.chain_output_to_input(atlas_affine=atlas_aff, native_affine=native_aff, atlas_to_native_ras=shift)
    cov = geometric_coverage((16, 16, 5), (16, 16, 16), out)
    assert cov[:, :, :9].all() and not cov[:, :, 10:].any()
    # a +4 mm native-to-atlas shift moves coverage up 4 slices; using the wrong direction moves it down
    n2a = np.eye(4)
    n2a[2, 3] = 4.0
    right = geometric_coverage((16, 16, 5), (16, 16, 16), au.chain_output_to_input(
        atlas_affine=atlas_aff, native_affine=native_aff, atlas_to_native_ras=np.linalg.inv(n2a)))
    wrong = geometric_coverage((16, 16, 5), (16, 16, 16), au.chain_output_to_input(
        atlas_affine=atlas_aff, native_affine=native_aff, atlas_to_native_ras=n2a))
    assert right[:, :, 4:13].all() and not right[:, :, :4].any()
    assert not np.array_equal(right, wrong)


def test_chain_with_rigid_and_downsampling_is_one_exact_4x4():
    atlas_aff = np.eye(4)
    native_aff = np.eye(4)
    rigid = au.rigid_voxel_map([0, 0, 0, 0, 0, 3.0], shape=(16, 16, 16), spacing_mm=(1, 1, 1))
    m = au.chain_output_to_input(atlas_affine=atlas_aff, native_affine=native_aff, atlas_to_native_ras=np.eye(4),
                                 rigid=rigid, scoring_to_atlas=au.downsample_to_1mm())
    assert m.shape == (4, 4) and np.allclose(m[3], [0, 0, 0, 1])
    assert np.allclose(m @ [1, 1, 1, 1], [2.5, 2.5, 5.5, 1])  # 2*1+0.5, then +3 mm rigid shift on z
