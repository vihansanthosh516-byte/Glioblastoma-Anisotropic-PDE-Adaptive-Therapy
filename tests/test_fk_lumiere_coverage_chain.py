"""Synthetic tests for the LUMIERE coverage chain. Known transforms, landmarks and physical-space resampling."""
import numpy as np
from scipy.ndimage import affine_transform

from src import fk_lumiere_coverage_chain as ch
from src.fk_longitudinal_prep import apply_rigid, moving_points_to_fixed
from src.track_b_coverage_scoring import common_support

ATLAS_SHAPE = (40, 44, 36)
ATLAS_AFF = np.array([[-1.0, 0, 0, 20], [0, 1, 0, -22], [0, 0, 1, -18], [0, 0, 0, 1]])   # LUMIERE-like flip


def _rot(deg, axis):
    a = np.deg2rad(deg)
    c, s = np.cos(a), np.sin(a)
    r = {2: [[c, -s, 0], [s, c, 0], [0, 0, 1]], 1: [[c, 0, s], [0, 1, 0], [-s, 0, c]]}[axis]
    m = np.eye(4)
    m[:3, :3] = r
    return m


def _native(spacing=(0.5, 0.5, 3.0), shape=(70, 80, 10), origin=(-15.0, -18.0, -12.0), rot=None):
    n = np.eye(4)
    n[:3, :3] = (np.eye(3) if rot is None else rot[:3, :3]) @ np.diag(spacing)
    n[:3, 3] = origin
    return n, shape


def _n2a(shift=(3.0, -2.0, 1.5), deg=7.0):
    t = _rot(deg, 2)
    t[:3, 3] = shift
    return t


def test_landmark_without_rigid_matches_physical_reasoning():
    N, nshape = _native(rot=_rot(-12, 1))
    T = _n2a()
    m = ch.visit_output_to_input(native_affine=N, tfm_native_to_atlas_ras=T, atlas_affine=ATLAS_AFF,
                                 atlas_shape=ATLAS_SHAPE)
    i = np.array([5, 7, 4, 1.0])
    atlas_mm = ATLAS_AFF @ np.r_[2 * i[:3], 1]                       # scoring voxel i -> atlas voxel 2i -> mm
    native_mm = np.linalg.solve(T, atlas_mm)                         # T maps native -> atlas, so invert
    native_vox = np.linalg.solve(N, native_mm)
    assert np.allclose(m @ i, native_vox)
    # and the forward direction sends that native voxel back to the same atlas point
    assert np.allclose(T @ N @ native_vox, atlas_mm)


def test_landmark_with_rigid_matches_apply_rigid_landmark_map():
    N, _ = _native()
    T = _n2a()
    params = [3.0, -4.0, 6.0, 2.0, -1.5, 3.0]
    m = ch.visit_output_to_input(native_affine=N, tfm_native_to_atlas_ras=T, atlas_affine=ATLAS_AFF,
                                 atlas_shape=ATLAS_SHAPE, rigid_params=params)
    i = np.array([6.0, 9.0, 5.0])
    fixed_vox = 2 * i
    # apply_rigid output (visit-2 grid) voxel x samples visit-k at R x; moving_points_to_fixed is its inverse
    moving_vox = np.linalg.solve(np.eye(3), (ch.rigid_voxel_map(params, shape=ATLAS_SHAPE, spacing_mm=(1, 1, 1))
                                             @ np.r_[fixed_vox, 1])[:3])
    assert np.allclose(moving_points_to_fixed(moving_vox[None], params, shape=np.array(ATLAS_SHAPE),
                                              spacing_mm=(1, 1, 1))[0], fixed_vox)
    native_vox = np.linalg.solve(N, np.linalg.solve(T, ATLAS_AFF @ np.r_[moving_vox, 1]))
    assert np.allclose(m @ np.r_[i, 1], native_vox)


def test_rigid_factor_reproduces_apply_rigid_resampling():
    params = [0.0, 0.0, 5.0, 3.0, 0.0, -2.0]
    vol = np.zeros(ATLAS_SHAPE)
    vol[10:20, 12:25, 8:15] = 1
    R = ch.rigid_voxel_map(params, shape=ATLAS_SHAPE, spacing_mm=(1, 1, 1))
    got = affine_transform(vol, R[:3, :3], offset=R[:3, 3], order=0)
    assert np.array_equal(got, apply_rigid(vol, params, spacing_mm=(1, 1, 1), order=0))


def _physical_coverage(m_world_to_native_vox, native_shape, atlas_shape, rigid=None):
    """Brute force: every scoring voxel centre, through physical space; covered if inside native centre hull."""
    out = np.zeros(ch.scoring_shape(atlas_shape), bool)
    R = np.eye(4) if rigid is None else rigid
    for idx in np.ndindex(out.shape):
        atlas_vox = R @ np.r_[2 * np.array(idx), 1.0]
        native = m_world_to_native_vox @ ATLAS_AFF @ atlas_vox
        x = native[:3]                    # covered = inside the hull of native voxel centres [0, n-1]
        out[idx] = bool(np.all(x >= -1e-6) and np.all(x <= np.array(native_shape) - 1 + 1e-6))
    return out


def test_coverage_agrees_with_physical_space_resampling_anisotropic_and_rotated():
    N, nshape = _native(rot=_rot(-10, 1))
    T = _n2a()
    for params in (None, [2.0, -3.0, 4.0, 1.5, 2.5, -2.0]):
        cov = ch.visit_coverage("v", native_shape=nshape, native_affine=N, tfm_native_to_atlas_ras=T,
                                atlas_affine=ATLAS_AFF, atlas_shape=ATLAS_SHAPE, rigid_params=params).mask
        R = None if params is None else ch.rigid_voxel_map(params, shape=ATLAS_SHAPE, spacing_mm=(1, 1, 1))
        ref = _physical_coverage(np.linalg.inv(N) @ np.linalg.inv(T), nshape, ATLAS_SHAPE, R)
        assert cov.any() and not cov.all()
        m = ch.visit_output_to_input(native_affine=N, tfm_native_to_atlas_ras=T, atlas_affine=ATLAS_AFF,
                                     atlas_shape=ATLAS_SHAPE, rigid_params=params)
        for idx in np.argwhere(cov != ref):       # only exact floating-point ties on the hull edge may differ
            x = (m @ np.r_[idx, 1.0])[:3]
            assert np.min(np.abs(np.r_[x, np.array(nshape) - 1 - x])) < 1e-6, (idx, x)
        assert (cov != ref).sum() <= 2


def test_wrong_direction_gives_different_coverage():
    N, nshape = _native()
    T = _n2a(shift=(6.0, 0.0, 4.0), deg=0.0)
    right = ch.visit_coverage("v", native_shape=nshape, native_affine=N, tfm_native_to_atlas_ras=T,
                              atlas_affine=ATLAS_AFF, atlas_shape=ATLAS_SHAPE).mask
    wrong = ch.visit_coverage("v", native_shape=nshape, native_affine=N, tfm_native_to_atlas_ras=np.linalg.inv(T),
                              atlas_affine=ATLAS_AFF, atlas_shape=ATLAS_SHAPE).mask
    assert not np.array_equal(right, wrong)


def test_boundary_full_and_empty_and_thin_slab():
    # native grid identical to atlas grid, identity transform -> full coverage
    full = ch.visit_coverage("v", native_shape=ATLAS_SHAPE, native_affine=ATLAS_AFF, tfm_native_to_atlas_ras=np.eye(4),
                             atlas_affine=ATLAS_AFF, atlas_shape=ATLAS_SHAPE).mask
    assert full.all()
    # native grid shifted far away -> empty coverage
    far = ATLAS_AFF.copy()
    far[:3, 3] += 500
    assert not ch.visit_coverage("v", native_shape=ATLAS_SHAPE, native_affine=far, tfm_native_to_atlas_ras=np.eye(4),
                                 atlas_affine=ATLAS_AFF, atlas_shape=ATLAS_SHAPE).mask.any()
    # 2 native slices of 6 mm: centres at atlas z 0 and 6 -> conservative hull rule covers atlas z 0..6,
    # i.e. scoring z 0..3 (atlas 0, 2, 4, 6); the outer half-voxel bands count as unobserved
    thin = ATLAS_AFF @ np.diag([1, 1, 6, 1])
    cov = ch.visit_coverage("v", native_shape=(40, 44, 2), native_affine=thin, tfm_native_to_atlas_ras=np.eye(4),
                            atlas_affine=ATLAS_AFF, atlas_shape=ATLAS_SHAPE).mask
    assert cov[:, :, :4].all() and not cov[:, :, 4:].any()


def test_common_support_over_three_visits_with_both_rigid_transforms():
    N, nshape = _native(shape=(90, 100, 14))
    covs = [ch.visit_coverage("visit1", native_shape=nshape, native_affine=N, tfm_native_to_atlas_ras=_n2a(),
                              atlas_affine=ATLAS_AFF, atlas_shape=ATLAS_SHAPE, rigid_params=[1, 0, 2, 2, 0, -1]),
            ch.visit_coverage("visit2", native_shape=nshape, native_affine=N, tfm_native_to_atlas_ras=_n2a(),
                              atlas_affine=ATLAS_AFF, atlas_shape=ATLAS_SHAPE),
            ch.visit_coverage("target", native_shape=nshape, native_affine=N, tfm_native_to_atlas_ras=_n2a(),
                              atlas_affine=ATLAS_AFF, atlas_shape=ATLAS_SHAPE, rigid_params=[0, 2, -1, -2, 1, 3])]
    sup = common_support(covs, ("visit1", "visit2", "target")).mask
    assert np.array_equal(sup, covs[0].mask & covs[1].mask & covs[2].mask)
    assert sup.sum() < covs[1].mask.sum()


def test_domain_uses_2i_rule():
    m = np.zeros(ATLAS_SHAPE, bool)
    m[4, 6, 8] = True
    m[5, 7, 9] = True
    d = ch.domain_on_scoring_grid(m)
    assert d.shape == ch.scoring_shape(ATLAS_SHAPE) and d[2, 3, 4] and d.sum() == 1
