"""LUMIERE scoring-grid -> native-voxel transform chain (amendment paper/trackB_lumiere_geometry_amendment_2026-11-10.md G4).

Direction of every factor (each maps an OUTPUT location to an INPUT location, as scipy.ndimage.affine_transform
expects):
  D     : scoring voxel i -> visit-2 atlas voxel 2i                         (G2)
  R_k   : visit-2 atlas voxel -> visit-k atlas voxel; the exact map apply_rigid samples with (identity for visit 2)
  A     : atlas voxel -> atlas RAS mm                                       (atlas header affine)
  inv(T): atlas RAS mm -> native RAS mm; stored ct1.tfm T maps native -> atlas (G1, brain-mask validated)
  inv(N): native RAS mm -> native voxel                                     (native header affine)
Composed: inv(N) @ inv(T) @ A @ R_k @ D.
"""
from __future__ import annotations

import numpy as np

from src.fk_lumiere_transform_audit import rigid_voxel_map
from src.track_b_coverage_scoring import coverage_from_geometry

SCORING_FACTOR = 2


def scoring_to_atlas() -> np.ndarray:
    return np.diag([float(SCORING_FACTOR)] * 3 + [1.0])


def scoring_shape(atlas_shape) -> tuple:
    return tuple(int((s + SCORING_FACTOR - 1) // SCORING_FACTOR) for s in atlas_shape)


def visit_output_to_input(*, native_affine, tfm_native_to_atlas_ras, atlas_affine, atlas_shape,
                          rigid_params=None, atlas_spacing_mm=(1.0, 1.0, 1.0)) -> np.ndarray:
    R = np.eye(4) if rigid_params is None else rigid_voxel_map(rigid_params, shape=atlas_shape,
                                                               spacing_mm=atlas_spacing_mm)
    m = (np.linalg.inv(np.asarray(native_affine, float)) @ np.linalg.inv(np.asarray(tfm_native_to_atlas_ras, float))
         @ np.asarray(atlas_affine, float) @ R @ scoring_to_atlas())
    if not np.allclose(m[3], [0, 0, 0, 1]) or not np.isfinite(m).all():
        raise ValueError("composed map is not a finite affine")
    return m


def visit_coverage(visit: str, *, native_shape, **kw):
    out_shape = scoring_shape(kw["atlas_shape"])
    return coverage_from_geometry(visit, tuple(int(s) for s in native_shape), out_shape, visit_output_to_input(**kw))


def domain_on_scoring_grid(atlas_mask_1mm) -> np.ndarray:
    """G2 nearest-neighbour rule: scoring voxel i takes the 1 mm value at index 2i."""
    return np.asarray(atlas_mask_1mm, bool)[::SCORING_FACTOR, ::SCORING_FACTOR, ::SCORING_FACTOR]
