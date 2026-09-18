#!/usr/bin/env python3
"""
Real-Patient DTI Integration
=============================
Loads patient-specific diffusion tensor imaging (DTI) volumes and builds a
3x3 (or 2x2 mid-slice) symmetric positive-definite diffusion tensor field
aligned to actual white-matter pathways (corpus callosum, association fibers),
replacing the synthetic diagonal-tract assumption in src/42_anisotropic_pde.py.

Tensor construction
-------------------
For each voxel with principal eigenvector v1(x) (the dominant white-matter
fiber direction) and eigenvalues l1 >= l2:

    D(x) = lam_parallel * v1 v1^T + lam_perp * (I - v1 v1^T)

where lam_parallel / lam_perp are scalars (default 0.013 / 0.0013 mm^2/day,
Swanson et al. 2003) optionally scaled by a patient-specific gene factor.

Supported inputs
----------------
- NIfTI/Analyze (.nii/.nii.gz) DTI tensor volumes via nibabel:
    * 4D (X, Y, Z, 6) lower-triangular symmetric tensor (DTI-TK / FSL)
    * 4D (X, Y, Z, 9) full 3x3 flattened row-major tensor
    * 5D (3, 3, X, Y, Z) explicit tensor field
- Pre-decomposed eigen-system: FA volume + 3 eigenvector volumes
  (.nii.gz with suffixes _fa, _v1, _v2 or stacked)
- Pre-baked .npz with keys: D_xx/D_xy/D_yy (2D) or fa/eigvecs (3D)

Optional DIPY tractography verification is provided (step 5 of the plan)
so the loaded eigenvectors can be visually checked against simulated invasion.

Dependencies: nibabel (required), dipy (optional, for tractography/resampling).
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.ndimage import gaussian_filter, zoom as nd_zoom

# Physical defaults (match src/42_anisotropic_pde.py)
D_PARALLEL_DEFAULT = 0.013    # mm^2/day  (along-tract)
D_PERPENDICULAR_DEFAULT = 0.0013   # mm^2/day (~D_parallel/10, cross-tract)
D_BASE = 0.0013              # mm^2/day isotropic baseline (gray matter)

# Lower-triangular Voigt ordering used by FSL / DTI-TK tensors:
# [Dxx, Dxy, Dxz, Dyy, Dyz, Dzz]
VOIGT_LOWER = [0, 1, 2, 4, 5, 8]


def _lazy_import_nibabel():
    try:
        import nibabel as nib
        return nib
    except ImportError as e:
        raise ImportError(
            "nibabel is required to load real DTI volumes. "
            "Install with: pip install nibabel"
        ) from e


def _lazy_import_dipy():
    try:
        import dipy
        return dipy
    except ImportError:
        return None


# --------------------------------------------------------------------------- #
# NIfTI / tensor field loading
# --------------------------------------------------------------------------- #
def load_tensor_volume(nifti_path: Path) -> np.ndarray:
    """Load a raw DTI tensor volume from NIfTI into a (3, 3, X, Y, Z) array."""
    nib = _lazy_import_nibabel()
    img = nib.load(str(nifti_path))
    data = np.asarray(img.get_fdata(dtype=np.float32))
    aff = img.affine

    if data.ndim == 5 and data.shape[0] == 3 and data.shape[1] == 3:
        return np.transpose(data, (0, 1, 2, 3, 4)).astype(np.float64)

    if data.ndim == 4:
        X, Y, Z, C = data.shape
        # Lower-triangular Voigt (6 components)
        if C == 6:
            full = np.zeros((3, 3, X, Y, Z), dtype=np.float64)
            comps = [data[..., i] for i in range(6)]
            full[0, 0] = comps[0]; full[0, 1] = comps[1]; full[0, 2] = comps[2]
            full[1, 0] = comps[1]; full[1, 1] = comps[3]; full[1, 2] = comps[4]
            full[2, 0] = comps[2]; full[2, 1] = comps[4]; full[2, 2] = comps[5]
            return full
        # Full 3x3 flattened row-major (9 components)
        if C == 9:
            full = np.transpose(data.reshape(X, Y, Z, 3, 3), (3, 4, 0, 1, 2))
            return full.astype(np.float64)

    raise ValueError(
        f"Unrecognized DTI tensor volume shape {data.shape}. Expected "
        "(X,Y,Z,6), (X,Y,Z,9) or (3,3,X,Y,Z)."
    )


def load_fa_and_eigvecs(
    fa_path: Path,
    v1_path: Optional[Path] = None,
    v2_path: Optional[Path] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load FA volume and eigenvector volumes (v1, v2). v3 = v1 x v2."""
    nib = _lazy_import_nibabel()
    fa = np.asarray(nib.load(str(fa_path)).get_fdata(dtype=np.float32))

    def _load_vec(p: Path) -> np.ndarray:
        d = np.asarray(nib.load(str(p)).get_fdata(dtype=np.float32))
        # Eigenvector volumes are commonly (X, Y, Z, 3) or (3, X, Y, Z)
        if d.ndim == 4 and d.shape[-1] == 3:
            return np.transpose(d, (3, 0, 1, 2))
        if d.ndim == 4 and d.shape[0] == 3:
            return d
        raise ValueError(f"Eigenvector volume {p} has unsupported shape {d.shape}")

    v1 = _load_vec(v1_path) if v1_path else None
    v2 = _load_vec(v2_path) if v2_path else None
    return fa, v1, v2


def eigvec_from_tensor(tensor: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Principal eigenvector v1(x) from a (3,3,X,Y,Z) tensor volume."""
    T = np.moveaxis(tensor, 0, -2)  # (3, 3, X, Y, Z) -> (3, X, Y, Z, 3)? -> need (X,Y,Z,3,3)
    # Build (X, Y, Z, 3, 3) by moving both 3x3 axes to the end
    T2 = np.moveaxis(tensor, [0, 1], [-2, -1])  # (X,Y,Z,3,3)
    flat = T2.reshape(-1, 3, 3)
    n = flat.shape[0]
    eigs = np.zeros((n, 3, 3), dtype=np.float64)
    vecs = np.zeros((n, 3, 3), dtype=np.float64)
    for i in range(n):
        w, v = np.linalg.eigh(flat[i])
        eigs[i] = w
        vecs[i] = v
    shape = T2.shape[:3]
    eigs = eigs.reshape(*shape, 3)     # (X,Y,Z,3) ascending eigenvalues
    vecs = vecs.reshape(*shape, 3, 3)  # (X,Y,Z,3,3) columns = eigenvectors
    # eigh returns ascending eigenvalues; v1 = eigenvector of largest eigenvalue
    l1 = eigs[..., 2]
    v1 = vecs[..., :, 2]
    v2 = vecs[..., :, 1]
    return l1, v1, v2


# --------------------------------------------------------------------------- #
# Resampling to model grid
# --------------------------------------------------------------------------- #
def resample_tensor_to_grid(
    tensor: np.ndarray,
    target_shape: Tuple[int, int, int],
    order: int = 1,
) -> np.ndarray:
    """Resample a (3,3,X,Y,Z) tensor volume to target_shape via component-wise zoom."""
    T = tensor
    src = T.shape[2:]
    zf = tuple(target_shape[i] / max(float(src[i]), 1.0) for i in range(3))
    out = np.zeros((3, 3, *target_shape), dtype=np.float64)
    for i in range(3):
        for j in range(3):
            out[i, j] = nd_zoom(T[i, j], zf, order=order)
    # Symmetrize after interpolation to kill small antisymmetric residuals
    out = 0.5 * (out + np.transpose(out, (1, 0, 2, 3, 4)))
    return out


def resample_volume(vol: np.ndarray, target_shape: Tuple[int, ...],
                    order: int = 1) -> np.ndarray:
    src = vol.shape
    if len(src) == 3:
        zf = tuple(target_shape[i] / max(float(src[i]), 1.0) for i in range(3))
        return nd_zoom(vol, zf, order=order)
    if len(src) == 4 and src[0] in (1, 3):
        zf = (1.0,) + tuple(target_shape[i] / max(float(src[i + 1]), 1.0) for i in range(3))
        return nd_zoom(vol, zf, order=order)
    raise ValueError(f"Unsupported volume shape for resampling: {src}")


def brain_mask_from_fa(fa: np.ndarray, fa_threshold: float = 0.1) -> np.ndarray:
    """Crude brain/white-matter mask from FA (typical brain FA > 0.1)."""
    mask = (fa > fa_threshold).astype(bool)
    # Clean speckle
    mask = nd_zoom_binmask(mask)
    return mask


def nd_zoom_binmask(mask: np.ndarray) -> np.ndarray:
    from scipy.ndimage import binary_opening, binary_closing
    struct = np.ones((3,) * mask.ndim, dtype=bool)
    return binary_closing(binary_opening(mask, structure=struct), structure=struct)


# --------------------------------------------------------------------------- #
# PatientTensorBuilder
# --------------------------------------------------------------------------- #
class PatientTensorBuilder:
    """Build a 2x2 (mid-slice) or 3x3 diffusion tensor field from real DTI.

    Parameters
    ----------
    tensor_field : np.ndarray, shape (3,3,X,Y,Z) or None
        Pre-loaded full DTI tensor volume. If None, fa_path + eigvec_path(s)
        must be supplied.
    fa_path : Path or None
        Fractional anisotropy NIfTI volume (used for brain masking).
    v1_path : Path or None
        Principal eigenvector NIfTI volume (X,Y,Z,3) or (3,X,Y,Z).
    gene_scale : float or np.ndarray
        Multiplicative patient-specific scaling on lam_parallel / lam_perp
        (e.g. derived from TARGET_GENES), preserving transcriptomic influence.
    """

    def __init__(
        self,
        tensor_volume: Optional[np.ndarray] = None,
        fa_path: Optional[Path] = None,
        v1_path: Optional[Path] = None,
        v2_path: Optional[Path] = None,
        lam_parallel: float = D_PARALLEL_DEFAULT,
        lam_perp: float = D_PERPENDICULAR_DEFAULT,
        d_base: float = D_BASE,
        gene_scale: float = 1.0,
        fa_threshold: float = 0.1,
        target_shape_3d: Optional[Tuple[int, int, int]] = None,
    ) -> None:
        self.tensor_volume = tensor_volume
        self.fa_path = fa_path
        self.v1_path = v1_path
        self.v2_path = v2_path
        self.lam_parallel = float(lam_parallel)
        self.lam_perp = float(lam_perp)
        self.d_base = float(d_base)
        self.gene_scale = gene_scale
        self.fa_threshold = fa_threshold
        self.target_shape_3d = target_shape_3d

        # Outputs (3D)
        self.D_3d: Optional[np.ndarray] = None   # (3,3,Nx,Ny,Nz)
        self.fa: Optional[np.ndarray] = None
        self.v1: Optional[np.ndarray] = None
        self.brain_mask: Optional[np.ndarray] = None
        # 2D mid-slice outputs (compatible with src/42 anisotropic solver)
        self.D_xx: Optional[np.ndarray] = None
        self.D_xy: Optional[np.ndarray] = None
        self.D_yy: Optional[np.ndarray] = None
        self.tract_mask: Optional[np.ndarray] = None
        self.theta_field: Optional[np.ndarray] = None
        self.lambda_1: Optional[np.ndarray] = None
        self.lambda_2: Optional[np.ndarray] = None

    # ------------------------------------------------------------------ #
    @classmethod
    def from_nifti(
        cls,
        tensor_path: Path,
        fa_path: Optional[Path] = None,
        target_shape_3d: Optional[Tuple[int, int, int]] = None,
        **kwargs,
    ) -> "PatientTensorBuilder":
        """Load a full DTI tensor NIfTI file and build the field."""
        tensor = load_tensor_volume(tensor_path)
        b = cls(tensor_volume=tensor, fa_path=fa_path,
                target_shape_3d=target_shape_3d, **kwargs)
        b.build_3d()
        return b

    @classmethod
    def from_fa_eigvecs(
        cls,
        fa_path: Path,
        v1_path: Path,
        v2_path: Optional[Path] = None,
        target_shape_3d: Optional[Tuple[int, int, int]] = None,
        **kwargs,
    ) -> "PatientTensorBuilder":
        """Build the tensor field from FA + eigenvector volumes."""
        fa, v1, v2 = load_fa_and_eigvecs(fa_path, v1_path, v2_path)
        b = cls(fa_path=fa_path, v1_path=v1_path, v2_path=v2_path,
                target_shape_3d=target_shape_3d, **kwargs)
        b.fa = fa
        b.v1 = v1
        b.v2 = v2
        b._build_tensor_from_eigvecs()
        return b

    # ------------------------------------------------------------------ #
    def build_3d(self) -> np.ndarray:
        """Resample / reconstruct the 3D tensor field and brain mask."""
        if self.tensor_volume is not None:
            if self.target_shape_3d is not None and \
                    self.tensor_volume.shape[2:] != self.target_shape_3d:
                self.tensor_volume = resample_tensor_to_grid(
                    self.tensor_volume, self.target_shape_3d
                )
            self.D_3d = self.tensor_volume
            # FA + v1 derived from the tensor
            l1, v1, v2 = eigvec_from_tensor(self.D_3d)
            self.lambda_1 = l1
            self.v1 = np.transpose(v1, (3, 0, 1, 2))  # (3,X,Y,Z)
        # FA for masking
        if self.fa_path is not None and self.fa is None:
            nib = _lazy_import_nibabel()
            fa = np.asarray(nib.load(str(self.fa_path)).get_fdata(dtype=np.float32))
            if self.target_shape_3d is not None and fa.shape != self.target_shape_3d:
                fa = resample_volume(fa, self.target_shape_3d, order=1)
            self.fa = fa
        # Mask
        if self.fa is not None:
            self.brain_mask = brain_mask_from_fa(self.fa, self.fa_threshold)
        else:
            self.brain_mask = np.ones(self.D_3d.shape[2:], dtype=bool)
        self._apply_mask_and_baseline()
        return self.D_3d

    def _build_tensor_from_eigvecs(self) -> np.ndarray:
        """Construct D = lam_par v1 v1^T + lam_perp (I - v1 v1^T)."""
        if self.v1 is None:
            raise ValueError("v1 eigenvector volume required to build tensor.")
        if self.target_shape_3d is not None and \
                self.v1.shape[1:] != self.target_shape_3d:
            self.v1 = resample_volume(self.v1, self.target_shape_3d, order=1)
            if self.v2 is not None:
                self.v2 = resample_volume(self.v2, self.target_shape_3d, order=1)
            if self.fa is not None:
                self.fa = resample_volume(self.fa, self.target_shape_3d, order=1)

        v1 = self.v1  # (3, X, Y, Z)
        # Normalize per-voxel
        nrm = np.sqrt((v1 ** 2).sum(axis=0))
        v1 = v1 / np.maximum(nrm[None], 1e-12)
        v1 = v1.reshape(self.v1.shape)
        X, Y, Z = self.v1.shape[1:]
        I = np.eye(3)[:, :, None, None, None]  # (3,3,1,1,1)
        # v1 v1^T -> (3,3,X,Y,Z)
        outer = v1[:, None] * v1[None, :]   # (3,3,X,Y,Z)
        lam_par = self.lam_parallel * self.gene_scale
        lam_perp = self.lam_perp * self.gene_scale
        D = lam_par * outer + lam_perp * (I - outer)
        self.D_3d = D

        if self.fa is not None:
            self.brain_mask = brain_mask_from_fa(self.fa, self.fa_threshold)
        else:
            self.brain_mask = np.ones((X, Y, Z), dtype=bool)
        self._apply_mask_and_baseline()
        return D

    def _apply_mask_and_baseline(self) -> None:
        """Outside-brain voxels -> isotropic D_BASE."""
        if self.D_3d is None or self.brain_mask is None:
            return
        mask = self.brain_mask
        I = np.eye(3)[:, :, None, None, None]
        D_base_iso = self.d_base * I * np.ones_like(self.D_3d)  # broadcast
        in_brain = mask[None, None, ...]
        # Where not brain, set to isotropic baseline
        self.D_3d = np.where(in_brain, self.D_3d, D_base_iso)
        # Symmetrize & SPD clamp
        self.D_3d = 0.5 * (self.D_3d + np.transpose(self.D_3d, (1, 0, 2, 3, 4)))
        # Tiny floor to guarantee positive-definiteness
        for i in range(3):
            self.D_3d[i, i] = np.maximum(self.D_3d[i, i], 1e-6)

    # ------------------------------------------------------------------ #
    def mid_slice_2d(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Extract the middle Z-slice returning (D_xx, D_xy, D_yy) for the
        2D anisotropic solver in src/42_anisotropic_pde.py."""
        if self.D_3d is None:
            self.build_3d()
        mid = self.D_3d.shape[4] // 2
        D2 = self.D_3d[:, :, :, :, mid]  # (3,3,X,Y)
        self.D_xx = D2[0, 0].astype(np.float64)
        self.D_xy = D2[0, 1].astype(np.float64)
        self.D_yy = D2[1, 1].astype(np.float64)
        # 2D tract mask = inside brain & anisotropic voxels
        if self.brain_mask is not None:
            b2 = self.brain_mask[:, :, mid]
        else:
            b2 = np.ones_like(self.D_xx, dtype=bool)
        # tensor -> 2x2 eigen system
        tr = self.D_xx + self.D_yy
        disc = np.sqrt(np.maximum((self.D_xx - self.D_yy) ** 2 + 4 * self.D_xy ** 2, 0))
        self.lambda_1 = 0.5 * (tr + disc)
        self.lambda_2 = 0.5 * (tr - disc)
        # theta = angle of principal eigenvector
        self.theta_field = 0.5 * np.arctan2(2 * self.D_xy, self.D_xx - self.D_yy)
        # tract mask = anisotropic in-brain voxels (lambda1 > lambda2)
        self.tract_mask = (b2 & (self.lambda_1 > self.lambda_2 * 1.05))
        return self.D_xx, self.D_xy, self.D_yy

    # ------------------------------------------------------------------ #
    def validate(self) -> Dict[str, float]:
        """SPD + symmetry validation over the 3D tensor field."""
        if self.D_3d is None:
            self.build_3d()
        T = self.D_3d
        sym_err = float(np.max(np.abs(T - np.transpose(T, (1, 0, 2, 3, 4)))))
        # min eigenvalue
        T2 = np.moveaxis(T, [0, 1], [-2, -1]).reshape(-1, 3, 3)
        eigs = np.linalg.eigvalsh(T2)
        min_eig = float(eigs.min())
        min_trace = float(np.trace(T, axis1=0, axis2=1).min())
        det = np.linalg.det(T2)
        min_det = float(det.min())
        metrics = {
            "symmetry_max_error": sym_err,
            "symmetry_pass": bool(sym_err < 1e-10),
            "min_eigenvalue": min_eig,
            "min_trace": min_trace,
            "min_determinant": min_det,
            "positive_definite_pass": bool(min_eig > 0),
        }
        print("[DTI] PatientTensorBuilder validation:")
        print(f"  symmetry max |D - D^T| = {sym_err:.3e} "
              f"{'PASS' if metrics['symmetry_pass'] else 'FAIL'}")
        print(f"  min eigenvalue           = {min_eig:.3e} "
              f"{'PASS' if metrics['positive_definite_pass'] else 'FAIL'}")
        return metrics

    # ------------------------------------------------------------------ #
    def verify_streamlines(
        self,
        seed: Optional[Tuple[int, int, int]] = None,
        n_seeds: int = 200,
        output_trk: Optional[Path] = None,
    ) -> Dict:
        """(Step 5) Run deterministic tractography on v1 to sanity-check the
        loaded eigenvectors against the simulated invasion front.

        Requires dipy. Returns a dict with the number of streamlines found
        and their mean length, and optionally saves a .trk track file.
        """
        dipy = _lazy_import_dipy()
        if dipy is None or self.v1 is None:
            return {"dipy_available": False, "streamlines": 0}
        from dipy.tracking import utils as dtu
        from dipy.tracking.local_tracking import LocalTracking
        from dipy.tracking.stopping_criterion import ThresholdStoppingCriterion
        from dipy.tracking.utils import random_seeds_from_mask
        try:
            from dipy.direction import DeterministicMaximumDirectionGetter
        except Exception:
            from dipy.tracking.localtracking import DeterministicMaximumDirectionGetter

        if self.fa is None:
            return {"dipy_available": True, "streamlines": 0, "error": "no FA volume"}
        crit = ThresholdStoppingCriterion(self.fa, self.fa_threshold)
        seed_mask = self.brain_mask
        if seed is not None:
            seed_mask = np.zeros_like(seed_mask)
            seed_mask[seed] = True

        seeds = random_seeds_from_mask(seed_mask, seeds_count=n_seeds,
                                       affine=np.eye(4), random_seed=42)
        dg = DeterministicMaximumDirectionGetter.from_pmf(
            np.moveaxis(self.v1, 0, -1).clip(min=0), max_angle=30.
        )
        tracker = LocalTracking(dg, crit, seeds, np.eye(4),
                                step_size=0.5, maxlen=200)
        streamlines = list(tracker)
        if output_trk is not None:
            from dipy.io.streamline import save_trk
            save_trk(str(output_trk), streamlines, np.eye(4),
                     self.fa.shape, squeeze_bbox=False)
        lengths = [len(s) for s in streamlines] if streamlines else [0]
        return {
            "dipy_available": True,
            "streamlines": len(streamlines),
            "mean_length_vox": float(np.mean(lengths)),
        }


# --------------------------------------------------------------------------- #
# Convenience: build a 2D tensor field from a real patient directory
# --------------------------------------------------------------------------- #
def build_patient_2d_tensor(
    tensor_path: Path,
    fa_path: Optional[Path] = None,
    target_grid_size: int = 100,
    gene_scale: float = 1.0,
    fa_threshold: float = 0.1,
) -> Dict[str, np.ndarray]:
    """Return a dict with D_xx, D_xy, D_yy, theta, tract_mask, lambda_1/2
    suitable for direct injection into src/42_anisotropic_pde.py's
    TensorFieldBuilder (overwrite D_xx/D_xy/D_yy and recompute)."""
    builder = PatientTensorBuilder.from_nifti(
        tensor_path,
        fa_path=fa_path,
        target_shape_3d=(target_grid_size,) * 3,
        gene_scale=gene_scale,
        fa_threshold=fa_threshold,
    )
    builder.mid_slice_2d()
    return {
        "D_xx": builder.D_xx,
        "D_xy": builder.D_xy,
        "D_yy": builder.D_yy,
        "theta": builder.theta_field,
        "tract_mask": builder.tract_mask,
        "lambda_1": builder.lambda_1,
        "lambda_2": builder.lambda_2,
    }


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Real-Patient DTI Loader (Proposal 1)")
    p.add_argument("--tensor", type=str, required=False,
                   help="NIfTI tensor volume (.nii/.nii.gz)")
    p.add_argument("--fa", type=str, required=False, help="FA NIfTI volume")
    p.add_argument("--grid-size", type=int, default=100)
    args = p.parse_args()
    if args.tensor:
        b = PatientTensorBuilder.from_nifti(Path(args.tensor),
                                            fa_path=Path(args.fa) if args.fa else None,
                                            target_shape_3d=(args.grid_size,) * 3)
        b.validate()
        D_xx, D_xy, D_yy = b.mid_slice_2d()
        print(f"2D mid-slice: D_xx {D_xx.shape}, tract voxels={b.tract_mask.sum()}")
    else:
        print("Provide --tensor <nii.gz> [--fa fa.nii.gz] [--grid-size 100]")
