"""Source-agnostic Track B preparation: visit triplets, origin-only support, MRI-only rigid registration QC.

Draft (2026-10-31). Nothing here reads target labels. Synthetic tests only; no patient data.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
from scipy import ndimage, optimize

# Frozen draft QC rule (set before any target mask is opened; change only by dated amendment).
REGISTRATION_NCC_THRESHOLD = 0.80
REGISTRATION_MAX_SHIFT_MM = 10.0
REGISTRATION_MAX_ROTATION_DEG = 10.0


@dataclass(frozen=True)
class Visit:
    patient_id: str
    week: int
    order: int  # same-week ordering suffix (week-000-1, week-000-2)
    has_sequences: bool
    has_label_map: bool
    has_support_mask: bool

    @property
    def usable(self) -> bool:
        return self.has_sequences and self.has_label_map and self.has_support_mask


def primary_triplets(visits: Sequence[Visit], *, min_gap_weeks: int, max_gap_weeks: int) -> dict:
    """First three consecutive usable post-operative visits per patient; keep if the v2->v3 gap is in the window.

    Uses only visit metadata (weeks and file presence). The window is applied to the first triplet only,
    so a patient is never searched for a later triplet that happens to fit.
    """
    by_patient: dict[str, list[Visit]] = {}
    for v in visits:
        by_patient.setdefault(v.patient_id, []).append(v)
    eligible, excluded = {}, {}
    for pid, vs in sorted(by_patient.items()):
        usable = sorted((v for v in vs if v.usable and v.week > 0), key=lambda v: (v.week, v.order))
        if len(usable) < 3:
            excluded[pid] = "fewer than 3 usable post-op visits"
            continue
        v1, v2, v3 = usable[:3]
        gap = v3.week - v2.week
        if min_gap_weeks <= gap <= max_gap_weeks:
            eligible[pid] = (v1, v2, v3)
        else:
            excluded[pid] = f"first-triplet gap {gap} weeks outside {min_gap_weeks}-{max_gap_weeks}"
    return {"eligible": eligible, "excluded": excluded}


def origin_support(origin_brain_mask, history_label_maps: Sequence = ()) -> np.ndarray:
    """Valid support at the forecast origin: the source-released brain mask of the origin visit.

    History label voxels outside that mask are added back, so a history core is never clipped. No target
    input exists in this signature, so target labels cannot define support.
    """
    support = np.asarray(origin_brain_mask, bool).copy()
    if not support.any():
        raise ValueError("origin brain mask is empty")
    for lab in history_label_maps:
        lab = np.asarray(lab)
        if lab.shape != support.shape:
            raise ValueError("history label map shape differs from origin brain mask")
        support |= lab > 0
    return support


def _ncc(a, b, mask) -> float:
    a, b = a[mask].astype(float), b[mask].astype(float)
    a, b = a - a.mean(), b - b.mean()
    den = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / den) if den > 0 else 0.0


def _rotation(angles_deg):
    ax, ay, az = np.deg2rad(angles_deg)
    rx = np.array([[1, 0, 0], [0, np.cos(ax), -np.sin(ax)], [0, np.sin(ax), np.cos(ax)]])
    ry = np.array([[np.cos(ay), 0, np.sin(ay)], [0, 1, 0], [-np.sin(ay), 0, np.cos(ay)]])
    rz = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
    return rz @ ry @ rx


def apply_rigid(volume, params, *, spacing_mm, order=1):
    """Resample ``volume`` under rigid params (3 rotations in degrees, 3 shifts in mm) about the grid centre."""
    vol = np.asarray(volume, float)
    sp = np.asarray(spacing_mm, float)
    rot = _rotation(params[:3])
    centre = (np.array(vol.shape) - 1) / 2.0
    shift_vox = np.asarray(params[3:], float) / sp
    # output voxel x samples input at R(x - c) + c + t (in voxel units, with anisotropic spacing handled)
    s, s_inv = np.diag(sp), np.diag(1.0 / sp)
    matrix = s_inv @ rot @ s
    offset = centre + shift_vox - matrix @ centre
    return ndimage.affine_transform(vol, matrix, offset=offset, order=order, mode="constant", cval=0.0)


def register_rigid_mri(fixed, moving, *, fixed_support, spacing_mm, downsample: int = 2) -> dict:
    """MRI-only rigid registration of ``moving`` onto ``fixed`` by maximising NCC inside ``fixed_support``.

    No label map is an input. Returns params, NCC before and after, and the frozen QC decision.
    """
    f = np.asarray(fixed, float)
    m = np.asarray(moving, float)
    sup = np.asarray(fixed_support, bool)
    if f.shape != m.shape or f.shape != sup.shape:
        raise ValueError("fixed, moving and support must share one grid")
    sl = tuple(slice(None, None, downsample) for _ in range(3))
    fs, ms, ss = f[sl], m[sl], sup[sl]
    sp_ds = tuple(float(s) * downsample for s in spacing_mm)

    def cost(p):
        return -_ncc(fs, apply_rigid(ms, p, spacing_mm=sp_ds), ss)

    res = optimize.minimize(cost, np.zeros(6), method="Powell",
                            options={"xtol": 1e-2, "ftol": 1e-4, "maxiter": 2000})
    coarse_success = bool(res.success)
    # Refine on the full grid from the coarse result.
    res = optimize.minimize(lambda p: -_ncc(f, apply_rigid(m, p, spacing_mm=spacing_mm), sup), res.x,
                            method="Powell", options={"xtol": 1e-3, "ftol": 1e-5, "maxiter": 2000})
    params = np.asarray(res.x, float)
    warped = apply_rigid(m, params, spacing_mm=spacing_mm)
    ncc_before, ncc_after = _ncc(f, m, sup), _ncc(f, warped, sup)
    qc = registration_qc(params, ncc_after)
    failed_stages = [name for name, success in (("coarse", coarse_success), ("full", bool(res.success)))
                     if not success]
    if failed_stages:
        qc["passed"] = False
        qc["reasons"].append("optimizer did not converge: " + ", ".join(failed_stages))
    return {"params": params, "warped": warped, "ncc_before": ncc_before, "ncc_after": ncc_after,
            "qc": qc}


def registration_qc(params, ncc_after) -> dict:
    """Frozen draft QC: pass only if NCC >= threshold and the transform is within the motion limits."""
    params = np.asarray(params, dtype=float)
    if params.shape != (6,) or not np.isfinite(params).all():
        return {"passed": False, "reasons": ["invalid or nonfinite transform parameters"],
                "max_rotation_deg": None, "shift_mm": None}
    rot = float(np.max(np.abs(params[:3])))
    shift = float(np.linalg.norm(params[3:]))
    reasons = []
    if not np.isfinite(ncc_after) or ncc_after < REGISTRATION_NCC_THRESHOLD:
        reasons.append(f"NCC {ncc_after:.3f} < {REGISTRATION_NCC_THRESHOLD}")
    if rot > REGISTRATION_MAX_ROTATION_DEG:
        reasons.append(f"rotation {rot:.1f} deg > {REGISTRATION_MAX_ROTATION_DEG}")
    if shift > REGISTRATION_MAX_SHIFT_MM:
        reasons.append(f"shift {shift:.1f} mm > {REGISTRATION_MAX_SHIFT_MM}")
    return {"passed": not reasons, "reasons": reasons, "max_rotation_deg": rot, "shift_mm": shift}


def exclusion_report(qc_by_patient: Mapping[str, Mapping]) -> dict:
    """Failure handling: a failed registration excludes the patient; report every exclusion and its reason."""
    failed = {pid: list(q["reasons"]) for pid, q in qc_by_patient.items() if not q["passed"]}
    return {"n_checked": len(qc_by_patient), "n_excluded": len(failed), "excluded": failed}
