#!/usr/bin/env python3
"""Forecast validation: anisotropic (DTI) vs isotropic tumour-growth models on real follow-up MRI.

The earlier version of this script started the PDE from a patient's tumour mask and scored
the result against that same mask, so it measured drift from the initial condition rather
than prediction. It also built the "anisotropic" operator as an isotropic Laplacian times a
scalar. This version is a forecast:

  scan 1 tumour mask  --(growth model, interval dt from the clinical sheet)-->  predicted scan 2
  score: Dice(predicted scan 2, observed scan 2 mask)

Data
  * MU-Glioma-Post (data/tcia/MU-Glioma-Post): longitudinal post-operative scans with tumour
    masks, all in the 240x240x155 1 mm SRI24 space. No DTI.
  * UCSF-PDGM (62 patients, one time point each, $UCSF_PDGM_DIR): raw DWI in native space
    plus scalar DTI maps and tumour masks in the same SRI24 space as MU-Glioma-Post.

Stage "atlas": fit a diffusion tensor to each UCSF-PDGM patient's raw DWI, check the gradient
sign convention with a fibre-coherence index, register native FA to the patient's SRI24 FA
(12-parameter affine, normalised cross-correlation), rotate the tensors into SRI24 space
(finite-strain reorientation) and average them over patients with tumour voxels excluded.
The result is a population DTI atlas at 2 mm.

Stage "forecast": for every MU-Glioma-Post patient with two scans and a positive interval,
solve du/dt = div(D grad u) + rho u (1 - u) on the 2 mm grid (brain domain, no-flux
boundary, full-tensor flux discretisation) from u0 = scan-1 mask and threshold u > 0.5.
Arms:
  anisotropic  D = d * T_atlas / MD_ref, optional anisotropy sharpening r (Jbabdi 2005)
  iso_same     isotropic, same local mean diffusivity and the same (rho, d) chosen for the
               anisotropic arm, so it differs from the anisotropic arm only in direction
  iso_homog    homogeneous isotropic D = d * I with its own fitted (rho, d)
  no_change    predicted scan 2 = scan-1 mask
(rho, d, r) are chosen per arm by 5-fold cross-validation over patients (mean Dice on the
training folds); every reported Dice is out-of-fold.

Stage "report": statistics and figure. --selftest checks the solver (mass conservation,
Fisher front speed, direction of anisotropic spread) and the Dice function.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path

import nibabel as nib
import numpy as np
import torch
from scipy import ndimage, optimize, stats

PROJECT_ROOT = Path(__file__).resolve().parent
OUT_DIR = PROJECT_ROOT / "output" / "forecast_validation"
UCSF_DIR = Path(os.environ.get("UCSF_PDGM_DIR", r"C:\Users\vihan\Downloads\ucsf-pdgm-grade23"))
MU_DIR = PROJECT_ROOT / "data" / "tcia" / "MU-Glioma-Post"
COHORT_JSON = PROJECT_ROOT / "output" / "mu_glioma_cohort.json"
ATLAS_NPZ = OUT_DIR / "dti_atlas_2mm.npz"
ATLAS_QC_JSON = OUT_DIR / "dti_atlas_qc.json"
GRID_JSON = OUT_DIR / "forecast_grid.json"
RESULTS_JSON = OUT_DIR / "forecast_results.json"
RESULTS_TSV = OUT_DIR / "forecast_results.tsv"
SELFTEST_JSON = OUT_DIR / "solver_selftest.json"
FIGURE = OUT_DIR / "forecast_results.png"
TAB, NL = "\t", "\n"
CACHE_DIR = OUT_DIR / "cache"                 # per-patient checkpoints (resume after interruption)

H = 2.0                       # mm, simulation / scoring grid
SRI_SHAPE = (240, 240, 155)
GRID2 = (120, 120, 77)        # 2 mm grid (last 1 mm slice dropped)
RHOS = [0.0, 0.01, 0.03, 0.1]            # 1/day
DS = [0.01, 0.03, 0.1, 0.3]              # mm^2/day, brain-median mean diffusivity
SHARPEN = [1.0, 10.0]                    # anisotropy sharpening factor r
U_VISIBLE = 0.5
MARGIN_VOX = 20                          # crop margin around scan-1 tumour (2 mm voxels)
ATLAS_MIN_COUNT = 10
N_FOLDS = 5
ATLAS_WORKERS = 2        # each atlas worker holds one ~0.5 GB DWI series in memory
FORECAST_WORKERS = 4     # 4 processes x 2 torch threads on 8 cores
torch.set_num_threads(os.cpu_count() or 4)


# ----------------------------------------------------------------------------- helpers
def to_2mm(vol: np.ndarray) -> np.ndarray:
    """Block-average a 240x240x155 1 mm volume onto the 120x120x77 2 mm grid."""
    v = np.asarray(vol, dtype=np.float32)[:240, :240, :154]
    return v.reshape(120, 2, 120, 2, 77, 2).mean(axis=(1, 3, 5))


def dice(a: np.ndarray, b: np.ndarray) -> float:
    a, b = a.astype(bool), b.astype(bool)
    s = a.sum() + b.sum()
    return float(2.0 * np.logical_and(a, b).sum() / s) if s > 0 else 1.0


def six_to_mat(d6: np.ndarray) -> np.ndarray:
    """(...,6) [xx,yy,zz,xy,xz,yz] -> (...,3,3)."""
    xx, yy, zz, xy, xz, yz = np.moveaxis(d6, -1, 0)
    return np.stack([np.stack([xx, xy, xz], -1), np.stack([xy, yy, yz], -1),
                     np.stack([xz, yz, zz], -1)], -2)


def mat_to_six(m: np.ndarray) -> np.ndarray:
    return np.stack([m[..., 0, 0], m[..., 1, 1], m[..., 2, 2], m[..., 0, 1], m[..., 0, 2],
                     m[..., 1, 2]], -1)


def fa_md(w: np.ndarray):
    md = w.mean(-1)
    num = np.sqrt(((w - md[..., None]) ** 2).sum(-1))
    den = np.sqrt((w ** 2).sum(-1))
    fa = np.sqrt(1.5) * np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    return fa, md


# ----------------------------------------------------------------------------- solver
class TensorFK:
    """Explicit full-tensor Fisher-Kolmogorov solver on a masked 3D grid.

    Flux form: face fluxes F_i = sum_j D_ij d_j u with the diagonal term as a compact face
    difference and the off-diagonal terms from averaged central differences (neighbours
    outside the domain replaced by the centre value). Faces touching a voxel outside the
    domain carry zero flux, so total mass is conserved exactly when rho = 0.
    """

    def __init__(self, d6: np.ndarray, mask: np.ndarray, h: float = H):
        self.h = h
        m = torch.as_tensor(mask.astype(np.float32))
        self.m = m
        D = torch.as_tensor(d6.astype(np.float32)) * m[..., None]
        comp = {k: D[..., i] for i, k in enumerate(["xx", "yy", "zz", "xy", "xz", "yz"])}
        self.Dmat = [[comp["xx"], comp["xy"], comp["xz"]],
                     [comp["xy"], comp["yy"], comp["yz"]],
                     [comp["xz"], comp["yz"], comp["zz"]]]
        self.face = []
        for ax in range(3):
            w = self._fa(m, ax) * self._fb(m, ax)
            self.face.append([0.5 * (self._fa(self.Dmat[ax][j], ax) + self._fb(self.Dmat[ax][j], ax)) * w
                              for j in range(3)])
        # Gershgorin bound on the largest eigenvalue -> stable explicit step
        self.has_off = bool(any(float(comp[k].abs().max()) > 0 for k in ("xy", "xz", "yz")))
        self.mnb = []
        for ax in range(3):
            n = m.shape[ax]
            z = torch.zeros_like(m.narrow(ax, 0, 1))
            self.mnb.append((torch.cat([m.narrow(ax, 1, n - 1), z], ax), torch.cat([z, m.narrow(ax, 0, n - 1)], ax)))
        rows = [sum(self.Dmat[i][j].abs() if i != j else self.Dmat[i][j] for j in range(3)) for i in range(3)]
        self.lam_max = float(torch.stack(rows).max().clamp_min(1e-12))

    @staticmethod
    def _fa(a, ax):   # value at cell i on the face between i and i+1
        return a.narrow(a.dim() - 3 + ax, 0, a.shape[a.dim() - 3 + ax] - 1)

    @staticmethod
    def _fb(a, ax):   # value at cell i+1
        return a.narrow(a.dim() - 3 + ax, 1, a.shape[a.dim() - 3 + ax] - 1)

    def _grad_c(self, u, ax):
        """Central difference along ax; out-of-domain neighbours replaced by centre value."""
        dim = u.dim() - 3 + ax
        n = u.shape[dim]
        up = torch.cat([u.narrow(dim, 1, n - 1), u.narrow(dim, n - 1, 1)], dim)
        un = torch.cat([u.narrow(dim, 0, 1), u.narrow(dim, 0, n - 1)], dim)
        mup, mun = self.mnb[ax]
        up = mup * up + (1 - mup) * u
        un = mun * un + (1 - mun) * u
        return (up - un) / (2 * self.h)

    def div_flux(self, u):
        g = [self._grad_c(u, ax) for ax in range(3)] if self.has_off else None
        out = torch.zeros_like(u)
        for ax in range(3):
            dim = u.dim() - 3 + ax
            du = (self._fb(u, ax) - self._fa(u, ax)) / self.h
            F = self.face[ax][ax] * du
            for j in range(3):
                if j != ax and g is not None:
                    F = F + self.face[ax][j] * 0.5 * (self._fa(g[j], ax) + self._fb(g[j], ax))
            z = torch.zeros_like(F.narrow(dim, 0, 1))
            Fp = torch.cat([z, F, z], dim)
            out = out + (Fp.narrow(dim, 1, Fp.shape[dim] - 1) - Fp.narrow(dim, 0, Fp.shape[dim] - 1)) / self.h
        return out

    def run(self, u0: np.ndarray, rhos, t_end: float, max_dt: float = 0.5):
        """Integrate a batch (one member per rho) to t_end; returns (B, X, Y, Z) numpy."""
        rho = torch.as_tensor(np.asarray(rhos, dtype=np.float32))[:, None, None, None]
        u = torch.as_tensor(u0.astype(np.float32))[None].repeat(len(rhos), 1, 1, 1) * self.m
        dt_max = min(max_dt, 0.9 * self.h ** 2 / (6.0 * self.lam_max),
                     0.05 / max(float(rho.max()), 1e-6))
        n = max(1, math.ceil(t_end / dt_max))
        dt = t_end / n
        with torch.no_grad():
            for _ in range(n):
                u = (u + dt * (self.div_flux(u) + rho * u * (1 - u))).clamp_(0.0, 1.0) * self.m
        return u.numpy(), n


def selftest() -> dict:
    res = {}
    rng = np.random.default_rng(0)
    # 1. mass conservation with a random SPD tensor field on an irregular domain
    shp = (28, 26, 24)
    zz, yy, xx = np.meshgrid(*[np.arange(s) for s in shp], indexing="ij")
    mask = ((xx - 12) ** 2 + (yy - 13) ** 2 + (zz - 14) ** 2) < 11 ** 2
    A = rng.normal(size=shp + (3, 3)) * 0.3
    T = A @ np.swapaxes(A, -1, -2) + 0.2 * np.eye(3)
    s = TensorFK(mat_to_six(T) * 0.3, mask)
    u0 = np.zeros(shp, np.float32)
    u0[10:18, 9:17, 8:16] = 1.0
    u0 *= mask
    ur = torch.as_tensor(rng.random(shp).astype(np.float32))[None] * s.m
    dv = s.div_flux(ur)
    # operator: face fluxes telescope, so the divergence sums to zero over the domain
    res["operator_sum_div_over_sum_abs_div"] = float(abs(dv.sum()) / dv.abs().sum())
    res["operator_conservation_pass"] = res["operator_sum_div_over_sum_abs_div"] < 1e-5
    # full run: clamping u >= 0 after off-diagonal undershoot at a sharp edge adds a little mass
    u, n = s.run(u0, [0.0], 100.0)
    res["run_mass_rel_change_with_clamp"] = float(abs(u.sum() - u0.sum()) / u0.sum())
    res["run_mass_pass"] = res["run_mass_rel_change_with_clamp"] < 1e-3
    # 2. Fisher front speed 2 sqrt(D rho), isotropic slab
    D, rho = 0.1, 0.05
    shp = (90, 4, 4)
    mask = np.ones(shp, bool)
    s = TensorFK(np.tile(np.array([D, D, D, 0, 0, 0], np.float32), shp + (1,)), mask)
    u0 = np.zeros(shp, np.float32)
    u0[:8] = 1.0
    pos = []
    for t in (200.0, 400.0):
        u, _ = s.run(u0, [rho], t)
        pos.append(float((u[0, :, 0, 0] > 0.5).sum() * H))
    speed = (pos[1] - pos[0]) / 200.0
    res["fisher_speed_numeric_mm_per_day"] = speed
    res["fisher_speed_analytic_mm_per_day"] = 2 * math.sqrt(D * rho)
    res["fisher_speed_rel_err"] = abs(speed / (2 * math.sqrt(D * rho)) - 1)
    res["fisher_speed_pass"] = res["fisher_speed_rel_err"] < 0.10
    # 3. anisotropic spread follows the principal direction (45 deg in the x-y plane)
    shp = (41, 41, 9)
    c = 1 / math.sqrt(2)
    R = np.array([[c, -c, 0], [c, c, 0], [0, 0, 1]])
    Tm = R @ np.diag([1.0, 0.1, 0.1]) @ R.T * 0.2
    s = TensorFK(np.tile(mat_to_six(Tm).astype(np.float32), shp + (1,)), np.ones(shp, bool))
    u0 = np.zeros(shp, np.float32)
    u0[19:22, 19:22, 3:6] = 1.0
    u, _ = s.run(u0, [0.0], 150.0)
    w = u[0].sum(2)
    X, Y = np.meshgrid(np.arange(41) - 20, np.arange(41) - 20, indexing="ij")
    C = np.array([[np.sum(w * X * X), np.sum(w * X * Y)], [np.sum(w * X * Y), np.sum(w * Y * Y)]]) / w.sum()
    ev, evec = np.linalg.eigh(C)
    ang = math.degrees(math.atan2(abs(evec[1, -1]), abs(evec[0, -1])))
    res["aniso_major_axis_deg"] = ang
    res["aniso_axis_ratio"] = float(math.sqrt(ev[-1] / ev[0]))
    res["aniso_pass"] = abs(ang - 45) < 5 and res["aniso_axis_ratio"] > 1.5
    # 4. Dice
    a = np.zeros((10, 10), bool); a[:5] = True
    b = np.zeros((10, 10), bool); b[2:7] = True
    res["dice_known_value"] = dice(a, b)                     # 2*30/(50+50) = 0.6
    res["dice_pass"] = abs(res["dice_known_value"] - 0.6) < 1e-12
    res["all_pass"] = all(v for k, v in res.items() if k.endswith("_pass"))
    return res


# ----------------------------------------------------------------------------- atlas stage
def ucsf_patients() -> list[str]:
    ids = []
    for d in sorted(UCSF_DIR.glob("UCSF-PDGM-*_nifti")):
        pid = d.name.replace("_nifti", "")
        if (d / f"{pid}_DTI_eddy_noreg.nii.gz").exists() and (d / f"{pid}_DTI_eddy_FA.nii.gz").exists():
            ids.append(pid)
    return ids


def fit_native_tensor(pid: str):
    d = UCSF_DIR / f"{pid}_nifti"
    img = nib.load(str(d / f"{pid}_DTI_eddy_noreg.nii.gz"))
    dwi = np.asarray(img.dataobj)                           # native dtype (int16), not float
    bval = np.loadtxt(UCSF_DIR / "UCSF-PDGM_DTI.bval").ravel()
    bvec = np.loadtxt(d / f"{pid}_DTI_eddy.eddy_rotated_bvecs")
    if np.linalg.det(img.affine[:3, :3]) > 0:     # FSL convention: bvecs in voxel frame, x flipped if det > 0
        bvec = bvec * np.array([-1, 1, 1])[:, None]
    b0 = bval < 50
    S0 = dwi[..., b0].astype(np.float32).mean(-1)
    from skimage.filters import threshold_otsu
    # Brain mask from the smoothed mean diffusion-weighted image: at b = 2000 tissue is bright
    # and CSF dark, so Otsu on b0 alone (bright CSF) can drop white matter.
    mean_dw = np.zeros(dwi.shape[:3], np.float32)
    for k in np.nonzero(~b0)[0]:
        mean_dw += dwi[..., k]
    mean_dw = ndimage.gaussian_filter(mean_dw / (~b0).sum(), 1.0)
    brain = mean_dw > threshold_otsu(mean_dw[mean_dw > 0])
    brain = ndimage.binary_fill_holes(ndimage.binary_closing(brain, iterations=3))
    lab, nlab = ndimage.label(brain)
    if nlab > 1:
        brain = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
    idx = np.nonzero(brain)
    S = dwi[idx][:, ~b0].astype(np.float32)
    g = bvec[:, ~b0].T
    B = bval[~b0][:, None] * np.stack([g[:, 0] ** 2, g[:, 1] ** 2, g[:, 2] ** 2, 2 * g[:, 0] * g[:, 1],
                                       2 * g[:, 0] * g[:, 2], 2 * g[:, 1] * g[:, 2]], 1)
    Y = -np.log(np.clip(S, 1, None) / np.clip(S0[idx], 1, None)[:, None])
    coef = Y @ np.linalg.pinv(B).T                          # mm^2/s, native voxel-axis frame
    D6 = np.zeros(S0.shape + (6,), np.float32)
    D6[idx] = coef
    vox = np.sqrt((img.affine[:3, :3] ** 2).sum(0))
    return D6, brain, img.affine, vox


def fibre_coherence(D6, brain, vox, flip) -> float:
    """Fraction of high-FA voxels whose neighbours one step along v1 (both ways) are also
    high-FA (after Schilling et al. 2019). Wrong gradient signs lower this index."""
    F = np.diag(flip).astype(np.float32)
    sel = brain.copy()
    T = six_to_mat(D6[sel])
    w, _ = np.linalg.eigh(T)
    fa, _ = fa_md(np.clip(w, 0, None))
    fa_vol = np.zeros(brain.shape, np.float32)
    fa_vol[sel] = fa
    hi = np.argwhere(fa_vol > 0.4)
    if len(hi) > 60000:
        hi = hi[np.random.default_rng(0).choice(len(hi), 60000, replace=False)]
    Th = F @ six_to_mat(D6[tuple(hi.T)]) @ F
    _, v = np.linalg.eigh(Th)
    v1 = v[..., -1]
    step = 1.5 * v1 / vox                                   # 1.5 mm along v1, in index units
    ok = []
    for sgn in (1, -1):
        p = (hi + sgn * step).T
        ok.append(ndimage.map_coordinates(fa_vol, p, order=0, mode="constant") > 0.4)
    return float(np.mean(ok[0] & ok[1]))


def register_affine(fixed: np.ndarray, moving: np.ndarray, A0: np.ndarray, t0: np.ndarray, fixed_step: float):
    """Find n = A s + t (fixed 1 mm index s -> moving index n) maximising NCC of FA.
    fixed is sampled on a grid with spacing fixed_step (1 mm units), cell centres at
    fixed_step * i + (fixed_step - 1) / 2."""
    sel = np.argwhere(fixed > 0.02)
    s = sel * fixed_step + (fixed_step - 1) / 2.0
    fv = fixed[tuple(sel.T)]
    fv = (fv - fv.mean()) / fv.std()

    def unpack(p):
        A = A0 @ (np.eye(3) + p[:9].reshape(3, 3))
        t = t0 + 10.0 * p[9:]
        return A, t

    def cost(p):
        A, t = unpack(p)
        n = s @ A.T + t
        mv = ndimage.map_coordinates(moving, n.T, order=1, mode="constant")
        sd = mv.std()
        return 1.0 if sd == 0 else -float(np.mean(fv * (mv - mv.mean()) / sd))

    p0 = np.zeros(12)
    r = optimize.minimize(cost, p0, method="Powell", options={"maxiter": 4000, "xtol": 1e-4, "ftol": 1e-6})
    A, t = unpack(r.x)
    return A, t, -float(r.fun)


def _rot(ax: int, deg: float) -> np.ndarray:
    c, s_ = math.cos(math.radians(deg)), math.sin(math.radians(deg))
    i, j = [k for k in range(3) if k != ax]
    R = np.eye(3)
    R[i, i], R[i, j], R[j, i], R[j, j] = c, -s_, s_, c
    return R


def rotation_init(fixed4, moving, A0, cs, cn):
    """Coarse search over rotations of the fixed grid about its centroid (pitch -40..40,
    roll and yaw -20..20 degrees) before the affine refinement; scans are often tilted."""
    sel = np.argwhere(fixed4 > 0.02)
    sp = sel * 4.0 + 1.5
    fv = fixed4[tuple(sel.T)]
    fv = (fv - fv.mean()) / fv.std()
    best = (-2.0, A0, cn - A0 @ cs)
    for p in range(-40, 41, 10):
        for r in range(-20, 21, 10):
            for y in range(-20, 21, 10):
                Rm = _rot(0, p) @ _rot(1, r) @ _rot(2, y)
                A = A0 @ Rm
                t = cn - A @ cs
                mv = ndimage.map_coordinates(moving, (sp @ A.T + t).T, order=1, mode="constant")
                if mv.std() == 0:
                    continue
                c = float(np.mean(fv * (mv - mv.mean()) / mv.std()))
                if c > best[0]:
                    best = (c, A, t)
    return best[1], best[2], best[0]


_THREAD_LIMIT = None


def _limit_threads(n: int) -> None:
    """Cap BLAS and torch threads inside a worker so workers do not oversubscribe the cores."""
    global _THREAD_LIMIT
    from threadpoolctl import threadpool_limits
    if _THREAD_LIMIT is None:
        _THREAD_LIMIT = threadpool_limits(n)
    torch.set_num_threads(n)


def atlas_patient(pid: str) -> dict:
    """Tensor fit, registration and reorientation for one UCSF-PDGM patient. Writes
    cache/atlas_<pid>.npz (valid voxel indices + SRI24 tensor) and returns the QC record."""
    cache = CACHE_DIR / f"atlas_{pid}.npz"
    if cache.exists():
        return json.loads(str(np.load(cache)["qc"]))
    _limit_threads(max(1, (os.cpu_count() or 8) // ATLAS_WORKERS))
    t0_ = time.time()
    s2 = np.stack(np.meshgrid(*[np.arange(n) for n in GRID2], indexing="ij"), -1).reshape(-1, 3)
    s1 = 2.0 * s2 + 0.5
    D6n, brain_n, aff_n, vox_n = fit_native_tensor(pid)
    flips = {"none": (1, 1, 1), "x": (-1, 1, 1), "y": (1, -1, 1), "z": (1, 1, -1)}
    fci = {k2: fibre_coherence(D6n, brain_n, vox_n, f) for k2, f in flips.items()}
    best_flip = max(fci, key=fci.get)
    Fm = np.diag(flips[best_flip]).astype(np.float32)
    Tn = Fm @ six_to_mat(D6n) @ Fm
    del D6n
    wn, _ = np.linalg.eigh(Tn[brain_n])
    fa_n = np.zeros(brain_n.shape, np.float32)
    fa_n[brain_n] = fa_md(np.clip(wn, 0, None))[0]

    d = UCSF_DIR / f"{pid}_nifti"
    fa_img = nib.load(str(d / f"{pid}_DTI_eddy_FA.nii.gz"))
    fa_s1 = np.clip(np.asarray(fa_img.dataobj, np.float32), 0, 1)
    seg = np.asarray(nib.load(str(d / f"{pid}_tumor_segmentation.nii.gz")).dataobj) > 0
    A0 = np.linalg.inv(aff_n[:3, :3]) @ fa_img.affine[:3, :3]
    cs = np.argwhere(fa_s1 > 0.05).mean(0)
    cn = np.argwhere(brain_n).mean(0)
    fa_s4 = fa_s1[:240, :240, :152].reshape(60, 4, 60, 4, 38, 4).mean(axis=(1, 3, 5))
    mov4 = ndimage.gaussian_filter(fa_n, 1.5)
    A0, t0, ncc_init = rotation_init(fa_s4, mov4, A0, cs, cn)
    A, t, _ = register_affine(fa_s4, mov4, A0, t0, 4.0)
    fa_s2 = to_2mm(fa_s1)
    A, t, ncc = register_affine(fa_s2, ndimage.gaussian_filter(fa_n, 0.7), A, t, 2.0)

    # resample tensor components and reorient (finite strain)
    n_idx = (s1 @ A.T + t).T
    T6 = mat_to_six(Tn)
    comps = np.stack([ndimage.map_coordinates(T6[..., i], n_idx, order=1, mode="constant")
                      for i in range(6)], -1).reshape(GRID2 + (6,))
    inside = ndimage.map_coordinates(brain_n.astype(np.float32), n_idx, order=1,
                                     mode="constant").reshape(GRID2) > 0.99
    J = np.diag(vox_n) @ A
    U, _, Vt = np.linalg.svd(J)
    R = U @ Vt
    Ts = R.T @ six_to_mat(comps) @ R
    ws, vs = np.linalg.eigh(Ts)
    fa_new, md_new = fa_md(np.clip(ws, 0, None))
    brain2 = (fa_s2 > 0.02) & inside
    fa_corr = float(np.corrcoef(fa_new[brain2], fa_s2[brain2])[0, 1])
    # corpus callosum check: mid-sagittal high-FA voxels should run left-right (axis 0)
    mid = np.zeros(GRID2, bool)
    mid[58:63, 40:80, 30:50] = True
    cc = mid & (fa_new > 0.5) & brain2
    cc_lr = float(np.mean(np.abs(vs[cc][:, 0, -1]) > 0.8)) if cc.sum() > 20 else None
    tum = ndimage.binary_dilation(to_2mm(seg) > 0, iterations=2)
    valid = brain2 & ~tum & (md_new > 0) & np.isfinite(ws).all(-1) & (ws[..., 0] >= 0)
    good = ncc > 0.5 and fa_corr > 0.5
    rec = {"patient_id": pid, "gradient_flip_by_coherence": best_flip, "fibre_coherence": fci,
           "registration_ncc_init_4mm": ncc_init, "registration_ncc_2mm": ncc,
           "fa_corr_vs_provided_sri24": fa_corr, "cc_fraction_left_right": cc_lr,
           "n_valid_voxels": int(valid.sum()), "included": bool(good), "seconds": round(time.time() - t0_, 1)}
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    tmp = CACHE_DIR / f"atlas_{pid}.tmp.npz"
    np.savez_compressed(tmp, idx=np.flatnonzero(valid).astype(np.int32),
                        d6=mat_to_six(Ts)[valid].astype(np.float32), qc=json.dumps(rec))
    tmp.replace(cache)
    return rec


def stage_atlas(limit: int | None = None) -> None:
    import multiprocessing as mp
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    ids = ucsf_patients()[:limit] if limit else ucsf_patients()
    qc = []
    with mp.get_context("spawn").Pool(ATLAS_WORKERS, maxtasksperchild=4) as pool:
        for k, rec in enumerate(pool.imap(atlas_patient, ids)):
            qc.append(rec)
            print(f"[{k + 1}/{len(ids)}] {rec['patient_id']} flip={rec['gradient_flip_by_coherence']} "
                  f"ncc={rec['registration_ncc_2mm']:.3f} fa_corr={rec['fa_corr_vs_provided_sri24']:.3f} "
                  f"cc_lr={rec['cc_fraction_left_right']} incl={rec['included']} {rec['seconds']}s", flush=True)
    acc = np.zeros((int(np.prod(GRID2)), 6), np.float64)
    cnt = np.zeros(int(np.prod(GRID2)), np.int32)
    for rec in qc:
        if rec["included"]:
            z = np.load(CACHE_DIR / f"atlas_{rec['patient_id']}.npz")
            acc[z["idx"]] += z["d6"]
            cnt[z["idx"]] += 1
    acc, cnt = acc.reshape(GRID2 + (6,)), cnt.reshape(GRID2)
    atlas = np.where(cnt[..., None] > 0, acc / np.maximum(cnt, 1)[..., None], 0).astype(np.float32)
    w, v = np.linalg.eigh(six_to_mat(atlas))
    fa, md = fa_md(np.clip(w, 0, None))
    inc = [q for q in qc if q["included"]]
    min_count = min(ATLAS_MIN_COUNT, max(1, len(inc)))      # < 10 only in --limit debug runs
    tissue = cnt >= min_count
    md_ref = float(np.median(md[tissue]))
    np.savez_compressed(ATLAS_NPZ, d6=atlas, count=cnt, fa=fa.astype(np.float32), md=md.astype(np.float32),
                        md_ref=md_ref, min_count=min_count)
    summary = {
        "n_patients_processed": len(qc), "n_included": len(inc),
        "gradient_flip_counts": {f: sum(q["gradient_flip_by_coherence"] == f for q in qc)
                                 for f in ("none", "x", "y", "z")},
        "registration_ncc_median": float(np.median([q["registration_ncc_2mm"] for q in qc])),
        "fa_corr_median": float(np.median([q["fa_corr_vs_provided_sri24"] for q in qc])),
        "cc_fraction_left_right_median": float(np.median([q["cc_fraction_left_right"] for q in inc
                                                          if q["cc_fraction_left_right"] is not None])),
        "atlas_voxels_count_ge_min": int(tissue.sum()), "atlas_min_count": min_count,
        "atlas_md_ref_mm2_per_s": md_ref,
        "atlas_fa_median_tissue": float(np.median(fa[tissue])),
        "atlas_fa_p90_tissue": float(np.percentile(fa[tissue], 90)),
        "inclusion_rule": "registration NCC > 0.5 and FA correlation with provided SRI24 FA > 0.5",
    }
    ATLAS_QC_JSON.write_text(json.dumps({"summary": summary, "patients": qc}, indent=2))
    print(json.dumps(summary, indent=2))


# ----------------------------------------------------------------------------- forecast stage
def load_atlas():
    z = np.load(ATLAS_NPZ)
    d6, cnt, md, md_ref = z["d6"], z["count"], z["md"], float(z["md_ref"])
    min_count = int(z["min_count"])
    w, v = np.linalg.eigh(six_to_mat(d6))
    w = np.clip(w, 0, None)
    tissue = (cnt >= min_count) & (md < 2.0 * md_ref)      # CSF (high MD) excluded
    return {"w": w / md_ref, "v": v, "md": md / md_ref, "tissue": tissue, "md_ref": md_ref}


def arm_tensor(atlas, arm: str, r: float, box) -> np.ndarray:
    """Dimensionless tensor field (brain-median mean diffusivity = 1) on the crop box."""
    w = atlas["w"][box]
    v = atlas["v"][box]
    if arm == "aniso":
        ws = w * np.array([1.0, 1.0, r], np.float32)
        ws = ws * (w.sum(-1, keepdims=True) / np.maximum(ws.sum(-1, keepdims=True), 1e-12))
        T = (v * ws[..., None, :]) @ np.swapaxes(v, -1, -2)
        return mat_to_six(T)
    if arm == "iso_same":
        md = atlas["md"][box]
        return np.stack([md, md, md, 0 * md, 0 * md, 0 * md], -1)
    if arm == "iso_homog":
        one = np.ones(w.shape[:3], np.float32)
        return np.stack([one, one, one, 0 * one, 0 * one, 0 * one], -1)
    raise ValueError(arm)


def forecast_pairs():
    cohort = json.loads(COHORT_JSON.read_text())
    pairs, excluded = [], []
    for p in cohort:
        tp = p["timepoints"]
        if len(tp) < 2:
            continue
        a, b = tp[0], tp[1]
        if a["day_from_diagnosis"] is None or b["day_from_diagnosis"] is None:
            excluded.append((p["patient_id"], "missing scan day"))
            continue
        dt = b["day_from_diagnosis"] - a["day_from_diagnosis"]
        if dt <= 0:
            excluded.append((p["patient_id"], f"non-positive interval {dt}"))
            continue
        pairs.append({"patient_id": p["patient_id"], "tp1": a["number"], "tp2": b["number"], "dt_days": dt})
    return pairs, excluded


def mu_paths(pid, tp):
    d = MU_DIR / pid / f"Timepoint_{tp}"
    return d / f"{pid}_Timepoint_{tp}_tumorMask.nii.gz", d / f"{pid}_Timepoint_{tp}_brain_t1n.nii.gz"


def forecast_patient(pair, atlas):
    m1p, t1p = mu_paths(pair["patient_id"], pair["tp1"])
    m2p, _ = mu_paths(pair["patient_id"], pair["tp2"])
    tum1 = to_2mm(np.asarray(nib.load(str(m1p)).dataobj) > 0) >= 0.5
    tum2 = to_2mm(np.asarray(nib.load(str(m2p)).dataobj) > 0) >= 0.5
    brain = to_2mm(np.asarray(nib.load(str(t1p)).dataobj) > 0) >= 0.5
    rec = {**pair, "v1_voxels": int(tum1.sum()), "v2_voxels": int(tum2.sum()),
           "grew": bool(tum2.sum() > tum1.sum())}
    if tum1.sum() == 0 or tum2.sum() == 0:
        rec["skip"] = "empty mask at 2 mm"
        return rec
    rec["dice"] = {"no_change": dice(tum1, tum2)}
    lo = np.maximum(np.argwhere(tum1).min(0) - MARGIN_VOX, 0)
    hi = np.minimum(np.argwhere(tum1).max(0) + MARGIN_VOX + 1, GRID2)
    box = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
    dom = (brain[box] & atlas["tissue"][box]) | tum1[box]
    u0 = tum1[box].astype(np.float32)
    n2 = int(tum2.sum())
    # shape-only secondary: distance-based uniform growth of scan 1 to the scan-2 volume
    dist = ndimage.distance_transform_edt(~tum1)
    dist[~brain] = np.inf
    rec["dice_volume_matched"] = {"uniform_dilation": dice(dist <= np.sort(dist.ravel())[n2 - 1], tum2)}
    rec["dice"]["grid"] = {}
    rec["dice_volume_matched"]["grid"] = {}
    runs = [("aniso", r) for r in SHARPEN] + [("iso_same", 1.0), ("iso_homog", 1.0)]
    steps = 0
    for arm, r in runs:
        T0 = arm_tensor(atlas, arm, r, box)
        T0[tum1[box] & ~atlas["tissue"][box]] = np.array([1, 1, 1, 0, 0, 0], np.float32)
        for dval in DS:
            solver = TensorFK(T0 * dval, dom)
            u, n = solver.run(u0, RHOS, pair["dt_days"])
            steps += n
            for i, rho in enumerate(RHOS):
                key = f"{arm}|r={r:g}|d={dval:g}|rho={rho:g}"
                full = np.zeros(GRID2, np.float32)
                full[box] = u[i]
                rec["dice"]["grid"][key] = dice(full > U_VISIBLE, tum2)
                # volume matched: top-n2 voxels by u (ties at u == 1 broken by distance to scan 1)
                score = full - 1e-6 * np.minimum(dist, 1e3)
                thr = np.sort(score.ravel())[-n2]
                rec["dice_volume_matched"]["grid"][key] = dice(score >= thr, tum2)
    rec["solver_steps"] = steps
    return rec


_ATLAS = None
_ATLAS_SHA = None


def _atlas_sha() -> str:
    import hashlib
    return hashlib.sha256(ATLAS_NPZ.read_bytes()).hexdigest()[:16]


def _forecast_init() -> None:
    global _ATLAS, _ATLAS_SHA
    _limit_threads(max(1, (os.cpu_count() or 8) // FORECAST_WORKERS))
    _ATLAS, _ATLAS_SHA = load_atlas(), _atlas_sha()


def forecast_worker(pair: dict) -> dict:
    """One patient; checkpointed to cache/forecast_<pid>.json, keyed by the atlas hash."""
    cache = CACHE_DIR / f"forecast_{pair['patient_id']}.json"
    if cache.exists():
        rec = json.loads(cache.read_text())
        if rec.get("atlas_sha") == _ATLAS_SHA and rec.get("dt_days") == pair["dt_days"]:
            return rec
    t0 = time.time()
    rec = forecast_patient(pair, _ATLAS)
    rec["seconds"] = round(time.time() - t0, 1)
    rec["atlas_sha"] = _ATLAS_SHA
    tmp = cache.with_suffix(".tmp")
    tmp.write_text(json.dumps(rec))
    tmp.replace(cache)
    return rec


def stage_forecast(limit: int | None = None) -> None:
    import multiprocessing as mp
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    pairs, excluded = forecast_pairs()
    if limit:
        pairs = pairs[:limit]
    recs = {}
    with mp.get_context("spawn").Pool(FORECAST_WORKERS, initializer=_forecast_init) as pool:
        for k, rec in enumerate(pool.imap_unordered(forecast_worker, pairs)):
            recs[rec["patient_id"]] = rec
            best = max(rec["dice"]["grid"].values()) if "dice" in rec and "grid" in rec["dice"] else None
            print(f"[{k + 1}/{len(pairs)}] {rec['patient_id']} dt={rec['dt_days']:.0f}d "
                  f"no_change={rec.get('dice', {}).get('no_change')} best_grid={best} {rec['seconds']}s", flush=True)
    ordered = [recs[p["patient_id"]] for p in pairs]       # deterministic order for the CV folds
    GRID_JSON.write_text(json.dumps({"rhos": RHOS, "ds_mm2_per_day": DS, "sharpen": SHARPEN,
                                     "u_visible": U_VISIBLE, "grid_mm": H, "atlas_sha": _atlas_sha(),
                                     "excluded": excluded, "patients": ordered}, indent=1))


# ----------------------------------------------------------------------------- report stage
def boot_ci(x, n=5000, seed=0):
    x = np.asarray(x)
    rng = np.random.default_rng(seed)
    m = x[rng.integers(0, len(x), (n, len(x)))].mean(1)
    return [float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))]


def cv_select(recs, table: str):
    """Out-of-fold Dice per arm. Parameters picked on training folds by mean Dice."""
    rng = np.random.default_rng(51)
    fold = rng.permutation(len(recs)) % N_FOLDS
    keys = list(recs[0][table]["grid"].keys())
    G = np.array([[r[table]["grid"][k] for k in keys] for r in recs])
    arm_of = np.array([k.split("|")[0] for k in keys])
    out = {a: np.zeros(len(recs)) for a in ("anisotropic", "iso_same", "iso_homog")}
    chosen = {a: [] for a in out}
    for f in range(N_FOLDS):
        tr, te = fold != f, fold == f
        mean_tr = G[tr].mean(0)
        ia = int(np.argmax(np.where(arm_of == "aniso", mean_tr, -np.inf)))
        ih = int(np.argmax(np.where(arm_of == "iso_homog", mean_tr, -np.inf)))
        _, _, d_a, rho_a = keys[ia].split("|")
        is_ = keys.index(f"iso_same|r=1|{d_a}|{rho_a}")
        for arm, i in (("anisotropic", ia), ("iso_same", is_), ("iso_homog", ih)):
            out[arm][te] = G[te, i]
            chosen[arm].append(keys[i])
    return out, chosen


def paired(a, b):
    d = np.asarray(a) - np.asarray(b)
    p = float(stats.wilcoxon(a, b).pvalue) if np.any(d != 0) else 1.0
    return {"mean_diff": float(d.mean()), "mean_diff_ci95": boot_ci(d), "median_diff": float(np.median(d)),
            "frac_first_better": float(np.mean(d > 0)), "wilcoxon_p": p, "n": int(len(d))}


def holm(pvals: dict) -> dict:
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m, out, run = len(items), {}, 0.0
    for i, (k, p) in enumerate(items):
        run = max(run, min(1.0, (m - i) * p))
        out[k] = run
    return out


def atlas_dispersion_qc() -> dict:
    """How much anisotropy the population average lost: FA of the mean tensor vs the mean of
    each patient's own FA, over atlas tissue voxels. A large drop means fibre orientations
    disagree across registered patients (affine registration only)."""
    z = np.load(ATLAS_NPZ)
    tissue = (z["count"] >= int(z["min_count"])).ravel()
    fsum = np.zeros(tissue.size)
    n = np.zeros(tissue.size)
    for f in sorted(CACHE_DIR.glob("atlas_*.npz")):
        q = np.load(f)
        if not json.loads(str(q["qc"]))["included"]:
            continue
        fa, _ = fa_md(np.clip(np.linalg.eigvalsh(six_to_mat(q["d6"])), 0, None))
        fsum[q["idx"]] += fa
        n[q["idx"]] += 1
    mean_pat = (fsum / np.maximum(n, 1))[tissue]
    fa_atlas = z["fa"].ravel()[tissue]
    return {"fa_of_mean_tensor_p50_p90_p99": [float(x) for x in np.percentile(fa_atlas, [50, 90, 99])],
            "mean_of_patient_fa_p50_p90_p99": [float(x) for x in np.percentile(mean_pat, [50, 90, 99])],
            "p90_ratio_atlas_over_patient": float(np.percentile(fa_atlas, 90) / np.percentile(mean_pat, 90))}


def stage_report() -> None:
    g = json.loads(GRID_JSON.read_text())
    recs = [r for r in g["patients"] if "skip" not in r]
    skipped = [(r["patient_id"], r["skip"]) for r in g["patients"] if "skip" in r]
    oof, chosen = cv_select(recs, "dice")
    oof["no_change"] = np.array([r["dice"]["no_change"] for r in recs])
    arms = ["anisotropic", "iso_same", "iso_homog", "no_change"]
    grew = np.array([r["grew"] for r in recs])

    def summ(x):
        return {"mean": float(np.mean(x)), "mean_ci95": boot_ci(x), "median": float(np.median(x)),
                "iqr": [float(np.percentile(x, 25)), float(np.percentile(x, 75))], "n": int(len(x))}

    primary = {a: summ(oof[a]) for a in arms}
    tests = {"anisotropic_vs_no_change": paired(oof["anisotropic"], oof["no_change"]),
             "iso_homog_vs_no_change": paired(oof["iso_homog"], oof["no_change"]),
             "anisotropic_vs_iso_same": paired(oof["anisotropic"], oof["iso_same"]),
             "anisotropic_vs_iso_homog": paired(oof["anisotropic"], oof["iso_homog"])}
    adj = holm({k: v["wilcoxon_p"] for k, v in tests.items()})
    for k in tests:
        tests[k]["wilcoxon_p_holm"] = adj[k]
    sub = {name: {a: summ(oof[a][mask]) for a in arms}
           for name, mask in (("grew", grew), ("shrank_or_same", ~grew)) if mask.sum() > 0}
    sub_tests = {name: {"anisotropic_vs_no_change": paired(oof["anisotropic"][mask], oof["no_change"][mask]),
                        "anisotropic_vs_iso_same": paired(oof["anisotropic"][mask], oof["iso_same"][mask])}
                 for name, mask in (("grew", grew), ("shrank_or_same", ~grew)) if mask.sum() > 5}

    # secondary, shape only (uses the scan-2 volume): growing tumours
    rg = [r for r in recs if r["grew"]]
    vm, vm_chosen = cv_select(rg, "dice_volume_matched")
    vm["uniform_dilation"] = np.array([r["dice_volume_matched"]["uniform_dilation"] for r in rg])
    secondary = {"n": len(rg), "arms": {a: summ(v) for a, v in vm.items()},
                 "anisotropic_vs_iso_same": paired(vm["anisotropic"], vm["iso_same"]),
                 "anisotropic_vs_uniform_dilation": paired(vm["anisotropic"], vm["uniform_dilation"]),
                 "chosen_parameters_per_fold": vm_chosen}
    qc = json.loads(ATLAS_QC_JSON.read_text())["summary"]
    st = json.loads(SELFTEST_JSON.read_text()) if SELFTEST_JSON.exists() else None
    best_arm = max(arms, key=lambda a: primary[a]["mean"])
    metrics = {
        "design": "forecast scan 2 from scan-1 tumour mask (MU-Glioma-Post, whole tumour = mask > 0, 2 mm grid); "
                  "parameters chosen per arm by 5-fold patient-level CV; all Dice out-of-fold",
        "data": {"cohort_json": str(COHORT_JSON.relative_to(PROJECT_ROOT)), "n_scored": len(recs),
                 "n_grew": int(grew.sum()), "n_shrank_or_same": int((~grew).sum()),
                 "excluded_pairs": g["excluded"], "skipped_patients": skipped,
                 "interval_days_median": float(np.median([r["dt_days"] for r in recs])),
                 "interval_days_range": [float(min(r["dt_days"] for r in recs)),
                                         float(max(r["dt_days"] for r in recs))]},
        "parameter_grid": {"rho_per_day": g["rhos"], "d_mm2_per_day": g["ds_mm2_per_day"],
                           "anisotropy_sharpening_r": g["sharpen"], "u_visible": g["u_visible"]},
        "primary_dice_out_of_fold": primary,
        "primary_tests": tests,
        "best_arm_by_mean_dice": best_arm,
        "subgroups": sub, "subgroup_tests": sub_tests,
        "chosen_parameters_per_fold": chosen,
        "secondary_volume_matched_growing_only": secondary,
        "dti_atlas_qc": qc,
        "dti_atlas_orientation_dispersion": atlas_dispersion_qc(),
        "solver_selftest": st,
        "limitations": [
            "MU-Glioma-Post has no DTI; the anisotropic arm uses a population DTI atlas built from 62 "
            "UCSF-PDGM patients registered to the same SRI24 space, not each patient's own tensors, and "
            "ignores mass-effect deformation.",
            "Averaging affinely registered tensors loses anisotropy (see dti_atlas_orientation_dispersion): "
            "fibre directions do not line up exactly across patients. The sharpening factor r = 10 was offered "
            "to compensate; cross-validation did not select it.",
            "The growth model has no treatment term; scan pairs are post-operative and under therapy.",
            "Whole-tumour masks (all labels) are forecast; sub-regions are not scored separately.",
            "The volume-matched secondary analysis uses the observed scan-2 volume and is a shape test, "
            "not a forecast.",
        ],
        "supersedes": ["output/cv_test_results.csv", "output/final_results_20260901_2213.csv",
                       "old run_improved_aniso.py (compared the simulation with its own initial mask)"],
    }
    RESULTS_JSON.write_text(json.dumps(metrics, indent=2))
    with RESULTS_TSV.open("w") as fh:
        cols = ["patient_id", "dt_days", "v1_voxels_2mm", "v2_voxels_2mm", "grew"] + [f"dice_{a}" for a in arms] + ["fold"]
        fh.write(TAB.join(cols) + NL)
        fold = np.random.default_rng(51).permutation(len(recs)) % N_FOLDS
        for i, r in enumerate(recs):
            row = [r["patient_id"], f"{r['dt_days']:g}", str(r["v1_voxels"]), str(r["v2_voxels"]), str(int(r["grew"]))]
            row += [f"{oof[a][i]:.4f}" for a in arms] + [str(fold[i])]
            fh.write(TAB.join(row) + NL)
    make_figure(oof, grew, vm)
    print(json.dumps({"primary": primary, "tests": tests}, indent=2))


def make_figure(oof, grew, vm) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    arms = ["anisotropic", "iso_same", "iso_homog", "no_change"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    ax[0].boxplot([oof[a] for a in arms], showmeans=True)
    ax[0].set_xticks(range(1, 5), arms)
    ax[0].set_ylabel("Dice vs observed scan 2 (out-of-fold)")
    ax[0].set_title(f"Forecast, all patients (n={len(grew)})")
    d = oof["anisotropic"] - oof["no_change"]
    d2 = oof["anisotropic"] - oof["iso_same"]
    ax[1].scatter(np.where(grew, 1, 0) + np.random.default_rng(0).uniform(-.15, .15, len(d)), d, s=10,
                  label="aniso - no_change")
    ax[1].scatter(np.where(grew, 3, 2) + np.random.default_rng(1).uniform(-.15, .15, len(d2)), d2, s=10,
                  label="aniso - iso_same")
    ax[1].axhline(0, color="k", lw=0.8)
    ax[1].set_xticks([0, 1, 2, 3], ["shrank", "grew", "shrank", "grew"])
    ax[1].set_ylabel("paired Dice difference")
    ax[1].legend(fontsize=8)
    ax[1].set_title("Paired differences by subgroup")
    va = ["anisotropic", "iso_same", "iso_homog", "uniform_dilation"]
    ax[2].boxplot([vm[a] for a in va], showmeans=True)
    ax[2].set_xticks(range(1, 5), ["aniso", "iso_same", "iso_homog", "dilation"])
    ax[2].set_title(f"Shape only, volume matched, growing (n={len(vm['anisotropic'])})")
    fig.tight_layout()
    fig.savefig(FIGURE, dpi=150)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["selftest", "atlas", "forecast", "report", "all"], default="all")
    ap.add_argument("--limit", type=int, default=None, help="process only the first N patients (debug)")
    a = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if a.stage in ("selftest", "all"):
        st = selftest()
        SELFTEST_JSON.write_text(json.dumps(st, indent=2))
        print(json.dumps(st, indent=2))
        if not st["all_pass"]:
            raise SystemExit("solver self-test failed")
    if a.stage in ("atlas", "all"):
        stage_atlas(a.limit)
    if a.stage in ("forecast", "all"):
        stage_forecast(a.limit)
    if a.stage in ("report", "all"):
        stage_report()


if __name__ == "__main__":
    main()
