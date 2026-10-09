"""Monotone (positivity-preserving) full-tensor Fisher-Kolmogorov solver (GRAND_PLAN #15; weak point W21).

Why: run_improved_aniso.TensorFK uses central differences for the off-diagonal tensor terms. That stencil can make u
negative under strong anisotropy, and the step then clamps u >= 0, which adds mass (2.8% in the script 97 stress case,
up to 42% in one aniso r = 10 production run, ledger 2q / 2w).

Method: non-negative directional splitting. At each voxel the tensor is written as
    D(x) = sum_k w_k(x) e_k e_k^T,   w_k >= 0,
over integer lattice offsets e_k. Then
    div(D grad u)(x) ~ sum_k [ c_k(x + e_k/2) (u(x + e_k) - u(x)) + c_k(x - e_k/2) (u(x - e_k) - u(x)) ] / h^2,
with c_k at the half point = mean of w_k at the two ends, set to 0 if either end is outside the domain (zero flux, so
mass is conserved exactly when rho = 0). All coefficients are >= 0, so with
    dt <= 0.9 h^2 / max_x sum_k (c_k+ + c_k-)   (enforced through lam_max for the inherited run())
an explicit Euler step is a convex combination plus the reaction and kill terms; with the existing limits
dt * rho <= 0.05 and dt * kill <= 0.05, u stays in [0, 1] without any clamp. The clamp in run() is kept only as a
bookkeeping check: it should record zero mass.

Two ways to get the weights (argument `split`):
  "selling" (default)  Selling's decomposition (Selling 1874; Conway & Sloane 1992, Proc R Soc A 436:55; used for
                       diffusion by Fehrenbach & Mirebeau 2014, J Math Imaging Vis 49:123, arXiv 1301.3925).
                       Selling's algorithm finds a D-obtuse superbase b0..b3 of Z^3 (b0+b1+b2+b3 = 0 and
                       <b_i, D b_j> <= 0 for i != j); then D = sum_{i<j} -<b_i, D b_j> e_ij e_ij^T with
                       e_ij = b_k x b_l ({i,j,k,l} = {0,1,2,3}). Six weights, all >= 0, exact for every SPD tensor.
                       Offsets get longer as anisotropy grows (cost O(log condition number)). A diagonal D gives the
                       three axes only, i.e. the 7-point operator of TensorFK.
  "nnls"               earlier version: NNLS on a fixed set of up to 37 offsets, smallest stencil first. Not exact
                       where D is not a non-negative combination of those offsets (kept for comparison).
The relative residual |D - sum w e e^T|_F / |D|_F per voxel is stored (self.rel_residual, self.residual_summary).
"""
from __future__ import annotations

import itertools
import math

import numpy as np
import torch
from scipy.optimize import nnls

import run_improved_aniso as ria


def _canonical(v):
    for c in v:
        if c != 0:
            return tuple(v) if c > 0 else tuple(-x for x in v)
    return None


def lattice_directions(wide: bool):
    """13 nearest directions (max |component| 1); with wide=True also the 24 primitive (2,1,0) and (2,1,1) types."""
    out = set()
    rng = range(-2, 3) if wide else range(-1, 2)
    for v in itertools.product(rng, repeat=3):
        if v == (0, 0, 0) or math.gcd(math.gcd(abs(v[0]), abs(v[1])), abs(v[2])) != 1:
            continue
        s = sorted(abs(x) for x in v)
        if max(s) == 2 and s not in ([0, 1, 2], [1, 1, 2]):
            continue
        out.add(_canonical(v))
    return sorted(out, key=lambda e: (max(map(abs, e)), sum(map(abs, e)), e))


def _design(dirs):
    # rows: xx, yy, zz, sqrt2*xy, sqrt2*xz, sqrt2*yz  (Frobenius weighting of the symmetric off-diagonals)
    r2 = math.sqrt(2.0)
    cols = [[e[0] * e[0], e[1] * e[1], e[2] * e[2], r2 * e[0] * e[1], r2 * e[0] * e[2], r2 * e[1] * e[2]] for e in dirs]
    return np.asarray(cols, dtype=np.float64).T


def _shift(a, e):
    """b[x] = a[x + e] on the last three dims, zero where x + e is outside the array."""
    out = torch.zeros_like(a)
    src, dst = [], []
    for ax, s in enumerate(e):
        n = a.shape[a.dim() - 3 + ax]
        if s >= 0:
            src.append(slice(s, n)); dst.append(slice(0, n - s))
        else:
            src.append(slice(0, n + s)); dst.append(slice(-s, n))
    lead = (slice(None),) * (a.dim() - 3)
    out[lead + tuple(dst)] = a[lead + tuple(src)]
    return out


_PAIRS = [(i, j, *[x for x in range(4) if x not in (i, j)]) for i in range(4) for j in range(i + 1, 4)]


def selling_superbase(D: np.ndarray, max_iter: int = 500) -> np.ndarray:
    """D: (N, 3, 3) SPD. Returns (N, 4, 3) integer superbases that are D-obtuse (Selling's algorithm, vectorised)."""
    b = np.zeros((D.shape[0], 4, 3), np.float64)
    b[:, 0, 0] = b[:, 1, 1] = b[:, 2, 2] = 1.0
    b[:, 3] = -1.0
    tol = 1e-12 * np.trace(D, axis1=1, axis2=2)
    for _ in range(max_iter):
        changed = False
        for i, j, k, l in _PAIRS:
            sel = np.einsum("ni,nij,nj->n", b[:, i], D, b[:, j]) > tol
            if sel.any():
                bi = b[sel, i].copy()
                b[sel, i] = -bi
                b[sel, k] += bi
                b[sel, l] += bi
                changed = True
        if not changed:
            return b
    raise RuntimeError("Selling's algorithm did not converge")


def selling_decompose(D: np.ndarray):
    """D: (N, 3, 3) SPD -> weights (N, 6) >= 0 and integer offsets (N, 6, 3) with D = sum_k w_k e_k e_k^T."""
    b = selling_superbase(D)
    w = [-np.einsum("ni,nij,nj->n", b[:, i], D, b[:, j]) for i, j, k, l in _PAIRS]
    off = [np.cross(b[:, k], b[:, l]) for i, j, k, l in _PAIRS]
    return np.clip(np.stack(w, 1), 0.0, None), np.rint(np.stack(off, 1)).astype(np.int64)


class TensorFKMonotone(ria.TensorFK):
    """Drop-in replacement for ria.TensorFK (same constructor and run()); only the diffusion operator differs."""

    def __init__(self, d6: np.ndarray, mask: np.ndarray, h: float = ria.H, tol: float = 1e-6,
                 split: str = "selling"):
        self.h = h
        self.split = split
        m_np = mask.astype(bool)
        self.m = torch.as_tensor(m_np.astype(np.float32))
        d6 = np.asarray(d6, np.float64) * m_np[..., None]
        r2 = math.sqrt(2.0)
        b_all = np.concatenate([d6[..., :3], r2 * d6[..., 3:]], axis=-1)
        if split == "selling":
            W, rel, n_wide = self._split_selling(d6, m_np)
        elif split == "nnls":
            W, rel, n_wide = self._split_nnls(m_np, b_all, tol)
        else:
            raise ValueError(split)
        self.rel_residual = rel
        inside = rel[m_np & (np.linalg.norm(b_all, axis=-1) > 0)]
        self.residual_summary = {
            "split": split, "n_offsets": len(self.dirs),
            "max_offset_length": float(max(np.linalg.norm(e) for e in self.dirs)) if self.dirs else 0.0,
            "n_voxels": int(inside.size), "n_voxels_wide_stencil": int(n_wide),
            "rel_residual_median": float(np.median(inside)) if inside.size else 0.0,
            "rel_residual_p95": float(np.quantile(inside, 0.95)) if inside.size else 0.0,
            "rel_residual_max": float(inside.max()) if inside.size else 0.0,
            "share_voxels_exact": float((inside <= tol).mean()) if inside.size else 1.0}
        self._build_coefficients(W)

    def _split_selling(self, d6, m_np):
        idx = np.nonzero(m_np & (np.abs(d6).sum(-1) > 0))
        s = d6[idx]
        D = np.stack([np.stack([s[:, 0], s[:, 3], s[:, 4]], -1),
                      np.stack([s[:, 3], s[:, 1], s[:, 5]], -1),
                      np.stack([s[:, 4], s[:, 5], s[:, 2]], -1)], 1)
        w, off = selling_decompose(D)
        # e and -e give the same term: make the first non-zero component positive
        lead = np.take_along_axis(off, (off != 0).argmax(-1)[..., None], -1)[..., 0]
        off = off * np.where(lead < 0, -1, 1)[..., None]
        keep = w > 0
        uniq, inv = np.unique(off[keep], axis=0, return_inverse=True)
        self.dirs = [tuple(int(x) for x in e) for e in uniq]
        W = np.zeros(m_np.shape + (len(self.dirs),), np.float64)
        vox = np.broadcast_to(np.arange(len(s))[:, None], keep.shape)[keep]
        np.add.at(W, (idx[0][vox], idx[1][vox], idx[2][vox], inv.ravel()), w[keep])
        of = off.astype(np.float64)
        rec = np.einsum("nk,nki,nkj->nij", w, of, of)
        rel = np.zeros(m_np.shape, np.float64)
        rel[idx] = np.linalg.norm((rec - D).reshape(len(s), -1), axis=1) / np.linalg.norm(D.reshape(len(s), -1), axis=1)
        n_wide = int(((np.abs(off) * keep[..., None]).max(axis=(1, 2)) >= 2).sum())
        return W, rel, n_wide

    def _split_nnls(self, m_np, b_all, tol):
        near, wide = lattice_directions(False), lattice_directions(True)
        wide_only = [e for e in wide if e not in near]
        self.dirs = near + wide_only
        # smallest stencil first (truncation error grows with |e|): axes -> + face diagonals -> 13 -> 37.
        # A diagonal D is then split on the axes only, i.e. the 7-point operator of TensorFK.
        levels = []
        for keep in (lambda e: sum(map(abs, e)) == 1, lambda e: max(map(abs, e)) == 1 and sum(map(abs, e)) <= 2,
                     lambda e: max(map(abs, e)) == 1, lambda e: True):
            cols = [k for k, e in enumerate(self.dirs) if keep(e)]
            levels.append((np.asarray(cols), _design([self.dirs[k] for k in cols])))
        W = np.zeros(m_np.shape + (len(self.dirs),), np.float64)
        rel = np.zeros(m_np.shape, np.float64)
        n_wide = 0
        cache = {}
        for idx in zip(*np.nonzero(m_np)):
            b = b_all[idx]
            nb = float(np.linalg.norm(b))
            if nb == 0.0:
                continue
            key = tuple(np.round(b, 12))
            if key not in cache:
                best = None
                for li, (cols, A) in enumerate(levels):
                    w, res = nnls(A, b)
                    if best is None or res < best[1] - 1e-15:
                        best = (cols, res, w, li)
                    if res / nb <= tol:
                        break
                cols, res, w, li = best
                full = np.zeros(len(self.dirs))
                full[cols] = w
                cache[key] = (full, res / nb, li == 3)
            w, r, used_wide = cache[key]
            W[idx] = w
            rel[idx] = r
            n_wide += int(used_wide)
        return W, rel, n_wide

    def _build_coefficients(self, W):
        Wt = torch.as_tensor(W.astype(np.float32))
        self.cp, self.cm = [], []
        diag = torch.zeros_like(self.m)
        for k, e in enumerate(self.dirs):
            wk = Wt[..., k]
            if float(wk.abs().max()) == 0.0:
                continue
            # bond x -> x + e is open only if both ends and the voxels the segment passes through are in the domain
            # (long Selling offsets must not jump across a sulcus or ventricle)
            path = self.m * _shift(self.m, e)
            n = max(map(abs, e))
            for t in range(1, n):
                path = path * _shift(self.m, tuple(int(np.floor(x * t / n + 0.5)) for x in e))
            cplus = 0.5 * (wk + _shift(wk, e)) * path
            cminus = _shift(cplus, tuple(-x for x in e))
            self.cp.append((e, cplus))
            self.cm.append((tuple(-x for x in e), cminus))
            diag = diag + cplus + cminus
        # inherited run(): dt_max = 0.9 h^2 / (6 lam_max)  ->  dt <= 0.9 h^2 / max(diag)
        self.lam_max = float(diag.max().clamp_min(1e-12)) / 6.0

    def div_flux(self, u):
        out = torch.zeros_like(u)
        for (e, c), (em, cm) in zip(self.cp, self.cm):
            out = out + c * (_shift(u, e) - u) + cm * (_shift(u, em) - u)
        return out / (self.h * self.h)
