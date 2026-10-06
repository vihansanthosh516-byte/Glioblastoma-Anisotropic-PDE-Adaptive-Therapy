#!/usr/bin/env python3
"""Script 101: why does DTI not help? Direction-alignment test (Phase 2.4; masterplan §86-88; plan A1.9, exploratory).

Hypotheses (masterplan §86), what this script can and cannot test:
  A  Macroscopic MRI tumour change is not strongly constrained by tract orientation.   -> TESTED here (alignment vs null)
  B  DTI resolution / tensor quality is insufficient.                                   -> atlas QC reported only
  C  The PDE already captures enough anisotropy without DTI.                             -> read from script 100
                                                                                         (iso_homog vs aniso contrasts)
  D  The DTI mapping (a population atlas, not the patient's own tensor) is too crude.   -> NOT testable with MU-Glioma
                                                                                         (no native DTI); stated as a limit

Test for A. For each primary pair, growth direction g = centroid(T \\ S) - centroid(S) (T next mask, S input mask, 2 mm,
core labels 1 and 3). Local fibre direction e = principal eigenvector of the mean atlas tensor over S dilated by 3 voxels
(6 mm). Alignment = |cos(g, e)|. For random directions in 3D, E|cos| = 0.5.
Null: pair each pair's growth direction with the tensor of ANOTHER pair (10,000 permutations): keeps the distribution of
growth directions and of tensor directions, breaks their link.
Pairs are used only if the new-growth region has >= 30 voxels and |g| >= 1 voxel (otherwise the direction is undefined).
Unit of analysis: patient (mean |cos| over the patient's pairs); bootstrap over patients.

Input : data/manifests/forecast_pairs_mu.csv, MU core masks, output/forecast_validation/dti_atlas_2mm.npz
Output: output/dti_alignment.json, output/dti_alignment_pairs.csv
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path as _Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
MAN = PROJECT_ROOT / "data" / "manifests"
sys.path.insert(0, str(PROJECT_ROOT))

import run_improved_aniso as ria  # noqa: E402
from src import patient_stats as ps  # noqa: E402
from src.run_manifest import write_run_manifest  # noqa: E402

MIN_GROWTH_VOX = 30
N_PERM = 10_000
SEED = ps.DEFAULT_SEED


def core(pid, tp):
    raw = np.isin(np.asarray(nib.load(str(ria.mu_paths(pid, tp)[0])).dataobj), (1, 3))
    return ria.to_2mm(raw) >= 0.5


def fibre_direction(atlas, S):
    """Principal eigenvector of the mean atlas tensor over S dilated by 3 voxels (tissue only); None if no tissue."""
    region = ndimage.binary_dilation(S, iterations=3) & atlas["tissue"]
    if region.sum() < 5:
        return None
    w, v = atlas["w"][region], atlas["v"][region]                  # (n, 3), (n, 3, 3)
    M = np.einsum("nij,nj,nkj->ik", v, w, v) / region.sum()
    ev, evec = np.linalg.eigh(M)
    return evec[:, -1], float(ev[-1] / max(ev.sum(), 1e-12))


def main():
    t0 = time.time()
    atlas = ria.load_atlas()
    pairs = pd.read_csv(MAN / "forecast_pairs_mu.csv")
    pairs = pairs[pairs["in_primary"]]
    rows = []
    for i, r in enumerate(pairs.itertuples(index=False)):
        S, T = core(r.patient_id, r.tp_in), core(r.patient_id, r.tp_out)
        G = T & ~S
        if S.sum() == 0 or G.sum() < MIN_GROWTH_VOX:
            rows.append({"patient_id": r.patient_id, "tp_in": r.tp_in, "tp_out": r.tp_out, "used": False,
                         "reason": "no seed" if S.sum() == 0 else f"growth region < {MIN_GROWTH_VOX} voxels"})
            continue
        g = np.argwhere(G).mean(0) - np.argwhere(S).mean(0)
        fd = fibre_direction(atlas, S)
        if np.linalg.norm(g) < 1.0 or fd is None:
            rows.append({"patient_id": r.patient_id, "tp_in": r.tp_in, "tp_out": r.tp_out, "used": False,
                         "reason": "growth direction undefined (centroid shift < 1 voxel)" if fd else "no tissue"})
            continue
        e, frac = fd
        rows.append({"patient_id": r.patient_id, "tp_in": r.tp_in, "tp_out": r.tp_out, "used": True,
                     "gx": g[0], "gy": g[1], "gz": g[2], "ex": e[0], "ey": e[1], "ez": e[2], "shift_vox": float(np.linalg.norm(g)),
                     "tensor_anisotropy_frac": frac, "abs_cos": float(abs(g @ e) / np.linalg.norm(g))})
        if i % 50 == 0:
            print(f"[{time.time() - t0:.0f}s] {i}/{len(pairs)}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "dti_alignment_pairs.csv", index=False)
    u = df[df["used"]].copy()
    pm = u.groupby("patient_id")["abs_cos"].mean()
    obs = float(pm.mean())
    # permutation null at the pair level: shuffle tensors across pairs, then average to patients as for the observed
    G = u[["gx", "gy", "gz"]].to_numpy()
    G = G / np.linalg.norm(G, axis=1, keepdims=True)
    E = u[["ex", "ey", "ez"]].to_numpy()
    rng = np.random.default_rng(SEED)
    pid_codes = pd.factorize(u["patient_id"])[0]
    cnt = np.bincount(pid_codes)
    null = np.empty(N_PERM)
    for k in range(N_PERM):
        c = np.abs((G * E[rng.permutation(len(E))]).sum(1))
        null[k] = (np.bincount(pid_codes, weights=c) / cnt).mean()
    lo, hi = ps.bootstrap_ci(pm.to_numpy(), n_boot=10_000)
    p = float((1 + (np.abs(null - null.mean()) >= abs(obs - null.mean())).sum()) / (N_PERM + 1))
    try:
        qc = json.loads((OUT / "forecast_validation" / "dti_atlas_qc.json").read_text())
    except Exception:
        qc = None
    res = {"script": "101_dti_alignment", "n_pairs_total": int(len(df)), "n_pairs_used": int(len(u)),
           "n_patients_used": int(pm.shape[0]), "excluded_reasons": df[~df["used"]]["reason"].value_counts().to_dict(),
           "mean_abs_cos_patient_level": obs, "ci95_patient_bootstrap": [lo, hi],
           "random_direction_expectation": 0.5, "permutation_null_mean": float(null.mean()),
           "permutation_null_ci95": [float(np.percentile(null, 2.5)), float(np.percentile(null, 97.5))],
           "permutation_p_two_sided": p,
           "median_centroid_shift_vox": float(u["shift_vox"].median()),
           "median_tensor_anisotropy_frac": float(u["tensor_anisotropy_frac"].median()),
           "hypothesis_A_reading": ("alignment indistinguishable from the permutation null supports A (tract orientation does not "
                                    "constrain macroscopic change); alignment clearly above null means orientation carries signal "
                                    "the PDE could use"),
           "hypothesis_B_atlas_qc": qc, "hypothesis_C": "see output/pde_manifest/results.json contrasts iso_homog vs aniso",
           "hypothesis_D": "not testable: MU-Glioma has no native DTI; the atlas is a population average",
           "caveats": ["centroid shift is a coarse growth direction; isotropic growth gives an undefined direction",
                       "tensor is the mean over a 6 mm shell around the input tumour of a population atlas, not the patient's own",
                       "2 mm grid; one voxel is large for small cores"]}
    (OUT / "dti_alignment.json").write_text(json.dumps(res, indent=1, default=float))
    write_run_manifest("dti_alignment_101", OUT / "dti_alignment.manifest.json", script="src/101_dti_alignment.py", seed=SEED,
                       config={"min_growth_vox": MIN_GROWTH_VOX, "n_perm": N_PERM}, inputs=[MAN / "forecast_pairs_mu.csv", ria.ATLAS_NPZ],
                       dataset="MU-Glioma-Post core", patient_split="n/a", primary_endpoint="n/a (exploratory)")
    print(json.dumps({k: v for k, v in res.items() if k not in ("hypothesis_B_atlas_qc",)}, indent=1, default=float))


if __name__ == "__main__":
    main()
