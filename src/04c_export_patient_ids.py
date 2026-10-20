"""
04c_export_patient_ids.py
=========================
Writes one patient ID per row of output/nn_X.npy (same row order as script 04).

Why: script 04 writes the patient-disjoint TEST split but not the patient IDs.
Scripts 05/07/08 need patient IDs to carve a patient-grouped VALIDATION split
from train in GBM_SPLIT_MODE=patient.

Additive only: does not change nn_X.npy, nn_y.npy or any split file.

Output:
    output/nn_patient_ids.npy   -- object array of patient ID strings (N_cells,)

Run:
    python src/04c_export_patient_ids.py
"""
import os
import numpy as np
import scanpy as sc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # project root
OUT_DIR = os.path.join(ROOT, "output")
H5_SUB = os.path.join(OUT_DIR, "02_adata_subsampled.h5ad")
NPY_X = os.path.join(OUT_DIR, "nn_X.npy")
NPY_PID = os.path.join(OUT_DIR, "nn_patient_ids.npy")

adata = sc.read_h5ad(H5_SUB, backed="r")  # backed: loads obs only, not the matrix
obs = adata.obs
assert "ID1" in obs.columns, "obs has no ID1 column"

pid = obs["ID1"].astype(str).to_numpy()
n_x = np.load(NPY_X, mmap_mode="r").shape[0]
assert len(pid) == n_x, f"row mismatch: obs={len(pid)} nn_X={n_x}"

np.save(NPY_PID, pid.astype(object), allow_pickle=True)
print(f"wrote {NPY_PID}: {len(pid)} cells, {len(np.unique(pid))} patients")
