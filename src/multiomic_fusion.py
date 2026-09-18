#!/usr/bin/env python3
"""
Multi-Omic & Epigenetic Feature Fusion (Proposal 2)
====================================================
Fuses transcriptomic (Neftel-state) fractions with DNA methylation,
copy-number alterations (CNV), and metabolomic flux features to predict
patient-specific proliferation (rho) and diffusivity (D) via cross-validated
elastic-net regression.

Feature vector layout (per patient)
-----------------------------------
    X = [ f_NPC, f_OPC, f_AC, f_MES,            # 4 Neftel-state fractions
          methylation_beta (M1 features),        e.g. MGMT promoter
          cnv_log2_ratio (M2 features),          e.g. EGFR / PDGFRA
          metabolic_flux (M3 features) ]         e.g. glycolysis / OXPHOS

Two independent ElasticNetCV models predict rho and D. The trained models are
pickled to output/multiomic_elasticnet.pkl; `predict_params_from_features`
loads the pickle and returns rho_mean, D_mean.

Storage
-------
A single tab-separated file `output/multiomic_features.tsv` holds the
per-patient feature matrix plus the ground-truth rho/D labels (when known, e.g.
from the inverse-estimation module 51). `build_multiomic_tsv(...)` synthesizes
this file when real omics are unavailable, enabling end-to-end CI on synthetic
cohorts.

Integration point: src/50_spatial_genomics_deconv.py map_fractions_to_parameters
is augmented with a `multiomic_predict` flag that delegates rho/D prediction
to this module.
"""
from __future__ import annotations

import pickle
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    import pandas as pd
    _HAVE_PANDAS = True
except ImportError:
    _HAVE_PANDAS = False
    pd = None  # type: ignore

# Forward-compat: sklearn is a declared dependency in Requirements.txt
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import KFold

# Physical ranges (kept in lockstep with src/50_spatial_genomics_deconv.py)
RHO_MIN, RHO_MAX = 0.005, 0.12
D_MIN, D_MAX = 0.01, 0.50

NEFTEL_STATES = ["NPC-like", "OPC-like", "AC-like", "MES-like"]

# Defaults for the synthetic multiomic schema
METHYLATION_FEATURES = ["MGMT_promoter_beta", "MGMT_body_beta",
                        "TERT_promoter_beta", "p16INK4a_beta"]
CNV_FEATURES = ["EGFR_log2ratio", "PDGFRA_log2ratio",
                "CDKN2A_log2ratio", "PTEN_log2ratio"]
METABOLIC_FEATURES = ["glycolysis_flux", "OXPHOS_flux",
                      "glutaminolysis_flux", "PPP_flux"]

MODEL_PATH = Path("output/multiomic_elasticnet.pkl")
FEATURES_TSV_PATH = Path("output/multiomic_features.tsv")


# --------------------------------------------------------------------------- #
# Feature schema / synthesis
# --------------------------------------------------------------------------- #
def feature_columns() -> List[str]:
    """Ordered column list for the feature matrix X (excluding patient_id)."""
    return (NEFTEL_STATES + METHYLATION_FEATURES
            + CNV_FEATURES + METABOLIC_FEATURES)


def feature_dim() -> int:
    """Dimension of the multi-omic feature vector (excludes patient id/labels)."""
    return len(feature_columns())


def build_multiomic_tsv(
    patient_ids: List[str],
    neftel_fractions: np.ndarray,
    rho_targets: Optional[np.ndarray] = None,
    D_targets: Optional[np.ndarray] = None,
    out_path: Path = FEATURES_TSV_PATH,
    seed: int = 42,
) -> Path:
    """Synthesize a multiomic_features.tsv (drops real omics where absent).

    Parameters
    ----------
    patient_ids : list of patient identifiers
    neftel_fractions : (N, 4) array of Neftel fractions
    rho_targets / D_targets : (N,) ground-truth labels (optional). When None,
        synthetic labels are generated from fractions + omic noise so the
        elastic-net training step has something to learn.
    out_path : output TSV path
    """
    if not _HAVE_PANDAS:
        raise ImportError("pandas is required to build multiomic_features.tsv")
    rng = np.random.default_rng(seed)
    n = len(patient_ids)
    if neftel_fractions.shape != (n, 4):
        raise ValueError(
            f"neftel_fractions must be shape ({n}, 4), got {neftel_fractions.shape}"
        )

    npc, opc, ac, mes = neftel_fractions.T
    # Methylation betas correlate with mesenchymal (MES) -> higher MGMT silencing
    meth = np.column_stack([
        np.clip(0.3 + 0.4 * mes + 0.1 * rng.standard_normal(n), 0, 1),
        np.clip(0.5 + 0.1 * rng.standard_normal(n), 0, 1),
        np.clip(0.2 + 0.3 * opc + 0.05 * rng.standard_normal(n), 0, 1),
        np.clip(0.6 - 0.2 * ac + 0.05 * rng.standard_normal(n), 0, 1),
    ])
    # CNV: EGFR amplified in proliferative NPC/OPC tumors
    cnv = np.column_stack([
        1.2 * (npc + opc) - 0.8 + 0.2 * rng.standard_normal(n),
        0.4 * mes - 0.2 + 0.15 * rng.standard_normal(n),
        -0.5 * opc + 0.1 * rng.standard_normal(n),
        -0.7 * npc - 0.3 * rng.standard_normal(n),
    ])
    # Metabolic: glycolysis high in proliferative core, OXPHOS shuffled
    metab = np.column_stack([
        2.0 * (npc + opc) + 0.3 * rng.standard_normal(n),
        1.0 + 0.4 * rng.standard_normal(n),
        0.5 * mes + 0.5 + 0.2 * rng.standard_normal(n),
        0.8 * ac + 0.2 * rng.standard_normal(n),
    ])

    X = np.column_stack([neftel_fractions, meth, cnv, metab])

    if rho_targets is None:
        rho_targets = np.clip(
            RHO_MIN + (RHO_MAX - RHO_MIN) * (npc + opc)
            + 0.01 * meth[:, 0]
            + 0.005 * (cnv[:, 0] + 1.0),
            RHO_MIN, RHO_MAX,
        )
    if D_targets is None:
        D_targets = np.clip(
            D_MIN + (D_MAX - D_MIN) * (ac + mes)
            + 0.01 * (1.0 - meth[:, 0])
            + 0.02 * metab[:, 0],
            D_MIN, D_MAX,
        )

    df = pd.DataFrame(X, columns=feature_columns())
    df.insert(0, "patient_id", patient_ids)
    df["rho"] = rho_targets
    df["D"] = D_targets

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"[multiomic] wrote {n} patient rows x {feature_dim()} features -> {out_path}")
    return out_path


def load_multiomic_tsv(path: Path = FEATURES_TSV_PATH) -> Tuple[np.ndarray, np.ndarray]:
    """Load multiomic_features.tsv into (X, y) for training.

    Returns
    -------
    X : (N, F) feature matrix (no patient id, no labels)
    y : (N, 2) matrix of [rho, D] labels
    """
    if not _HAVE_PANDAS:
        raise ImportError("pandas is required to load multiomic_features.tsv")
    df = pd.read_csv(path, sep="\t")
    cols = feature_columns()
    X = df[cols].to_numpy(dtype=np.float64)
    y = df[["rho", "D"]].to_numpy(dtype=np.float64)
    return X, y


def get_patient_vector(patient_id: str, path: Path = FEATURES_TSV_PATH) -> np.ndarray:
    """Load one patient's multi-omic vector X (length = feature_dim())."""
    if not _HAVE_PANDAS:
        raise ImportError("pandas is required to read multi-omic features")
    df = pd.read_csv(path, sep="\t")
    row = df[df["patient_id"] == patient_id]
    if len(row) == 0:
        raise KeyError(f"patient {patient_id} not in {path}")
    return row[feature_columns()].to_numpy(dtype=np.float64).squeeze()


# --------------------------------------------------------------------------- #
# Elastic-Net training & inference
# --------------------------------------------------------------------------- #
def train_elasticnet(
    X: np.ndarray,
    y: np.ndarray,
    n_splits: int = 5,
    n_jobs: int = 1,
    random_state: int = 42,
    model_path: Path = MODEL_PATH,
) -> Dict:
    """Train two ElasticNetCV models (rho, D) with cross-validated alpha & l1_ratio."""
    # Clip negative rho/D from noise-floor minus regularization
    y_clipped = np.clip(y, [RHO_MIN, D_MIN], [RHO_MAX, D_MAX])
    rho_y, D_y = y_clipped[:, 0], y_clipped[:, 1]

    kf = KFold(n_splits=min(n_splits, len(y)), shuffle=True,
               random_state=random_state)

    print("[multiomic] training elastic-net (rho)...")
    rho_model = ElasticNetCV(
        l1_ratio=[0.1, 0.5, 0.7, 0.9, 1.0],
        eps=1e-3, cv=kf, n_jobs=n_jobs,
        random_state=random_state,
    )
    rho_model.fit(X, rho_y)
    print(f"  rho  alpha={rho_model.alpha_:.5f}  l1_ratio={rho_model.l1_ratio_:.2f}  R2={rho_model.score(X, rho_y):.3f}")

    print("[multiomic] training elastic-net (D)...")
    D_model = ElasticNetCV(
        l1_ratio=[0.1, 0.5, 0.7, 0.9, 1.0],
        eps=1e-3, cv=kf, n_jobs=n_jobs,
        random_state=random_state,
    )
    D_model.fit(X, D_y)
    print(f"  D    alpha={D_model.alpha_:.5f}  l1_ratio={D_model.l1_ratio_:.2f}  R2={D_model.score(X, D_y):.3f}")

    bundle = {
        "rho_model": rho_model,
        "D_model": D_model,
        "feature_columns": feature_columns(),
        "rho_range": [RHO_MIN, RHO_MAX],
        "D_range": [D_MIN, D_MAX],
    }
    model_path.parent.mkdir(parents=True, exist_ok=True)
    with open(model_path, "wb") as f:
        pickle.dump(bundle, f)
    print(f"[multiomic] saved trained models -> {model_path}")
    return {
        "rho_r2": float(rho_model.score(X, rho_y)),
        "D_r2": float(D_model.score(X, D_y)),
        "rho_alpha": float(rho_model.alpha_),
        "D_alpha": float(D_model.alpha_),
        "model_path": str(model_path),
    }


def predict_params_from_features(
    X: np.ndarray,
    model_path: Path = MODEL_PATH,
) -> Tuple[np.ndarray, np.ndarray]:
    """Inference-time entry point used by src/50_spatial_genomics_deconv.py.

    Returns (rho_mean, D_mean) per row of X, clipped to physiological bounds.
    """
    with open(model_path, "rb") as f:
        bundle = pickle.load(f)
    rho = bundle["rho_model"].predict(np.asarray(X, dtype=np.float64))
    D = bundle["D_model"].predict(np.asarray(X, dtype=np.float64))
    rho = np.clip(rho, bundle["rho_range"][0], bundle["rho_range"][1])
    D = np.clip(D, bundle["D_range"][0], bundle["D_range"][1])
    return rho.astype(np.float32), D.astype(np.float32)


# --------------------------------------------------------------------------- #
# Single-pixel convenience used by deconvolution module
# --------------------------------------------------------------------------- #
def predict_one(
    neftel_fraction: np.ndarray,
    patient_id: Optional[str] = None,
    multiomic_path: Path = FEATURES_TSV_PATH,
    model_path: Path = MODEL_PATH,
) -> Tuple[float, float]:
    """Predict (rho, D) for one spot/voxel given the local Neftel fraction.

    If a patient_id is supplied, the patient's methylation/CNV/metabolic features
    are loaded once and combined with the spot-level Neftel vector.
    """
    if not Path(model_path).exists():
        # Graceful fallback: train on the spot-level fraction model
        if not Path(multiomic_path).exists():
            build_multiomic_tsv(
                patient_ids=[f"P{i}" for i in range(8)],
                neftel_fractions=np.random.default_rng(0).dirichlet(
                    np.ones(4), size=8
                ),
            )
        train_elasticnet(*load_multiomic_tsv(multiomic_path), model_path=model_path)

    if patient_id is not None and Path(multiomic_path).exists():
        pt = get_patient_vector(patient_id, multiomic_path)
        meth = pt[4:4 + len(METHYLATION_FEATURES)]
        cnv = pt[4 + len(METHYLATION_FEATURES):4 + len(METHYLATION_FEATURES) + len(CNV_FEATURES)]
        metab = pt[4 + len(METHYLATION_FEATURES) + len(CNV_FEATURES):]
    else:
        # Spot-level omic guess: neutral methylation, no CNV, mean metabolism
        meth = np.full(len(METHYLATION_FEATURES), 0.5)
        cnv = np.zeros(len(CNV_FEATURES))
        metab = np.full(len(METABOLIC_FEATURES), 1.0)
    # Ensure Neftel vector is 1D length 4
    nf = np.asarray(neftel_fraction, dtype=np.float64).reshape(4)
    if abs(nf.sum() - 1.0) > 1e-3:
        nf = nf / max(nf.sum(), 1e-12)
    X = np.concatenate([nf, meth, cnv, metab]).reshape(1, -1)
    rho, D = predict_params_from_features(X, model_path=model_path)
    return float(rho[0]), float(D[0])


# --------------------------------------------------------------------------- #
# CLI: train / verify the pipeline end-to-end on synthetic cohort
# --------------------------------------------------------------------------- #
def _main() -> None:
    import argparse
    p = argparse.ArgumentParser(description="Multi-omic feature fusion (Proposal 2)")
    p.add_argument("--n-patients", type=int, default=8)
    p.add_argument("--features-tsv", type=str, default=str(FEATURES_TSV_PATH))
    p.add_argument("--model-path", type=str, default=str(MODEL_PATH))
    args = p.parse_args()

    rng = np.random.default_rng(0)
    pids = [f"P{i:02d}" for i in range(args.n_patients)]
    fracs = rng.dirichlet(np.ones(4), size=args.n_patients)
    build_multiomic_tsv(pids, fracs, out_path=Path(args.features_tsv))

    X, y = load_multiomic_tsv(Path(args.features_tsv))
    metrics = train_elasticnet(X, y, model_path=Path(args.model_path))
    print(f"\n[multiomic] rho R2 = {metrics['rho_r2']:.3f}")
    print(f"[multiomic] D   R2 = {metrics['D_r2']:.3f}")

    rho_pred, D_pred = predict_params_from_features(X[:3], model_path=Path(args.model_path))
    print(f"[multiomic] sample predictions:")
    for i, (r, dd) in enumerate(zip(rho_pred, D_pred)):
        print(f"  P{i:02d}: rho={r:.4f}/day  D={dd:.4f} mm^2/day")


if __name__ == "__main__":
    _main()
