"""PREDICT-GBM data access and evaluation (Amendment 6 A6.6, Amendment 7).

Evaluation functions reproduce BrainLesion/PredictGBM `predict_gbm/evaluation/evaluate.py` and `metrics.py`
(commit db3dd1b, Apache-2.0): standard plan = Euclidean distance <= 15 voxels (1 mm) from the pre-op core
(labels 1, 3), masked to the brain (holes filled); model plan = top-k voxels of the model map inside the brain at the
standard plan's voxel count, with the distance-fade fallback when the top-k contains non-positive scores;
coverage = |recurrence AND plan| / |recurrence| (1.0 if no recurrence). Enhancing recurrence = label 3;
"all" recurrence = labels 1, 2, 3 (cavity label 4 ignored).

TEST LOCK (A6.6): recurrence masks of the test patients (LUMIERE, RHUH; lists from the PredictGBM repo) are refused
unless configs/h2_frozen.yaml exists AND the git tag h2-frozen exists. Development = TUM patients.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path as _Path

import nibabel as nib
import numpy as np
from scipy.ndimage import binary_fill_holes, distance_transform_edt

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
DATA = PROJECT_ROOT / "data" / "external" / "predict_gbm" / "predict_gbm"
LISTS = PROJECT_ROOT / "data" / "external" / "PredictGBM_code" / "predict_gbm" / "data" / "datasets"
FROZEN = PROJECT_ROOT / "configs" / "h2_frozen.yaml"
CTV_MARGIN = 15
PUBLISHED_MODELS = ["unet", "gliodil", "sbtc", "pinngbm", "lmi", "gliomap"]   # sbtc = LOTI in the paper


def test_ids() -> set[str]:
    import os
    if os.environ.get("GBM_PG_TEST_IDS"):
        return set(json.loads(_Path(os.environ["GBM_PG_TEST_IDS"]).read_text()))
    ids = set()
    for f in ("lumiere.json", "rhuh.json"):
        ids |= {p["patient_id"] for p in json.loads((LISTS / f).read_text())["patients"]}
    return ids


def all_ids() -> list[str]:
    import os
    if os.environ.get("GBM_PG_PACK"):
        return sorted({k.split("|")[0] for k in np.load(os.environ["GBM_PG_PACK"]).files})
    return sorted(p.name for p in DATA.iterdir() if p.is_dir())


def dev_ids() -> list[str]:
    t = test_ids()
    return [i for i in all_ids() if i not in t and i.startswith("tum_gbm")]


def is_frozen() -> bool:
    if not FROZEN.exists():
        return False
    r = subprocess.run(["git", "tag", "--list", "h2-frozen"], cwd=PROJECT_ROOT, capture_output=True, text=True)
    return r.stdout.strip() == "h2-frozen"


_PACK = None   # GBM_PG_PACK = npz from src/113_pack_predictgbm_for_kaggle.py (Kaggle runs)


def _load(pid: str, name: str) -> np.ndarray:
    import os
    if os.environ.get("GBM_PG_PACK"):
        global _PACK
        if _PACK is None:
            _PACK = np.load(os.environ["GBM_PG_PACK"])
        key = f"{pid}|{name}"
        if key + "|bits" in _PACK.files:
            return np.unpackbits(_PACK[key + "|bits"])[: 240 * 240 * 155].reshape(240, 240, 155)
        return _PACK[key]
    return np.asarray(nib.load(str(DATA / pid / name)).dataobj)


def seg_labels(a: np.ndarray) -> np.ndarray:
    """As PredictGBM load_segmentation: round to the nearest integer (some files store 1 as 0.996)."""
    return np.rint(np.asarray(a, np.float64)).astype(np.int16)


def load_case(pid: str, with_recurrence: bool = False) -> dict:
    """Pre-op inputs (always allowed) and, if asked, the recurrence (locked for test patients until the freeze)."""
    c = {"pid": pid, "seg": seg_labels(_load(pid, "tumor_seg.nii.gz")),
         "wm": _load(pid, "wm_pbmap.nii.gz").astype(np.float32), "gm": _load(pid, "gm_pbmap.nii.gz").astype(np.float32),
         "csf": _load(pid, "csf_pbmap.nii.gz").astype(np.float32),
         "brain": binary_fill_holes(_load(pid, "t1c_bet_mask.nii.gz") > 0)}
    if with_recurrence:
        if pid in test_ids() and not is_frozen():
            raise PermissionError(f"{pid} is a test patient; recurrence locked until configs/h2_frozen.yaml + tag h2-frozen")
        c["rec"] = seg_labels(_load(pid, "recurrence_preop.nii.gz"))
    return c


def core(seg: np.ndarray) -> np.ndarray:
    return (seg == 1) | (seg == 3)


def standard_plan(seg: np.ndarray, brain: np.ndarray, margin: int = CTV_MARGIN) -> np.ndarray:
    plan = distance_transform_edt(~core(seg)) <= margin
    return plan & brain


def distance_fade(binary: np.ndarray) -> np.ndarray:
    inside = binary.astype(bool)
    if not inside.any():
        return np.zeros(binary.shape, np.float32)
    if inside.all():
        return np.ones(binary.shape, np.float32)
    d = distance_transform_edt(~inside)
    fade = 1.0 - d / d.max() if d.max() > 0 else np.zeros_like(d)
    fade[inside] = 1.0
    return fade.astype(np.float32)


def topk_plan(scores: np.ndarray, k: int, mask: np.ndarray) -> np.ndarray:
    out = np.zeros(scores.shape, bool)
    cand = np.flatnonzero(mask)
    if k <= 0 or cand.size == 0:
        return out
    if k >= cand.size:
        out.flat[cand] = True
        return out
    s = np.nan_to_num(scores.astype(np.float32).flat[cand], nan=-np.inf, posinf=np.finfo(np.float32).max, neginf=-np.inf)
    order = np.argsort(s, kind="stable")
    sel = s[order[-k:]]
    if not (np.any(sel <= 0) or np.any(~np.isfinite(sel))):
        out.flat[cand[order[-k:]]] = True
        return out
    fade = distance_fade((scores > 0) & mask).flat[cand]
    order = np.argsort(np.nan_to_num(fade, nan=-np.inf), kind="stable")
    out.flat[cand[order[-k:]]] = True
    return out


def resample_to(a: np.ndarray, shape) -> np.ndarray:
    """As ants.resample_image(use_voxels=True, interp_type=0 = linear): same origin, new spacing = old * n_old / n_new,
    so new index i samples old index i * n_old / n_new."""
    if a.shape == tuple(shape):
        return a
    from scipy.ndimage import map_coordinates
    grid = np.meshgrid(*[np.arange(n) * (o / n) for n, o in zip(shape, a.shape)], indexing="ij")
    return map_coordinates(a.astype(np.float32), grid, order=1, mode="nearest")


def model_plan(pred: np.ndarray, seg: np.ndarray, brain: np.ndarray, k: int | None = None) -> np.ndarray:
    """k = standard-plan voxel count; pass it when known to skip recomputing the standard plan."""
    p = resample_to(pred.astype(np.float32), seg.shape)
    if np.isin(np.unique(p), [0, 1]).all():
        p = distance_fade(p > 0)
    p = np.clip(p, 0.0, 1.0)
    return topk_plan(p, int(standard_plan(seg, brain).sum()) if k is None else int(k), brain)


def rec_enhancing(rec: np.ndarray) -> np.ndarray:
    return rec == 3


def rec_all(rec: np.ndarray) -> np.ndarray:
    return np.isin(rec, (1, 2, 3))


def coverage(rec_mask: np.ndarray, plan: np.ndarray) -> float:
    n = rec_mask.sum()
    return 1.0 if n == 0 else float((rec_mask & plan).sum() / n)
