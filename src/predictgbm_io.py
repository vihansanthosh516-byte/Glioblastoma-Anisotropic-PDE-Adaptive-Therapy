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
    ids = set()
    for f in ("lumiere.json", "rhuh.json"):
        ids |= {p["patient_id"] for p in json.loads((LISTS / f).read_text())["patients"]}
    return ids


def all_ids() -> list[str]:
    return sorted(p.name for p in DATA.iterdir() if p.is_dir())


def dev_ids() -> list[str]:
    t = test_ids()
    return [i for i in all_ids() if i not in t and i.startswith("tum_gbm")]


def is_frozen() -> bool:
    if not FROZEN.exists():
        return False
    r = subprocess.run(["git", "tag", "--list", "h2-frozen"], cwd=PROJECT_ROOT, capture_output=True, text=True)
    return r.stdout.strip() == "h2-frozen"


def _load(pid: str, name: str) -> np.ndarray:
    return np.asarray(nib.load(str(DATA / pid / name)).dataobj)


def load_case(pid: str, with_recurrence: bool = False) -> dict:
    """Pre-op inputs (always allowed) and, if asked, the recurrence (locked for test patients until the freeze)."""
    c = {"pid": pid, "seg": _load(pid, "tumor_seg.nii.gz").astype(np.int16),
         "wm": _load(pid, "wm_pbmap.nii.gz").astype(np.float32), "gm": _load(pid, "gm_pbmap.nii.gz").astype(np.float32),
         "csf": _load(pid, "csf_pbmap.nii.gz").astype(np.float32),
         "brain": binary_fill_holes(_load(pid, "t1c_bet_mask.nii.gz") > 0)}
    if with_recurrence:
        if pid in test_ids() and not is_frozen():
            raise PermissionError(f"{pid} is a test patient; recurrence locked until configs/h2_frozen.yaml + tag h2-frozen")
        c["rec"] = _load(pid, "recurrence_preop.nii.gz").astype(np.int16)
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


def model_plan(pred: np.ndarray, seg: np.ndarray, brain: np.ndarray) -> np.ndarray:
    p = pred.astype(np.float32)
    if np.isin(np.unique(p), [0, 1]).all():
        p = distance_fade(p > 0)
    p = np.clip(p, 0.0, 1.0)
    return topk_plan(p, int(standard_plan(seg, brain).sum()), brain)


def rec_enhancing(rec: np.ndarray) -> np.ndarray:
    return rec == 3


def rec_all(rec: np.ndarray) -> np.ndarray:
    return np.isin(rec, (1, 2, 3))


def coverage(rec_mask: np.ndarray, plan: np.ndarray) -> float:
    n = rec_mask.sum()
    return 1.0 if n == 0 else float((rec_mask & plan).sum() / n)
