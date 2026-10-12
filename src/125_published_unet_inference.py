#!/usr/bin/env python3
"""Script 125: run the released PREDICT-GBM U-Net (rebuild in unet_published.py) leak-free (Amendment 9, A9.5).

--stage validate : choose the rebuild variant (norm x activation x seg encoding) that best reproduces the RELEASED
                   out-of-fold predictions (unet_pred.nii.gz) on 8 TUM dev patients, using the 3 checkpoints of each
                   patient's held-out fold. No recurrence is read. Output: output/predictgbm/h3/unet_validate.json
--stage dev      : per TUM dev patient, mean logit of the 3 checkpoints whose fold lists the patient as TEST (A9.5b).
                   Output (outside git): <UCSF_H3_DIR>/unet_dev_logits.npz
--stage ucsf     : per UCSF H-3 patient (script 124 inputs), mean logit of all 15 checkpoints (none saw UCSF).
                   Output (outside git): <UCSF_H3_DIR>/unet_ucsf_logits.npz
Inputs exactly as the authors' container (main.py / transform.py, read 2026-10-12, not run): channel order csf, wm, gm,
tumour seg (raw labels), T1c, FLAIR, ADC (missing -> zeros, as the container does); 240 x 240 x 155 symmetric
zero-padded to 256 x 256 x 160; output = mean logit of the models, Gaussian smoothing sigma 1, times the hole-filled
brain mask. Weights: torch.load(weights_only=True), outside the repo.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
from pathlib import Path as _Path

import numpy as np
import torch

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "predictgbm" / "h3"
H3 = _Path(os.environ.get("UCSF_H3_DIR", r"C:\Users\vihan\data\ucsf_h3"))
WDIR = _Path(os.environ.get("UNET_WEIGHTS", r"C:\Users\vihan\data\unet_weights"))
SPLITS = PROJECT_ROOT / "data" / "external" / "unet_splits" / "ensemble_data_splits.json"
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import predictgbm_io as pg  # noqa: E402
from unet_published import PublishedUNet  # noqa: E402

torch.set_num_threads(os.cpu_count() or 8)
W = {p.stem: p for p in WDIR.rglob("ens_*_crv_*.tar")}
_NETS = {}


def test_fold():
    sp = json.loads(SPLITS.read_text())["ensemble_members"]["ensemble_1"]
    return {e["patient_id"]: k for k in range(5) for e in sp[f"fold_{k}"]["test"]}


def net(name, norm, act):
    key = (name, norm, act)
    if key not in _NETS:
        n = PublishedUNet(7, norm, act)
        n.load_state_dict(torch.load(W[name], map_location="cpu", weights_only=True))
        _NETS.clear() if len(_NETS) > 15 else None
        _NETS[key] = n.eval()
    return _NETS[key]


def stack(seg, wm, gm, csf, t1c, flair, segmode):
    seg = np.rint(seg).astype(np.float32)
    s = {"raw": seg, "div3": seg / 3.0, "div4": seg / 4.0}[segmode]
    x = np.zeros((7, 256, 256, 160), np.float32)
    x[:6, 8:248, 8:248, 2:157] = np.stack([csf, wm, gm, s, t1c, flair])[..., :155]
    return torch.from_numpy(x)[None]


def predict(x, names, norm, act, brain=None):
    from scipy.ndimage import binary_fill_holes, gaussian_filter
    with torch.no_grad():
        out = np.mean([net(n, norm, act)(x)[0, 0].numpy() for n in names], 0)
    out = gaussian_filter(out, 1.0)[8:248, 8:248, 2:157]
    return out if brain is None else out * binary_fill_holes(brain)


def pg_inputs(pid, segmode):
    g = lambda n: np.asarray(pg._load(pid, n), np.float32)   # noqa: E731
    return stack(g("tumor_seg.nii.gz"), g("wm_pbmap.nii.gz"), g("gm_pbmap.nii.gz"), g("csf_pbmap.nii.gz"),
                 g("t1c_bet_normalized.nii.gz"), g("flair_bet_normalized.nii.gz"), segmode)


def pgbrain(pid):
    return np.asarray(pg._load(pid, "t1c_bet_mask.nii.gz")) > 0


def fold_names(k):
    return [f"ens_{e}_crv_{k + 1}" for e in (1, 2, 3)]


def validate():
    from scipy.stats import spearmanr
    fold = test_fold()
    ids = sorted(pg.dev_ids())[:8]
    t0 = time.time()

    def score(pred, pid):
        ref = np.asarray(pg._load(pid, "unet_pred.nii.gz"), np.float32)
        m = np.asarray(pg._load(pid, "t1c_bet_mask.nii.gz")) > 0
        return round(float(spearmanr(pred[m][::20], ref[m][::20])[0]), 4), round(float(np.abs(pred[m] - ref[m]).mean()), 4)

    res = {}
    for norm, act, seg in itertools.product(("instance",), ("leaky_e2", "relu"), ("raw",)):
        res[f"{norm}|{act}|{seg}"] = [score(predict(pg_inputs(p, seg), [f"ens_1_crv_{fold[p] + 1}"], norm, act), p) for p in ids[:2]]
        print(norm, act, seg, res[f"{norm}|{act}|{seg}"], f"{time.time() - t0:.0f}s", flush=True)
    best = max(res, key=lambda v: np.mean([a for a, _ in res[v]]))
    norm, act, seg = best.split("|")
    conf = {p: score(predict(pg_inputs(p, seg), fold_names(fold[p]), norm, act, pgbrain(p)), p) for p in ids}
    out = {"script": "125_published_unet_inference", "stage": "validate", "screen_member1_2patients": res, "best": best,
           "confirm_3members_8patients": conf, "confirm_mean_spearman": float(np.mean([a for a, _ in conf.values()]))}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "unet_validate.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out["confirm_3members_8patients"]), "BEST", best, out["confirm_mean_spearman"])


def variant():
    return json.loads((OUT / "unet_validate.json").read_text())["best"].split("|")


def dev():
    norm, act, seg = variant()
    fold = test_fold()
    logits = {}
    for i, pid in enumerate(sorted(pg.dev_ids())):
        logits[pid] = predict(pg_inputs(pid, seg), fold_names(fold[pid]), norm, act, pgbrain(pid)).astype(np.float16)
        print(i, pid, flush=True)
    np.savez_compressed(H3 / "unet_dev_logits.npz", **logits)


def ucsf():
    norm, act, seg = variant()
    logits = {}
    names = sorted(W)
    for f in sorted((H3 / "inputs").glob("*.npz")):
        z = np.load(f)
        t1c, fl, b = z["t1c"], z["flair"], z["brain"]

        def nrm(a):   # BrainLes PercentileNormalizer as in the authors' package: 0.1 / 99.9 percentiles of the whole
            a = a * b    # skull-stripped volume (zeros included), mapped to [0, 1] and clipped (checked on released maps)
            lo, hi = np.percentile(a, [0.1, 99.9])
            return np.clip((a - lo) / max(hi - lo, 1e-6), 0, 1) * b

        x = stack(z["seg"], z["wm"], z["gm"], z["csf"], nrm(t1c), nrm(fl), seg)
        logits[f.stem] = predict(x, names, norm, act, b).astype(np.float16)
        print(f.stem, flush=True)
    np.savez_compressed(H3 / "unet_ucsf_logits.npz", **logits)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["validate", "dev", "ucsf"], required=True)
    {"validate": validate, "dev": dev, "ucsf": ucsf}[ap.parse_args().stage]()
