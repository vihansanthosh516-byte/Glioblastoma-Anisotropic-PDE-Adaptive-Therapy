#!/usr/bin/env python3
"""Script 115: H-2 learned arm (GRAND_PLAN 12 B2) and hybrid arm (B3): small 3D U-Net on PREDICT-GBM development patients.

Development (TUM) only, 5-fold cross-validation by patient (fixed seed): each patient is scored by a model that never
saw it, so the dev coverage is out-of-fold, not in-sample. Test patients are never loaded (dev_ids; loader lock).
Inputs at 2 mm (block mean of 1 mm maps): core (labels 1, 3), edema (label 2), wm, gm, csf, brain, distance from the core
(mm / 50). Hybrid arm adds a growth-model channel: isotropic tissue model M1 density at lambda 2 mm (the script 112
selected cell), grown to the standard-plan volume. Target: all recurrence (labels 1, 2, 3) inside the brain.
MRI arms (unet_mri, hybrid_mri; user option 1, 2026-10-10) add the pre-op T1c and FLAIR at 2 mm from the private
pack of script 116 (GBM_PG_MRI); same training settings, so the MRI effect is isolated.
Loss: BCE + soft Dice. Output: logits, upsampled to 1 mm, plan = predictgbm_io.model_plan (scores mapped to 0-1 by
a sigmoid so the benchmark clip keeps the ranking). Coverage of enhancing and all recurrence per patient.
Output: output/predictgbm/unet_dev.json, unet_dev_pairs.csv
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
from pathlib import Path as _Path

import numpy as np

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "predictgbm"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import predictgbm_io as pg  # noqa: E402
from scipy.ndimage import distance_transform_edt  # noqa: E402

_spec = importlib.util.spec_from_file_location("s112", PROJECT_ROOT / "src" / "112_predictgbm_models_dev.py")
s112 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s112)

SEED, K = 20261009, 5
SHAPE = (120, 120, 77)


def m1_density(c, core2, std):
    """Script 112 M1 at lambda 2 mm (its selected cell), full 2 mm grid."""
    from solver_monotone import TensorFKMonotone
    wm, gm, brain2 = s112.to2(c["wm"]), s112.to2(c["gm"]), s112.to2(c["brain"]) >= 0.5
    lo = np.maximum(np.argwhere(core2).min(0) - s112.BOX_MM // 2, 0)
    hi = np.minimum(np.argwhere(core2).max(0) + s112.BOX_MM // 2 + 1, core2.shape)
    box = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
    scale = (wm + gm / 10.0)[box]
    scale[core2[box]] = np.maximum(scale[core2[box]], 1.0)
    dom = (brain2[box] & (scale > 0.05)) | core2[box]
    D = 2.0 * 2.0 * s112.RHO
    fk = TensorFKMonotone(s112.tensor(None, scale * D, box), dom, h=2.0)
    u, t = core2[box].astype(np.float32), 0.0
    while t < s112.T_MAX and (u >= s112.U_FRONT).sum() < std.sum() / 8.0:
        u = fk.run(u, [s112.RHO], 2.0)[0][0]
        t += 2.0
    full = np.zeros(core2.shape, np.float32)
    full[box] = u
    return full


_MRI = None


def mri2(pid):
    global _MRI
    if _MRI is None:
        _MRI = np.load(os.environ["GBM_PG_MRI"])
    return [_MRI[f"{pid}|{m}"].astype(np.float32) for m in ("t1c_bet_normalized.nii.gz", "flair_bet_normalized.nii.gz")]


def build(pid, hybrid, mri=False):
    c = pg.load_case(pid, with_recurrence=True)
    core1 = pg.core(c["seg"])
    std = pg.standard_plan(c["seg"], c["brain"])
    core2 = s112.to2(core1) >= 0.5
    dist = distance_transform_edt(~core1)[:240, :240, :154]
    dist2 = dist.reshape(120, 2, 120, 2, 77, 2).min((1, 3, 5)) / 50.0
    x = [core2, s112.to2(c["seg"] == 2), s112.to2(c["wm"]), s112.to2(c["gm"]), s112.to2(c["csf"]),
         s112.to2(c["brain"]), np.minimum(dist2, 3.0)]
    if hybrid:
        x.append(m1_density(c, core2, std) if core2.any() else np.zeros(SHAPE, np.float32))
    if mri:
        x += mri2(pid)
    y = s112.to2(pg.rec_all(c["rec"]) & c["brain"]) >= 0.5
    return {"x": np.stack(x).astype(np.float32), "y": y.astype(np.float32),
            "brain2": (s112.to2(c["brain"]) >= 0.5).astype(np.float32)}


def make_net(cin):
    import torch.nn as nn

    def blk(a, b):
        return nn.Sequential(nn.Conv3d(a, b, 3, padding=1), nn.InstanceNorm3d(b), nn.LeakyReLU(0.1),
                             nn.Conv3d(b, b, 3, padding=1), nn.InstanceNorm3d(b), nn.LeakyReLU(0.1))

    class UNet(nn.Module):
        def __init__(self):
            super().__init__()
            f = [16, 32, 64, 128]
            self.e = nn.ModuleList([blk(cin, f[0]), blk(f[0], f[1]), blk(f[1], f[2]), blk(f[2], f[3])])
            self.up = nn.ModuleList([nn.ConvTranspose3d(f[i + 1], f[i], 2, 2) for i in range(3)])
            self.d = nn.ModuleList([blk(2 * f[i], f[i]) for i in range(3)])
            self.pool = nn.MaxPool3d(2)
            self.head = nn.Conv3d(f[0], 1, 1)

        def forward(self, x):
            s = []
            for i, e in enumerate(self.e):
                x = e(x if i == 0 else self.pool(x))
                s.append(x)
            for i in (2, 1, 0):
                x = self.d[i](torch_cat([self.up[i](x), s[i]]))
            return self.head(x)

    return UNet()


def torch_cat(xs):
    import torch
    return torch.cat(xs, 1)


PAD = (120, 120, 80)   # divisible by 8


def pad(a):
    out = np.zeros(a.shape[:-3] + PAD, np.float32)
    out[..., :77] = a
    return out


def train_fold(data, tr, epochs, dev):
    import torch
    torch.manual_seed(SEED)
    net = make_net(data[tr[0]]["x"].shape[0]).to(dev)
    opt = torch.optim.AdamW(net.parameters(), 1e-3, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs)
    rng = np.random.default_rng(SEED)
    for ep in range(epochs):
        net.train()
        tot = 0.0
        for i in rng.permutation(tr):
            d = data[i]
            x, y, b = (torch.from_numpy(pad(v))[None].to(dev) for v in (d["x"], d["y"], d["brain2"]))
            if rng.random() < 0.5:   # left-right flip augmentation
                x, y, b = x.flip(2), y.flip(1), b.flip(1)
            logit = net(x)[:, 0]
            p = torch.sigmoid(logit) * b
            bce = torch.nn.functional.binary_cross_entropy_with_logits(logit, y, weight=b)
            dice = 1 - (2 * (p * y).sum() + 1) / (p.sum() + y.sum() + 1)
            loss = bce + dice
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += float(loss)
        sched.step()
        if ep % 10 == 0 or ep == epochs - 1:
            print(f"  epoch {ep} loss {tot / len(tr):.4f}", flush=True)
    return net


def run(arm, epochs):
    import torch
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    hybrid, mri = arm.startswith("hybrid"), arm.endswith("_mri")
    ids = pg.dev_ids()
    t0 = time.time()
    data = {}
    for pid in ids:
        data[pid] = build(pid, hybrid, mri)
    print(f"built {len(data)} in {time.time() - t0:.0f}s on {dev}", flush=True)
    order = np.random.default_rng(SEED).permutation(ids)
    folds = [list(order[k::K]) for k in range(K)]
    rows = []
    for k in range(K):
        tr = [p for j in range(K) if j != k for p in folds[j]]
        net = train_fold(data, tr, epochs, dev)
        net.eval()
        for pid in folds[k]:
            with torch.no_grad():
                logit = net(torch.from_numpy(pad(data[pid]["x"]))[None].to(dev))[0, 0, :, :, :77].cpu().numpy()
            c = pg.load_case(pid, with_recurrence=True)
            std = pg.standard_plan(c["seg"], c["brain"])
            score = 1.0 / (1.0 + np.exp(-s112.up1(logit)))
            plan = pg.model_plan(score, c["seg"], c["brain"], k=int(std.sum()))
            re, ra = pg.rec_enhancing(c["rec"]), pg.rec_all(c["rec"])
            rows.append({"pid": pid, "fold": k, "n_rec_enh": int(re.sum()), "standard_enh": pg.coverage(re, std),
                         "standard_all": pg.coverage(ra, std), f"{arm}_enh": pg.coverage(re, plan),
                         f"{arm}_all": pg.coverage(ra, plan)})
        print(f"fold {k} done [{(time.time() - t0) / 60:.1f} min]", flush=True)
    return rows


def main():
    import pandas as pd
    from src import patient_stats as ps
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", default="unet,hybrid")
    ap.add_argument("--epochs", type=int, default=60)
    ap.add_argument("--merge-prev", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    df = None
    for arm in a.arms.split(","):
        d = pd.DataFrame(run(arm, a.epochs))
        df = d if df is None else df.merge(d[["pid", f"{arm}_enh", f"{arm}_all"]], on="pid")
    if a.merge_prev and (OUT / "unet_dev_pairs.csv").exists():   # add earlier arms (same folds, same seed)
        prev = pd.read_csv(OUT / "unet_dev_pairs.csv")
        keep = [c for c in prev.columns if c.endswith(("_enh", "_all")) and not c.startswith("standard") and c not in df]
        df = df.merge(prev[["pid"] + keep], on="pid")
    tag = "_mri" if "mri" in a.arms else ""
    df.to_csv(OUT / f"unet_dev{tag}_pairs.csv", index=False)
    use = df[df["n_rec_enh"] > 0]
    res = {"script": "115_predictgbm_unet_dev", "split": "development (TUM), 5-fold out-of-fold", "seed": SEED,
           "epochs": a.epochs, "n_patients": int(len(use)), "standard_enh_mean": float(use["standard_enh"].mean()),
           "arms": {}}
    for arm in [c[:-4] for c in df.columns if c.endswith("_enh") and c not in ("standard_enh",)]:
        res["arms"][arm] = {"enh_mean": float(use[f"{arm}_enh"].mean()), "all_mean": float(use[f"{arm}_all"].mean()),
                            "delta_vs_standard_enh": ps.summarize_delta((use[f"{arm}_enh"] - use["standard_enh"]).to_numpy())}
    for a2, b2 in (("unet_mri", "unet"), ("hybrid_mri", "unet_mri")):
        if f"{a2}_enh" in df.columns and f"{b2}_enh" in df.columns:
            res[f"{a2}_minus_{b2}_enh"] = ps.summarize_delta((use[f"{a2}_enh"] - use[f"{b2}_enh"]).to_numpy())
    if "hybrid" in res["arms"] and "unet" in res["arms"]:
        res["hybrid_minus_unet_enh"] = ps.summarize_delta((use["hybrid_enh"] - use["unet_enh"]).to_numpy())
    res["weakest_points"] = ["Development patients; architecture and epochs set before the run, not tuned, but only TUM.",
                             "2 mm inputs; no MRI intensities (segmentation and tissue maps only)."]
    (OUT / f"unet_dev{tag}.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=float)[:1500])


if __name__ == "__main__":
    main()
