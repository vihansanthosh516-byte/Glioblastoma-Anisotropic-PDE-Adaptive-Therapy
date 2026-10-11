#!/usr/bin/env python3
"""Script 122: H-2 FINAL model (frozen in configs/h2_frozen.yaml): train once, predict the locked test patients.

Final model by the GRAND_PLAN 15 selection rule (highest out-of-fold dev enhancing coverage at equal volume; on the
137 common patients: hybrid_rank_mu 0.8529 > ensemble GliODIL+U-Net 0.8500 > hybrid_mri_rank_mu 0.8519; ledger 2aw).
Arm hybrid_rank_mu of script 115: inputs core, edema, wm, gm, csf, brain, distance, M1 density; BCE + Dice + soft top-k
ranking loss; training data = all 142 TUM development patients + the 334 MU pairs (script 121 pack); 15 epochs (the
validated setting), seed 20261009, same network and optimiser as script 115.
Test patients (LUMIERE, RHUH): ONLY pre-op inputs are read (no recurrence; none is in the Kaggle pack). Output: 2 mm
logits per test patient. Scoring is script 123, which needs the freeze (loader lock).
Output: output/predictgbm/final/test_logits_2mm.npz, final_train.json
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path as _Path

import numpy as np

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "predictgbm" / "final"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import predictgbm_io as pg  # noqa: E402
from scipy.ndimage import distance_transform_edt  # noqa: E402

_s = importlib.util.spec_from_file_location("s115", PROJECT_ROOT / "src" / "115_predictgbm_unet_dev.py")
s115 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(s115)
s112 = s115.s112

ARM, EPOCHS = "hybrid_rank_mu", 15


def build_inputs_only(pid):
    """Same channels as s115.build(hybrid=True, mri=False), from pre-op inputs only (no recurrence read)."""
    c = pg.load_case(pid, with_recurrence=False)
    core1 = pg.core(c["seg"])
    std = pg.standard_plan(c["seg"], c["brain"])
    core2 = s112.to2(core1) >= 0.5
    dist = distance_transform_edt(~core1)[:240, :240, :154] if core1.any() else np.full((240, 240, 154), 150.0)
    dist2 = dist.reshape(120, 2, 120, 2, 77, 2).min((1, 3, 5)) / 50.0
    x = [core2, s112.to2(c["seg"] == 2), s112.to2(c["wm"]), s112.to2(c["gm"]), s112.to2(c["csf"]),
         s112.to2(c["brain"]), np.minimum(dist2, 3.0)]
    x.append(s115.m1_density(c, core2, std) if core2.any() else np.zeros(s115.SHAPE, np.float32))
    return np.stack(x).astype(np.float32)


def main():
    import torch
    OUT.mkdir(parents=True, exist_ok=True)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    t0 = time.time()
    data = {pid: s115.build(pid, True, False) for pid in pg.dev_ids()}
    mu = [f"mu:{k}" for k in s115._mu()["keys"]]
    for k in mu:
        data[k] = s115.build_mu(k[3:], True, False)
    print(f"built {len(data)} training samples in {time.time() - t0:.0f}s", flush=True)
    net = s115.train_fold(data, list(pg.dev_ids()) + mu, EPOCHS, dev, rank=True)
    net.eval()
    test = sorted(pg.test_ids() & set(pg.all_ids()))
    logits = {}
    for pid in test:
        x = build_inputs_only(pid)
        with torch.no_grad():
            logits[pid] = net(torch.from_numpy(s115.pad(x))[None].to(dev))[0, 0, :, :, :77].cpu().numpy().astype(np.float16)
    np.savez_compressed(OUT / "test_logits_2mm.npz", **logits)
    torch.save(net.state_dict(), OUT / "final_model_state.pt")
    rec = {"script": "122_predictgbm_final_train_predict", "arm": ARM, "epochs": EPOCHS, "seed": s115.SEED,
           "n_train_tum": len(pg.dev_ids()), "n_train_mu": len(mu), "n_test_predicted": len(logits),
           "recurrence_read_for_test": False, "minutes": round((time.time() - t0) / 60, 1)}
    (OUT / "final_train.json").write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
