#!/usr/bin/env python3
"""Script 136: H-4 step 2 retrain (Amendment 13, A13.2) -- Kaggle.

Trains arm hybrid_mri_rank_mu once on the pooled development set:
  TUM (pg.dev_ids(), full, not out-of-fold) + MU pairs (GBM_MU_PACK) + the script 135 pack
  (GBM_H4_RETRAIN_PACK: LUMIERE 61 + RHUH 40 + UCSF 59, time1 targets) = 636 samples.
Identical settings to script 132 (script 115 train_fold, 15 epochs, seed 20261009). No test data of any kind
is read or predicted; the plan rule stays G10_1.0 and the freeze happens later (A13.6).
Output: output/predictgbm/h4_retrain/h4_retrain_model_state.pt, h4_retrain_train.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import time
from pathlib import Path as _Path

import numpy as np

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "predictgbm" / "h4_retrain"
import sys
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import predictgbm_io as pg  # noqa: E402

_s = importlib.util.spec_from_file_location("s115", PROJECT_ROOT / "src" / "115_predictgbm_unet_dev.py")
s115 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(s115)

EPOCHS = 15


def main():
    import torch
    OUT.mkdir(parents=True, exist_ok=True)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    t0 = time.time()
    data = s115.DiskData(os.environ["GBM_DISK_CACHE"]) if os.environ.get("GBM_DISK_CACHE") else {}
    for pid in pg.dev_ids():
        data[pid] = s115.build(pid, True, True)
    mu = [f"mu:{k}" for k in s115._mu()["keys"]]
    for k in mu:
        data[k] = s115.build_mu(k[3:], True, True)
    z = np.load(os.environ["GBM_H4_RETRAIN_PACK"])
    new_ids = sorted({k.split("|")[0] for k in z.files})
    for pid in new_ids:
        data[pid] = {"x": z[f"{pid}|x"].astype(np.float16), "y": z[f"{pid}|y"].astype(np.float32),
                     "ye": z[f"{pid}|ye"].astype(np.float32), "k2": int(z[f"{pid}|k2"]),
                     "brain2": z[f"{pid}|brain2"].astype(np.float32)}
    tr = list(pg.dev_ids()) + mu + new_ids
    print(f"built {len(data)} samples in {time.time() - t0:.0f}s ({len(pg.dev_ids())} TUM, {len(mu)} MU, {len(new_ids)} new)", flush=True)
    net = s115.train_fold(data, tr, EPOCHS, dev, rank=True)
    sp = OUT / "h4_retrain_model_state.pt"
    torch.save(net.state_dict(), sp)
    import hashlib
    rec = {"script": "136_h4_retrain_train", "arm": "hybrid_mri_rank_mu", "epochs": EPOCHS, "seed": s115.SEED,
           "n_train_tum": len(pg.dev_ids()), "n_train_mu": len(mu), "n_train_new": len(new_ids),
           "n_train_total": len(tr), "amendment": "A13.2", "test_data_read": False,
           "model_state_sha256": hashlib.sha256(sp.read_bytes()).hexdigest(),
           "minutes": round((time.time() - t0) / 60, 1)}
    (OUT / "h4_retrain_train.json").write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
