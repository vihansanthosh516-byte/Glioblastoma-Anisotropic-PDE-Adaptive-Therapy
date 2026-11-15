"""Read-only post hoc check: regression to the mean by input-size tertile on the spent outside sets. Writes nothing."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(sys.argv[1])
for name in ["track_b_v8_ucsd_external_20261114", "track_b_v8_ucsf_alptdg_20261114"]:
    p = pd.read_csv(ROOT / "output" / name / "predictions.csv")
    if "used_before" in p: p = p[~p.used_before]
    ch = p.y_out - p.persistence
    q = pd.qcut(p.persistence, 3, labels=["small", "mid", "large"])
    print(name, "n", len(p), "| signed mean change by input tertile", ch.groupby(q, observed=True).mean().round(3).to_dict(),
          "| median |change|", ch.abs().groupby(q, observed=True).median().round(3).to_dict(),
          "| corr(y_in, change)", round(float(np.corrcoef(p.persistence, ch)[0, 1]), 3),
          "| mean y_in", round(p.persistence.mean(), 3), "sd y_in", round(p.persistence.std(), 3))
    nz = p.y_out > 0
    print("   excluding zero-outcome pairs: corr", round(float(np.corrcoef(p.persistence[nz], ch[nz])[0, 1]), 3), "n", int(nz.sum()))
