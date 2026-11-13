"""Controller model set for the Track C v2 research run, read from the config.

Decision D1 (paper/trackc_amendment_unblock_2026-11-12.md): the controller uses the
config's 5-level synthetic_q_signal_grid x the 3x3x3 (a, b, d) grid = 135 models, so
every true world q level (0.6-0.9) is inside the controller's model set. The older
81-model grid in src/212_trackc_signal_fisher_grid.py stays unchanged as the record
of the historical diagnostics that used it.
"""
from __future__ import annotations

import json
from itertools import product
from pathlib import Path

from trackc_integrated_signal_model import SignalModelParameters

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = PROJECT_ROOT / "configs" / "trackABC_integrated_v2.json"


def _config(config=None) -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8")) if config is None else config


def research_q_grid(config=None) -> tuple:
    return tuple(float(q) for q in _config(config)["signal_model"]["synthetic_q_signal_grid"])


def research_model_grid(config=None) -> tuple:
    """All (q, a, b, d) controller models in stable order."""
    c = _config(config)
    coeff = tuple(float(v) for v in c["signal_model"]["a_b_d_grid"])
    models = tuple(
        SignalModelParameters(q=q, a=a, b=b, d=d)
        for q, a, b, d in product(research_q_grid(c), coeff, coeff, coeff)
    )
    maximum = int(c["controller_comparison"]["candidate_model_set_maximum"])
    if len(models) > maximum:
        raise ValueError(f"model set has {len(models)} models; config maximum is {maximum}")
    sampling = {float(q) for q in c["signal_model"]["primary_q_signal_sampling_levels"]}
    if not sampling <= set(research_q_grid(c)):
        raise ValueError("a sampled world q level is outside the controller q grid")
    return models
