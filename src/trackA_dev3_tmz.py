"""Track A TMZ development falsification 3: origin-level confounding.

Plan: paper/trackA_dev3_tmz_plan_2026-11-23.md (frozen before run). Generated data only.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.preprocessing import StandardScaler

from .trackA_dev2 import (PENALTY, W, admit_direction, augmented, balanced_fit, blind_direction,
                          origin_mae, summarize, within_direction)
from .trackA_tmz_model import Cultures, _fit_shared, _groups, _shuffle_links

DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "output/trackA_dev3_tmz/result.json"
SEED = 202611233
REPLICATES = 20
CONFOUND_GENE = 4  # zero-based column index
SCENARIOS = {  # (train confound, eval confound, effect)
    "confounded_transfer": (0.8, 0.0, W),
    "confounded_both": (0.8, 0.8, W),
    "strong_confound": (1.5, 0.0, W),
    "plain_transfer": (0.0, 0.0, W),
    "null": (0.0, 0.0, 0 * W),
}


def cultures(rng, prefix, count, *, effect, confound):
    sizes = rng.integers(1, 5, size=count)
    groups = np.repeat(np.arange(count), sizes)
    common, noise = rng.normal(size=(count, 8)), rng.normal(size=(len(groups), 8))
    z = rng.normal(size=count)
    x = 0.6 * common[groups] + noise
    x[:, CONFOUND_GENE] += 0.8 * z[groups]
    y = (0.6 * common @ effect + confound * z)[groups] + noise @ effect
    y += rng.normal(scale=0.3, size=count)[groups] + rng.normal(scale=0.3, size=len(groups))
    return Cultures(x, y, [f"{prefix}-{g:03d}" for g in groups],
                    [f"{prefix}-c{i:04d}" for i in range(len(groups))])


def replicate(rng, scenario):
    train_c, eval_c, effect = SCENARIOS[scenario]
    train = cultures(rng, "t", 40, effect=effect, confound=train_c)
    ev = cultures(rng, "e", 100, effect=effect, confound=eval_c)
    _, inv, counts = _groups(train.origins)
    scale = StandardScaler().fit(train.expression, sample_weight=1.0 / (len(counts) * counts[inv]))
    x, ex, y, o = scale.transform(train.expression), scale.transform(ev.expression), train.viability, train.origins
    shuffled = _shuffle_links(o, ["assay"] * len(o), int(rng.integers(1 << 31)))
    preds = {"origin_balanced_ridge": balanced_fit(x, o, y).predict(ex)}
    admitted, cosine = admit_direction(x, y, o, np.random.default_rng(int(rng.integers(1 << 31))))
    d_true = within_direction(x, y, o)
    preds["replicated_direction"] = augmented(x, o, y, d_true, ex) if admitted else preds["origin_balanced_ridge"]
    preds["direction_unchecked"] = augmented(x, o, y, d_true, ex)
    preds["shuffled_link_direction"] = augmented(x, o, y, within_direction(x, y, shuffled), ex)
    preds["outcome_blind_direction"] = augmented(x, o, y, blind_direction(x, o), ex)
    intercept, slope = _fit_shared(x, y, o, o, PENALTY, 1.0)
    preds["one_slope_contrast"] = intercept + ex @ slope
    return {"admitted": bool(admitted), "cosine": cosine,
            "confound_weight_true_direction": float(d_true[CONFOUND_GENE]),
            "origin_mae": {a: origin_mae(p, ev.viability, ev.origins) for a, p in preds.items()}}


def run():
    result = {}
    for scen, seq in zip(SCENARIOS, np.random.SeedSequence(SEED).spawn(len(SCENARIOS)), strict=True):
        reps = [replicate(np.random.default_rng(s), scen) for s in seq.spawn(REPLICATES)]
        diff = lambda arm: summarize([r["origin_mae"][arm] - r["origin_mae"]["replicated_direction"] for r in reps])
        result[scen] = {
            "mean_origin_mae": {a: float(np.mean([r["origin_mae"][a] for r in reps])) for a in reps[0]["origin_mae"]},
            "base_minus_candidate": diff("origin_balanced_ridge"),
            "shuffled_minus_candidate": diff("shuffled_link_direction"),
            "admitted_fraction": float(np.mean([r["admitted"] for r in reps])),
            "replicates": reps,
        }
    ok = lambda s, k: result[s][k]["ci95"][0] > 0
    rules = {
        "confounded_transfer_beats_base": ok("confounded_transfer", "base_minus_candidate"),
        "confounded_transfer_beats_shuffle": ok("confounded_transfer", "shuffled_minus_candidate"),
        "strong_confound_beats_base": ok("strong_confound", "base_minus_candidate"),
        "strong_confound_beats_shuffle": ok("strong_confound", "shuffled_minus_candidate"),
        "confounded_both_safe": result["confounded_both"]["base_minus_candidate"]["mean"] >= -0.01,
        "null_safe": result["null"]["base_minus_candidate"]["mean"] >= -0.01,
    }
    return {"scope": "generated development data only", "plan": "paper/trackA_dev3_tmz_plan_2026-11-23.md",
            "seed": SEED, "replicates": REPLICATES, "rules": rules, "scenarios": result}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.output.with_suffix(".partial")
    with tmp.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=1, allow_nan=False)
        stream.write("\n")
    tmp.rename(args.output)
    print(json.dumps(result["rules"], indent=1))


if __name__ == "__main__":
    main()
