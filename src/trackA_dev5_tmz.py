"""Track A TMZ development falsification 5: add Mundlak and within-only comparators.

Plan: paper/trackA_dev5_tmz_plan_2026-11-23.md (frozen before run). Generated data only.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.preprocessing import StandardScaler

from .trackA_dev2 import balanced_fit, origin_mae, summarize, within_direction
from .trackA_dev4_tmz import SCENARIOS, admit_contrast, contrast_predict, cultures
from .trackA_tmz_model import _groups, _means, _shuffle_links

DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "output/trackA_dev5_tmz/result.json"
SEED = 202611235
REPLICATES = 20
COMPARATORS = ("origin_balanced_ridge", "mundlak_ridge", "within_only", "shuffled_link_contrast")


def with_means(x, origins):
    _, inv, counts = _groups(origins)
    return np.column_stack((x, _means(x, inv, counts)[inv]))


def within_only(x, y, origins, eval_x):
    beta = within_direction(x, y, origins)
    if beta is None:
        return None
    _, inv, counts = _groups(origins)
    score, mean_y = _means(x @ beta, inv, counts), _means(y, inv, counts)
    slope, intercept = np.polyfit(score, mean_y, 1)
    return intercept + slope * (eval_x @ beta)


def replicate(rng, scenario):
    between, within, tc, ec, size = SCENARIOS[scenario]
    train = cultures(rng, "t", 40, between=between, within=within, confound=tc, max_size=size)
    ev = cultures(rng, "e", 100, between=between, within=within, confound=ec, max_size=size)
    _, inv, counts = _groups(train.origins)
    scale = StandardScaler().fit(train.expression, sample_weight=1.0 / (len(counts) * counts[inv]))
    x, ex, y, o = scale.transform(train.expression), scale.transform(ev.expression), train.viability, train.origins
    shuffled = _shuffle_links(o, ["assay"] * len(o), int(rng.integers(1 << 31)))
    base = balanced_fit(x, o, y).predict(ex)
    preds = {"origin_balanced_ridge": base}
    admitted, gain = admit_contrast(x, y, o, np.random.default_rng(int(rng.integers(1 << 31))))
    preds["contrast_unchecked"] = contrast_predict(x, y, o, o, ex)
    preds["gated_contrast"] = preds["contrast_unchecked"] if admitted else base
    preds["shuffled_link_contrast"] = contrast_predict(x, y, o, shuffled, ex)
    preds["mundlak_ridge"] = balanced_fit(with_means(x, o), o, y).predict(with_means(ex, ev.origins))
    wo = within_only(x, y, o, ex)
    preds["within_only"] = base if wo is None else wo
    return {"admitted": bool(admitted), "gate_gain": gain,
            "origin_mae": {a: origin_mae(p, ev.viability, ev.origins) for a, p in preds.items()}}


def run():
    result = {}
    for scen, seq in zip(SCENARIOS, np.random.SeedSequence(SEED).spawn(len(SCENARIOS)), strict=True):
        reps = [replicate(np.random.default_rng(s), scen) for s in seq.spawn(REPLICATES)]
        means = {a: float(np.mean([r["origin_mae"][a] for r in reps])) for a in reps[0]["origin_mae"]}
        best = min(COMPARATORS, key=means.get)
        result[scen] = {
            "mean_origin_mae": means, "strongest_comparator": best,
            "strongest_minus_candidate": summarize(
                [r["origin_mae"][best] - r["origin_mae"]["gated_contrast"] for r in reps]),
            "admitted_fraction": float(np.mean([r["admitted"] for r in reps])),
            "replicates": reps,
        }
    rules = {}
    for s in ("confounded_transfer", "strong_confound", "unequal_confound"):
        rules[f"{s}_beats_strongest"] = result[s]["strongest_minus_candidate"]["ci95"][0] > 0
    for s in ("reversal", "confounded_both", "null"):
        r = result[s]
        rules[f"{s}_safe"] = r["mean_origin_mae"]["gated_contrast"] <= r["mean_origin_mae"][r["strongest_comparator"]] + 0.01
    return {"scope": "generated development data only", "plan": "paper/trackA_dev5_tmz_plan_2026-11-23.md",
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
