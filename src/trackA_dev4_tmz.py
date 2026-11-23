"""Track A TMZ development falsification 4: gated one-slope contrast.

Plan: paper/trackA_dev4_tmz_plan_2026-11-23.md (frozen before run). Generated data only.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.preprocessing import StandardScaler

from .trackA_dev2 import (PENALTY, W, admit_direction, augmented, balanced_fit, origin_mae,
                          summarize, within_direction)
from .trackA_tmz_model import Cultures, _fit_shared, _groups, _shuffle_links

DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "output/trackA_dev4_tmz/result.json"
SEED = 202611234
REPLICATES = 20
CONFOUND_GENE = 4
SCENARIOS = {  # between, within, train confound, eval confound, max sisters
    "reversal": (-W, W, 0.0, 0.0, 4),
    "confounded_transfer": (W, W, 0.8, 0.0, 4),
    "strong_confound": (W, W, 1.5, 0.0, 4),
    "confounded_both": (W, W, 0.8, 0.8, 4),
    "plain_transfer": (W, W, 0.0, 0.0, 4),
    "unequal_confound": (W, W, 0.8, 0.0, 8),
    "null": (0 * W, 0 * W, 0.0, 0.0, 4),
}


def cultures(rng, prefix, count, *, between, within, confound, max_size):
    sizes = rng.integers(1, max_size + 1, size=count)
    groups = np.repeat(np.arange(count), sizes)
    common, noise = rng.normal(size=(count, 8)), rng.normal(size=(len(groups), 8))
    z = rng.normal(size=count)
    x = 0.6 * common[groups] + noise
    x[:, CONFOUND_GENE] += 0.8 * z[groups]
    y = (0.6 * common @ between + confound * z)[groups] + noise @ within
    y += rng.normal(scale=0.3, size=count)[groups] + rng.normal(scale=0.3, size=len(groups))
    return Cultures(x, y, [f"{prefix}-{g:03d}" for g in groups],
                    [f"{prefix}-c{i:04d}" for i in range(len(groups))])


def contrast_predict(x, y, origins, links, eval_x):
    intercept, slope = _fit_shared(x, y, origins, links, PENALTY, 1.0)
    return intercept + eval_x @ slope


def admit_contrast(x, y, origins, rng):
    names = rng.permutation(np.unique(origins))
    half = np.isin(origins, names[: len(names) // 2])
    gain = 0.0
    for fit, test in ((half, ~half), (~half, half)):
        _, _, sizes = _groups(origins[fit])
        if not np.any(sizes > 1):
            return False, None
        base = balanced_fit(x[fit], origins[fit], y[fit]).predict(x[test])
        con = contrast_predict(x[fit], y[fit], origins[fit], origins[fit], x[test])
        gain += origin_mae(base, y[test], origins[test]) - origin_mae(con, y[test], origins[test])
    return gain > 0, float(gain)


def replicate(rng, scenario):
    between, within, tc, ec, size = SCENARIOS[scenario]
    train = cultures(rng, "t", 40, between=between, within=within, confound=tc, max_size=size)
    ev = cultures(rng, "e", 100, between=between, within=within, confound=ec, max_size=size)
    _, inv, counts = _groups(train.origins)
    scale = StandardScaler().fit(train.expression, sample_weight=1.0 / (len(counts) * counts[inv]))
    x, ex, y, o = scale.transform(train.expression), scale.transform(ev.expression), train.viability, train.origins
    shuffled = _shuffle_links(o, ["assay"] * len(o), int(rng.integers(1 << 31)))
    preds = {"origin_balanced_ridge": balanced_fit(x, o, y).predict(ex)}
    admitted, gain = admit_contrast(x, y, o, np.random.default_rng(int(rng.integers(1 << 31))))
    preds["contrast_unchecked"] = contrast_predict(x, y, o, o, ex)
    preds["gated_contrast"] = preds["contrast_unchecked"] if admitted else preds["origin_balanced_ridge"]
    preds["shuffled_link_contrast"] = contrast_predict(x, y, o, shuffled, ex)
    ok, _ = admit_direction(x, y, o, np.random.default_rng(int(rng.integers(1 << 31))))
    preds["replicated_direction"] = (augmented(x, o, y, within_direction(x, y, o), ex)
                                     if ok else preds["origin_balanced_ridge"])
    return {"admitted": bool(admitted), "gate_gain": gain,
            "origin_mae": {a: origin_mae(p, ev.viability, ev.origins) for a, p in preds.items()}}


def run():
    result = {}
    for scen, seq in zip(SCENARIOS, np.random.SeedSequence(SEED).spawn(len(SCENARIOS)), strict=True):
        reps = [replicate(np.random.default_rng(s), scen) for s in seq.spawn(REPLICATES)]
        diff = lambda arm: summarize([r["origin_mae"][arm] - r["origin_mae"]["gated_contrast"] for r in reps])
        result[scen] = {
            "mean_origin_mae": {a: float(np.mean([r["origin_mae"][a] for r in reps])) for a in reps[0]["origin_mae"]},
            "base_minus_candidate": diff("origin_balanced_ridge"),
            "shuffled_minus_candidate": diff("shuffled_link_contrast"),
            "admitted_fraction": float(np.mean([r["admitted"] for r in reps])),
            "replicates": reps,
        }
    ok = lambda s, k: result[s][k]["ci95"][0] > 0
    safe = lambda s: result[s]["base_minus_candidate"]["mean"] >= -0.01
    rules = {}
    for s in ("confounded_transfer", "strong_confound", "unequal_confound"):
        rules[f"{s}_beats_base"] = ok(s, "base_minus_candidate")
        rules[f"{s}_beats_shuffle"] = ok(s, "shuffled_minus_candidate")
    for s in ("reversal", "confounded_both", "null"):
        rules[f"{s}_safe"] = safe(s)
    return {"scope": "generated development data only", "plan": "paper/trackA_dev4_tmz_plan_2026-11-23.md",
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
