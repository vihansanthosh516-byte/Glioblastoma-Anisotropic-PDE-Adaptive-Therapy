"""Track A TMZ development falsification 6: disjoint-halves selector over established methods.

Plan: paper/trackA_dev6_tmz_plan_2026-11-23.md (frozen before run). Generated data only.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.preprocessing import StandardScaler

from .trackA_dev2 import balanced_fit, origin_mae
from .trackA_dev4_tmz import SCENARIOS, contrast_predict, cultures
from .trackA_dev5_tmz import with_means, within_only
from .trackA_tmz_model import _groups, _shuffle_links

DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "output/trackA_dev6_tmz/result.json"
SEED = 202611236
REPLICATES = 20
COMPARATORS = ("within_only", "origin_balanced_ridge", "mundlak_ridge", "shuffled_link_contrast")


def methods(x, y, o, ex, eo):
    base = balanced_fit(x, o, y).predict(ex)
    wo = within_only(x, y, o, ex)
    return {
        "within_only": base if wo is None else wo,
        "origin_balanced_ridge": base,
        "mundlak_ridge": balanced_fit(with_means(x, o), o, y).predict(with_means(ex, eo)),
    }


def select(x, y, o, rng):
    names = rng.permutation(np.unique(o))
    half = np.isin(o, names[: len(names) // 2])
    totals = {}
    for fit, test in ((half, ~half), (~half, half)):
        for arm, p in methods(x[fit], y[fit], o[fit], x[test], o[test]).items():
            totals[arm] = totals.get(arm, 0.0) + origin_mae(p, y[test], o[test])
    return min(totals, key=totals.get)


def replicate(rng, scenario):
    between, within, tc, ec, size = SCENARIOS[scenario]
    train = cultures(rng, "t", 40, between=between, within=within, confound=tc, max_size=size)
    ev = cultures(rng, "e", 100, between=between, within=within, confound=ec, max_size=size)
    _, inv, counts = _groups(train.origins)
    scale = StandardScaler().fit(train.expression, sample_weight=1.0 / (len(counts) * counts[inv]))
    x, ex, y, o = scale.transform(train.expression), scale.transform(ev.expression), train.viability, train.origins
    shuffled = _shuffle_links(o, ["assay"] * len(o), int(rng.integers(1 << 31)))
    preds = methods(x, y, o, ex, ev.origins)
    preds["shuffled_link_contrast"] = contrast_predict(x, y, o, shuffled, ex)
    picked = select(x, y, o, np.random.default_rng(int(rng.integers(1 << 31))))
    preds["halves_selector"] = preds[picked]
    return {"picked": picked, "origin_mae": {a: origin_mae(p, ev.viability, ev.origins) for a, p in preds.items()}}


def run():
    result = {}
    for scen, seq in zip(SCENARIOS, np.random.SeedSequence(SEED).spawn(len(SCENARIOS)), strict=True):
        reps = [replicate(np.random.default_rng(s), scen) for s in seq.spawn(REPLICATES)]
        means = {a: float(np.mean([r["origin_mae"][a] for r in reps])) for a in reps[0]["origin_mae"]}
        best = min(means[a] for a in COMPARATORS)
        result[scen] = {"mean_origin_mae": means,
                        "regret": {a: means[a] - best for a in (*COMPARATORS, "halves_selector")},
                        "picked": {a: sum(r["picked"] == a for r in reps) for a in COMPARATORS[:3]},
                        "replicates": reps}
    max_regret = {a: max(v["regret"][a] for v in result.values()) for a in (*COMPARATORS, "halves_selector")}
    rules = {"selector_lowest_max_regret": all(max_regret["halves_selector"] < max_regret[a] for a in COMPARATORS)}
    return {"scope": "generated development data only", "plan": "paper/trackA_dev6_tmz_plan_2026-11-23.md",
            "seed": SEED, "replicates": REPLICATES, "rules": rules, "max_regret": max_regret, "scenarios": result}


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
    print(json.dumps({"rules": result["rules"], "max_regret": result["max_regret"]}, indent=1))


if __name__ == "__main__":
    main()
