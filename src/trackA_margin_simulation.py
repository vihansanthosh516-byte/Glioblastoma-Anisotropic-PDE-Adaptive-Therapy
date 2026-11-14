"""Run the prospectively specified Track A synthetic mechanics screen."""

import json
from pathlib import Path

import numpy as np

from trackA_margin_pairs import fit_margin_pairs


def simulate(rng, donors, scenario):
    offsets = rng.normal(size=(donors, 1, 24))
    core = offsets + rng.normal(scale=0.55, size=(donors, 4, 24))
    margin = offsets + rng.normal(scale=0.55, size=core.shape)
    if scenario == "shared_plus_specific":
        margin[..., :2] += [0.6, -0.6]
    if scenario != "null":
        margin[:, 0, 2:4] += [3, -3]
    return core, margin


def concordance(core, margin):
    return float(((margin > core) + 0.5 * (margin == core)).mean())


def rank_arm(train, test, genes, *, weakest_state=True):
    model = fit_margin_pairs(*train, genes, 1, weakest_state=weakest_state)
    return {
        "concordance": concordance(*(model.score(values) for values in test)),
        "selected_pair": [
            [genes[first], genes[second]]
            for first, second in zip(model.first, model.second)
        ],
    }


def run():
    rng = np.random.default_rng(7461)
    genes = [f"G{i:02}" for i in range(24)]
    scenarios = {}
    for scenario in ("shared_plus_specific", "specific_only", "null"):
        train = simulate(rng, 96, scenario)
        test = simulate(rng, 256, scenario)
        pooled_train = tuple(x.mean(axis=1, keepdims=True) for x in train)
        pooled_test = tuple(x.mean(axis=1, keepdims=True) for x in test)
        centroid = (pooled_train[1] - pooled_train[0]).mean(axis=(0, 1))
        scenarios[scenario] = {
            "minimum_state": rank_arm(train, test, genes),
            "mean_state": rank_arm(train, test, genes, weakest_state=False),
            "pooled_tsp": rank_arm(pooled_train, pooled_test, genes),
            "centroid": {
                "concordance": concordance(*(x[:, 0] @ centroid for x in pooled_test)),
            },
        }

    shared = scenarios["shared_plus_specific"]
    passed = all(
        shared["minimum_state"]["concordance"] > shared[arm]["concordance"]
        for arm in ("mean_state", "pooled_tsp")
    ) and all(0.40 <= arm["concordance"] <= 0.60 for arm in scenarios["null"].values())
    return {"label": "simulation", "scenarios": scenarios, "passed": passed}


if __name__ == "__main__":
    destination = Path(__file__).resolve().parents[1] / "output/trackA_margin_simulation/result.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("x", encoding="utf-8") as stream:
        json.dump(run(), stream, indent=2)
        stream.write("\n")
