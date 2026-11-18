"""Run the frozen Track A support threshold simulation."""
import json
from pathlib import Path
import numpy as np


def votes(rng, patients, probability):
    core, edge = rng.choice([-1, 1], size=(2, patients, 30))
    response = rng.random(patients) < probability
    core[response, :5] = -1
    edge[response, :5] = 1
    return edge - core


def score(train, test, method):
    mean = train.mean(axis=0)
    positive = (train > 0).mean(axis=0)
    negative = (train < 0).mean(axis=0)
    direction = np.where(mean >= 0, 1, -1)
    if method != "mean":
        direction = np.where(positive > negative, 1, np.where(negative > positive, -1, direction))
    support = np.where(direction > 0, positive, negative)
    oriented_mean = mean * direction
    indices = np.arange(train.shape[1])
    if method == "mean":
        order = np.lexsort((indices, -oriented_mean))
        selected = order[oriented_mean[order] > 0][:5]
    else:
        order = np.lexsort((indices, -oriented_mean, -support))
        cutoff = 1 if method == "veto" else .70
        selected = order[support[order] >= cutoff][:5]
    differences = test[:, selected] @ direction[selected]
    return float(np.mean((differences > 0) + .5 * (differences == 0))), not len(selected)


def main():
    rng = np.random.default_rng(20261117)
    methods = ("support", "veto", "mean")
    results = {}
    for probability in (0, .45, .75):
        runs = {method: [] for method in methods}
        for _ in range(500):
            train, test = votes(rng, 10, probability), votes(rng, 100, probability)
            for method in methods:
                runs[method].append(score(train, test, method))
        results[str(probability)] = {
            method: dict(zip(("mean_concordance", "empty_fraction"), np.mean(runs[method], axis=0).tolist()))
            for method in methods
        }
    passed = abs(results["0"]["support"]["mean_concordance"] - .5) <= .02 and all(
        results[regime]["support"]["mean_concordance"] - results[regime][method]["mean_concordance"] >= .02
        for regime in ("0.45", "0.75") for method in ("veto", "mean")
    )
    output = Path(__file__).resolve().parents[1] / "output/trackA_support_simulation_20261117/result.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"label": "simulation", "passed": passed, "regimes": results}, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
