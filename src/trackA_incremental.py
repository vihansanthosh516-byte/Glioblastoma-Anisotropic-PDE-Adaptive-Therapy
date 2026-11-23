"""Run the frozen Track A nested-feature simulation."""
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression


def main():
    rng = np.random.default_rng(20261127)
    regimes = {
        "spatial_signal": ([.35, .35, .80], [.35, .35, .80]),
        "expression_only": ([.35, .35, 0], [.35, .35, 0]),
        "null": ([0, 0, 0], [0, 0, 0]),
        "artifact": ([0, 0, .80], [0, 0, -.80]),
    }
    results = {}
    names = ("composition", "plus_expression", "plus_spatial")
    for regime, (train_shift, test_shift) in regimes.items():
        scores = {name: [] for name in names}
        for _ in range(100):
            train = rng.normal(size=(80, 3)) + train_shift
            test = rng.normal(size=(200, 3)) + test_shift
            labels = np.repeat([1, 0], len(train))
            for size, name in enumerate(names, 1):
                features = train[:, :size]
                model = LogisticRegression(C=1, fit_intercept=False)
                model.fit(np.concatenate([features, -features]), labels)
                margins = model.decision_function(test[:, :size])
                scores[name].append(float(np.mean((margins > 0) + .5 * (margins == 0))))
        results[regime] = {name: float(np.mean(values)) for name, values in scores.items()}

    def gap(regime):
        row = results[regime]
        return row["plus_spatial"] - max(row[name] for name in names[:-1])

    checks = {
        "spatial_gain": gap("spatial_signal") >= .05,
        "expression_tolerance": gap("expression_only") >= -.02,
        "null": abs(results["null"]["plus_spatial"] - .5) <= .03,
        "artifact": abs(results["artifact"]["plus_spatial"] - .5) <= .03,
    }
    output = Path(__file__).resolve().parents[1] / "output/trackA_incremental/result.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"label": "simulation", "passed": all(checks.values()),
                                  "checks": checks, "concordance": results}, indent=2) + "\n",
                      encoding="utf-8")


if __name__ == "__main__":
    main()
