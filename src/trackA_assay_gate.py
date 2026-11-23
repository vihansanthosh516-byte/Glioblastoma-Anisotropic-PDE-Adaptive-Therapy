"""Run the frozen Track A assay-control simulation."""
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression


def main():
    rng = np.random.default_rng(20261128)
    regimes = {
        "spatial_signal": ([.35, .35, .80], [.35, .35, .80]),
        "expression_only": ([.35, .35, 0], [.35, .35, 0]),
        "null": ([0, 0, 0], [0, 0, 0]),
        "artifact": ([0, 0, .80], [0, 0, -.80]),
    }
    results, inclusion = {}, {}
    names = ("composition", "plus_expression", "plus_spatial")
    regimes["mixed"] = ([.35, .35, 1.60], [.35, .35, 0])
    for regime, (train_shift, test_shift) in regimes.items():
        full_count = 0
        scores = {name: [] for name in names}
        for _ in range(100):
            train = rng.normal(size=(80, 3)) + train_shift
            test = rng.normal(size=(200, 3)) + test_shift
            control = rng.normal(loc=.80 if regime in {"artifact", "mixed"} else 0, size=40)
            full = abs(control.mean()) <= .40
            full_count += full
            labels = np.repeat([1, 0], len(train))
            for size, name in enumerate(names, 1):
                features = train[:, :size]
                model = LogisticRegression(C=1, fit_intercept=False)
                model.fit(np.concatenate([features, -features]), labels)
                margins = model.decision_function(test[:, :size])
                scores[name].append(float(np.mean((margins > 0) + .5 * (margins == 0))))
            scores.setdefault("gated", []).append(scores["plus_spatial" if full else "plus_expression"][-1])
        inclusion[regime] = float(full_count / 100)
        results[regime] = {name: float(np.mean(values)) for name, values in scores.items()}

    def gap(regime):
        row = results[regime]
        return row["gated"] - max(row[name] for name in names[:-1])

    checks = {
        "spatial_gain": gap("spatial_signal") >= .05,
        "expression_tolerance": gap("expression_only") >= -.02,
        "null": abs(results["null"]["gated"] - .5) <= .03,
        "artifact": abs(results["artifact"]["gated"] - .5) <= .03,
    }
    checks["mixed"] = results["mixed"]["gated"] - results["mixed"]["plus_expression"] >= -.02
    output = Path(__file__).resolve().parents[1] / "output/trackA_assay_gate/result.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"label": "simulation", "passed": all(checks.values()),
                                  "checks": checks, "concordance": results, "spatial_inclusion": inclusion}, indent=2) + "\n",
                      encoding="utf-8")


if __name__ == "__main__":
    main()
