"""Stress test the frozen Track A assay-control gate."""
import json
from itertools import product
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression


def main():
    rng = np.random.default_rng(20261129)
    results = {}
    for size, sensitivity in product((8, 40), (0, .5, 1)):
        regimes = {}
        for regime in ("spatial_signal", "artifact"):
            rows = []
            artifact = regime == "artifact"
            for _ in range(100):
                train = rng.normal(size=(80, 3)) + ([0, 0, .8] if artifact else [.35, .35, .8])
                test = rng.normal(size=(200, 3)) + ([0, 0, -.8] if artifact else [.35, .35, .8])
                controls = rng.normal(.8 * sensitivity if artifact else 0, size=size)
                include = abs(controls.mean()) <= .4
                scores = []
                for width in (2, 3):
                    x = train[:, :width]
                    model = LogisticRegression(C=1, fit_intercept=False).fit(
                        np.concatenate([x, -x]), np.repeat([1, 0], len(x)))
                    margin = model.decision_function(test[:, :width])
                    scores.append(np.mean((margin > 0) + .5 * (margin == 0)))
                rows.append([*scores, scores[int(include)], float(include)])
            regimes[regime] = dict(zip(("omit", "include", "gate", "inclusion_rate"),
                                      np.mean(rows, axis=0).tolist()))
        signal, artifact = regimes.values()
        passed = signal["gate"] - signal["omit"] >= .05 and abs(artifact["gate"] - .5) <= .03
        results[f"{size}_{sensitivity}"] = {"passed": passed, **regimes}
    out = Path(__file__).resolve().parents[1] / "output/trackA_assay_transport/result.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"label": "simulation", "passed": all(r["passed"] for r in results.values()),
                               "settings": results}, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
