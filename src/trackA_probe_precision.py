"""Frozen paired classifier comparison with independent assay controls."""
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression


def concordance(train, test):
    labels = np.repeat([1, 0], len(train))
    model = LogisticRegression(C=1, fit_intercept=False).fit(
        np.concatenate([train, -train]), labels)
    margin = model.decision_function(test)
    return float(np.mean((margin > 0) + .5 * (margin == 0)))


def main():
    rng = np.random.default_rng(20261202)
    regimes = {"clean": (.8, 0, 0), "mixed": (.8, .8, -.8),
               "artifact": (0, .8, -.8), "null": (0, 0, 0),
               "blind_mixed": (.8, .8, -.8), "half_mixed": (.8, .8, -.8)}
    results = {}
    for name, (biology, dev_offset, test_offset) in regimes.items():
        shared = .35 if biology else 0
        rows = []
        sensitivity = {"blind_mixed": 0, "half_mixed": .5}.get(name, 1)
        for _ in range(100):
            train = rng.normal(size=(80, 3)) + [shared, shared, biology + dev_offset]
            test = rng.normal(size=(200, 3)) + [shared, shared, biology + test_offset]
            dev_control = rng.normal(sensitivity * dev_offset, size=40).mean()
            test_control = rng.normal(sensitivity * test_offset, size=40).mean()
            omit = concordance(train[:, :2], test[:, :2])
            full = concordance(train, test)
            corrected = concordance(train - [0, 0, dev_control], test - [0, 0, test_control])
            probes = rng.normal(sensitivity, size=(2, 160))
            qualified = bool(np.all(probes.mean(axis=1) - 1.645 / np.sqrt(160) > .60))
            rows.append([omit, full, omit if abs(dev_control) > .4 else full,
                         corrected, corrected if qualified else omit, qualified])
        results[name] = dict(zip(("omit", "full", "train_gate", "corrected", "qualified", "qualification_rate"),
                                 np.mean(rows, axis=0).tolist()))

    def gap(name):
        row = results[name]
        return row["qualified"] - max(v for k, v in row.items() if k not in ("qualified", "qualification_rate"))

    checks = {"mixed": gap("mixed") >= .05, "clean": gap("clean") >= -.02,
              **{k: abs(results[k]["qualified"] - .5) <= .03 for k in ("artifact", "null")}}
    for name in ("blind_mixed", "half_mixed"):
        row = results[name]
        checks[name] = row["qualified"] - row["omit"] >= -.02 and row["qualification_rate"] <= .10

    output = Path(__file__).resolve().parents[1] / "output/trackA_probe_precision/result.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"label": "simulation", "passed": all(checks.values()),
                                  "checks": checks, "concordance": results}, indent=2) + "\n",
                      encoding="utf-8")


if __name__ == "__main__":
    main()
