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
    rng = np.random.default_rng(20261130)
    regimes = {"clean": (.8, 0, 0), "mixed": (.8, .8, -.8),
               "artifact": (0, .8, -.8), "null": (0, 0, 0),
               "blind_mixed": (.8, .8, -.8)}
    results = {}
    for name, (biology, dev_offset, test_offset) in regimes.items():
        shared = .35 if biology else 0
        rows = []
        for _ in range(100):
            train = rng.normal(size=(80, 3)) + [shared, shared, biology + dev_offset]
            test = rng.normal(size=(200, 3)) + [shared, shared, biology + test_offset]
            dev_control = rng.normal(0 if name == "blind_mixed" else dev_offset, size=40).mean()
            test_control = rng.normal(0 if name == "blind_mixed" else test_offset, size=40).mean()
            omit = concordance(train[:, :2], test[:, :2])
            full = concordance(train, test)
            corrected = concordance(train - [0, 0, dev_control], test - [0, 0, test_control])
            rows.append([omit, full, omit if abs(dev_control) > .4 else full, corrected])
        results[name] = dict(zip(("omit", "full", "train_gate", "corrected"),
                                 np.mean(rows, axis=0).tolist()))

    def gap(name):
        row = results[name]
        return row["corrected"] - max(v for k, v in row.items() if k != "corrected")

    checks = {"mixed": gap("mixed") >= .05, "clean": gap("clean") >= -.02,
              **{k: abs(results[k]["corrected"] - .5) <= .03 for k in ("artifact", "null")}}
    broad = results["blind_mixed"]["corrected"] - results["blind_mixed"]["omit"] >= .05
    output = Path(__file__).resolve().parents[1] / "output/trackA_crossover/result.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"label": "simulation", "primary_pass": all(checks.values()),
                                  "broad_pass": all(checks.values()) and broad,
                                  "checks": checks, "blind_probe_pass": broad,
                                  "concordance": results}, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
