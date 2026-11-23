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
    rng = np.random.default_rng(20261204)
    regimes = {"clean": (.8, 0, 0), "mixed": (.8, .8, -.8),
               "artifact": (0, .8, -.8), "null": (0, 0, 0),
               "blind_mixed": (.8, .8, -.8), "half_mixed": (.8, .8, -.8), "mismatched_probe": (.8, .8, -.8)}
    results = {}
    for name, (biology, dev_offset, test_offset) in regimes.items():
        shared = .35 if biology else 0
        rows = []
        sensitivity = {"blind_mixed": 0, "half_mixed": .5, "mismatched_probe": 0}.get(name, 1)
        for _ in range(100):
            train = rng.normal(size=(80, 3)) + [shared, shared, biology + dev_offset]
            test = rng.normal(size=(200, 3)) + [shared, shared, biology + test_offset]
            dev_control = rng.normal(sensitivity * dev_offset, size=40).mean()
            test_control = rng.normal(sensitivity * test_offset, size=40).mean()
            omit = concordance(train[:, :2], test[:, :2])
            full = concordance(train, test)
            corrected = concordance(train - [0, 0, dev_control], test - [0, 0, test_control])
            probes = rng.normal(1 if name == "mismatched_probe" else sensitivity, size=(2, 160))
            qualified = bool(np.all(probes.mean(axis=1) - 1.645 / np.sqrt(160) > .60))
            rows.append([omit, full, omit if abs(dev_control) > .4 else full,
                         corrected, corrected if qualified else omit, qualified])
        results[name] = dict(zip(("omit", "full", "train_gate", "corrected", "qualified", "qualification_rate"),
                                 np.mean(rows, axis=0).tolist()))

    def checks_for(model, qualification):
        def gap(name, arms):
            row = results[name]
            return row[model] - max(row[k] for k in arms)

        baselines = ("omit", "full", "train_gate", "corrected")
        checks = {
            "mixed": gap("mixed", baselines[:-1]) >= .05 and gap("mixed", ("corrected",)) >= -.02,
            "clean": gap("clean", baselines) >= -.02,
            **{k: abs(results[k][model] - .5) <= .03 for k in ("artifact", "null")},
        }
        for name in ("blind_mixed", "half_mixed"):
            checks[name] = gap(name, ("omit",)) >= -.02 and results[name][qualification] <= .10
        checks["blind_advantage"] = gap("blind_mixed", ("corrected",)) >= .05
        checks["mismatch"] = gap("mismatched_probe", ("omit",)) >= -.02
        return {"passed": all(checks.values()), "checks": checks}

    for name, row in results.items():
        aligned = name != "mismatched_probe"
        row["audit_rule"] = row["qualified"] if aligned else row["omit"]
        row["audit_qualification_rate"] = row["qualification_rate"] if aligned else 0
    evaluation = {"prior_gate": checks_for("qualified", "qualification_rate"),
                  "audit_rule": checks_for("audit_rule", "audit_qualification_rate")}
    path = Path(__file__).resolve().parents[1] / "output/trackA_probe_mismatch/result.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"label": "simulation", "idealized_alignment_flag": True,
                                "evaluation": evaluation, "concordance": results}, indent=2) + "\n",
                    encoding="utf-8")


if __name__ == "__main__":
    main()
