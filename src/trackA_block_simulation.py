"""Frozen tissue block voting simulation."""
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression


def counts(rng, n, regime, depths):
    log = (rng.normal(size=(n, 1, 1, 24))
           + rng.normal(scale=.5, size=(n, 2, 1, 24))
           + rng.normal(scale=.6, size=(n, 2, 4, 24)))
    if regime.startswith("diffuse"):
        log[:, 1, :, :6] += 1
    elif regime == "localized":
        log[np.arange(n), 1, rng.integers(4, size=n), :6] += 1
    abundance = np.exp(log)
    return rng.multinomial(depths, abundance / abundance.sum(axis=-1, keepdims=True))


def contrasts(x):
    pooled = x.sum(axis=2)
    votes = np.sign(x[..., ::2] - x[..., 1::2]).mean(axis=2)
    pooled_votes = np.sign(pooled[..., ::2] - pooled[..., 1::2])
    block_log = np.log1p(10000 * x / x.sum(axis=-1, keepdims=True)).mean(axis=2)
    pooled_log = np.log1p(10000 * pooled / pooled.sum(axis=-1, keepdims=True))
    return {k: v[:, 1] - v[:, 0] for k, v in {
        "soft": votes, "hard": np.sign(votes), "pooled_vote": pooled_votes,
        "block_logistic": block_log, "pooled_logistic": pooled_log,
    }.items()}


def concordance(train, test, method):
    if method.endswith("logistic"):
        x = np.concatenate([train, -train])
        y = np.repeat([1, 0], len(train))
        model = LogisticRegression(C=1, fit_intercept=False).fit(x, y)
        scores = model.decision_function(test)
    else:
        mean = train.mean(axis=0)
        selected = np.flatnonzero(mean)
        chosen = selected[np.argsort(-abs(mean[selected]), kind="stable")[:5]]
        scores = test[:, chosen] @ np.sign(mean[chosen])
    return float(np.mean((scores > 0) + .5 * (scores == 0)))


def main():
    rng = np.random.default_rng(20261118)
    results = {}
    for regime in ("diffuse_unequal", "localized", "null", "diffuse_equal"):
        depths = [500, 500, 500, 500 if regime == "diffuse_equal" else 10000]
        runs = []
        for _ in range(100):
            train = contrasts(counts(rng, 24, regime, depths))
            test = contrasts(counts(rng, 80, regime, depths))
            runs.append([concordance(train[k], test[k], k) for k in train])
        results[regime] = dict(zip(train, np.mean(runs, axis=0).tolist()))

    def gap(regime):
        scores = results[regime]
        return scores["soft"] - max(v for k, v in scores.items() if k != "soft")

    checks = {
        "diffuse_superiority": gap("diffuse_unequal") >= .02,
        "null": abs(results["null"]["soft"] - .5) <= .02,
        "localized": gap("localized") >= -.05,
        "equal_depth": gap("diffuse_equal") >= -.05,
    }
    output = Path(__file__).resolve().parents[1] / "output/trackA_block_simulation_20261117/result.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"label": "simulation", "passed": all(checks.values()),
                                  "checks": checks, "mean_concordance": results}, indent=2) + "\n",
                      encoding="utf-8")


if __name__ == "__main__":
    main()
