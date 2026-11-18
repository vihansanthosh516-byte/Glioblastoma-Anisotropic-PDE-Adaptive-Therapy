import json
from pathlib import Path
import numpy as np
from trackA_block_simulation import counts, contrasts, concordance, LogisticRegression


def draw(rng, size, regime):
    depths = [500, 500, 500, 10000]
    if regime in ("diffuse_unequal", "null"):
        return counts(rng, size, regime, depths)
    log = (rng.normal(size=(size, 1, 1, 24))
           + rng.normal(scale=.5, size=(size, 2, 1, 24))
           + rng.normal(scale=.6, size=(size, 2, 4, 24)))
    shift = np.tile([1, -1], 3)
    if regime == "diffuse_reversal":
        log[:, 1, :, :6] += shift
    else:
        log[np.arange(size), 1, rng.integers(4, size=size), :6] += shift
    abundance = np.exp(log)
    return rng.multinomial(depths, abundance / abundance.sum(axis=-1, keepdims=True))


def ensemble(train, test):
    x = train["block_logistic"]
    model = LogisticRegression(C=1, fit_intercept=False).fit(
        np.concatenate([x, -x]), np.repeat([1, 0], len(x)))
    train_logit = model.decision_function(x)
    test_logit = model.decision_function(test["block_logistic"])
    mean = train["soft"].mean(axis=0)
    idx = np.flatnonzero(mean)
    idx = idx[np.argsort(-abs(mean[idx]), kind="stable")[:5]]
    weights = np.sign(mean[idx])
    train_vote = train["soft"][:, idx] @ weights
    test_vote = test["soft"][:, idx] @ weights
    scores = test_logit / (np.sqrt(np.mean(train_logit**2)) or 1)
    scores += test_vote / (np.sqrt(np.mean(train_vote**2)) or 1)
    return float(np.mean((scores > 0) + .5 * (scores == 0)))


def main():
    rng = np.random.default_rng(20261119)
    results = {}
    for regime in ("diffuse_unequal", "null", "diffuse_reversal", "localized_reversal"):
        rows = []
        for _ in range(100):
            train = contrasts(draw(rng, 24, regime))
            test = contrasts(draw(rng, 80, regime))
            rows.append([concordance(train[k], test[k], k) for k in train] + [ensemble(train, test)])
        results[regime] = dict(zip([*train, "ensemble"], np.mean(rows, axis=0).tolist()))
    gaps = {k: d["ensemble"] - max(v for m, v in d.items() if m != "ensemble")
            for k, d in results.items()}
    checks = {"null": abs(results["null"]["ensemble"] - .5) <= .02,
              "shared_shift": gaps["diffuse_unequal"] >= -.02,
              "diffuse_reversal": gaps["diffuse_reversal"] >= .01,
              "localized_reversal": gaps["localized_reversal"] >= .01}
    path = Path(__file__).resolve().parents[1] / "output/trackA_order_simulation/result.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"label": "simulation", "passed": all(checks.values()), "checks": checks,
                                "concordance": results, "gaps": gaps}, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
