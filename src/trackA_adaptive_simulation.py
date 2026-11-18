import json
from collections import Counter
from pathlib import Path

import numpy as np
from trackA_order_simulation import draw, contrasts, concordance, ensemble, LogisticRegression


def margins(train, test):
    x = train["block_logistic"]
    model = LogisticRegression(C=1, fit_intercept=False).fit(
        np.concatenate([x, -x]), np.repeat([1, 0], len(x)))
    mean = train["soft"].mean(axis=0)
    idx = np.flatnonzero(mean)
    idx = idx[np.argsort(-abs(mean[idx]), kind="stable")[:5]]
    signs = np.sign(mean[idx])
    fit = np.column_stack([model.decision_function(x), train["soft"][:, idx] @ signs])
    predict = np.column_stack([model.decision_function(test["block_logistic"]), test["soft"][:, idx] @ signs])
    scales = np.sqrt(np.mean(fit**2, axis=0))
    return predict / np.where(scales == 0, 1, scales)


def score(margin):
    return float(np.mean((margin > 0) + .5 * (margin == 0)))


def adaptive(train, test):
    indices = np.arange(len(train["soft"]))
    held = np.empty((len(indices), 2))
    for fold in np.array_split(indices, 4):
        keep = ~np.isin(indices, fold)
        held[fold] = margins({k: v[keep] for k, v in train.items()},
                             {k: v[fold] for k, v in train.items()})
    weights = (0, .125, .25, .5)
    weight = max(weights, key=lambda w: score(held @ [1-w, w]) - .04*w)
    return score(margins(train, test) @ [1-weight, weight]), weight


def main():
    rng = np.random.default_rng(20261120)
    results, selected = {}, {}
    for regime in ("diffuse_unequal", "null", "diffuse_reversal", "localized_reversal"):
        rows, weights = [], Counter()
        for _ in range(100):
            train = contrasts(draw(rng, 24, regime))
            test = contrasts(draw(rng, 80, regime))
            value, weight = adaptive(train, test)
            weights[weight] += 1
            rows.append([concordance(train[k], test[k], k) for k in train] +
                        [ensemble(train, test), value])
        results[regime] = dict(zip([*train, "fixed", "adaptive"], np.mean(rows, axis=0).tolist()))
        selected[regime] = dict(weights)
    checks = {"null": abs(results["null"]["adaptive"] - .5) <= .02,
              "shared_vs_fixed": results["diffuse_unequal"]["adaptive"] -
              results["diffuse_unequal"]["fixed"] >= .05}
    for regime, scores in results.items():
        if regime != "null":
            checks[regime] = scores["adaptive"] - max(v for k, v in scores.items() if k != "adaptive") >= -.02
    out = Path(__file__).resolve().parents[1] / "output/trackA_adaptive_simulation/result.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"label": "simulation", "passed": all(checks.values()), "checks": checks,
                               "concordance": results, "selected_weights": selected}, indent=2) + "\n")


if __name__ == "__main__":
    main()
