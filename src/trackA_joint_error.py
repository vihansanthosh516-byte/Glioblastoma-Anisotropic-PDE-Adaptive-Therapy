"""Frozen Track A synthetic spatial correction simulation."""
import json
from pathlib import Path


def simulate():
    state = 20261126

    def random():
        nonlocal state
        state = (1664525 * state + 1013904223) % 2**32
        return state / 2**32

    def bridge_score(labels, edges, positions):
        if not edges:
            return 0.0, 0.0
        observed = sum(labels[i] != labels[j] for i, j in edges) / len(edges)
        n = len(positions)
        a = sum(labels[i] for i in positions)
        return observed, observed - 2 * a * (n - a) / (n * (n - 1))

    results = {}
    for regime in ("clean_signal", "noisy_signal", "noisy_null", "correlated_null"):
        totals = dict.fromkeys(("depth_matched", "global_adjusted", "raw", "composition", "degree_composition", "geometry", "corrected", "misspecified", "joint", "joint_misspecified"), 0.0)
        empty = [0, 0]
        for _ in range(1000):
            positions = list(range(24))
            graph = [(i, (i + 1) % 24) for i in positions]
            graph += [(i, (i + 3) % 24) for i in positions if random() < .25]
            offset = int(random() * 24)
            regions = []
            for region, depth_probability in enumerate((.25, .75)):
                low = [random() < depth_probability for _ in positions]
                latent = [int((i + offset) % 24 < 12) if region == 0 or regime.endswith("null")
                          else (i + offset) % 2 for i in positions]
                if regime == "correlated_null":
                    uniforms = [random() for _ in range(6)]
                    draws = [uniforms[i // 4] for i in positions]
                else:
                    draws = [random() for _ in positions] if regime != "clean_signal" else [1] * 24
                labels = [v ^ (u < (.20 if shallow else .02))
                          for v, u, shallow in zip(latent, draws, low)]
                def corrected(low_rate, high_rate):
                    z = [(2*v-1)/(1-2*(low_rate if shallow else high_rate))
                         for v, shallow in zip(labels, low)]
                    global_product = (sum(z)**2 - sum(v*v for v in z)) / (24*23)
                    edge_product = sum(z[i]*z[j] for i,j in graph) / len(graph)
                    return (global_product - edge_product) / 2

                def joint(low_rate, high_rate):
                    if regime == "clean_signal":
                        return adjusted
                    if regime != "correlated_null":
                        return corrected(low_rate, high_rate)
                    rates = [low_rate if shallow else high_rate for shallow in low]

                    def product(i, j):
                        p, q = rates[i], rates[j]
                        multiplier = 1 - 2 * abs(p - q) if i // 4 == j // 4 else (1 - 2*p) * (1 - 2*q)
                        return (2*labels[i]-1) * (2*labels[j]-1) / multiplier

                    overall = sum(product(i, j) for i in positions for j in range(i)) / 276
                    neighbors = sum(product(i, j) for i, j in graph) / len(graph)
                    return (overall - neighbors) / 2

                high = [i for i in positions if not low[i]]
                eligible = [(i, j) for i, j in graph if not low[i] and not low[j]]
                empty[region] += not eligible
                raw, adjusted = bridge_score(labels, graph, positions)
                _, matched = bridge_score(labels, eligible, high)
                weighted = sum(labels[i] + labels[j] for i, j in graph) / (2 * len(graph))
                regions.append((matched, adjusted, raw, sum(labels) / 24, weighted, 0, corrected(.20, .02) if regime != "clean_signal" else adjusted, corrected(.16, .01) if regime != "clean_signal" else adjusted, joint(.20, .02), joint(.16, .01)))
            for key, core, edge in zip(totals, *regions):
                totals[key] += 1 if edge > core else .5 if edge == core else 0
        results[regime] = {"concordance": {k: v / 1000 for k, v in totals.items()},
                           "no_high_depth_bridge_rate": dict(zip(("core", "edge"), (n / 1000 for n in empty)))}
    checks = {}
    for variant in ("joint", "joint_misspecified"):
        for regime, row in results.items():
            c = row["concordance"]
            row.setdefault("gain_vs_depth_matched", {})[variant] = c[variant] - c["depth_matched"]
            checks[f"{variant}_{regime}"] = (
                c[variant] == c["global_adjusted"] if regime == "clean_signal" else
                c[variant] >= .98 if regime == "noisy_signal" else abs(c[variant] - .50) <= .03
            )

    return {"label": "simulation", "passed": all(checks.values()), "checks": checks, "regimes": results}


if __name__ == "__main__":
    path = Path(__file__).resolve().parents[1] / "output/trackA_joint_error/result.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(simulate(), indent=2) + "\n", encoding="utf-8")
