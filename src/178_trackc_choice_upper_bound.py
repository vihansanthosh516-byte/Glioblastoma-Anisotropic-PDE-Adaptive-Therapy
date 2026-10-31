#!/usr/bin/env python3
"""Post-hoc, clairvoyant policy-choice bound for development screen 177.

The selector sees both policies' realized replicate-averaged outcomes for each
scenario. It cannot be implemented prospectively and does not test a policy.
No sealed results or patient data are read. The comparison covers exactly two
fixed policies chosen for context, not all possible adaptive schedules.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / 'output' / 'trackc_dev_plasticity_sensitivity' / 'dev_sensitivity_results.json'
OUTPUT = ROOT / 'output' / 'trackc_policy_choice_upper_bound' / 'development_bound.json'
POLICIES = ('front_loaded_q1', 'wide_containment')
EXPECTED_SOURCE_SHA256 = 'e0246ecd49e329fbc319057459f62fa47e2940b9d3f648514c72bbc645824e24'


def lower_tenth(values: np.ndarray) -> float:
    ordered = np.sort(np.asarray(values, dtype=float))
    count = 0.1 * ordered.size
    whole = int(np.floor(count))
    value = ordered[:whole].sum()
    if whole < ordered.size:
        value += (count - whole) * ordered[whole]
    return float(value / count)


def main() -> None:
    if OUTPUT.exists():
        raise SystemExit(f'Refusing to overwrite {OUTPUT}')
    raw = SOURCE.read_bytes()
    if hashlib.sha256(raw).hexdigest() != EXPECTED_SOURCE_SHA256:
        raise ValueError('Development screen hash differs from the recorded result')
    source = json.loads(raw)
    n, reps = source['n_scenarios'], source['noise_replicates']
    if n != 36 or reps != 2 or len(source['per_scenario_outcomes']) != 9:
        raise ValueError('Unexpected development-screen dimensions')
    cells = {}
    for cell, policies in source['per_scenario_outcomes'].items():
        means = []
        for policy in POLICIES:
            outcomes = np.asarray(policies[policy]['ttp_days_replicate_by_scenario'], dtype=float)
            if outcomes.shape != (reps, n) or not np.isfinite(outcomes).all():
                raise ValueError(f'Invalid outcomes in {cell}: {policy}')
            means.append(outcomes.mean(axis=0))
        front, containment = means
        oracle = np.maximum(front, containment)
        front_tail, containment_tail = lower_tenth(front), lower_tenth(containment)
        cells[cell] = {
            'front_mean_days': float(front.mean()),
            'containment_mean_days': float(containment.mean()),
            'oracle_mean_days': float(oracle.mean()),
            'oracle_mean_gain_over_better_fixed_days': float(oracle.mean() - max(front.mean(), containment.mean())),
            'front_worst_tenth_days': front_tail,
            'containment_worst_tenth_days': containment_tail,
            'oracle_worst_tenth_days': lower_tenth(oracle),
            'oracle_tail_gain_over_better_fixed_days': float(lower_tenth(oracle) - max(front_tail, containment_tail)),
            'front_wins': int(np.sum(front > containment)),
            'containment_wins': int(np.sum(containment > front)),
            'ties': int(np.sum(front == containment)),
        }
    result = {
        'label': 'SIMULATION; post-hoc clairvoyant two-policy scenario-choice upper bound; not implementable',
        'source_sha256': hashlib.sha256(raw).hexdigest(),
        'policies': list(POLICIES), 'n_scenarios': n, 'n_replicates': reps,
        'selection': 'per scenario, maximum of the two realized replicate-average outcomes',
        'scope': 'upper bound only for choosing one of these two fixed policies per scenario, not for adaptive switching or a broader policy class',
        'results': cells,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(f'Wrote {OUTPUT}; {len(cells)} cells; development only')


if __name__ == '__main__':
    main()
