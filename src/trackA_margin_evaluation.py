"""Evaluate frozen Track A arm scores on shared, matched donor pairs."""

from collections.abc import Mapping, Sequence
from math import comb, isfinite


def _concordance(pair):
    core, margin = pair
    if not isfinite(core) or not isfinite(margin):
        raise ValueError("Every arm must give finite scores for both regions")
    return 1.0 if margin > core else 0.0 if margin < core else 0.5


def _sign_p(wins, losses):
    discordant = wins + losses
    return sum(comb(discordant, k) for k in range(wins, discordant + 1)) / (1 << discordant)


def _holm(raw):
    adjusted = {}
    previous = 0.0
    for index, arm in enumerate(sorted(raw, key=raw.get)):
        previous = max(previous, min(1.0, (len(raw) - index) * raw[arm]))
        adjusted[arm] = previous
    return adjusted


def minimum_possible_adjusted_p(donor_count: int, comparison_count: int) -> float:
    """Best attainable first Holm bound, assuming all donors favor the candidate."""
    if donor_count < 1 or comparison_count < 1:
        raise ValueError("Need at least one donor and comparator")
    return min(1.0, comparison_count / (1 << donor_count))


def evaluate(
    scores: Mapping[str, Mapping[str, tuple[float, float]]],
    candidate: str,
    comparators: Sequence[str],
    *,
    minimum_donors: int = 6,
    alpha: float = 0.05,
) -> dict:
    """Score each arm from mappings of donor ID to (core score, margin score).

    The caller supplies the comparator list from the frozen manifest. Donor and
    study eligibility must be verified before generating the input scores.
    """
    arms = [candidate, *comparators]
    if not comparators or len(set(arms)) != len(arms) or set(scores) != set(arms):
        raise ValueError("Score arms must match the candidate and frozen comparator list")
    if minimum_donors < 1 or not 0 < alpha < 1:
        raise ValueError("Invalid donor minimum or significance threshold")
    donors = sorted(scores[candidate])
    if len(donors) < minimum_donors or any(not donor.strip() for donor in donors):
        raise ValueError("Donor IDs are missing or the donor minimum was missed")
    if any(set(scores[arm]) != set(donors) for arm in arms):
        raise ValueError("All arms must score exactly the same donors")

    concordances = {
        arm: {donor: _concordance(scores[arm][donor]) for donor in donors}
        for arm in arms
    }
    means = {arm: sum(values.values()) / len(donors) for arm, values in concordances.items()}
    comparisons = {}
    for arm in comparators:
        differences = {
            donor: concordances[candidate][donor] - concordances[arm][donor]
            for donor in donors
        }
        wins = sum(value > 0 for value in differences.values())
        losses = sum(value < 0 for value in differences.values())
        comparisons[arm] = {
            "mean_difference": sum(differences.values()) / len(donors),
            "wins": wins,
            "losses": losses,
            "ties": len(donors) - wins - losses,
            "sign_p": _sign_p(wins, losses),
            "per_donor_difference": differences,
        }
    adjusted = _holm({arm: result["sign_p"] for arm, result in comparisons.items()})
    for arm, result in comparisons.items():
        result["holm_p"] = adjusted[arm]
    return {
        "donor_count": len(donors),
        "alpha": alpha,
        "concordance": means,
        "per_donor_concordance": concordances,
        "comparisons": comparisons,
        "minimum_possible_first_holm_p": minimum_possible_adjusted_p(len(donors), len(comparators)),
        "passed": all(
            result["mean_difference"] > 0 and result["holm_p"] <= alpha
            for result in comparisons.values()
        ),
    }
