"""Preflight development-world seeds before any benchmark execution."""
from __future__ import annotations

from numbers import Integral

import trackc_seed_registry as registry

SCENARIO_ROLES = ("world", "selection_scenario", "calibration_scenario", "evaluation_scenario", "calibration_signal", "evaluation_signal")
POLICY_ROLES = ("selection_policy", "q_training", "calibration_policy", "evaluation_policy", "schedule_mixture")


def validate_development_schedule(worlds, *, entries=None, purpose="stage2_dev"):
    """Validate a proposed seed schedule; this does not reserve or spend seeds.

    Each world is a mapping with a unique world_id and the exact scalar seed
    roles listed above. All numeric seeds are distinct across worlds and streams.
    purpose: reserved block that scenario seeds must lie in. The trackc_v2 purposes also
    reserve a policy block, so policy seeds are then checked against that block.
    """
    policy_purpose = purpose if str(purpose).startswith("trackc_v2_") else None
    worlds = tuple(worlds)
    if not worlds:
        raise ValueError("at least one world is required")
    entries = registry.ENTRIES if entries is None else entries
    seen_ids, seen_seeds = set(), set()
    roles = set(SCENARIO_ROLES + POLICY_ROLES)
    for world in worlds:
        if not isinstance(world, dict) or set(world) != roles | {"world_id"}:
            raise ValueError("each world must have exactly the declared seed roles")
        world_id = world["world_id"]
        if not isinstance(world_id, str) or not world_id.strip() or world_id != world_id.strip() or world_id in seen_ids:
            raise ValueError("world IDs must be unique canonical nonempty strings")
        seen_ids.add(world_id)
        for role in SCENARIO_ROLES + POLICY_ROLES:
            seed = world[role]
            if isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0:
                raise ValueError(f"{role} seed must be a nonnegative integer")
            seed = int(seed)
            if seed in seen_seeds:
                raise ValueError(f"seed {seed} reused across worlds or streams")
            seen_seeds.add(seed)
            if role in SCENARIO_ROLES:
                problems = registry.check_free(scenario=(seed, seed), purpose=purpose, entries=entries)
            else:
                problems = registry.check_free(policy=(seed, seed), purpose=policy_purpose, entries=entries)
            if problems:
                raise ValueError(f"{role} seed {seed}: {'; '.join(problems)}")
    return {"world_count": len(worlds), "seed_count": len(seen_seeds), "status": "validated_not_reserved"}
