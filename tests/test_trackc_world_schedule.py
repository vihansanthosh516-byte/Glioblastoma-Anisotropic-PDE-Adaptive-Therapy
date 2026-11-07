import pytest

import trackc_world_schedule as schedule
import trackc_seed_registry as registry


def _entries():
    return [registry._entry("reserved", "development", "reserved", [(1, 1000)], [], purpose="stage2_dev")]


def _world(name, start):
    values = {role: start + i for i, role in enumerate(schedule.SCENARIO_ROLES)}
    values.update({role: 1001 + start + i for i, role in enumerate(schedule.POLICY_ROLES)})
    return {"world_id": name, **values}


def test_validates_but_does_not_reserve_two_worlds():
    entries = _entries()
    assert schedule.validate_development_schedule([_world("a", 1), _world("b", 20)], entries=entries) == {
        "world_count": 2, "seed_count": 22, "status": "validated_not_reserved"
    }
    assert len(entries) == 1


def test_rejects_reuse_and_registry_conflicts():
    first, second = _world("a", 1), _world("b", 20)
    second["evaluation_policy"] = first["world"]
    with pytest.raises(ValueError, match="reused"):
        schedule.validate_development_schedule([first, second], entries=_entries())
    entries = _entries() + [registry._entry("old", "used", "dev_used", [(1, 1)], [])]
    with pytest.raises(ValueError, match="overlaps"):
        schedule.validate_development_schedule([first], entries=entries)


def test_rejects_unknown_roles_and_outside_development_block():
    world = _world("a", 1)
    world["world"] = 10000
    with pytest.raises(ValueError, match="reserved block"):
        schedule.validate_development_schedule([world], entries=_entries())
    world = _world("a", 1)
    world["extra"] = 77
    with pytest.raises(ValueError, match="exactly"):
        schedule.validate_development_schedule([world], entries=_entries())
