# Track C development-world seed capacity (2026-11-07)

**Read-only preflight; no seeds reserved or used.** At commit `e7d2cc4`, the registry declares the Stage 2 development scenario block as 20261800–20262099. Checking each integer with `check_free(scenario=(seed, seed), purpose="stage2_dev")` finds 89 registry-free values. The proposed world schedule requires six distinct scenario-class seeds per world: world generation, three patient cohorts, and two signal generations. Thus the current block supports at most `floor(89 / 6) = 14` complete worlds under that schedule.

This is an upper bound, not a selected sample size or a precision calculation. Five additional policy-class seeds per world must also pass the cross-stream registry checks, and all assigned numbers must be distinct. Each predeclared response family, signal anchor, and missing-scan setting needs its own planned allocation. The registry cannot detect undocumented prior use, and checking does not atomically reserve a seed. If a defensible precision target requires more worlds, the solution is a dated, disjoint seed-block amendment and a resource plan before outcomes are viewed; do not lower the target based on results.

Reproduction, from the repository root using Python 3.14:

```powershell
$env:PYTHONPATH='src'
py -3.14 -c 'import trackc_seed_registry as r; lo,hi=r.STAGE2_DEV_SCENARIO_BLOCK; free=[x for x in range(lo,hi+1) if not r.check_free(scenario=(x,x),purpose="stage2_dev")]; print("block",lo,hi,"free",len(free),"world_upper_bound",len(free)//6)'
```

The count depends on the registry version at the time of checking. The outcome protocol, seed reservation, world-count decision, and multi-world scoring remain open. No policy outcome, patient data, or sealed seed was accessed.
