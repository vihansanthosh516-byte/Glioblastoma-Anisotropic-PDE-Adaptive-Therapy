import sys
sys.path.insert(0, "src")
from resistance_env import make_resistance_env


def run_fixed(env, policy):
    obs, _ = env.reset()
    R = 0.0
    for step in range(env.max_steps):
        day = step + 1
        if policy == "always_combo":
            a = 3
        elif policy == "always_rest":
            a = 0
        elif policy == "stupp":
            a = 3 if 20 <= day < 50 else (1 if (day % 28) < 5 and day >= 50 else 0)
        elif policy == "holiday":
            a = 3 if 20 <= day < 80 else 0
        obs, r, term, trunc, _ = env.step(a)
        R += r
        if term or trunc:
            break
    rf = env.solver.u_r.sum() / max(env.solver.u.sum(), 1e-12)
    dose = sum(1 for t in env.trajectory if t["action"] > 0)
    vol = env.solver.u.sum() * env.solver.dx ** 3
    print(f"{policy:14s}: R={R:9.1f}  vol={vol:9.1f}  res={rf:.3f}  dose={dose}")


env = make_resistance_env(0.025, 0.0012, (16, 16, 16))
for p in ["always_rest", "always_combo", "stupp", "holiday"]:
    run_fixed(env, p)