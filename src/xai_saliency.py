#!/usr/bin/env python3
"""
Explainable AI (XAI) Saliency Maps for RL (Proposal 3)
======================================================
Produces policy-saliency maps over the spatial tumour field so clinicians
can see where the RL agent "sees" risky growth that triggers a dose.

Mathematics
-----------
For a policy pi(.|s), the saliency of the 3D tumour density field u(x,y,z)
at a given decision step is:

    S(u) = | d / d u  log prob(a* | s(u)) |

where a* = argmax_a pi(a|s) is the chosen action, and s(u) is the
environment observation constructed from the tumour tensor plus biosensor
readings. We compute `S(u)` via `torch.autograd.grad(log_prob, u_tensor)`,
then take absolute value, normalize, optionally Gaussian-smooth, and
overlay the heatmap on the tumour density slice (2D) or render a 3D volume
(Plotly isosurface).

A small differentiable adapter lets any of these environments feed the
saliency computation:

    - src/58_rl_adaptive_steering.py   ->  GbmTherapyEnv  (Discrete policy)
    - src/rl/chronotherapy_env.py      ->  ChronotherapyEnv (Box policy)

For the discrete-action envs returning a Categorical policy, the
log-probability of the chosen action is differentiable through both the
policy network parameters and the *input tensor* u that constructs the
scalar observation `[norm_vol, u_max, ...]`.

Usage
-----
    from src.xai_saliency import compute_saliency
    sal = compute_saliency(env, policy, target="action")  # returns (Z,Y,X) saliency
    saliency_to_html(env, policy, day=10, patient_id="PAT_01",
                     output="output/saliency_day_10.html")
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple

import numpy as np

# Optional torch / plotly deps; the helpers degrade gracefully
try:
    import torch
    import torch.nn.functional as F
    _HAVE_TORCH = True
except ImportError:
    torch = None  # type: ignore
    _HAVE_TORCH = False

try:
    import scipy.ndimage as _ndi
except ImportError:
    _ndi = None  # type: ignore

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAVE_MPL = True
except ImportError:
    plt = None  # type: ignore
    _HAVE_MPL = False

# The torch / gymnasium RL stack live in src/58_rl_adaptive_steering.py and
# src/rl/chronotherapy_env.py. We resolve them lazily so importing this
# module does not blow up when those deps are absent.
REPO_ROOT = Path(__file__).resolve().parent.parent


def _load_module_by_path(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _import_rl_module() -> Optional[Any]:
    """Load src/58_rl_adaptive_steering.py (filename starts with digit)."""
    p = REPO_ROOT / "src" / "58_rl_adaptive_steering.py"
    if not p.exists():
        return None
    try:
        return _load_module_by_path("rl_adaptive_steering", p)
    except Exception as e:  # pragma: no cover - import guarded
        print(f"[xai_saliency] could not import rl_adaptive_steering: {e}")
        return None


# --------------------------------------------------------------------------- #
# Differentiable observation adapter
# --------------------------------------------------------------------------- #
def _make_obs_from_u(u_tensor: "torch.Tensor", env_attrs: Dict[str, float]) -> "torch.Tensor":
    """Reconstruct the 5-D observation of GbmTherapyEnv from u + toxicities.

    Matches src/58_rl_adaptive_steering.py FastPDESolver.get_observation, but
    is *differentiable* in u_tensor.
    """
    vol = u_tensor.sum() * env_attrs.get("dx_cubed", 1.0)
    u_max = u_tensor.amax()
    day_frac = env_attrs.get("day_frac", 0.0)
    norm_vol = torch.clamp(vol / max(env_attrs.get("initial_volume", 1e-6),
                                    1e-6), 0.0, 1.0)
    chemo = torch.tensor(float(env_attrs.get("chemo_tox", 0.0)))
    rad = torch.tensor(float(env_attrs.get("rad_tox", 0.0)))
    return torch.stack([
        norm_vol, torch.clamp(u_max, 0.0, 1.0),
        torch.tensor(day_frac),
        torch.clamp(chemo, 0.0, 1.0),
        torch.clamp(rad, 0.0, 1.0),
    ])


# --------------------------------------------------------------------------- #
# Core saliency computation
# --------------------------------------------------------------------------- #
def compute_saliency(
    env: Any,
    policy: Any,
    target: str = "action",
    smooth_sigma: float = 1.0,
) -> np.ndarray:
    """Return a (Z, Y, X) saliency map over the env's tumour density field.

    Parameters
    ----------
    env : GbmTherapyEnv (or compatible). Must expose `.solver.u` as the 3D
        tumour tensor on a numpy array. The current env state determines the
        saliency, i.e. you should `step()` it to the desired day first.
    policy : a torch.nn.Module whose forward(x) returns a distribution with
        a `.log_prob(action)` and `.probs` attribute (Categorical / Normal).
    target : "action" (chosen action) or "value" (max prob as a value-ish proxy)
    """
    if not _HAVE_TORCH:
        raise ImportError("torch is required for saliency computation")

    u_np = getattr(getattr(env, "solver", env), "u", None)
    if u_np is None:
        raise AttributeError("env.solver.u (tumour field) not found")
    u_gpu = torch.tensor(np.asarray(u_np, dtype=np.float32),
                         requires_grad=True)

    solver = getattr(env, "solver", env)
    attrs = dict(
        dx_cubed=getattr(solver, "dx", 1.0) ** 3,
        initial_volume=float(getattr(solver, "initial_volume", 1.0) or 1e-6),
        day_frac=float(getattr(solver, "step_count", 0)) /
                  max(getattr(env, "max_steps", 90), 1),
        chemo_tox=float(getattr(solver, "chemo_tox", 0.0)),
        rad_tox=float(getattr(solver, "rad_tox", 0.0)),
    )
    obs_tensor = _make_obs_from_u(u_gpu, attrs).unsqueeze(0)

    dist = policy(obs_tensor)
    if target == "action":
        action = dist.probs.argmax(dim=-1, keepdim=True)
        log_prob = dist.log_prob(action.squeeze(-1)).sum()
    else:
        log_prob = dist.probs.max(dim=-1)[0].sum() + 1e-6
    grads = torch.autograd.grad(
        log_prob, u_gpu, retain_graph=False, create_graph=False
    )[0]
    sal = grads.detach().cpu().numpy()
    sal = np.abs(sal).astype(np.float32)
    if sal.max() > 1e-12:
        sal = sal / sal.max()
    if smooth_sigma > 0 and _ndi is not None:
        sal = _ndi.gaussian_filter(sal, sigma=smooth_sigma)
        sal = sal / max(sal.max(), 1e-12)
    return sal


# --------------------------------------------------------------------------- #
# Convenience: run the environment up to a given day, then capture saliency
# --------------------------------------------------------------------------- #
def rollout_saliency(
    env: Any,
    policy: Any,
    days: int = 90,
    capture_interval: int = 5,
    target: str = "action",
) -> Dict[int, Dict[str, np.ndarray]]:
    """Roll the env forward, capturing (tumour, action, saliency) every
    ``capture_interval`` steps."""
    try:
        env.reset()
    except Exception:
        pass
    captures: Dict[int, Dict[str, np.ndarray]] = {}
    for step in range(days):
        obs = get_obs(env)
        action = policy_action(policy, obs)
        try:
            env.step(int(action))
        except Exception:
            break
        if step % max(capture_interval, 1) == 0 or int(action) > 0:
            sal = compute_saliency(env, policy, target=target)
            u = np.asarray(getattr(env.solver, "u", None), dtype=np.float32)
            captures[step] = {
                "tumour": u,
                "saliency": sal,
                "action": int(action),
                "dose_on": int(action) > 0,
            }
    return captures


def get_obs(env: Any) -> np.ndarray:
    obs = env.solver.get_observation()
    return np.asarray(obs, dtype=np.float32)


def policy_action(policy: Any, obs: np.ndarray) -> int:
    if _HAVE_TORCH:
        with torch.no_grad():
            t = torch.as_tensor(obs, dtype=torch.float32).unsqueeze(0)
            dist = policy(t)
            return int(dist.probs.argmax(dim=-1).item())
    return 0


# --------------------------------------------------------------------------- #
# Visualization: 2D slice + 3D isosurface (Plotly HTML)
# --------------------------------------------------------------------------- #
def render_saliency_2d(
    tumour: np.ndarray,
    saliency: np.ndarray,
    action: int,
    day: int,
):
    """Return a matplotlib Figure overlaying saliency on the mid-z tumour slice."""
    if not _HAVE_MPL:
        raise ImportError("matplotlib required for 2D saliency rendering")
    if tumour.ndim == 3:
        sl = tumour[tumour.shape[0] // 2]
        sal_sl = saliency[saliency.shape[0] // 2]
    else:
        sl, sal_sl = tumour, saliency
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].imshow(sl, cmap="hot", origin="lower")
    axes[0].set_title(f"Day {day}: Tumour density")
    axes[1].imshow(sl, cmap="gray", origin="lower", alpha=0.4)
    axes[1].imshow(sal_sl, cmap="jet", origin="lower", alpha=0.7)
    axes[1].set_title(f"Day {day}: Saliency (action={action})")
    for a in axes:
        a.axis("off")
    plt.tight_layout()
    return fig


def render_saliency_3d_html(
    tumour: np.ndarray,
    saliency: np.ndarray,
    action: int,
    day: int,
    output_path: Path,
) -> Path:
    """Write a Plotly HTML viewer with two isosurfaces: tumour + saliency."""
    try:
        import plotly.graph_objects as go
    except ImportError:
        if not _HAVE_MPL:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(
                f"<html><body>Day {day} saliency: neither plotly nor "
                f"matplotlib available; max={float(saliency.max()):.3f}</body></html>"
            )
            return output_path
        fig = render_saliency_2d(tumour, saliency, action, day)
        png = str(output_path).replace(".html", ".png")
        fig.savefig(png, dpi=200)
        plt.close(fig)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            f"<html><body><h3>Day {day} Saliency (action={action})</h3>"
            f"<img src='{output_path.stem}.png' /></body></html>"
        )
        return output_path
    Z, Y, X = np.mgrid[0:tumour.shape[0], 0:tumour.shape[1], 0:tumour.shape[2]]
    sal_thr = max(0.3, saliency.max() * 0.6)
    sal_mask = saliency >= sal_thr
    tmask = tumour >= max(0.05, tumour.max() * 0.4)
    if not tmask.any():
        tmask[...] = True
    if not sal_mask.any():
        sal_mask[...] = True
    fig = go.Figure(data=[
        go.Isosurface(
            x=X[tmask].ravel(), y=Y[tmask].ravel(), z=Z[tmask].ravel(),
            value=tumour[tmask].ravel(),
            isomin=float(tumour[tmask].min()),
            isomax=float(tumour[tmask].max()),
            surface_count=1, opacity=0.3, colorscale="Greys",
            name="Tumour",
        ),
        go.Isosurface(
            x=X[sal_mask].ravel(), y=Y[sal_mask].ravel(), z=Z[sal_mask].ravel(),
            value=saliency[sal_mask].ravel(),
            isomin=float(sal_thr),
            isomax=float(saliency[sal_mask].max() + 1e-9),
            surface_count=1, opacity=0.85, colorscale="Jet",
            name="Saliency",
        ),
    ])
    fig.update_layout(
        title=f"Day {day}: Saliency map (action={action})",
        scene=dict(xaxis_title="x", yaxis_title="y", zaxis_title="z"),
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_path))
    return output_path


# --------------------------------------------------------------------------- #
# End-to-end saliency report (steps 4 & 5 combined)
# --------------------------------------------------------------------------- #
def saliency_to_html(
    env: Any,
    policy: Any,
    day: int,
    patient_id: str,
    output: Path = Path("output/saliency_day_00.html"),
) -> Path:
    """Drive the env to ``day``, compute saliency, write a self-contained HTML."""
    obs = get_obs(env)
    for _ in range(day):
        try:
            obs, _, done, _, _ = env.step(policy_action(policy, obs))
        except Exception:
            break
        if done:
            break
    sal = compute_saliency(env, policy, target="action")
    u = np.asarray(env.solver.u, dtype=np.float32)
    action = int(policy_action(policy, obs))
    return render_saliency_3d_html(u, sal, action, day, Path(output))


def generate_xai_report(
    env: Any,
    policy: Any,
    n_days: int = 90,
    capture_interval: int = 10,
    patient_id: str = "PAT_01",
    output_dir: Path = Path("output/xai_reports"),
) -> Dict:
    """Generates the per-day dose-on saliency HTML pages plus a compact
    summary report consumed by the HIL UI (Phase 14 iframe).

    Returns a manifest of generated HTML pages.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    captures = rollout_saliency(env, policy, days=n_days,
                                capture_interval=capture_interval,
                                target="action")
    pages = []
    html_index = ["<html><head><title>XAI Saliency Report</title></head><body>",
                  f"<h2>XAI Saliency Report - {patient_id}</h2>",
                  "<ul>"]
    for day, cap in sorted(captures.items()):
        if not cap["dose_on"]:
            continue
        day_html = output_dir / f"saliency_day_{day:02d}.html"
        render_saliency_3d_html(cap["tumour"], cap["saliency"],
                                cap["action"], day, day_html)
        pages.append(day_html)
        html_index.append(
            f'<li><a href="{day_html.name}">Day {day} '
            f'(action={cap["action"]}, dose-on)</a></li>'
        )
    html_index.append("</ul></body></html>")
    index = output_dir / f"xai_report_{patient_id}.html"
    index.write_text("\n".join(html_index))
    print(f"[xai] generated {len(pages)} saliency pages + index -> {index}")
    return {
        "patient_id": patient_id,
        "n_days": n_days,
        "n_dose_on": len(pages),
        "index_html": str(index),
        "pages": [str(p) for p in pages],
    }


# --------------------------------------------------------------------------- #
# Minimal self-test using the existing PolicyNetwork architecture
# --------------------------------------------------------------------------- #
def _self_test() -> None:
    if not _HAVE_TORCH:
        print("[xai] torch not available; skipping self-test")
        return
    rl = _import_rl_module()
    if rl is None:
        print("[xai] rl module not importable; skipping self-test")
        return
    try:
        env = rl.GbmTherapyEnv(training=True)
        env.reset()
        policy = rl.PolicyNetwork()
        sal = compute_saliency(env, policy, target="action")
        u = np.asarray(env.solver.u, dtype=np.float32)
        assert sal.shape == u.shape, (sal.shape, u.shape)
        print(f"[xai] self-test OK: saliency {sal.shape} "
              f"max={sal.max():.3f} sum={sal.sum():.3f}")
    except Exception as e:
        print(f"[xai] self-test ERROR: {type(e).__name__}: {e}")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="XAI Saliency Maps for RL (Proposal 3)")
    p.add_argument("--day", type=int, default=10, help="Decision day to inspect")
    p.add_argument("--patient-id", type=str, default="PAT_01")
    p.add_argument("--output", type=str,
                   default="output/saliency_day_10.html")
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        _self_test()
        sys.exit(0)
    rl = _import_rl_module()
    if rl is None:
        print("[xai] rl module not available; run with --self-test"); sys.exit(2)
    env = rl.GbmTherapyEnv(training=False)
    env.reset()
    policy = rl.PolicyNetwork()
    saliency_to_html(env, policy, args.day, args.patient_id,
                     Path(args.output))
    print(f"[xai] wrote HTML -> {args.output}")
