#!/usr/bin/env python3
"""
visualization/view_3d_saliency.py
=================================
Loads a saved policy checkpoint + a patient environment, runs the env to a
chosen decision day, and writes a self-contained Plotly HTML saliency page
(``saliency_day_XX.html``) that can be embedded by the HIL UI iframe.

This is the Visualization entrypoint referenced by Proposal 3, Step 4.

Usage:
    python visualization/view_3d_saliency.py \\
        --checkpoint output/phase5_rl_policy.pth \\
        --patient-id PAT_01 \\
        --day 10 \\
        --output output/saliency_day_10.html

If no checkpoint is supplied, a freshly-initialized PolicyNetwork from
src/58_rl_adaptive_steering.py is used so the script always produces an
artifact (useful for smoke-testing the XAI pipeline).
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
# Make `src` importable both as 'src.*' and via filename modules
SRC = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Visualize RL policy saliency maps (Proposal 3)"
    )
    parser.add_argument(
        "--checkpoint", type=str, default=None,
        help="Path to a torch state_dict for PolicyNetwork "
             "(from src/58_rl_adaptive_steering.py). When omitted a fresh "
             "PolicyNetwork is used."
    )
    parser.add_argument("--patient-id", type=str, default="PAT_01")
    parser.add_argument("--day", type=int, default=10)
    parser.add_argument("--output", type=str,
                        default="output/saliency_day_10.html")
    parser.add_argument("--full-report", action="store_true",
                        help="Generate the full multi-day XAI report instead "
                             "of a single-day page.")
    parser.add_argument("--n-days", type=int, default=90)
    parser.add_argument("--capture-interval", type=int, default=10)
    args = parser.parse_args()

    try:
        import torch
    except ImportError:
        print("[view_3d_saliency] torch is required to compute saliency.")
        sys.exit(2)

    sys.path.insert(0, str(SRC))
    rl = _load("rl_adaptive_steering", SRC / "58_rl_adaptive_steering.py")
    saliency_mod = _load("xai_saliency", SRC / "xai_saliency.py")

    env = rl.GbmTherapyEnv(training=False)
    env.reset()
    policy = rl.PolicyNetwork()
    if args.checkpoint and Path(args.checkpoint).exists():
        try:
            state = torch.load(args.checkpoint, map_location="cpu")
            if isinstance(state, dict) and "model_state" in state:
                state = state["model_state"]
            policy.load_state_dict(state)
            print(f"[view_3d_saliency] loaded checkpoint -> {args.checkpoint}")
        except Exception as e:
            print(f"[view_3d_saliency] WARN checkpoint load failed: {e}")

    out_path = Path(args.output)
    if args.full_report:
        report = saliency_mod.generate_xai_report(
            env, policy, n_days=args.n_days,
            capture_interval=args.capture_interval,
            patient_id=args.patient_id,
            output_dir=out_path.parent / f"xai_reports_{args.patient_id}",
        )
        print(f"[view_3d_saliency] report index -> {report['index_html']}"
              f" ({report['n_dose_on']} dose-on days captured)")
        return

    saliency_mod.saliency_to_html(env, policy, args.day, args.patient_id, out_path)
    print(f"[view_3d_saliency] wrote -> {out_path}")


if __name__ == "__main__":
    main()
