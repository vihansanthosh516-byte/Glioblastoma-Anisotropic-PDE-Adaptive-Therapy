#!/usr/bin/env python3
"""
Phase 15: Virtual Clinical Trial - 1000-Patient In-Silico Evaluation
=====================================================================
Large-scale virtual cohort evaluation of the trained circadian-aware PPO policy
against standard-of-care (Stupp protocol) and adaptive therapy baselines.

Requires GPU for fast FNO inference and parallel patient simulation.
"""
import os
import sys
import json
import numpy as np
import torch
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

from rl.chronotherapy_env import ChronotherapyEnv
from stable_baselines3 import PPO
from stable_baselines3.common.vec_env import DummyVecEnv, VecNormalize

OUTPUT_DIR = Path("/kaggle/working/output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_trained_policy(model_path: str, vecnorm_path: str):
    """Load trained PPO policy with VecNormalize."""
    # Create eval env
    def make_env():
        env = ChronotherapyEnv(
            grid_size=32,
            dt_hours=2.0,
            max_episode_hours=48,
            circadian=True,
            seed=42
        )
        return env
    
    eval_env = DummyVecEnv([make_env])
    eval_env = VecNormalize.load(vecnorm_path, eval_env)
    eval_env.training = False
    eval_env.norm_reward = False
    
    model = PPO.load(model_path, env=eval_env)
    return model, eval_env


def run_virtual_trial(
    model,
    eval_env,
    n_patients: int = 1000,
    max_episode_hours: int = 168,  # 7 days
    dt_hours: float = 2.0,
) -> Dict:
    """
    Run virtual clinical trial with 1000 patients.
    
    Compares:
    - PPO policy (trained)
    - Standard Stupp protocol (control)
    - Adaptive threshold therapy (baseline)
    """
    print(f"[Phase 15] Starting virtual trial with {n_patients} patients...")
    print(f"[Phase 15] Episode length: {max_episode_hours}h, dt={2.0}h")
    
    results = {
        "ppo": {"volumes": [], "rewards": [], "toxicity": [], "clearance": 0, "survival": []},
        "stupp": {"volumes": [], "rewards": [], "toxicity": [], "clearance": 0, "survival": []},
        "adaptive": {"volumes": [], "rewards": [], "toxicity": [], "clearance": 0, "survival": []},
    }
    
    # Get base env for manual control
    base_env = eval_env.envs[0].unwrapped
    
    for patient_idx in range(n_patients):
        if patient_idx % 100 == 0:
            print(f"[Phase 15] Progress: {patient_idx}/{n_patients} patients")
        
        # --- PPO Policy ---
        obs = eval_env.reset()
        ep_reward = 0
        terminated = False
        
        while not terminated:
            action, _ = model.predict(obs, deterministic=True)
            obs, reward, done, info = eval_env.step(action)
            terminated = done[0] if isinstance(done, np.ndarray) else done
        
        # Get final volume from base env
        ppo_volume = eval_env.envs[0].unwrapped.u.sum()
        results["ppo"]["volumes"].append(ppo_volume)
        results["ppo"]["clearance"] += 1 if ppo_volume < 10 else 0
        
        # --- Stupp Protocol (Standard of Care) ---
        # Simplified: fixed schedule TMZ + RT
        base_env.reset()
        base_env.drug_concs = {"TMZ": 0.0, "Inhibitor": 0.0, "Radiation": 0.0}
        base_env.last_dose_time = {"TMZ": -100, "Inhibitor": -100, "Radiation": -100}
        base_env.cum_toxicity = {"TMZ": 0.0, "Inhibitor": 0.0, "Radiation": 0.0}
        
        stupp_volume = 0
        for step in range(int(max_episode_hours / 2.0)):
            # Stupp: TMZ days 1-5, 15-19, 29-33... RT days 1-30
            day = step * 2
            in_tmz = (day % 28) < 5
            in_rt = day < 30
            
            dose_tmz = 1.0 if in_tmz else 0.0
            dose_rad = 1.0 if in_rt else 0.0
            
            action = np.array([float(dose_tmz), 0.0, float(dose_rad), 0.0], dtype=np.float32)
            obs, reward, terminated, truncated, _ = base_env.envs[0].step(action)
            if terminated or truncated:
                break
        
        stupp_volume = base_env.unwrapped.u.sum()
        results["stupp"]["volumes"].append(stupp_volume)
        results["stupp"]["clearance"] += 1 if stupp_volume < 10 else 0
        
        # --- Adaptive Threshold Therapy ---
        base_env.reset()
        base_env.drug_concs = {"TMZ": 0.0, "Inhibitor": 0.0, "Radiation": 0.0}
        base_env.last_dose_time = {"TMZ": -100, "Inhibitor": -100, "Radiation": -100}
        base_env.cum_toxicity = {"TMZ": 0.0, "Inhibitor": 0.0, "Radiation": 0.0}
        
        vol_prev = base_env.unwrapped.u.sum()
        for step in range(int(max_episode_hours / 2.0)):
            vol_curr = base_env.unwrapped.u.sum()
            
            # Threshold: dose if volume > 120% of baseline
            if vol_curr > vol_prev * 1.2:
                action = np.array([1.0, 0.0, 1.0, 0.0], dtype=np.float32)
            else:
                action = np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float32)
            
            obs, reward, terminated, truncated, _ = base_env.envs[0].step(action)
            vol_prev = vol_curr
            if terminated or truncated:
                break
        
        adaptive_volume = base_env.unwrapped.u.sum()
        results["adaptive"]["volumes"].append(adaptive_volume)
        results["adaptive"]["clearance"] += 1 if adaptive_volume < 10 else 0
        
        # Progress
        if (patient_idx + 1) % 100 == 0:
            print(f"  Patient {patient_idx+1}: PPO={ppo_volume:.0f}, Stupp={stupp_volume:.0f}, Adaptive={adaptive_volume:.0f}")
    
    return results


def compute_statistics(results: Dict) -> Dict:
    """Compute statistical summary of trial results."""
    stats = {}
    for arm in ["ppo", "stupp", "adaptive"]:
        vols = np.array(results[arm]["volumes"])
        stats[arm] = {
            "n": len(vols),
            "mean_volume": float(np.mean(vols)),
            "std_volume": float(np.std(vols)),
            "median_volume": float(np.median(vols)),
            "clearance_rate": results[arm]["clearance"] / len(vols) * 100,
            "volume_p25": float(np.percentile(vols, 25)),
            "volume_p75": float(np.percentile(vols, 75)),
        }
    
    # Pairwise comparisons
    from scipy import stats as scipy_stats
    stats["comparisons"] = {
        "ppo_vs_stupp": {
            "ttest_p": float(scipy_stats.ttest_rel(results["ppo"]["volumes"], results["stupp"]["volumes"]).pvalue),
            "wilcoxon_p": float(scipy_stats.wilcoxon(results["ppo"]["volumes"], results["stupp"]["volumes"]).pvalue),
            "cohens_d": float((np.mean(results["stupp"]["volumes"]) - np.mean(results["ppo"]["volumes"])) / 
                            np.std(results["stupp"]["volumes"] + results["ppo"]["volumes"])),
        },
        "ppo_vs_adaptive": {
            "ttest_p": float(scipy_stats.ttest_rel(results["ppo"]["volumes"], results["adaptive"]["volumes"]).pvalue),
            "wilcoxon_p": float(scipy_stats.wilcoxon(results["ppo"]["volumes"], results["adaptive"]["volumes"]).pvalue),
            "cohens_d": float((np.mean(results["adaptive"]["volumes"]) - np.mean(results["ppo"]["volumes"])) / 
                            np.std(results["adaptive"]["volumes"] + results["ppo"]["volumes"])),
        },
    }
    
    return stats


def main():
    print("=" * 70)
    print("Phase 15: Virtual Clinical Trial - 1000-Patient Cohort")
    print("=" * 70)
    
    # Paths
    model_path = "/kaggle/working/gbm-repo/output/phase13_ppo_chronotherapy/ppo_chronotherapy_final.zip"
    vecnorm_path = "/kaggle/working/gbm-repo/output/phase13_ppo_chronotherapy/vecnormalize.pkl"
    
    if not os.path.exists(model_path):
        print(f"[ERROR] Model not found at {model_path}")
        print("[INFO] Run Phase 13 training first to generate model")
        return
    
    # Load model
    print(f"[Phase 15] Loading model from {model_path}")
    model, eval_env = load_trained_policy(model_path, "/kaggle/working/gbm-repo/" + vecnorm_path)
    
    # Run trial
    results = run_virtual_trial(model, eval_env, n_patients=1000, max_episode_hours=168)
    
    # Statistics
    stats = compute_statistics(results)
    
    # Save results
    output = {
        "timestamp": datetime.now().isoformat(),
        "n_patients": 1000,
        "max_episode_hours": 168,
        "results": results,
        "statistics": stats,
    }
    
    output_path = OUTPUT_DIR / "phase15_virtual_trial_results.json"
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    
    # Print summary
    print("\n" + "=" * 70)
    print("PHASE 15 VIRTUAL TRIAL RESULTS")
    print("=" * 70)
    
    for arm in ["ppo", "stupp", "adaptive"]:
        s = stats[arm]
        print(f"\n{arm.upper()}:")
        print(f"  Mean Volume: {s['mean_volume']:.1f} ± {s['std_volume']:.1f} mm³")
        print(f"  Median Volume: {s['median_volume']:.1f} mm³")
        print(f"  Clearance Rate: {s['clearance_rate']:.1f}%")
        print(f"  IQR: [{s['volume_p25']:.1f}, {s['volume_p75']:.1f}] mm³")
    
    print(f"\nPPO vs Stupp: p={stats['comparisons']['ppo_vs_stupp']['ttest_p']:.4f}, "
          f"d={stats['comparisons']['ppo_vs_stupp']['cohens_d']:.2f}")
    print(f"PPO vs Adaptive: p={stats['comparisons']['ppo_vs_adaptive']['ttest_p']:.4f}, "
          f"d={stats['comparisons']['ppo_vs_adaptive']['cohens_d']:.2f}")
    
    print(f"\n[SUCCESS] Results saved to {output_path}")


if __name__ == "__main__":
    main()