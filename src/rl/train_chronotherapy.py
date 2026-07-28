#!/usr/bin/env python3
"""
Phase 13: Circadian-Aware RL Training (PPO/SAC)
================================================
Train PPO agent on ChronotherapyEnv for circadian-aware adaptive chronotherapy.

Key Features:
- Continuous action space: [TMZ_rate, RT_rate, Adjuvant_rate, Phase_target]
- Circadian-aware reward shaping: timing bonus for circadian alignment
- Subpopulation resistance penalty: spatial entropy minimization
- Curriculum learning: progressive episode length
"""
import os
import sys
import numpy as np
import torch
import gymnasium as gym
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from stable_baselines3 import PPO
from stable_baselines3.common.vec_env import DummyVecEnv, VecNormalize
from stable_baselines3.common.callbacks import EvalCallback, CheckpointCallback
from stable_baselines3.common.monitor import Monitor
from stable_baselines3.common.logger import configure

from rl.chronotherapy_env import ChronotherapyEnv


def make_env(grid_size=32, dt_hours=2.0, max_episode_hours=48, circadian=True, seed=42, rank=0):
    """Factory for creating environments."""
    def _init():
        env = ChronotherapyEnv(
            grid_size=grid_size,
            dt_hours=dt_hours,
            max_episode_hours=max_episode_hours,
            circadian=circadian,
            seed=seed + rank
        )
        env = Monitor(env)  # Log episode stats
        return env
    return _init


def train_ppo(
    total_timesteps: int = 200_000,
    n_envs: int = 4,
    grid_size: int = 32,
    dt_hours: float = 2.0,
    max_episode_hours: int = 48,
    circadian: bool = True,
    learning_rate: float = 3e-4,
    n_steps: int = 2048,
    batch_size: int = 64,
    n_epochs: int = 10,
    gamma: float = 0.99,
    gae_lambda: float = 0.95,
    clip_range: float = 0.2,
    ent_coef: float = 0.01,
    vf_coef: float = 0.5,
    max_grad_norm: float = 0.5,
    seed: int = 42,
    save_path: str = "output/phase13_ppo_chronotherapy",
    log_dir: str = "output/logs/phase13_ppo",
):
    """
    Train PPO agent on circadian-aware chronotherapy environment.
    
    Args:
        total_timesteps: Total training timesteps
        n_envs: Number of parallel environments
        grid_size: Tumor grid resolution
        dt_hours: Hours per RL step
        max_episode_hours: Max hours per episode
        circadian: Enable circadian rhythms
        learning_rate: PPO learning rate
        n_steps: Steps per update per env
        batch_size: Minibatch size
        n_epochs: PPO epochs per update
        gamma: Discount factor
        gae_lambda: GAE lambda
        clip_range: PPO clip range
        ent_coef: Entropy coefficient
        vf_coef: Value function coefficient
        max_grad_norm: Gradient clipping
        seed: Random seed
        save_path: Model save path
        log_dir: TensorBoard log directory
    """
    print("=" * 70)
    print("Phase 13: Circadian-Aware PPO Training for Chronotherapy")
    print("=" * 70)
    print(f"Total timesteps: {total_timesteps:,}")
    print(f"Parallel envs: {n_envs}")
    print(f"Grid size: {grid_size}^3, dt={dt_hours}h, max_ep={max_episode_hours}h")
    print(f"Circadian: {circadian}")
    print(f"Learning rate: {learning_rate}, n_steps: {n_steps}")
    print(f"Batch: {batch_size}, epochs: {n_epochs}")
    print(f"Seed: {seed}")
    print("=" * 70)
    
    # Create directories
    os.makedirs(save_path, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    
    # Set seeds
    np.random.seed(seed)
    torch.manual_seed(seed)
    
    # Create vectorized environments
    env_fns = [make_env(grid_size, dt_hours, max_episode_hours, circadian, seed, i) 
               for i in range(n_envs)]
    vec_env = DummyVecEnv(env_fns)
    
    # Normalize observations and rewards
    vec_env = VecNormalize(
        vec_env, 
        norm_obs=True, 
        norm_reward=True, 
        clip_obs=10.0,
        gamma=gamma
    )
    
    # Create model
    model = PPO(
        "MlpPolicy",
        vec_env,
        learning_rate=learning_rate,
        n_steps=n_steps,
        batch_size=batch_size,
        n_epochs=n_epochs,
        gamma=gamma,
        gae_lambda=gae_lambda,
        clip_range=clip_range,
        ent_coef=ent_coef,
        vf_coef=vf_coef,
        max_grad_norm=max_grad_norm,
        verbose=1,
        seed=seed,
        tensorboard_log=log_dir,
        device="auto",
    )
    
    # Callbacks
    eval_env = DummyVecEnv([make_env(grid_size, dt_hours, max_episode_hours, circadian, seed + 100, 0)])
    eval_env = VecNormalize(eval_env, norm_obs=True, norm_reward=False, training=False)
    # Sync normalization stats
    eval_env.obs_rms = vec_env.obs_rms
    eval_env.ret_rms = vec_env.ret_rms
    
    eval_callback = EvalCallback(
        eval_env,
        best_model_save_path=save_path,
        log_path=log_dir,
        eval_freq=max(10000 // n_envs, 1),
        deterministic=True,
        render=False,
        n_eval_episodes=5,
    )
    
    checkpoint_callback = CheckpointCallback(
        save_freq=max(50000 // n_envs, 1),
        save_path=save_path,
        name_prefix="ppo_chrono",
        save_replay_buffer=False,
        save_vecnormalize=True,
    )
    
    print("\n[Training] Starting PPO training...")
    print(f"[Training] Total timesteps: {total_timesteps:,}")
    print(f"[Training] Log dir: {log_dir}")
    print(f"[Training] Save path: {save_path}")
    
    model.learn(
        total_timesteps=total_timesteps,
        callback=[eval_callback, checkpoint_callback],
        progress_bar=True,
    )
    
    # Save final model and normalization stats
    model.save(os.path.join(save_path, "ppo_chronotherapy_final"))
    vec_env.save(os.path.join(save_path, "vecnormalize.pkl"))
    
    print(f"\n[SUCCESS] Training complete! Model saved to {save_path}")
    
    return model, vec_env


def evaluate_model(
    model_path: str,
    vecnorm_path: str,
    n_episodes: int = 10,
    grid_size: int = 32,
    dt_hours: float = 2.0,
    max_episode_hours: int = 48,
    circadian: bool = True,
    seed: int = 42,
):
    """Evaluate trained model."""
    print(f"\n[Evaluation] Loading model from {model_path}")
    
    # Create eval env
    eval_env = DummyVecEnv([make_env(grid_size, 2.0, 48, circadian, seed + 200, 0)])
    eval_env = VecNormalize.load(vecnorm_path, eval_env)
    eval_env.training = False
    eval_env.norm_reward = False
    
    model = PPO.load(model_path, env=eval_env)
    
    episode_rewards = []
    episode_lengths = []
    volumes = []
    
    for ep in range(n_episodes):
        obs = eval_env.reset()
        done = False
        ep_reward = 0
        ep_len = 0
        
        while not done:
            action, _ = model.predict(obs, deterministic=True)
            obs, reward, done, info = eval_env.step(action)
            ep_reward += reward[0] if isinstance(reward, np.ndarray) else reward
            ep_len += 1
        
        episode_rewards.append(ep_reward)
        episode_lengths.append(ep_len)
        # Get volume from inner env
        volumes.append(eval_env.envs[0].unwrapped.unwrapped.u.sum() if hasattr(eval_env.envs[0], 'unwrapped') else 0)
        
        print(f"  Episode {ep+1}: Reward={ep_reward:.2f}, Length={ep_len}")
    
    print(f"\n[Evaluation] Mean Reward: {np.mean(episode_rewards):.2f} ± {np.std(episode_rewards):.2f}")
    print(f"[Evaluation] Mean Length: {np.mean(episode_lengths):.1f} ± {np.std(episode_lengths):.1f}")
    
    return episode_rewards, episode_lengths


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Phase 13: Circadian-Aware PPO Training")
    parser.add_argument("--timesteps", type=int, default=200_000)
    parser.add_argument("--n-envs", type=int, default=4)
    parser.add_argument("--grid-size", type=int, default=32)
    parser.add_argument("--dt-hours", type=float, default=2.0)
    parser.add_argument("--max-ep-hours", type=int, default=48)
    parser.add_argument("--no-circadian", action="store_true")
    parser.add_argument("--lr", type=float, default=3e-4)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--eval-only", action="store_true")
    parser.add_argument("--model-path", type=str, default="")
    
    args = parser.parse_args()
    
    if args.eval_only:
        if not args.model_path:
            print("Error: --model-path required for eval-only mode")
            exit(1)
        evaluate_model(
            model_path=args.model_path,
            vecnorm_path=os.path.join(os.path.dirname(args.model_path), "vecnormalize.pkl"),
            n_episodes=10,
            grid_size=args.grid_size,
            dt_hours=args.dt_hours,
            max_episode_hours=args.max_ep_hours,
            circadian=not args.no_circadian,
            seed=args.seed,
        )
    else:
        train_ppo(
            total_timesteps=args.timesteps,
            n_envs=args.n_envs,
            grid_size=args.grid_size,
            dt_hours=args.dt_hours,
            max_episode_hours=args.max_ep_hours,
            circadian=not args.no_circadian,
            learning_rate=args.lr,
            seed=args.seed,
        )