import os
import subprocess
import sys

repo_dir = "/kaggle/working/gbm-repo"
if not os.path.exists(repo_dir):
    print("Cloning repository...")
    subprocess.run(["git", "clone", "https://github.com/vihansanthosh/gbm-pde-phase1.git", repo_dir], check=True)

os.chdir(repo_dir)

print("============================================================")
print("RUNNING DIRECTLY: Phase 5 RL Adaptive Steering")
print("============================================================")
subprocess.run([sys.executable, "src/58_rl_adaptive_steering.py"], check=True)
