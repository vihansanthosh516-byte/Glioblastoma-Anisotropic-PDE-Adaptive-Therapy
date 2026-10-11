# Kaggle GPU runner: script 115 for H-3 candidates (Amendment 9), development TUM 5-fold out-of-fold.
import glob, os, shutil, subprocess, sys, time
T0 = time.time()
ARMS = "hybrid_mri_rank_mu_s0,hybrid_mri_rank_enhw_mu_s0"
PACK = glob.glob("/kaggle/input/**/pg_inputs.npz", recursive=True)[0]
DS = os.path.dirname(PACK)
REPO = "/kaggle/working/repo"
code = glob.glob("/kaggle/input/**/src/115_predictgbm_unet_dev.py", recursive=True)[0]
shutil.copytree(os.path.dirname(os.path.dirname(code)), REPO, dirs_exist_ok=True)
for f in glob.glob(os.path.join(DS, "*")):
    if os.path.isdir(os.path.join(f, "src")):
        pass
env = dict(os.environ, GBM_PG_PACK=PACK, GBM_PG_TEST_IDS=os.path.join(DS, "test_ids.json"),
           GBM_MU_PACK=glob.glob("/kaggle/input/**/mu_train_2mm.npz", recursive=True)[0],
           GBM_PG_MRI=glob.glob("/kaggle/input/**/pg_mri2mm.npz", recursive=True)[0], PYTHONUNBUFFERED="1")
subprocess.run(["nvidia-smi"])
r = subprocess.run([sys.executable, "-u", "src/115_predictgbm_unet_dev.py", "--arms", ARMS, "--epochs", "60", "--tag", "_h3A"], cwd=REPO, env=env)
for p in glob.glob(os.path.join(REPO, "output", "predictgbm", "*h3*")) + glob.glob(os.path.join(REPO, "output", "predictgbm", "oof_logits_*.npz")):
    shutil.copy2(p, "/kaggle/working/")
print(f"exit {r.returncode}; done in {(time.time() - T0) / 3600:.2f} h", flush=True)
