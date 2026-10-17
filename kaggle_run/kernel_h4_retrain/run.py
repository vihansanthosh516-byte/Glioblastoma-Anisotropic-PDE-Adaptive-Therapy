# Script 136 on Kaggle: H-4 step 2 retrain (Amendment 13 A13.2), arm hybrid_mri_rank_mu on the pooled set.
# Inputs: predictgbm-inputs-private, gbm-code-private, predictgbm-mri2mm-private, mu-train-2mm-private,
#         h4-retrain-pack-private (kaggle datasets create -p kaggle_run/pack_h4 --private)
import glob, os, shutil, subprocess, sys, time
T0 = time.time()
f = lambda n: glob.glob(f"/kaggle/input/**/{n}", recursive=True)[0]
PACK = f("pg_inputs.npz")
code = f("src/136_h4_retrain_train.py")
REPO = "/kaggle/working/repo"
shutil.copytree(os.path.dirname(os.path.dirname(code)), REPO, dirs_exist_ok=True)
env = dict(os.environ, GBM_PG_PACK=PACK, GBM_PG_TEST_IDS=os.path.join(os.path.dirname(PACK), "test_ids.json"),
           GBM_MU_PACK=f("mu_train_2mm.npz"), GBM_PG_MRI=f("pg_mri2mm.npz"),
           GBM_H4_RETRAIN_PACK=f("h4_retrain_pack.npz"), GBM_DISK_CACHE="/kaggle/working/cache",
           PYTHONUNBUFFERED="1")
subprocess.run(["nvidia-smi"])
r = subprocess.run([sys.executable, "-u", "src/136_h4_retrain_train.py"], cwd=REPO, env=env)
for p in glob.glob(os.path.join(REPO, "output", "predictgbm", "h4_retrain", "*")):
    shutil.copy2(p, "/kaggle/working/")
print(f"exit {r.returncode}; done in {(time.time() - T0) / 3600:.2f} h", flush=True)
