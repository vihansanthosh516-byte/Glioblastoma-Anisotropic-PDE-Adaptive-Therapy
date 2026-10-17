# Script 132 on Kaggle: H-3 FINAL model (Amendment 11, C1 = hybrid_mri_rank_mu). Train once, predict UCSF inputs only.
# Inputs: predictgbm-inputs-private, gbm-code-private, predictgbm-mri2mm-private, mu-train-2mm-private, ucsf-h3-inputs-private
import glob, os, shutil, subprocess, sys, time
T0 = time.time()
f = lambda n: glob.glob(f"/kaggle/input/**/{n}", recursive=True)[0]
PACK = f("pg_inputs.npz")
code = f("src/132_h3_final_train_predict.py")
REPO = "/kaggle/working/repo"
shutil.copytree(os.path.dirname(os.path.dirname(code)), REPO, dirs_exist_ok=True)
env = dict(os.environ, GBM_PG_PACK=PACK, GBM_PG_TEST_IDS=os.path.join(os.path.dirname(PACK), "test_ids.json"),
           GBM_MU_PACK=f("mu_train_2mm.npz"), GBM_PG_MRI=f("pg_mri2mm.npz"), GBM_UCSF_PACK=f("ucsf_h3_inputs_2mm.npz"),
           PYTHONUNBUFFERED="1")
subprocess.run(["nvidia-smi"])
r = subprocess.run([sys.executable, "-u", "src/132_h3_final_train_predict.py"], cwd=REPO, env=env)
for p in glob.glob(os.path.join(REPO, "output", "predictgbm", "h3_final", "*")):
    shutil.copy2(p, "/kaggle/working/")
print(f"exit {r.returncode}; done in {(time.time() - T0) / 3600:.2f} h", flush=True)
