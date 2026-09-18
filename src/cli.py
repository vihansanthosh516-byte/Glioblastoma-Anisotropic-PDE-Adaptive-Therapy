"""CLI entry point for the Glioblastoma Anisotropic PDE & Adaptive Therapy Pipeline."""

import os
import sys
import subprocess

# Project root from env, falling back to the repository root (three levels up
# from this file: src/cli.py -> src/ -> <project root>). This works on any
# platform without hardcoded user paths.
PROJECT_ROOT = os.getenv("GBM_PROJECT_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..")
)

SCRIPTS = {
    "1": "src/01_load_and_filter.py",
    "2": "src/02_preprocess_umap_de.py",
    "3": "src/03_finalize_de_and_export_for_nn.py",
    "4": "src/04_export_for_attention_model.py",
    "5": "src/05_attention_gated_network.py",
    "6": "src/06_method1_classical_baseline.py",
    "7": "src/07_method2_transformer.py",
    "8": "src/08_method3_hybrid.py",
    "9": "src/09_benchmark_comparison.py",
    "10": "src/10_cvae_pretrain.py",
    "11": "src/11_cvae_extract_latent.py",
    "12": "src/12_gat_build_graph.py",
    "13": "src/13_gat_train.py",
    "14": "src/14_cgat_evaluate.py",
    "15": "src/15_baseline_scvi.py",
    "16": "src/16_baseline_nmf.py",
    "16c": "src/55_cgat_train.py",
    "17": "src/17_gradient_diagnostic.py",
    "18": "src/18_csgt_framework.py",
    "19": "src/19_phenotypic_velocity.py",
    "20": "src/20_fokker_planck_solver.py",
    "21": "src/21_drift_diffusion_analysis.py",
    "22": "src/22_saddle_point_proof.py",
    "23": "src/23_transfer_entropy_engine.py",
    "24": "src/24_causal_grn_builder.py",
    "25": "src/25_pid_analysis.py",
    "26": "src/26_grn_validation.py",
    "27": "src/27_aba_lattice.py",
    "28": "src/28_fisher_kolmogorov_pde.py",
    "29": "src/29_invasion_simulator.py",
    "30": "src/30_aba_analysis.py",
    "31": "src/31_virtual_knockout_engine.py",
    "32": "src/32_combinatorial_screen.py",
    "33": "src/33_therapeutic_index.py",
    "34": "src/34_drug_gating_report.py",
    "35": "src/35_ivygap_clinical_ingest.py",
    "36": "src/36_survival_analysis.py",
    "37": "src/37_clinical_validation_report.py",
    "38": "src/38_real_cohort_ingest.py",
    "39": "src/39_penalized_survival.py",
    "40": "src/40_spatial_recurrence_mapper.py",
    "41": "src/41_dose_response_model.py",
    "42": "src/42_anisotropic_pde.py",
    "43": "src/43_stromal_feedback.py",
    "44": "src/44_adaptive_therapy.py",
    "45": "src/45_validation_synthesis.py",
    "46": "src/46_sensitivity_analysis.py",
    "47": "src/47_optimal_control.py",
    "48": "src/48_3d_extension.py",
    "49": "src/49_interactive_3d_dashboard.py",
    "50": "src/50_clinical_cdss_app.py",
    "51": "src/51_inverse_parameter_estimation.py",
    "52": "src/52_robust_mpc_controller.py",
    "53": "src/53_spatial_metrics.py",    "54": "src/54_cgat_build_graph.py",    "55": "src/55_cgat_train.py",    "56": "src/56_multi_scale_framework.py",
}


def run_script(script_name: str) -> int:
    """Run a single pipeline script."""
    script_path = os.path.join(PROJECT_ROOT, SCRIPTS[script_name])
    if not os.path.exists(script_path):
        print(f"[CLI] ERROR: script '{script_name}' -> {script_path} not found", file=sys.stderr)
        return 1
    env = os.environ.copy()
    env["GBM_PROJECT_ROOT"] = PROJECT_ROOT
    result = subprocess.run(
        [sys.executable, script_path],
        env=env,
        cwd=PROJECT_ROOT,
    )
    return result.returncode


def list_scripts():
    """List available scripts."""
    print("Available pipeline scripts:")
    print("=" * 60)
    for key, path in sorted(SCRIPTS.items(), key=lambda x: int(x[0]) if x[0].isdigit() else 999):
        label = f"{key}: {os.path.basename(path)}"
        print(f"  {label}")
    print("=" * 60)


def main():
    """Main CLI entry point."""
    if len(sys.argv) < 2:
        list_scripts()
        print("\nUsage: python -m src.cli <script_number> | --list")
        print("       python -m src.cli --run-all       (run all scripts sequentially)")
        print("       GBM_PROJECT_ROOT=/path python -m src.cli <script_number>")
        sys.exit(1)

    arg = sys.argv[1]

    if arg == "--list" or arg == "-l":
        list_scripts()
        return

    if arg == "--run-all":
        print("Running all scripts sequentially...")
        print("=" * 60)
        failed = 0
        total = len(SCRIPTS)
        for key in sorted(SCRIPTS.keys(), key=lambda x: int(x) if x.isdigit() else 999):
            script_name = SCRIPTS[key]
            rc = run_script(key)
            if rc != 0:
                print(f"[CLI] FAILED at script {key} ({os.path.basename(script_path) if False else script_name})")
                failed += 1
            else:
                print(f"[CLI] OK   script {key}")
        print("=" * 60)
        print(f"Completed: {total - failed}/{total} succeeded; {failed} failed")
        if failed > 0:
            sys.exit(1)
        return

    # Run specific script
    script_key = arg
    if script_key not in SCRIPTS:
        print(f"[CLI] ERROR: unknown script key '{script_key}'. Use --list to see options.", file=sys.stderr)
        sys.exit(1)

    rc = run_script(script_key)
    sys.exit(rc)


if __name__ == "__main__":
    main()