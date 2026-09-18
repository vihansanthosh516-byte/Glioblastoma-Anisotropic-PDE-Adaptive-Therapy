#!/usr/bin/env bash
# =============================================================================
# docker-entrypoint.sh - GBM Digital Twin benchmark orchestrator (Proposal 5)
# =============================================================================
# Entry modes (set CMD or first arg):
#   --benchmark      Run the full Track B + Track C pipeline (default)
#   --track-b        Run only Month 7 -> 10 cohort (Track B)
#   --track-c        Run only the Digital Twin Reactor phases (Track C)
#   --train-multiomic  Train the ElasticNet multi-omic fusion models
#   --uq-ensemble    Train/render the FNO ensemble UQ trajectory for /data
#   --tests          Run the pytest suite (tests/test_*.py)
#   --serve          Start the interactive 3D dashboard on :7860
#   <cmd>            Execute arbitrary command inside /app
#
# Data/output conventions:
#   /data    -> raw DICOM / NIfTI DTI / BraTS patient dirs (read-only mount)
#   /output  -> generated artifacts (JSON/PNG/HTML/npz mount)
# =============================================================================
set -euo pipefail

MODE="${1:-}"
[ -n "${MODE}" ] && shift || true

APP_DIR="/app"
OUTPUT_DIR="${OUTPUT_DIR:-/output}"
DATA_DIR="${DATA_DIR:-/data}"
PYTHON="${PYTHON:-python}"

mkdir -p "$OUTPUT_DIR"

run_track_b() {
    echo "############### TRACK B (Month 7 -> 10 PDE cohort) ###############"
    bash "$APP_DIR/run_all.sh"
}

run_track_c() {
    echo "############### TRACK C (Digital Twin Reactor) ###############"
    local scripts=(
        50_spatial_genomics_deconv.py
        51_inverse_parameter_estimation.py
        52_robust_mpc_controller.py
        53_spatial_metrics.py
        58_rl_adaptive_steering.py
        59_sensitivity_analysis.py
        60_baselines_and_ablation.py
        61_rl_convergence_diagnostics.py
        62_biomarker_bootstrap_stability.py
        63_reward_sensitivity.py
        64_virtual_cohort_simulation.py
        65_generate_final_report.py
    )
    for s in "${scripts[@]}"; do
        echo "########## Track C - $s ##########"
        "$PYTHON" "$APP_DIR/src/$s" --test 2>/dev/null || \
        "$PYTHON" "$APP_DIR/src/$s" || true
    done
}

run_train_multiomic() {
    echo "############### Multi-omic elastic-net training (Proposal 2) ###############"
    "$PYTHON" "$APP_DIR/src/multiomic_fusion.py" --n-patients 8 \
        --features-tsv "$OUTPUT_DIR/multiomic_features.tsv" \
        --model-path "$OUTPUT_DIR/multiomic_elasticnet.pkl"
}

run_uq_ensemble() {
    echo "############### FNO ensemble UQ (Proposal 4) ###############"
    "$PYTHON" "$APP_DIR/src/uq_fno_ensemble.py" \
        --model-dir "$OUTPUT_DIR/fno_ensemble" \
        --horizon-days 30 --M 200 --patient-id UQ_DEMO \
        --output-dir "$OUTPUT_DIR/uq"
}

run_tests() {
    echo "############### Pytest suite ###############"
    for t in inverse_estimation robust_mpc spatial_metrics multiomic_fusion; do
        echo "---- $t ----"
        "$PYTHON" "$APP_DIR/tests/test_${t}.py" || true
    done
}

run_serve() {
    echo "############### Interactive dashboard on :7860 ###############"
    exec "$PYTHON" "$APP_DIR/src/49_interactive_3d_dashboard.py" --port 7860 \
        --output-dir "$OUTPUT_DIR" || \
    exec "$PYTHON" "$APP_DIR/visualization/view_3d_time_slider.py" \
        --input-dir "$OUTPUT_DIR/time_series" --output "$OUTPUT_DIR/dashboard.html"
}

case "$MODE" in
    --benchmark)
        run_track_b
        run_train_multiomic
        run_track_c
        run_uq_ensemble
        run_tests
        ;&  # fall-through to generate the final visualization
    --track-b)
        run_track_b
        ;;
    --track-c)
        run_track_c
        ;;
    --train-multiomic)
        run_train_multiomic
        ;;
    --uq-ensemble)
        run_uq_ensemble
        ;;
    --tests)
        run_tests
        ;;
    --serve)
        run_serve
        ;;
    --help|help|-h)
        sed -n '1,30p' "$0"
        ;;
    "")
        run_track_b
        run_train_multiomic
        run_track_c
        run_uq_ensemble
        run_tests
        ;;
    *)
        exec "$@"
        ;;
esac

echo "=============================================================="
echo " Benchmark complete.  Artifacts in $OUTPUT_DIR"
echo "=============================================================="
ls -la "$OUTPUT_DIR" || true
