# =============================================================================
# Dockerfile - GBM Digital Twin & RL Adaptive Therapy Framework (Proposal 5)
# =============================================================================
# Reproducible, stateless container that runs the full Track B + C end-to-end
# pipeline (aniso PDE -> inverse est -> robust MPC -> 3D DTI -> RL -> SA ->
# baselines -> virtual cohort -> report) and emits JSON/PNG/HTML artifacts to
# /output.  Mount raw DICOM/DTI/BraTS data at /data; results land in /output.
#
# Build:
#   docker build -t gbm-digital-twin:$(git rev-parse --short HEAD) .
#   docker build -t gbm-digital-twin:latest -t gbm-digital-twin:$(git rev-parse --short HEAD) .
#
# Run full benchmark:
#   docker compose run --rm benchmark
#   docker run --rm -v "$(pwd)/data:/data" -v "$(pwd)/output:/output" \
#     gbm-digital-twin:latest
#
# Run an interactive shell inside the container:
#   docker run --rm -it -v "$(pwd)/data:/data" -v "$(pwd)/output:/output" \
#     --entrypoint /bin/bash gbm-digital-twin:latest
#
# Spinning up the CPU base image keeps the image portable; for GPU work
# swap FROM to pytorch/pytorch:2.2.0-cuda12.1-cudnn8-runtime and rebuild.
# =============================================================================
FROM python:3.12-slim

LABEL org.opencontainers.image.title="GBM-Digital-Twin"
LABEL org.opencontainers.image.description="Biophysical Glioblastoma Digital Twin & RL Adaptive Therapy Framework"
LABEL org.opencontainers.image.licenses="MIT"

# ---------- System libraries (ffmpeg for viz, GL for matplotlib headless) ---
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        ffmpeg \
        libglib2.0-0 \
        libgl1-mesa-glx \
        build-essential \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# ---------- Python dependencies -------------------------------------------------
# Install pinned core deps first so this layer is cacheable when only the
# repo code changes. Includes the new deps from Proposals 1-4:
#   nibabel + dipy     (real DTI ingestion)
#   scikit-learn + h5py (multi-omic elastic-net + HDF5)
#   torch (CPU)        (FNO rollout, RL, XAI saliency, UQ ensemble)
#   gymnasium + plotly  (closed-loop env + HIL/XAI dashboards)
COPY Requirements.txt /app/Requirements.txt
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r Requirements.txt && \
    pip install --no-cache-dir \
        numpy scipy matplotlib pillow pandas \
        gymnasium \
        plotly \
        pymc pytensor \
        "--extra-index-url" "https://download.pytorch.org/whl/cpu" torch

# ---------- Repo code ----------------------------------------------------------
COPY . /app

# ---------- Stateless mount points ---------------------------------------------
# /data  : raw DICOM / NIfTI DTI / BraTS patient directories (read-only)
# /output: generated JSON/PNG/HTML/npz artifacts
RUN mkdir -p /data /output

# Make the helper scripts runnable
RUN chmod +x /app/docker-entrypoint.sh /app/run_all.sh 2>/dev/null || true

ENV PYTHONUNBUFFERED=1 \
    PYTHONPATH=/app:/app/src \
    MPLBACKEND=Agg

EXPOSE 7860

# ---------- Entrypoint ---------------------------------------------------------
# Default entrypoint orchestrates the full benchmark chain.
# Override with `--entrypoint /bin/bash` for debugging, or pass additional
# CLI args that are appended to the benchmark command.
ENTRYPOINT ["/app/docker-entrypoint.sh"]
CMD ["--benchmark"]
