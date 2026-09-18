# =============================================================================
# Makefile - GBM Digital Twin convenience targets (Proposal 5, step 8)
# =============================================================================
# Mirrors Proposal 5 of the roadmap:
#   - build/run/test the containerized benchmark suite
#   - tag images with the git commit SHA for versioning & CI
#   - run the in-container pytest suite
#   - shell into a debugging container
#
# Requires docker (and docker compose) on PATH.  Git is used for SHA tagging.
# =============================================================================
.DEFAULT_GOAL := help

PY        ?= python
IMAGE     ?= gbm-digital-twin
COMPOSE   ?= docker compose
SHA       := $(shell git rev-parse --short HEAD 2>/dev/null || echo dev)
TAG       ?= latest

.PHONY: help build build-sha run test shell clean prune multiomic uq

help:  ## Show this help
	@echo "GBM Digital Twin Docker benchmark targets:"
	@echo "  make build          Build $(IMAGE):latest image"
	@echo "  make build-sha      Build and tag with the git commit SHA"
	@echo "  make run            Run the full benchmark (Track B + C + UQ + tests)"
	@echo "  make test           Run the pytest suite inside the container"
	@echo "  make shell          Open an interactive shell mounted to ./data ./output"
	@echo "  make multiomic      Train the ElasticNet multi-omic models"
	@echo "  make uq             Render the FNO-ensemble UQ trajectory"
	@echo "  make clean          Remove ./output/* (keeps data mount untouched)"
	@echo "  make prune          Remove dangling docker images"
	@echo ""
	@echo "Image tag in use: $(IMAGE):$(TAG)  (SHA=$(SHA))"

build:  ## Build the Docker image (latest tag)
	$(if $(filter $(shell uname),Linux),docker build -t $(IMAGE):latest .,\
	docker build -t $(IMAGE):latest .)

build-sha: build  ## Also tag with the git commit SHA
	@docker tag $(IMAGE):latest $(IMAGE):$(SHA)
	@echo "Tagged: $(IMAGE):$(SHA)"

run: build  ## Run the full containerized benchmark
	$(COMPOSE) run --rm benchmark

test: build  ## Run the pytest suite inside the container
	$(COMPOSE) run --rm --entrypoint /app/docker-entrypoint.sh benchmark --tests

shell: build  ## Open interactive bash with ./data and ./output mounted
	$(COMPOSE) run --rm ash

multiomic: build  ## Train ElasticNet multi-omic models (Proposal 2)
	$(COMPOSE) run --rm --entrypoint /app/docker-entrypoint.sh benchmark --train-multiomic

uq: build  ## Render FNO-ensemble UQ trajectory (Proposal 4)
	$(COMPOSE) run --rm --entrypoint /app/docker-entrypoint.sh benchmark --uq-ensemble

clean:  ## Remove generated artifacts in ./output
	@rm -rf output/* 2>/dev/null || true

prune:  ## Remove dangling docker images / build cache
	@docker image prune -f || true
