#!/bin/bash
#SBATCH --job-name=metabolic_pipeline
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=metabolic_pipeline.%j.out

set -euo pipefail

PROJECT_DIR="/path/to/Metagenomic Metabolic Profiling"
CONFIG_PATH="${PROJECT_DIR}/config/config.cluster.yaml"

cd "$PROJECT_DIR"

ENV_NAME="$(tr -d '\r\n' < .conda-env)"

conda run -n "$ENV_NAME" python main.py --config "$CONFIG_PATH"
