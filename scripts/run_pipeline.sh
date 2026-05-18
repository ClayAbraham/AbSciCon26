#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${PROJECT_ROOT}"

ENV_NAME="$(tr -d '\r\n' < "${PROJECT_ROOT}/.conda-env")"

if [[ -z "${ENV_NAME}" ]]; then
  echo "Missing conda environment name in .conda-env" >&2
  exit 1
fi

conda run -n "${ENV_NAME}" python "${PROJECT_ROOT}/main.py" "$@"
