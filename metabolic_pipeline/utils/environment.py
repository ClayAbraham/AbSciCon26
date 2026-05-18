"""Conda environment helpers."""

from __future__ import annotations

import os
from pathlib import Path


def project_root() -> Path:
    """Return the repository root."""
    return Path(__file__).resolve().parents[2]


def load_project_conda_env_name() -> str | None:
    """Load the expected conda environment name from the project marker file."""
    env_file = project_root() / ".conda-env"
    if not env_file.exists():
        return None
    value = env_file.read_text(encoding="utf-8").strip()
    return value or None


def current_conda_env_name() -> str | None:
    """Detect the active conda environment name from shell variables."""
    explicit = os.environ.get("CONDA_DEFAULT_ENV")
    if explicit:
        return explicit

    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        return Path(prefix).name

    return None


def ensure_expected_conda_env(required_env: str | None) -> None:
    """Require that the active conda environment matches the configured one."""
    if not required_env:
        return

    current_env = current_conda_env_name()
    if current_env == required_env:
        return

    if current_env is None:
        raise EnvironmentError(
            "This project must run inside the conda environment "
            f"'{required_env}', but no active conda environment was detected."
        )

    raise EnvironmentError(
        "This project must run inside the conda environment "
        f"'{required_env}', but the active environment is '{current_env}'."
    )
