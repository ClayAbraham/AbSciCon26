"""Version capture utilities."""

from __future__ import annotations

import platform
from pathlib import Path

from metabolic_pipeline.utils.environment import current_conda_env_name
from metabolic_pipeline.utils.subprocess_utils import run_command


def _safe_tool_version(command: list[str]) -> str:
    try:
        result = run_command(command)
        first_line = result.stdout.strip().splitlines()[0] if result.stdout.strip() else "unknown"
        return first_line
    except Exception:
        return "unavailable"


def collect_versions(
    prodigal_executable: str,
    kofamscan_executable: str,
    include_eggnog: bool,
) -> dict[str, str]:
    """Collect tool and runtime versions."""
    versions = {
        "python": platform.python_version(),
        "conda_env": current_conda_env_name() or "unavailable",
        "prodigal": _safe_tool_version([prodigal_executable, "-v"]),
        "kofamscan": _safe_tool_version([kofamscan_executable, "--version"]),
    }
    if include_eggnog:
        versions["eggnog_mapper"] = _safe_tool_version(["emapper.py", "--version"])
    return versions


def write_versions(path: Path, versions: dict[str, str]) -> None:
    """Write versions to text file."""
    lines = [f"{key}\t{value}" for key, value in versions.items()]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
