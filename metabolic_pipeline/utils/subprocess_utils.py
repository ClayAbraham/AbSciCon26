"""Subprocess helpers."""

from __future__ import annotations

import subprocess
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


@dataclass
class CommandResult:
    """Captured subprocess output."""

    returncode: int
    stdout: str
    stderr: str


def run_command(
    command: Sequence[str],
    cwd: Path | None = None,
    dry_run: bool = False,
) -> CommandResult:
    """Run a shell command with captured output and clear failures."""
    logging.getLogger("metabolic_pipeline").info("CMD: %s", " ".join(command))
    if dry_run:
        return CommandResult(returncode=0, stdout="DRY_RUN", stderr="")
    proc = subprocess.run(
        list(command),
        cwd=str(cwd) if cwd else None,
        text=True,
        capture_output=True,
        check=False,
    )
    result = CommandResult(returncode=proc.returncode, stdout=proc.stdout, stderr=proc.stderr)
    if result.returncode != 0:
        joined = " ".join(command)
        raise RuntimeError(
            f"Subprocess failed (code {result.returncode}): {joined}\n"
            f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    return result
