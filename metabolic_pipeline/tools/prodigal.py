"""Prodigal wrapper."""

from __future__ import annotations

import shutil
from pathlib import Path

from metabolic_pipeline.utils.subprocess_utils import run_command


def _validate_prodigal_executable(executable: str) -> None:
    exe_path = shutil.which(executable) if Path(executable).name == executable else executable
    if exe_path is None and not Path(executable).exists():
        raise FileNotFoundError(f"Missing Prodigal executable: {executable}")


def run_prodigal_for_sample(
    contig_fasta: Path,
    sample_id: str,
    output_dir: Path,
    executable: str,
    mode: str,
) -> tuple[Path, Path]:
    """Run Prodigal on one contig FASTA and return protein FASTA plus GFF paths."""
    _validate_prodigal_executable(executable)
    output_dir.mkdir(parents=True, exist_ok=True)

    proteins_path = output_dir / f"{sample_id}.faa"
    gff_path = output_dir / f"{sample_id}.gff"

    command = [
        executable,
        "-i",
        str(contig_fasta),
        "-a",
        str(proteins_path),
        "-o",
        str(gff_path),
        "-p",
        mode,
        "-q",
    ]
    run_command(command)

    if not proteins_path.exists():
        raise FileNotFoundError(f"Missing Prodigal protein FASTA: {proteins_path}")
    if not gff_path.exists():
        raise FileNotFoundError(f"Missing Prodigal GFF output: {gff_path}")
    return proteins_path, gff_path
