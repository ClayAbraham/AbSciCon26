"""KOfamScan wrapper."""

from __future__ import annotations

import shutil
from pathlib import Path

from metabolic_pipeline.utils.subprocess_utils import run_command


def _validate_kofam_paths(
    executable: str,
    profiles_dir: Path | None,
    ko_list: Path | None,
) -> None:
    exe_path = shutil.which(executable) if Path(executable).name == executable else executable
    if exe_path is None and not Path(executable).exists():
        raise FileNotFoundError(f"Missing KOfam executable: {executable}")
    if profiles_dir is None or not profiles_dir.exists():
        raise FileNotFoundError(f"Missing KOfam profiles directory: {profiles_dir}")
    if ko_list is None or not ko_list.exists():
        raise FileNotFoundError(f"Missing KOfam ko_list file: {ko_list}")


def run_kofamscan_for_sample(
    faa_path: Path,
    sample_id: str,
    output_dir: Path,
    executable: str,
    profiles_dir: Path,
    ko_list: Path,
    threads: int,
) -> Path:
    """Run KOfamScan on a single sample FASTA."""
    _validate_kofam_paths(executable, profiles_dir, ko_list)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{sample_id}.kofam.tsv"
    command = [
        executable,
        "--cpu",
        str(threads),
        "-p",
        str(profiles_dir),
        "-k",
        str(ko_list),
        "-o",
        str(out_path),
        str(faa_path),
    ]
    run_command(command)
    if not out_path.exists():
        raise FileNotFoundError(f"Missing annotation output: {out_path}")
    return out_path

