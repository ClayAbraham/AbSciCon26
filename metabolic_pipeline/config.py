"""Configuration handling for the metabolic profiling pipeline."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class Config:
    """Runtime configuration."""

    input_path: Path
    output_dir: Path
    threads: int
    run_kofam: bool
    prodigal_executable: str
    prodigal_mode: str
    kofamscan_executable: str
    kofam_profiles_dir: Path | None
    kofam_ko_list: Path | None
    parse_eggnog: bool
    eggnog_annotation_paths: list[Path]
    pathway_definition_file: Path
    overwrite: bool
    required_conda_env: str | None

    def to_serializable_dict(self) -> dict[str, Any]:
        """Serialize config to primitives for snapshot writing."""
        out = asdict(self)
        for key, value in out.items():
            if isinstance(value, Path):
                out[key] = str(value)
            if isinstance(value, list):
                out[key] = [str(v) if isinstance(v, Path) else v for v in value]
        return out


def _load_yaml(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError("Config YAML root must be a mapping.")
    return data


def _to_path_or_none(value: Any) -> Path | None:
    if value is None or value == "":
        return None
    return Path(str(value))


def _to_path_list(value: Any) -> list[Path]:
    if value is None:
        return []
    if isinstance(value, (str, Path)):
        return [Path(str(value))]
    if isinstance(value, list):
        return [Path(str(v)) for v in value]
    raise ValueError("eggnog_annotation_paths must be a path or list of paths.")


def build_config(
    cli_input_path: Path | None,
    cli_output_dir: Path | None,
    cli_config_path: Path | None,
    cli_threads: int | None,
    cli_pathways: Path | None,
    cli_run_kofam: bool,
    cli_parse_eggnog: bool,
    cli_eggnog_paths: list[Path] | None,
    cli_overwrite: bool,
) -> Config:
    """Merge config YAML with CLI overrides and validate."""
    raw = _load_yaml(cli_config_path)

    input_path = cli_input_path or _to_path_or_none(raw.get("input_path"))
    output_dir = cli_output_dir or _to_path_or_none(raw.get("output_dir"))
    threads = cli_threads if cli_threads is not None else int(raw.get("threads", 1))
    run_kofam = bool(raw.get("run_kofam", False)) or cli_run_kofam
    parse_eggnog = bool(raw.get("parse_eggnog", False)) or cli_parse_eggnog
    overwrite = bool(raw.get("overwrite", False)) or cli_overwrite
    required_conda_env_raw = raw.get("required_conda_env")
    required_conda_env = (
        str(required_conda_env_raw).strip() if required_conda_env_raw not in (None, "") else None
    )

    pathway_definition_file = cli_pathways or _to_path_or_none(raw.get("pathway_definition_file"))
    prodigal_executable = str(raw.get("prodigal_executable", "prodigal"))
    prodigal_mode = str(raw.get("prodigal_mode", "meta")).strip() or "meta"
    kofamscan_executable = str(raw.get("kofamscan_executable", "exec_annotation"))
    kofam_profiles_dir = _to_path_or_none(raw.get("kofam_profiles_dir"))
    kofam_ko_list = _to_path_or_none(raw.get("kofam_ko_list"))

    eggnog_annotation_paths = _to_path_list(raw.get("eggnog_annotation_paths"))
    if cli_eggnog_paths:
        eggnog_annotation_paths = cli_eggnog_paths

    if input_path is None:
        raise ValueError("Missing required setting: input_path (or --input).")
    if output_dir is None:
        raise ValueError("Missing required setting: output_dir (or --output).")
    if pathway_definition_file is None:
        raise ValueError("Missing required setting: pathway_definition_file (or --pathways).")
    if threads < 1:
        raise ValueError("threads must be >= 1.")
    if prodigal_mode not in {"single", "meta"}:
        raise ValueError("prodigal_mode must be 'single' or 'meta'.")

    if run_kofam and (kofam_profiles_dir is None or kofam_ko_list is None):
        raise ValueError(
            "run_kofam is enabled but kofam_profiles_dir or kofam_ko_list is missing."
        )

    if parse_eggnog and not eggnog_annotation_paths:
        # Allowed: parse later from directory in runner if empty and input is dir.
        eggnog_annotation_paths = []

    return Config(
        input_path=input_path,
        output_dir=output_dir,
        threads=threads,
        run_kofam=run_kofam,
        prodigal_executable=prodigal_executable,
        prodigal_mode=prodigal_mode,
        kofamscan_executable=kofamscan_executable,
        kofam_profiles_dir=kofam_profiles_dir,
        kofam_ko_list=kofam_ko_list,
        parse_eggnog=parse_eggnog,
        eggnog_annotation_paths=eggnog_annotation_paths,
        pathway_definition_file=pathway_definition_file,
        overwrite=overwrite,
        required_conda_env=required_conda_env,
    )
