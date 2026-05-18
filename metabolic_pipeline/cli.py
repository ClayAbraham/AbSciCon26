"""Command line interface for the metabolic profiling pipeline."""

from __future__ import annotations

import argparse
from pathlib import Path

from metabolic_pipeline.config import Config, build_config
from metabolic_pipeline.utils.environment import ensure_expected_conda_env, load_project_conda_env_name
from metabolic_pipeline.runner import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser."""
    parser = argparse.ArgumentParser(
        prog="metabolic-pipeline",
        description="Predict proteins from contig FASTA files, assign KOs, and score pathways.",
    )
    parser.add_argument("--input", dest="input_path", type=Path, required=False)
    parser.add_argument("--output", dest="output_dir", type=Path, required=False)
    parser.add_argument("--config", dest="config_path", type=Path, required=False)
    parser.add_argument("--threads", dest="threads", type=int, required=False)
    parser.add_argument("--pathways", dest="pathway_definition_file", type=Path, required=False)
    parser.add_argument("--run-kofam", dest="run_kofam", action="store_true")
    parser.add_argument("--parse-eggnog", dest="parse_eggnog", nargs="?", const=True, default=None)
    parser.add_argument("--overwrite", dest="overwrite", action="store_true")
    return parser


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    return build_parser().parse_args()


def _resolve_parse_eggnog_value(value: object) -> tuple[bool, list[Path] | None]:
    if value is None:
        return False, None
    if value is True:
        return True, None
    return True, [Path(str(value))]


def main() -> None:
    """CLI entry point."""
    ensure_expected_conda_env(load_project_conda_env_name())
    args = parse_args()
    parse_eggnog_flag, eggnog_paths = _resolve_parse_eggnog_value(args.parse_eggnog)
    config: Config = build_config(
        cli_input_path=args.input_path,
        cli_output_dir=args.output_dir,
        cli_config_path=args.config_path,
        cli_threads=args.threads,
        cli_pathways=args.pathway_definition_file,
        cli_run_kofam=args.run_kofam,
        cli_parse_eggnog=parse_eggnog_flag,
        cli_eggnog_paths=eggnog_paths,
        cli_overwrite=args.overwrite,
    )
    run_pipeline(config)
