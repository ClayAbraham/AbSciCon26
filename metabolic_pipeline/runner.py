"""Workflow orchestration."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from metabolic_pipeline.config import Config
from metabolic_pipeline.ko_normalizer import ko_summary, normalize_kos
from metabolic_pipeline.merge_annotations import merge_annotations
from metabolic_pipeline.parsers.eggnog_parser import parse_eggnog_file
from metabolic_pipeline.parsers.kofam_parser import parse_kofam_output
from metabolic_pipeline.parsers.prodigal_parser import parse_prodigal_outputs
from metabolic_pipeline.pathway.definitions import load_pathway_definitions
from metabolic_pipeline.pathway.scorer import score_pathways
from metabolic_pipeline.reports import write_reports
from metabolic_pipeline.tools.kofamscan import run_kofamscan_for_sample
from metabolic_pipeline.tools.prodigal import run_prodigal_for_sample
from metabolic_pipeline.utils.environment import ensure_expected_conda_env
from metabolic_pipeline.utils.logging_utils import setup_logger
from metabolic_pipeline.utils.versions import collect_versions, write_versions


def _discover_contig_files(input_path: Path) -> list[Path]:
    if not input_path.exists():
        raise FileNotFoundError(f"Missing input path: {input_path}")
    if input_path.is_file():
        if input_path.suffix.lower() not in {".fasta", ".fa", ".fna"}:
            raise ValueError(f"Input file is not .fasta/.fa/.fna: {input_path}")
        return [input_path]

    candidates = sorted(
        list(input_path.glob("*.fasta")) + list(input_path.glob("*.fa")) + list(input_path.glob("*.fna"))
    )
    if not candidates:
        raise FileNotFoundError(f"No .fasta, .fa, or .fna files found in directory: {input_path}")
    return candidates


def _discover_eggnog_files(paths: list[Path]) -> list[Path]:
    discovered: list[Path] = []
    for p in paths:
        if p.is_file():
            discovered.append(p)
        elif p.is_dir():
            discovered.extend(sorted(p.glob("*.annotations")))
            discovered.extend(sorted(p.glob("*.tsv")))
    return discovered


def _sample_to_file_map(input_files: list[Path]) -> dict[str, Path]:
    mapping: dict[str, Path] = {}
    for path in input_files:
        sid = path.stem
        if sid in mapping:
            raise ValueError(f"Duplicate sample_id from file stems: {sid}")
        mapping[sid] = path
    return mapping


def _require_writable_output(output_dir: Path, overwrite: bool) -> None:
    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise FileExistsError(
            f"Output directory is not empty: {output_dir}. Use --overwrite to continue."
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    probe = output_dir / ".write_test"
    try:
        probe.write_text("ok", encoding="utf-8")
        probe.unlink()
    except Exception as exc:
        raise PermissionError(f"Write permission problem in output directory: {output_dir}") from exc


def run_pipeline(config: Config) -> None:
    """Run full workflow."""
    ensure_expected_conda_env(config.required_conda_env)
    _require_writable_output(config.output_dir, overwrite=config.overwrite)
    logger = setup_logger(config.output_dir)
    logger.info("Starting pipeline run.")
    logger.info("Input path: %s", config.input_path)

    contig_files = _discover_contig_files(config.input_path)
    sample_file_map = _sample_to_file_map(contig_files)
    logger.info("Discovered %d contig FASTA files.", len(contig_files))

    prodigal_dir = config.output_dir / "prodigal"
    prodigal_dir.mkdir(parents=True, exist_ok=True)

    fasta_frames: list[pd.DataFrame] = []
    protein_file_map: dict[str, Path] = {}
    for sample_id, contig_path in sample_file_map.items():
        logger.info("Running Prodigal for sample: %s", sample_id)
        protein_faa, prodigal_gff = run_prodigal_for_sample(
            contig_fasta=contig_path,
            sample_id=sample_id,
            output_dir=prodigal_dir,
            executable=config.prodigal_executable,
            mode=config.prodigal_mode,
        )
        protein_file_map[sample_id] = protein_faa
        fasta_frames.append(parse_prodigal_outputs(protein_faa, prodigal_gff, sample_id=sample_id))

    fasta_df = pd.concat(fasta_frames, ignore_index=True) if fasta_frames else pd.DataFrame()
    logger.info("Parsed Prodigal protein metadata rows: %d", len(fasta_df))

    kofam_dir = config.output_dir / "kofam"
    kofam_dir.mkdir(parents=True, exist_ok=True)
    kofam_frames: list[pd.DataFrame] = []
    for sample_id, _contig_path in sample_file_map.items():
        logger.info("Processing sample: %s", sample_id)
        if config.run_kofam:
            logger.info("Running KOfamScan for %s", sample_id)
            output_path = run_kofamscan_for_sample(
                faa_path=protein_file_map[sample_id],
                sample_id=sample_id,
                output_dir=kofam_dir,
                executable=config.kofamscan_executable,
                profiles_dir=config.kofam_profiles_dir or Path(),
                ko_list=config.kofam_ko_list or Path(),
                threads=config.threads,
            )
        else:
            output_path = kofam_dir / f"{sample_id}.kofam.tsv"
            if not output_path.exists():
                raise FileNotFoundError(
                    f"run_kofam is disabled and expected KOfam output is missing: {output_path}"
                )
        frame = parse_kofam_output(output_path, sample_id=sample_id)
        kofam_frames.append(frame)

    kofam_df = pd.concat(kofam_frames, ignore_index=True) if kofam_frames else pd.DataFrame()
    logger.info("KOfam parsed rows: %d", len(kofam_df))

    eggnog_df: pd.DataFrame | None = None
    if config.parse_eggnog:
        discovered = _discover_eggnog_files(config.eggnog_annotation_paths)
        if not discovered:
            raise FileNotFoundError(
                "parse_eggnog is enabled but no eggNOG annotation files were found."
            )
        frames: list[pd.DataFrame] = []
        for path in discovered:
            sid = path.stem.split(".")[0]
            frames.append(parse_eggnog_file(path, sample_id=sid))
        eggnog_df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        logger.info("eggNOG parsed rows: %d", len(eggnog_df))

    merged = merge_annotations(fasta_df=fasta_df, kofam_df=kofam_df, eggnog_df=eggnog_df)
    normalized = normalize_kos(merged)
    ko_sum = ko_summary(normalized)
    logger.info("Normalized KO rows: %d", len(normalized))

    pathway_defs = load_pathway_definitions(config.pathway_definition_file)
    all_sample_ids = sorted(fasta_df["sample_id"].dropna().unique().tolist())
    path_summary, path_long = score_pathways(normalized, pathway_defs, sample_ids=all_sample_ids)
    logger.info("Scored %d pathways across %d samples.", len(pathway_defs), merged["sample_id"].nunique())

    write_reports(
        output_dir=config.output_dir,
        annotations_merged=merged,
        ko_summary=ko_sum,
        pathway_summary=path_summary,
        pathway_long=path_long,
    )
    logger.info("Report writing complete.")

    snapshot_path = config.output_dir / "config_snapshot.yaml"
    snapshot_path.write_text(
        yaml.safe_dump(config.to_serializable_dict(), sort_keys=True),
        encoding="utf-8",
    )

    versions = collect_versions(
        config.prodigal_executable,
        config.kofamscan_executable,
        include_eggnog=config.parse_eggnog,
    )
    write_versions(config.output_dir / "versions.txt", versions)
    logger.info("Pipeline complete.")
