from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .analysis import (
    compute_feature_statistics,
    compute_group_means,
    compute_group_tests,
    compute_pca,
    compute_sample_correlation,
    compute_sample_distances,
    compute_sample_statistics,
    filter_matrix_by_prevalence,
)
from .config import ProjectConfig, load_config
from .io_utils import list_workbook_sheets, load_metadata, load_sheet, normalize_sheet_to_matrix
from .plotting import (
    plot_correlation_heatmap,
    plot_heatmap,
    plot_pca,
    plot_sample_totals,
    plot_top_features_bar,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Excel-driven metabolic profiling analysis.")
    parser.add_argument(
        "--config",
        default="Metabolic_Profiling_v2/config.example.yaml",
        help="Path to the YAML config file.",
    )
    return parser.parse_args()


def _write_table(frame: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, sep="\t")


def _build_group_series(
    config: ProjectConfig,
    metadata: pd.DataFrame | None,
) -> pd.Series | None:
    if (
        metadata is None
        or not config.sample_id_column
        or not config.group_column
        or config.sample_id_column not in metadata.columns
        or config.group_column not in metadata.columns
    ):
        return None

    groups = metadata[[config.sample_id_column, config.group_column]].dropna().copy()
    groups[config.sample_id_column] = groups[config.sample_id_column].astype(str).str.strip()
    groups[config.group_column] = groups[config.group_column].astype(str).str.strip()
    return groups.drop_duplicates(subset=[config.sample_id_column]).set_index(config.sample_id_column)[
        config.group_column
    ]


def run_from_config(config: ProjectConfig) -> int:
    workbook_path = config.resolve_path(config.workbook_path)
    output_dir = config.resolve_path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not workbook_path.exists():
        raise FileNotFoundError(
            f"Workbook not found: {workbook_path}. Place it in the input folder or update the config."
        )

    workbook_overview = list_workbook_sheets(workbook_path)
    _write_table(workbook_overview, output_dir / "workbook_sheet_overview.tsv")

    metadata = load_metadata(config, workbook_path)
    group_series = _build_group_series(config, metadata)

    dataset_summary_rows: list[dict[str, object]] = []

    for sheet_config in config.sheets:
        raw_sheet = load_sheet(workbook_path, sheet_config.sheet_name)
        matrix = normalize_sheet_to_matrix(raw_sheet, sheet_config)
        matrix = filter_matrix_by_prevalence(matrix, config.analysis.min_prevalence_fraction)

        dataset_dir = output_dir / sheet_config.dataset_name
        dataset_dir.mkdir(parents=True, exist_ok=True)

        _write_table(matrix, dataset_dir / "normalized_matrix.tsv")

        feature_stats = compute_feature_statistics(matrix)
        sample_stats = compute_sample_statistics(matrix)
        correlation = compute_sample_correlation(matrix, method=config.analysis.correlation_method)
        distances = compute_sample_distances(matrix)

        _write_table(feature_stats, dataset_dir / "feature_statistics.tsv")
        _write_table(sample_stats, dataset_dir / "sample_statistics.tsv")
        _write_table(correlation, dataset_dir / "sample_correlation.tsv")
        _write_table(distances, dataset_dir / "sample_distances.tsv")

        plot_top_features_bar(
            feature_stats=feature_stats,
            output_path=dataset_dir / "top_features_bar.png",
            title=f"{sheet_config.dataset_name}: top pathways/features",
            top_n=config.analysis.top_n_bar,
            metric="mean",
        )
        plot_heatmap(
            matrix=matrix,
            feature_stats=feature_stats,
            output_path=dataset_dir / "top_features_heatmap.png",
            title=f"{sheet_config.dataset_name}: top feature heatmap",
            top_n=config.analysis.top_n_heatmap,
        )
        plot_sample_totals(
            sample_stats=sample_stats,
            output_path=dataset_dir / "sample_totals.png",
            title=f"{sheet_config.dataset_name}: total signal by sample",
        )
        plot_correlation_heatmap(
            correlation=correlation,
            output_path=dataset_dir / "sample_correlation_heatmap.png",
            title=f"{sheet_config.dataset_name}: sample correlation",
        )

        pca_result = compute_pca(matrix, n_components=config.analysis.pca_components)
        if pca_result is not None:
            pca_frame, explained = pca_result
            _write_table(pca_frame, dataset_dir / "sample_pca.tsv")
            plot_pca(
                coordinates=pca_frame,
                explained_variance=explained,
                output_path=dataset_dir / "sample_pca.png",
                title=f"{sheet_config.dataset_name}: PCA of samples",
                group_series=group_series,
            )

        if (
            metadata is not None
            and config.sample_id_column
            and config.group_column
        ):
            group_means = compute_group_means(
                matrix=matrix,
                metadata=metadata,
                sample_id_column=config.sample_id_column,
                group_column=config.group_column,
            )
            if group_means is not None:
                _write_table(group_means, dataset_dir / "group_means.tsv")

            group_tests = compute_group_tests(
                matrix=matrix,
                metadata=metadata,
                sample_id_column=config.sample_id_column,
                group_column=config.group_column,
            )
            if group_tests is not None:
                _write_table(group_tests, dataset_dir / "group_tests.tsv")

        dataset_summary_rows.append(
            {
                "dataset_name": sheet_config.dataset_name,
                "sheet_name": sheet_config.sheet_name,
                "dataset_type": sheet_config.dataset_type,
                "n_features": int(matrix.shape[0]),
                "n_samples": int(matrix.shape[1]),
                "top_feature": str(feature_stats.index[0]) if not feature_stats.empty else "",
                "top_feature_mean": float(feature_stats.iloc[0]["mean"]) if not feature_stats.empty else 0.0,
            }
        )

    summary_frame = pd.DataFrame(dataset_summary_rows)
    _write_table(summary_frame, output_dir / "analysis_summary.tsv")
    return 0


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    return run_from_config(config)
