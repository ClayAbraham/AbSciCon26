"""Report writing helpers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ANNOTATIONS_COLUMNS = [
    "sample_id",
    "protein_id",
    "contig_id",
    "ko_id",
    "ec_number",
    "product",
    "preferred_name",
    "source_tool",
    "score",
    "evalue",
    "passed_threshold",
    "raw_annotation",
]

KO_SUMMARY_COLUMNS = ["sample_id", "ko_id", "protein_count", "contig_count"]

PATHWAY_SUMMARY_COLUMNS = [
    "sample_id",
    "pathway_id",
    "pathway_name",
    "required_total",
    "required_found",
    "required_missing",
    "optional_found",
    "marker_found",
    "completeness_fraction",
    "completeness_percent",
    "status_label",
    "confidence_label",
]

PATHWAY_LONG_COLUMNS = ["sample_id", "pathway_id", "pathway_name", "ko_id", "ko_role", "present"]


def _ensure_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        if col not in out.columns:
            out[col] = None
    return out[columns]


def write_reports(
    output_dir: Path,
    annotations_merged: pd.DataFrame,
    ko_summary: pd.DataFrame,
    pathway_summary: pd.DataFrame,
    pathway_long: pd.DataFrame,
) -> None:
    """Write all required TSV report files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    ann = _ensure_columns(annotations_merged, ANNOTATIONS_COLUMNS).sort_values(
        by=["sample_id", "protein_id", "source_tool", "ko_id"], na_position="last"
    )
    ko = _ensure_columns(ko_summary, KO_SUMMARY_COLUMNS).sort_values(by=["sample_id", "ko_id"])
    psum = _ensure_columns(pathway_summary, PATHWAY_SUMMARY_COLUMNS).sort_values(
        by=["sample_id", "pathway_id"]
    )
    plong = _ensure_columns(pathway_long, PATHWAY_LONG_COLUMNS).sort_values(
        by=["sample_id", "pathway_id", "ko_id", "ko_role"]
    )

    ann.to_csv(output_dir / "annotations_merged.tsv", sep="\t", index=False)
    ko.to_csv(output_dir / "ko_summary.tsv", sep="\t", index=False)
    psum.to_csv(output_dir / "pathway_summary.tsv", sep="\t", index=False)
    plong.to_csv(output_dir / "pathway_long.tsv", sep="\t", index=False)

