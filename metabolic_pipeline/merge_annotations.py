"""Annotation merge logic."""

from __future__ import annotations

import pandas as pd


MASTER_COLUMNS = [
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


def merge_annotations(
    fasta_df: pd.DataFrame,
    kofam_df: pd.DataFrame,
    eggnog_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Merge FASTA metadata with KOfam and optional eggNOG rows."""
    fasta_core = fasta_df[["sample_id", "protein_id", "contig_id"]].copy()

    merged_frames: list[pd.DataFrame] = []

    if not kofam_df.empty:
        k = fasta_core.merge(kofam_df, on=["sample_id", "protein_id"], how="right")
        for col in ["ec_number", "product", "preferred_name"]:
            k[col] = None
        merged_frames.append(k)

    if eggnog_df is not None and not eggnog_df.empty:
        e = fasta_core.merge(eggnog_df, on=["sample_id", "protein_id"], how="right")
        e["raw_annotation"] = None
        e["passed_threshold"] = None
        merged_frames.append(e)

    if merged_frames:
        merged = pd.concat(merged_frames, ignore_index=True, sort=False)
    else:
        merged = fasta_core.copy()
        for col in [
            "ko_id",
            "ec_number",
            "product",
            "preferred_name",
            "source_tool",
            "score",
            "evalue",
            "passed_threshold",
            "raw_annotation",
        ]:
            merged[col] = None

    for col in MASTER_COLUMNS:
        if col not in merged.columns:
            merged[col] = None
    return merged[MASTER_COLUMNS].copy()

