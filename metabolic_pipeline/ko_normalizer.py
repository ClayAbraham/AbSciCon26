"""KO normalization and summary generation."""

from __future__ import annotations

import re

import pandas as pd

KO_RE = re.compile(r"^K\d{5}$")
KO_LIKE_RE = re.compile(r"^K\d+$")


def normalize_kos(annotation_df: pd.DataFrame) -> pd.DataFrame:
    """Normalize KO values to one valid KO per row."""
    if annotation_df.empty:
        return annotation_df.copy()
    rows: list[dict[str, object]] = []
    for _, row in annotation_df.iterrows():
        raw_ko = row.get("ko_id")
        if raw_ko is None or (isinstance(raw_ko, float) and pd.isna(raw_ko)):
            continue
        parts = re.split(r"[,\s;|]+", str(raw_ko).strip())
        for part in parts:
            if not part:
                continue
            if KO_LIKE_RE.match(part) and not KO_RE.match(part):
                raise ValueError(f"Invalid KO format: {part}")
            if not KO_RE.match(part):
                continue
            new_row = row.to_dict()
            new_row["ko_id"] = part
            rows.append(new_row)
    if not rows:
        return pd.DataFrame(columns=annotation_df.columns)
    out = pd.DataFrame(rows)
    out = out.drop_duplicates(subset=["sample_id", "protein_id", "ko_id", "source_tool"])
    return out.reset_index(drop=True)


def ko_summary(normalized_df: pd.DataFrame) -> pd.DataFrame:
    """Create KO summary table by sample."""
    if normalized_df.empty:
        return pd.DataFrame(columns=["sample_id", "ko_id", "protein_count", "contig_count"])
    grouped = (
        normalized_df.groupby(["sample_id", "ko_id"], dropna=False)
        .agg(
            protein_count=("protein_id", "nunique"),
            contig_count=("contig_id", lambda s: s.dropna().nunique()),
        )
        .reset_index()
    )
    return grouped[["sample_id", "ko_id", "protein_count", "contig_count"]]
