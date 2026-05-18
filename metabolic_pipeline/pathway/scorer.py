"""Pathway scoring logic."""

from __future__ import annotations

import pandas as pd

from metabolic_pipeline.pathway.definitions import PathwayDefinition


def _status_label(required_total: int, required_found: int, marker_found: int) -> str:
    if required_total == 0 and marker_found > 0:
        return "marker_only"
    if required_total == 0:
        return "not_detected"
    completeness = required_found / required_total
    if required_found == required_total:
        return "complete"
    if completeness >= 0.8:
        return "near_complete"
    if required_found > 0:
        return "partial"
    if marker_found > 0:
        return "marker_only"
    return "not_detected"


def _confidence_label(status: str, marker_found: int) -> str:
    if status == "complete" and marker_found > 0:
        return "high"
    if status == "near_complete":
        return "medium"
    if status == "partial" and marker_found > 0:
        return "medium"
    if status in {"partial", "marker_only"}:
        return "low"
    return "none"


def score_pathways(
    normalized_df: pd.DataFrame,
    pathway_definitions: list[PathwayDefinition],
    sample_ids: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score pathways by sample and create long KO presence table."""
    if sample_ids is None:
        sample_ids = (
            sorted(normalized_df["sample_id"].dropna().unique().tolist()) if not normalized_df.empty else []
        )

    summary_rows: list[dict[str, object]] = []
    long_rows: list[dict[str, object]] = []

    sample_kos: dict[str, set[str]] = {}
    for sid in sample_ids:
        kos = set(normalized_df.loc[normalized_df["sample_id"] == sid, "ko_id"].dropna().astype(str).tolist())
        sample_kos[sid] = kos

    for sid in sample_ids:
        kos = sample_kos[sid]
        for definition in pathway_definitions:
            required = set(definition.required_kos)
            optional = set(definition.optional_kos)
            markers = set(definition.marker_kos)

            req_found = len(required.intersection(kos))
            req_total = len(required)
            req_missing = req_total - req_found
            opt_found = len(optional.intersection(kos))
            marker_found = len(markers.intersection(kos))
            completeness_fraction = (req_found / req_total) if req_total else 0.0
            completeness_percent = completeness_fraction * 100.0
            status = _status_label(req_total, req_found, marker_found)
            confidence = _confidence_label(status, marker_found)

            summary_rows.append(
                {
                    "sample_id": sid,
                    "pathway_id": definition.pathway_id,
                    "pathway_name": definition.pathway_name,
                    "required_total": req_total,
                    "required_found": req_found,
                    "required_missing": req_missing,
                    "optional_found": opt_found,
                    "marker_found": marker_found,
                    "completeness_fraction": completeness_fraction,
                    "completeness_percent": completeness_percent,
                    "status_label": status,
                    "confidence_label": confidence,
                }
            )

            role_map = {
                "required": definition.required_kos,
                "optional": definition.optional_kos,
                "marker": definition.marker_kos,
            }
            for role, ko_list in role_map.items():
                for ko in ko_list:
                    long_rows.append(
                        {
                            "sample_id": sid,
                            "pathway_id": definition.pathway_id,
                            "pathway_name": definition.pathway_name,
                            "ko_id": ko,
                            "ko_role": role,
                            "present": 1 if ko in kos else 0,
                        }
                    )

    summary_df = pd.DataFrame(summary_rows)
    long_df = pd.DataFrame(long_rows)
    if summary_df.empty:
        summary_df = pd.DataFrame(
            columns=[
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
        )
    if long_df.empty:
        long_df = pd.DataFrame(
            columns=["sample_id", "pathway_id", "pathway_name", "ko_id", "ko_role", "present"]
        )
    return summary_df, long_df
