from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class AnalysisConfig:
    top_n_bar: int = 15
    top_n_heatmap: int = 25
    min_prevalence_fraction: float = 0.0
    correlation_method: str = "pearson"
    pca_components: int = 2


@dataclass
class SheetConfig:
    dataset_name: str
    sheet_name: str
    dataset_type: str
    feature_column: str
    sample_id: str | None = None
    sample_column: str | None = None
    value_column: str | None = None
    combine_feature_columns: bool = False
    include_columns: list[str] = field(default_factory=list)
    exclude_columns: list[str] = field(default_factory=list)


@dataclass
class ProjectConfig:
    config_path: Path
    workbook_path: str
    output_dir: str
    sheets: list[SheetConfig]
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)
    metadata_sheet: str | None = None
    sample_id_column: str | None = None
    group_column: str | None = None

    @property
    def base_dir(self) -> Path:
        return self.config_path.parent

    def resolve_path(self, value: str) -> Path:
        path = Path(value)
        if path.is_absolute():
            return path
        return (self.base_dir / path).resolve()


def _load_sheet_config(data: dict[str, Any]) -> SheetConfig:
    return SheetConfig(
        dataset_name=str(data["dataset_name"]),
        sheet_name=str(data["sheet_name"]),
        dataset_type=str(data["dataset_type"]),
        feature_column=str(data["feature_column"]),
        sample_id=str(data["sample_id"]) if data.get("sample_id") else None,
        sample_column=str(data["sample_column"]) if data.get("sample_column") else None,
        value_column=str(data["value_column"]) if data.get("value_column") else None,
        combine_feature_columns=bool(data.get("combine_feature_columns", False)),
        include_columns=[str(item) for item in data.get("include_columns", [])],
        exclude_columns=[str(item) for item in data.get("exclude_columns", [])],
    )


def load_config(config_path: str | Path) -> ProjectConfig:
    config_path = Path(config_path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}

    sheets = [_load_sheet_config(item) for item in raw.get("sheets", [])]
    if not sheets:
        raise ValueError("Config must include at least one sheet definition under 'sheets'.")

    analysis_raw = raw.get("analysis", {})
    analysis = AnalysisConfig(
        top_n_bar=int(analysis_raw.get("top_n_bar", 15)),
        top_n_heatmap=int(analysis_raw.get("top_n_heatmap", 25)),
        min_prevalence_fraction=float(analysis_raw.get("min_prevalence_fraction", 0.0)),
        correlation_method=str(analysis_raw.get("correlation_method", "pearson")),
        pca_components=int(analysis_raw.get("pca_components", 2)),
    )

    return ProjectConfig(
        config_path=config_path,
        workbook_path=str(raw["workbook_path"]),
        output_dir=str(raw.get("output_dir", "output")),
        sheets=sheets,
        analysis=analysis,
        metadata_sheet=str(raw["metadata_sheet"]) if raw.get("metadata_sheet") else None,
        sample_id_column=str(raw["sample_id_column"]) if raw.get("sample_id_column") else None,
        group_column=str(raw["group_column"]) if raw.get("group_column") else None,
    )
