from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import ProjectConfig, SheetConfig


def list_workbook_sheets(workbook_path: Path) -> pd.DataFrame:
    excel_file = pd.ExcelFile(workbook_path)
    rows: list[dict[str, object]] = []
    for sheet_name in excel_file.sheet_names:
        frame = excel_file.parse(sheet_name=sheet_name)
        rows.append(
            {
                "sheet_name": sheet_name,
                "n_rows": int(frame.shape[0]),
                "n_columns": int(frame.shape[1]),
                "columns": ", ".join(str(column) for column in frame.columns),
            }
        )
    return pd.DataFrame(rows)


def load_sheet(workbook_path: Path, sheet_name: str) -> pd.DataFrame:
    return pd.read_excel(workbook_path, sheet_name=sheet_name)


def load_metadata(config: ProjectConfig, workbook_path: Path) -> pd.DataFrame | None:
    if not config.metadata_sheet:
        return None
    return load_sheet(workbook_path, config.metadata_sheet)


def _combine_feature_columns(clean: pd.DataFrame, sheet_config: SheetConfig) -> pd.DataFrame:
    if not sheet_config.combine_feature_columns:
        return clean

    if sheet_config.feature_column not in clean.columns:
        return clean

    feature_position = list(clean.columns).index(sheet_config.feature_column)
    reserved_columns = {
        sheet_config.sample_column,
        sheet_config.value_column,
        *sheet_config.include_columns,
        *sheet_config.exclude_columns,
    }
    candidate_columns = []
    for column in list(clean.columns)[feature_position:]:
        if column in reserved_columns:
            continue
        candidate_columns.append(column)

    if not candidate_columns:
        return clean

    feature_frame = clean[candidate_columns].fillna("")
    clean[sheet_config.feature_column] = (
        feature_frame.astype(str)
        .apply(lambda row: " ".join(item.strip() for item in row if str(item).strip()), axis=1)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )
    return clean


def normalize_sheet_to_matrix(frame: pd.DataFrame, sheet_config: SheetConfig) -> pd.DataFrame:
    clean = frame.copy()
    clean.columns = [str(column).strip() for column in clean.columns]
    clean = _combine_feature_columns(clean, sheet_config)

    if sheet_config.feature_column not in clean.columns:
        raise ValueError(
            f"Feature column '{sheet_config.feature_column}' was not found in sheet "
            f"'{sheet_config.sheet_name}'. Available columns: {list(clean.columns)}"
        )

    if sheet_config.sample_column and sheet_config.value_column:
        missing = [
            column
            for column in (sheet_config.sample_column, sheet_config.value_column)
            if column not in clean.columns
        ]
        if missing:
            raise ValueError(
                f"Missing required columns {missing} in sheet '{sheet_config.sheet_name}'."
            )

        matrix = (
            clean[[sheet_config.feature_column, sheet_config.sample_column, sheet_config.value_column]]
            .dropna(subset=[sheet_config.feature_column, sheet_config.sample_column])
            .assign(
                **{
                    sheet_config.feature_column: lambda df: df[sheet_config.feature_column]
                    .astype(str)
                    .str.strip(),
                    sheet_config.sample_column: lambda df: df[sheet_config.sample_column]
                    .astype(str)
                    .str.strip(),
                    sheet_config.value_column: lambda df: pd.to_numeric(
                        df[sheet_config.value_column], errors="coerce"
                    ).fillna(0.0),
                }
            )
            .pivot_table(
                index=sheet_config.feature_column,
                columns=sheet_config.sample_column,
                values=sheet_config.value_column,
                aggfunc="sum",
                fill_value=0.0,
            )
        )
        return matrix.sort_index()

    candidate_columns = [
        column
        for column in clean.columns
        if column != sheet_config.feature_column and column not in set(sheet_config.exclude_columns)
    ]
    if sheet_config.include_columns:
        candidate_columns = [
            column for column in candidate_columns if column in set(sheet_config.include_columns)
        ]

    if not candidate_columns:
        raise ValueError(
            f"No candidate sample columns were found for sheet '{sheet_config.sheet_name}'."
        )

    working = clean[[sheet_config.feature_column, *candidate_columns]].copy()
    working[sheet_config.feature_column] = working[sheet_config.feature_column].astype(str).str.strip()
    for column in candidate_columns:
        working[column] = pd.to_numeric(working[column], errors="coerce")

    numeric_columns = [column for column in candidate_columns if working[column].notna().any()]
    if not numeric_columns:
        raise ValueError(
            f"No numeric sample columns were detected in sheet '{sheet_config.sheet_name}'."
        )

    matrix = (
        working[[sheet_config.feature_column, *numeric_columns]]
        .dropna(subset=[sheet_config.feature_column])
        .groupby(sheet_config.feature_column, dropna=False)[numeric_columns]
        .sum(min_count=1)
        .fillna(0.0)
    )
    matrix.columns = [str(column).strip() for column in matrix.columns]
    if sheet_config.sample_id and len(matrix.columns) == 1:
        matrix.columns = [sheet_config.sample_id]
    return matrix.sort_index()
