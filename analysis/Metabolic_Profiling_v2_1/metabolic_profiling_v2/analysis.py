from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f_oneway, kruskal
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def filter_matrix_by_prevalence(matrix: pd.DataFrame, minimum_fraction: float) -> pd.DataFrame:
    if minimum_fraction <= 0:
        return matrix
    prevalence = (matrix > 0).sum(axis=1) / max(matrix.shape[1], 1)
    return matrix.loc[prevalence >= minimum_fraction].copy()


def compute_feature_statistics(matrix: pd.DataFrame) -> pd.DataFrame:
    prevalence_count = (matrix > 0).sum(axis=1)
    prevalence_fraction = prevalence_count / max(matrix.shape[1], 1)
    mean_values = matrix.mean(axis=1)
    std_values = matrix.std(axis=1, ddof=1).fillna(0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        coefficient_of_variation = (std_values / mean_values.replace(0, np.nan)).replace(
            [np.inf, -np.inf], np.nan
        )

    stats = pd.DataFrame(
        {
            "sum": matrix.sum(axis=1),
            "mean": mean_values,
            "median": matrix.median(axis=1),
            "std": std_values,
            "min": matrix.min(axis=1),
            "max": matrix.max(axis=1),
            "prevalence_count": prevalence_count,
            "prevalence_fraction": prevalence_fraction,
            "coefficient_of_variation": coefficient_of_variation.fillna(0.0),
        }
    )
    return stats.sort_values(["mean", "sum"], ascending=[False, False])


def compute_sample_statistics(matrix: pd.DataFrame) -> pd.DataFrame:
    stats = pd.DataFrame(
        {
            "total_signal": matrix.sum(axis=0),
            "detected_features": (matrix > 0).sum(axis=0),
            "mean_feature_value": matrix.mean(axis=0),
            "median_feature_value": matrix.median(axis=0),
            "std_feature_value": matrix.std(axis=0, ddof=1).fillna(0.0),
        }
    )
    return stats.sort_values("total_signal", ascending=False)


def compute_sample_correlation(matrix: pd.DataFrame, method: str = "pearson") -> pd.DataFrame:
    return matrix.corr(method=method)


def compute_sample_distances(matrix: pd.DataFrame) -> pd.DataFrame:
    sample_frame = matrix.transpose()
    if sample_frame.shape[0] < 2:
        return pd.DataFrame(index=sample_frame.index, columns=sample_frame.index, data=0.0)
    distances = squareform(pdist(sample_frame.values, metric="euclidean"))
    return pd.DataFrame(distances, index=sample_frame.index, columns=sample_frame.index)


def compute_pca(matrix: pd.DataFrame, n_components: int = 2) -> tuple[pd.DataFrame, list[float]] | None:
    sample_frame = matrix.transpose()
    if sample_frame.shape[0] < 2 or sample_frame.shape[1] < 2:
        return None

    variable_columns = sample_frame.loc[:, sample_frame.std(axis=0, ddof=0) > 0]
    if variable_columns.shape[1] < 2:
        return None

    scaler = StandardScaler()
    scaled = scaler.fit_transform(variable_columns)
    component_count = min(n_components, scaled.shape[0], scaled.shape[1])
    if component_count < 2:
        return None

    pca = PCA(n_components=component_count)
    coords = pca.fit_transform(scaled)
    coord_frame = pd.DataFrame(
        coords,
        index=sample_frame.index,
        columns=[f"PC{i + 1}" for i in range(component_count)],
    )
    explained = [float(value) for value in pca.explained_variance_ratio_]
    return coord_frame, explained


def compute_group_means(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_id_column: str,
    group_column: str,
) -> pd.DataFrame | None:
    if sample_id_column not in metadata.columns or group_column not in metadata.columns:
        return None

    groups = metadata[[sample_id_column, group_column]].dropna().copy()
    groups[sample_id_column] = groups[sample_id_column].astype(str).str.strip()
    groups[group_column] = groups[group_column].astype(str).str.strip()
    sample_frame = matrix.transpose().copy()
    sample_frame.index = sample_frame.index.astype(str)
    merged = sample_frame.merge(groups, left_index=True, right_on=sample_id_column, how="inner")
    if merged.empty:
        return None

    feature_columns = [column for column in merged.columns if column not in {sample_id_column, group_column}]
    group_means = merged.groupby(group_column)[feature_columns].mean().transpose()
    return group_means


def compute_group_tests(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_id_column: str,
    group_column: str,
) -> pd.DataFrame | None:
    if sample_id_column not in metadata.columns or group_column not in metadata.columns:
        return None

    groups = metadata[[sample_id_column, group_column]].dropna().copy()
    groups[sample_id_column] = groups[sample_id_column].astype(str).str.strip()
    groups[group_column] = groups[group_column].astype(str).str.strip()

    sample_frame = matrix.transpose().copy()
    sample_frame.index = sample_frame.index.astype(str)
    merged = sample_frame.merge(groups, left_index=True, right_on=sample_id_column, how="inner")
    if merged.empty:
        return None

    feature_columns = [column for column in merged.columns if column not in {sample_id_column, group_column}]
    results: list[dict[str, float | str | int]] = []

    for feature in feature_columns:
        grouped_values = []
        group_sizes = []
        for _, values in merged.groupby(group_column)[feature]:
            clean_values = values.dropna().astype(float).tolist()
            if clean_values:
                grouped_values.append(clean_values)
                group_sizes.append(len(clean_values))

        if len(grouped_values) < 2:
            continue

        enough_for_anova = all(size >= 2 for size in group_sizes)
        anova_p = math.nan
        if enough_for_anova:
            _, anova_p = f_oneway(*grouped_values)

        try:
            _, kruskal_p = kruskal(*grouped_values)
        except ValueError:
            kruskal_p = math.nan

        results.append(
            {
                "feature": str(feature),
                "n_groups": len(grouped_values),
                "group_sizes": ",".join(str(size) for size in group_sizes),
                "anova_pvalue": float(anova_p) if not math.isnan(anova_p) else math.nan,
                "kruskal_pvalue": float(kruskal_p) if not math.isnan(kruskal_p) else math.nan,
            }
        )

    if not results:
        return None

    result_frame = pd.DataFrame(results).set_index("feature")
    return result_frame.sort_values(["kruskal_pvalue", "anova_pvalue"], na_position="last")
