from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


sns.set_theme(style="whitegrid")


def _prepare_output(output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)


def plot_top_features_bar(
    feature_stats: pd.DataFrame,
    output_path: Path,
    title: str,
    top_n: int,
    metric: str = "mean",
) -> None:
    top_frame = feature_stats.head(top_n).iloc[::-1].reset_index()
    feature_label = top_frame.columns[0]
    _prepare_output(output_path)
    plt.figure(figsize=(10, max(6, top_n * 0.35)))
    sns.barplot(x=metric, y=feature_label, data=top_frame, hue=feature_label, dodge=False, palette="crest")
    legend = plt.gca().get_legend()
    if legend is not None:
        legend.remove()
    plt.title(title)
    plt.xlabel(metric.replace("_", " ").title())
    plt.ylabel("Pathway / Feature")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_heatmap(
    matrix: pd.DataFrame,
    feature_stats: pd.DataFrame,
    output_path: Path,
    title: str,
    top_n: int,
) -> None:
    top_features = feature_stats.head(top_n).index
    heatmap_frame = matrix.loc[top_features]
    if heatmap_frame.empty:
        return
    _prepare_output(output_path)
    plt.figure(figsize=(max(8, heatmap_frame.shape[1] * 1.2), max(8, top_n * 0.35)))
    sns.heatmap(heatmap_frame, cmap="mako", linewidths=0.2)
    plt.title(title)
    plt.xlabel("Sample")
    plt.ylabel("Pathway / Feature")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_sample_totals(sample_stats: pd.DataFrame, output_path: Path, title: str) -> None:
    frame = sample_stats.reset_index().rename(columns={"index": "sample_id"})
    _prepare_output(output_path)
    plt.figure(figsize=(10, max(5, len(frame) * 0.35)))
    sns.barplot(x="total_signal", y="sample_id", data=frame, palette="flare")
    plt.title(title)
    plt.xlabel("Total Signal")
    plt.ylabel("Sample")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_correlation_heatmap(correlation: pd.DataFrame, output_path: Path, title: str) -> None:
    if correlation.empty:
        return
    _prepare_output(output_path)
    plt.figure(figsize=(max(7, correlation.shape[1]), max(6, correlation.shape[0] * 0.8)))
    sns.heatmap(correlation, cmap="vlag", center=0, annot=False, linewidths=0.2)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_pca(
    coordinates: pd.DataFrame,
    explained_variance: list[float],
    output_path: Path,
    title: str,
    group_series: pd.Series | None = None,
) -> None:
    if coordinates.empty or coordinates.shape[1] < 2:
        return

    _prepare_output(output_path)
    plot_frame = coordinates.copy()
    plot_frame["sample_id"] = plot_frame.index
    if group_series is not None:
        aligned = group_series.reindex(plot_frame["sample_id"]).fillna("Unassigned")
        plot_frame["group"] = aligned.values

    plt.figure(figsize=(8, 6))
    if "group" in plot_frame.columns:
        sns.scatterplot(
            data=plot_frame,
            x="PC1",
            y="PC2",
            hue="group",
            s=90,
            palette="tab10",
        )
    else:
        sns.scatterplot(data=plot_frame, x="PC1", y="PC2", s=90, color="#2c7fb8")

    for _, row in plot_frame.iterrows():
        plt.text(row["PC1"], row["PC2"], str(row["sample_id"]), fontsize=9, ha="left", va="bottom")

    x_label = f"PC1 ({explained_variance[0] * 100:.1f}% var)"
    y_label = f"PC2 ({explained_variance[1] * 100:.1f}% var)"
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
