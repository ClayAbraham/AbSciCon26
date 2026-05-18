from __future__ import annotations

import argparse
import csv
import html
import math
import textwrap
from pathlib import Path


OVERVIEW_PATHWAYS = {
    "map01100",
    "map01110",
    "map01120",
    "map01200",
    "map01210",
    "map01212",
    "map01220",
    "map01230",
    "map01232",
    "map01240",
    "map01250",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render 1A pathway charts as SVG files.")
    parser.add_argument("--input-dir", required=True, help="Directory containing extracted TSV files.")
    parser.add_argument("--output-dir", required=True, help="Directory for charts.")
    parser.add_argument("--sample-id", required=True, help="Sample identifier used in chart titles and filenames.")
    return parser.parse_args()


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def nice_ticks(max_value: float, count: int = 6) -> list[float]:
    if max_value <= 0:
        return [0, 1]
    raw_step = max_value / max(count, 1)
    magnitude = 10 ** math.floor(math.log10(raw_step))
    step = magnitude
    for candidate in (1, 2, 5, 10):
        test_step = candidate * magnitude
        if test_step >= raw_step:
            step = test_step
            break
    tick_max = math.ceil(max_value / step) * step
    ticks = [0.0]
    current = step
    while current <= tick_max + (step * 0.001):
        ticks.append(float(current))
        current += step
    return ticks


def wrap_label(label: str, width: int = 34, max_lines: int = 2) -> list[str]:
    wrapped = textwrap.wrap(label, width=width) or [label]
    if len(wrapped) > max_lines:
        kept = wrapped[: max_lines]
        kept[-1] = kept[-1][: max(0, width - 3)].rstrip() + "..."
        return kept
    return wrapped


def svg_header(width: int, height: int) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
    )


def text_svg(x: float, y: float, text: str, **attrs: str | float) -> str:
    attr_text = " ".join(f'{key}="{html.escape(str(value), quote=True)}"' for key, value in attrs.items())
    return f'<text x="{x:.2f}" y="{y:.2f}" {attr_text}>{html.escape(text)}</text>'


def line_svg(x1: float, y1: float, x2: float, y2: float, **attrs: str | float) -> str:
    attr_text = " ".join(f'{key}="{html.escape(str(value), quote=True)}"' for key, value in attrs.items())
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" {attr_text} />'


def rect_svg(x: float, y: float, width: float, height: float, **attrs: str | float) -> str:
    attr_text = " ".join(f'{key}="{html.escape(str(value), quote=True)}"' for key, value in attrs.items())
    return (
        f'<rect x="{x:.2f}" y="{y:.2f}" width="{width:.2f}" height="{height:.2f}" '
        f'{attr_text} />'
    )


def circle_svg(cx: float, cy: float, r: float, **attrs: str | float) -> str:
    attr_text = " ".join(f'{key}="{html.escape(str(value), quote=True)}"' for key, value in attrs.items())
    return f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r:.2f}" {attr_text} />'


def render_horizontal_bar_chart(
    records: list[dict[str, object]],
    label_key: str,
    value_key: str,
    output_path: Path,
    title: str,
    subtitle: str,
    x_axis_title: str,
    number_format: str,
    axis_max: float | None = None,
) -> None:
    bar_height = 30
    bar_gap = 14
    count = len(records)
    width = 1700
    height = 180 + count * (bar_height + bar_gap) + 100
    left_margin = 520
    right_margin = 140
    top_margin = 140
    bottom_margin = 90
    plot_width = width - left_margin - right_margin
    plot_height = height - top_margin - bottom_margin

    max_value = max(float(record[value_key]) for record in records) if records else 1.0
    axis_limit = axis_max if axis_max is not None else max_value * 1.1
    ticks = nice_ticks(axis_limit, count=6)
    axis_limit = max(ticks)

    parts = [svg_header(width, height)]
    parts.append(rect_svg(0, 0, width, height, fill="white"))
    parts.append(text_svg(width / 2, 48, title, **{"text-anchor": "middle", "font-size": 22, "font-weight": "700", "font-family": "Arial"}))
    parts.append(
        text_svg(
            width / 2,
            78,
            subtitle,
            **{"text-anchor": "middle", "font-size": 13, "fill": "#555", "font-family": "Arial"},
        )
    )

    x0 = left_margin
    y0 = top_margin
    x1 = left_margin + plot_width
    y1 = top_margin + plot_height

    for tick in ticks:
        x = x0 + (tick / axis_limit) * plot_width
        parts.append(line_svg(x, y0, x, y1, stroke="#e3e3e3", **{"stroke-width": 1}))
        parts.append(
            text_svg(
                x,
                y1 + 28,
                format(tick, number_format),
                **{"text-anchor": "middle", "font-size": 12, "fill": "#333", "font-family": "Arial"},
            )
        )

    parts.append(line_svg(x0, y0, x0, y1, stroke="#222", **{"stroke-width": 1.4}))
    parts.append(line_svg(x0, y1, x1, y1, stroke="#222", **{"stroke-width": 1.4}))
    parts.append(
        text_svg(
            (x0 + x1) / 2,
            height - 24,
            x_axis_title,
            **{"text-anchor": "middle", "font-size": 15, "font-weight": "700", "font-family": "Arial"},
        )
    )

    for index, record in enumerate(records):
        value = float(record[value_key])
        label = str(record[label_key])
        bar_y = y0 + index * (bar_height + bar_gap)
        bar_w = 0 if axis_limit == 0 else (value / axis_limit) * plot_width

        parts.append(
            rect_svg(
                x0,
                bar_y,
                bar_w,
                bar_height,
                fill="#2c6aa0",
                stroke="#1f4e77",
                **{"stroke-width": 1},
            )
        )

        wrapped = wrap_label(label, width=34, max_lines=2)
        if len(wrapped) == 1:
            label_start_y = bar_y + bar_height / 2 + 5
        else:
            label_start_y = bar_y + bar_height / 2 - 5
        for line_index, line in enumerate(wrapped):
            parts.append(
                text_svg(
                    x0 - 12,
                    label_start_y + line_index * 15,
                    line,
                    **{"text-anchor": "end", "font-size": 13, "fill": "#111", "font-family": "Arial"},
                )
            )

        value_x = min(x0 + bar_w + 8, x1 - 5)
        parts.append(
            text_svg(
                value_x,
                bar_y + bar_height / 2 + 5,
                format(value, number_format),
                **{"text-anchor": "start", "font-size": 12, "fill": "#111", "font-family": "Arial"},
            )
        )

    parts.append("</svg>")
    write_text(output_path, "\n".join(parts))


def render_scatter_chart(records: list[dict[str, object]], output_path: Path, sample_id: str) -> None:
    width = 1400
    height = 900
    left_margin = 110
    right_margin = 70
    top_margin = 120
    bottom_margin = 90
    plot_width = width - left_margin - right_margin
    plot_height = height - top_margin - bottom_margin

    max_x = max(float(record["observed_kos"]) for record in records) * 1.05
    max_total = max(float(record["total_kos"]) for record in records)
    x_ticks = nice_ticks(max_x, count=6)
    y_ticks = list(range(0, 101, 10))

    parts = [svg_header(width, height)]
    parts.append(rect_svg(0, 0, width, height, fill="white"))
    parts.append(text_svg(width / 2, 48, f"Sample {sample_id} Pathway Completeness vs Observed KOs", **{"text-anchor": "middle", "font-size": 22, "font-weight": "700", "font-family": "Arial"}))
    parts.append(
        text_svg(
            width / 2,
            78,
            "Point size scales with total KOs; high-completeness pathways with larger denominators are more robust",
            **{"text-anchor": "middle", "font-size": 13, "fill": "#555", "font-family": "Arial"},
        )
    )

    x0 = left_margin
    y0 = top_margin
    x1 = left_margin + plot_width
    y1 = top_margin + plot_height

    for tick in x_ticks:
        x = x0 + (tick / max(x_ticks)) * plot_width
        parts.append(line_svg(x, y0, x, y1, stroke="#e5e5e5", **{"stroke-width": 1}))
        parts.append(
            text_svg(
                x,
                y1 + 28,
                f"{tick:.0f}",
                **{"text-anchor": "middle", "font-size": 12, "fill": "#333", "font-family": "Arial"},
            )
        )

    for tick in y_ticks:
        y = y1 - (tick / 100.0) * plot_height
        parts.append(line_svg(x0, y, x1, y, stroke="#e5e5e5", **{"stroke-width": 1}))
        parts.append(
            text_svg(
                x0 - 10,
                y + 4,
                str(tick),
                **{"text-anchor": "end", "font-size": 12, "fill": "#333", "font-family": "Arial"},
            )
        )

    parts.append(line_svg(x0, y0, x0, y1, stroke="#222", **{"stroke-width": 1.4}))
    parts.append(line_svg(x0, y1, x1, y1, stroke="#222", **{"stroke-width": 1.4}))
    parts.append(text_svg((x0 + x1) / 2, height - 24, "Observed KOs", **{"text-anchor": "middle", "font-size": 15, "font-weight": "700", "font-family": "Arial"}))
    parts.append(
        f'<text x="36" y="{(y0 + y1) / 2:.2f}" transform="rotate(-90 36 {(y0 + y1) / 2:.2f})" '
        f'text-anchor="middle" font-size="15" font-weight="700" font-family="Arial">Completeness (%)</text>'
    )

    for record in records:
        observed = float(record["observed_kos"])
        completeness = float(record["completeness_pct"])
        total = float(record["total_kos"])
        x = x0 + (observed / max(x_ticks)) * plot_width
        y = y1 - (completeness / 100.0) * plot_height
        radius = 3.5 + 8.5 * math.sqrt(total / max_total)
        parts.append(circle_svg(x, y, radius, fill="#d25a28", opacity="0.55", stroke="#9d3e17", **{"stroke-width": 1}))

    parts.append("</svg>")
    write_text(output_path, "\n".join(parts))


def render_html_index(output_dir: Path, sample_id: str) -> None:
    chart_files = [
        (f"{sample_id}_top_pathway_abundance.svg", "Top specific pathways by KO count"),
        (f"{sample_id}_top_pathway_completeness.svg", "Highest pathway completeness with denominator filter"),
        (f"{sample_id}_top_observed_kos.svg", "Pathways with highest observed KOs"),
        (f"{sample_id}_top_ko_counts.svg", "Top KO assignments"),
        (f"{sample_id}_completeness_vs_observed_kos.svg", "Completeness versus observed KOs"),
    ]
    blocks = []
    for filename, caption in chart_files:
        blocks.append(
            f'<section><h2>{html.escape(caption)}</h2>'
            f'<img src="{html.escape(filename)}" alt="{html.escape(caption)}" style="max-width: 100%; border: 1px solid #ddd;" /></section>'
        )
    html_doc = (
        f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>Sample {html.escape(sample_id)} Charts</title>"
        "<style>body{font-family:Arial,sans-serif;margin:32px;line-height:1.4}section{margin-bottom:40px}</style>"
        f"</head><body><h1>Sample {html.escape(sample_id)} Metabolic Profiling Charts</h1>"
        + "".join(blocks)
        + "</body></html>"
    )
    write_text(output_dir / f"{sample_id}_chart_index.html", html_doc)


def main() -> int:
    args = parse_args()
    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    sample_id = args.sample_id
    output_dir.mkdir(parents=True, exist_ok=True)

    annotated = load_tsv(input_dir / "annotated_pathways.tsv")
    completeness = load_tsv(input_dir / "pathway_completeness.tsv")
    ko_counts = load_tsv(input_dir / "ko_counts.tsv")

    for row in annotated:
        row["ko_count"] = float(row["ko_count"])
    for row in completeness:
        row["observed_kos"] = float(row["observed_kos"])
        row["total_kos"] = float(row["total_kos"])
        row["completeness_pct"] = float(row["completeness_pct"])
    for row in ko_counts:
        row["count"] = float(row["count"])

    specific_annotated = [row for row in annotated if row["pathway_id"] not in OVERVIEW_PATHWAYS]
    specific_completeness = [row for row in completeness if row["pathway_id"] not in OVERVIEW_PATHWAYS]

    top_abundance = sorted(specific_annotated, key=lambda row: float(row["ko_count"]), reverse=True)[:15]
    top_completeness = sorted(
        [row for row in specific_completeness if float(row["total_kos"]) >= 20],
        key=lambda row: (float(row["completeness_pct"]), float(row["observed_kos"])),
        reverse=True,
    )[:15]
    top_observed = sorted(
        specific_completeness,
        key=lambda row: (float(row["observed_kos"]), float(row["completeness_pct"])),
        reverse=True,
    )[:15]
    top_ko = sorted(ko_counts, key=lambda row: float(row["count"]), reverse=True)[:20]

    render_horizontal_bar_chart(
        records=top_abundance,
        label_key="pathway_name",
        value_key="ko_count",
        output_path=output_dir / f"{sample_id}_top_pathway_abundance.svg",
        title=f"Sample {sample_id} Top Specific Pathways by KO Count",
        subtitle="KEGG overview maps excluded so broad summary bins do not dominate the ranking",
        x_axis_title="KO Count",
        number_format=",.0f",
    )
    render_horizontal_bar_chart(
        records=top_completeness,
        label_key="pathway_name",
        value_key="completeness_pct",
        output_path=output_dir / f"{sample_id}_top_pathway_completeness.svg",
        title=f"Sample {sample_id} Highest Pathway Completeness",
        subtitle="Pathways filtered to total KOs >= 20 so completeness percentages rest on a meaningful denominator",
        x_axis_title="Completeness (%)",
        number_format=",.1f",
        axis_max=100.0,
    )
    render_horizontal_bar_chart(
        records=top_observed,
        label_key="pathway_name",
        value_key="observed_kos",
        output_path=output_dir / f"{sample_id}_top_observed_kos.svg",
        title=f"Sample {sample_id} Pathways with Highest Observed KOs",
        subtitle="KEGG overview maps excluded to focus on specific pathways instead of summary categories",
        x_axis_title="Observed KOs",
        number_format=",.0f",
    )
    render_horizontal_bar_chart(
        records=top_ko,
        label_key="ko_id",
        value_key="count",
        output_path=output_dir / f"{sample_id}_top_ko_counts.svg",
        title=f"Sample {sample_id} Top KO Assignments",
        subtitle="Numeric ranking from the ko_counts worksheet using true KO hit counts",
        x_axis_title="KO Count",
        number_format=",.0f",
    )
    render_scatter_chart(completeness, output_dir / f"{sample_id}_completeness_vs_observed_kos.svg", sample_id)
    render_html_index(output_dir, sample_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
