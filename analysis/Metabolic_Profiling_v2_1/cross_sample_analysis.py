from __future__ import annotations

import argparse
import csv
import html
import math
import re
import statistics
from pathlib import Path
from urllib import error, request


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
    parser = argparse.ArgumentParser(description="Cross-sample metabolic profiling analysis.")
    parser.add_argument("--input-root", required=True, help="Directory containing sample *_charts folders.")
    parser.add_argument("--output-dir", required=True, help="Directory for combined outputs.")
    parser.add_argument(
        "--min-pathway-abundance-count",
        type=float,
        default=10.0,
        help="Per-sample pathway abundance values below this threshold are treated as noise and set to zero.",
    )
    parser.add_argument(
        "--min-ko-count",
        type=float,
        default=20.0,
        help="Per-sample KO count values below this threshold are treated as noise and set to zero.",
    )
    parser.add_argument(
        "--min-observed-kos",
        type=float,
        default=3.0,
        help="Pathway completeness values with fewer observed KOs than this are treated as low-signal.",
    )
    parser.add_argument(
        "--min-total-kos",
        type=float,
        default=20.0,
        help="Pathway completeness values with total KOs below this threshold are treated as low-signal.",
    )
    return parser.parse_args()


def natural_sample_key(value: str) -> tuple[int, str]:
    match = re.match(r"(\d+)(.*)", value)
    if match:
        return int(match.group(1)), match.group(2)
    return 10**9, value


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def std(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    return statistics.stdev(values)


def pearson_correlation(values_a: list[float], values_b: list[float]) -> float:
    if len(values_a) != len(values_b) or not values_a:
        return 0.0
    mean_a = mean(values_a)
    mean_b = mean(values_b)
    centered_a = [value - mean_a for value in values_a]
    centered_b = [value - mean_b for value in values_b]
    denom_a = math.sqrt(sum(value * value for value in centered_a))
    denom_b = math.sqrt(sum(value * value for value in centered_b))
    if denom_a == 0.0 and denom_b == 0.0:
        return 1.0
    if denom_a == 0.0 or denom_b == 0.0:
        return 0.0
    numerator = sum(a * b for a, b in zip(centered_a, centered_b))
    return numerator / (denom_a * denom_b)


def bray_curtis(values_a: list[float], values_b: list[float]) -> float:
    numerator = sum(abs(a - b) for a, b in zip(values_a, values_b))
    denominator = sum(a + b for a, b in zip(values_a, values_b))
    if denominator == 0.0:
        return 0.0
    return numerator / denominator


def jaccard_presence(values_a: list[float], values_b: list[float]) -> float:
    present_a = [value > 0 for value in values_a]
    present_b = [value > 0 for value in values_b]
    intersection = sum(1 for a, b in zip(present_a, present_b) if a and b)
    union = sum(1 for a, b in zip(present_a, present_b) if a or b)
    if union == 0:
        return 1.0
    return intersection / union


def mean_absolute_difference(values_a: list[float], values_b: list[float]) -> float:
    return mean([abs(a - b) for a, b in zip(values_a, values_b)])


def log10_1p(value: float) -> float:
    return math.log10(1.0 + value)


def normalize_profile(values: list[float]) -> list[float]:
    total = sum(values)
    if total <= 0.0:
        return [0.0 for _ in values]
    return [value / total for value in values]


def kl_divergence(prob_a: list[float], prob_b: list[float]) -> float:
    total = 0.0
    for value_a, value_b in zip(prob_a, prob_b):
        if value_a <= 0.0:
            continue
        if value_b <= 0.0:
            continue
        total += value_a * math.log2(value_a / value_b)
    return total


def jensen_shannon_divergence(values_a: list[float], values_b: list[float]) -> float:
    prob_a = normalize_profile(values_a)
    prob_b = normalize_profile(values_b)
    if sum(prob_a) == 0.0 and sum(prob_b) == 0.0:
        return 0.0
    midpoint = [(value_a + value_b) / 2.0 for value_a, value_b in zip(prob_a, prob_b)]
    return 0.5 * kl_divergence(prob_a, midpoint) + 0.5 * kl_divergence(prob_b, midpoint)


def zscore_row(values: list[float]) -> list[float]:
    row_mean = mean(values)
    row_std = std(values)
    if row_std == 0.0:
        return [0.0 for _ in values]
    return [(value - row_mean) / row_std for value in values]


def nice_ticks(max_value: float, count: int = 6) -> list[float]:
    if max_value <= 0:
        return [0.0, 1.0]
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


def wrap_label(label: str, width: int = 34, max_lines: int = 2) -> list[str]:
    words = label.split()
    if not words:
        return [label]
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = word if not current else current + " " + word
        if len(candidate) <= width:
            current = candidate
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        lines[-1] = lines[-1][: max(0, width - 3)].rstrip() + "..."
    return lines


def interpolate_stops(ratio: float, stops: list[tuple[float, tuple[int, int, int]]]) -> str:
    ratio = max(0.0, min(1.0, ratio))

    for index in range(len(stops) - 1):
        left_pos, left_rgb = stops[index]
        right_pos, right_rgb = stops[index + 1]
        if ratio <= right_pos:
            local = 0.0 if right_pos == left_pos else (ratio - left_pos) / (right_pos - left_pos)
            r = int(left_rgb[0] + local * (right_rgb[0] - left_rgb[0]))
            g = int(left_rgb[1] + local * (right_rgb[1] - left_rgb[1]))
            b = int(left_rgb[2] + local * (right_rgb[2] - left_rgb[2]))
            return f"rgb({r},{g},{b})"

    r, g, b = stops[-1][1]
    return f"rgb({r},{g},{b})"


def diverging_color(value: float, minimum: float = -2.5, maximum: float = 2.5) -> str:
    if maximum <= minimum:
        return "rgb(245,247,240)"

    clamped = max(minimum, min(maximum, value))
    if minimum < 0 < maximum:
        if clamped <= 0:
            ratio = 0.0 if minimum == 0 else (clamped - minimum) / (0 - minimum)
            return interpolate_stops(
                ratio,
                [
                    (0.00, (19, 35, 84)),
                    (0.55, (0, 119, 182)),
                    (1.00, (245, 247, 240)),
                ],
            )
        ratio = 0.0 if maximum == 0 else clamped / maximum
        return interpolate_stops(
            ratio,
            [
                (0.00, (245, 247, 240)),
                (0.45, (255, 214, 10)),
                (1.00, (249, 65, 68)),
            ],
        )

    ratio = (clamped - minimum) / (maximum - minimum)
    return interpolate_stops(
        ratio,
        [
            (0.00, (19, 35, 84)),
            (0.55, (0, 119, 182)),
            (1.00, (245, 247, 240)),
        ],
    )


def sequential_color(value: float, minimum: float, maximum: float) -> str:
    if maximum <= minimum:
        ratio = 0.0
    else:
        ratio = (value - minimum) / (maximum - minimum)
    ratio = max(0.0, min(1.0, ratio))
    r = int(247 - ratio * (247 - 40))
    g = int(252 - ratio * (252 - 120))
    b = int(240 - ratio * (240 - 142))
    return f"rgb({r},{g},{b})"


def vibrant_completeness_color(value: float, minimum: float, maximum: float) -> str:
    if maximum <= minimum:
        ratio = 0.0
    else:
        ratio = (value - minimum) / (maximum - minimum)
    return interpolate_stops(
        ratio,
        [
            (0.00, (19, 35, 84)),
            (0.25, (0, 119, 182)),
            (0.50, (0, 191, 99)),
            (0.75, (255, 214, 10)),
            (1.00, (249, 65, 68)),
        ],
    )


def heatmap_color(value: float, scheme: str, minimum: float, maximum: float) -> str:
    if scheme in {"diverging", "vibrant_diverging"}:
        return diverging_color(value, minimum, maximum)
    if scheme in {"vibrant_completeness", "vibrant_sequential"}:
        return vibrant_completeness_color(value, minimum, maximum)
    return sequential_color(value, minimum, maximum)


def parse_kegg_ko_entries(payload: str) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    for block in payload.split("///"):
        lines = [line.rstrip() for line in block.splitlines() if line.strip()]
        if not lines:
            continue

        current_field = ""
        fields: dict[str, list[str]] = {}
        for line in lines:
            field = line[:12].strip()
            value = line[12:].strip()
            if field:
                current_field = field
                fields.setdefault(field, []).append(value)
            elif current_field:
                fields.setdefault(current_field, []).append(value)

        entry_value = fields.get("ENTRY", [""])[0]
        ko_id = entry_value.split()[0]
        if not ko_id:
            continue

        symbol = " ".join(fields.get("SYMBOL", [])).strip()
        name = " ".join(fields.get("NAME", [])).strip()
        if symbol and name:
            label = f"{symbol}; {name}"
        elif name:
            label = name
        elif symbol:
            label = symbol
        else:
            label = ko_id

        records[ko_id] = {
            "ko_id": ko_id,
            "symbol": symbol,
            "name": name,
            "label": label,
        }
    return records


def load_kegg_ko_cache(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    rows = load_tsv(path)
    cache: dict[str, dict[str, str]] = {}
    for row in rows:
        ko_id = row.get("ko_id", "").strip()
        if ko_id:
            cache[ko_id] = {
                "ko_id": ko_id,
                "symbol": row.get("symbol", "").strip(),
                "name": row.get("name", "").strip(),
                "label": row.get("label", ko_id).strip() or ko_id,
            }
    return cache


def write_kegg_ko_cache(path: Path, cache: dict[str, dict[str, str]]) -> None:
    rows = [cache[ko_id] for ko_id in sorted(cache)]
    write_tsv(path, ["ko_id", "symbol", "name", "label"], rows)


def fetch_kegg_ko_labels(
    ko_ids: list[str],
    cache_path: Path,
    batch_size: int = 10,
) -> dict[str, dict[str, str]]:
    cache = load_kegg_ko_cache(cache_path)
    missing = [ko_id for ko_id in ko_ids if ko_id not in cache]

    for start in range(0, len(missing), batch_size):
        batch = missing[start : start + batch_size]
        if not batch:
            continue
        url = "https://rest.kegg.jp/get/" + "+".join(f"ko:{ko_id}" for ko_id in batch)
        with request.urlopen(url, timeout=30) as response:
            payload = response.read().decode("utf-8", errors="replace")
        cache.update(parse_kegg_ko_entries(payload))

    write_kegg_ko_cache(cache_path, cache)
    return {ko_id: cache[ko_id] for ko_id in ko_ids if ko_id in cache}


def pairwise_distance(vectors: list[list[float]]) -> dict[tuple[int, int], float]:
    distances: dict[tuple[int, int], float] = {}
    for index_a in range(len(vectors)):
        for index_b in range(index_a + 1, len(vectors)):
            corr = pearson_correlation(vectors[index_a], vectors[index_b])
            distances[(index_a, index_b)] = 1.0 - corr
    return distances


def average_cluster_distance(
    members_a: list[int],
    members_b: list[int],
    distances: dict[tuple[int, int], float],
) -> float:
    values: list[float] = []
    for member_a in members_a:
        for member_b in members_b:
            if member_a == member_b:
                continue
            key = (member_a, member_b) if member_a < member_b else (member_b, member_a)
            values.append(distances.get(key, 0.0))
    return mean(values) if values else 0.0


def cluster_order(vectors: list[list[float]]) -> list[int]:
    if len(vectors) <= 1:
        return list(range(len(vectors)))

    distances = pairwise_distance(vectors)
    clusters = [
        {
            "members": [index],
            "order": [index],
            "min_member": index,
        }
        for index in range(len(vectors))
    ]

    while len(clusters) > 1:
        best_pair: tuple[int, int] | None = None
        best_score: tuple[float, int, int] | None = None

        for left_index in range(len(clusters)):
            for right_index in range(left_index + 1, len(clusters)):
                left = clusters[left_index]
                right = clusters[right_index]
                distance = average_cluster_distance(left["members"], right["members"], distances)
                score = (round(distance, 12), min(left["min_member"], right["min_member"]), max(left["min_member"], right["min_member"]))
                if best_score is None or score < best_score:
                    best_score = score
                    best_pair = (left_index, right_index)

        assert best_pair is not None
        left_index, right_index = best_pair
        left = clusters[left_index]
        right = clusters[right_index]
        merged = {
            "members": left["members"] + right["members"],
            "order": left["order"] + right["order"],
            "min_member": min(left["min_member"], right["min_member"]),
        }
        new_clusters = []
        for index, cluster in enumerate(clusters):
            if index not in {left_index, right_index}:
                new_clusters.append(cluster)
        new_clusters.append(merged)
        clusters = new_clusters

    return clusters[0]["order"]


def transpose(matrix: list[list[float]]) -> list[list[float]]:
    if not matrix:
        return []
    return [list(column) for column in zip(*matrix)]


def format_label(feature_id: str, feature_name: str | None) -> str:
    if feature_name:
        return feature_name
    return feature_id


def render_clustered_heatmap(
    output_path: Path,
    title: str,
    subtitle: str,
    sample_ids: list[str],
    feature_ids: list[str],
    feature_labels: dict[str, str],
    matrix: list[list[float]],
    scheme: str,
    legend_title: str,
    cluster_columns: bool = False,
) -> tuple[list[str], list[str]]:
    if not matrix:
        raise ValueError(f"No data available for heatmap {title}")

    sample_order_indices = cluster_order(transpose(matrix)) if cluster_columns else list(range(len(sample_ids)))
    row_order_indices = cluster_order(matrix)

    ordered_samples = [sample_ids[index] for index in sample_order_indices]
    ordered_features = [feature_ids[index] for index in row_order_indices]
    ordered_matrix = [
        [matrix[row_index][column_index] for column_index in sample_order_indices]
        for row_index in row_order_indices
    ]

    cell_width = 56
    cell_height = 18
    width = 260 + len(ordered_samples) * cell_width + 460
    height = 190 + len(ordered_features) * cell_height + 120
    left_margin = 420
    top_margin = 150
    legend_x = left_margin + len(ordered_samples) * cell_width + 60
    legend_y = top_margin

    if scheme in {"diverging", "vibrant_diverging"}:
        legend_min = -2.5
        legend_max = 2.5
        legend_ticks = [-2, -1, 0, 1, 2]
    else:
        legend_min = 0.0
        legend_max = 100.0
        legend_ticks = [0, 20, 40, 60, 80, 100]

    parts = [svg_header(width, height)]
    parts.append(rect_svg(0, 0, width, height, fill="white"))
    parts.append(text_svg(width / 2, 42, title, **{"text-anchor": "middle", "font-size": 22, "font-weight": "700", "font-family": "Arial"}))
    parts.append(text_svg(width / 2, 70, subtitle, **{"text-anchor": "middle", "font-size": 13, "fill": "#555", "font-family": "Arial"}))

    for sample_index, sample_id in enumerate(ordered_samples):
        x = left_margin + sample_index * cell_width + cell_width / 2
        y = top_margin - 18
        parts.append(
            f'<text x="{x:.2f}" y="{y:.2f}" transform="rotate(-45 {x:.2f} {y:.2f})" '
            f'text-anchor="end" font-size="12" font-family="Arial">{html.escape(sample_id)}</text>'
        )

    for row_index, feature_id in enumerate(ordered_features):
        label = feature_labels[feature_id]
        wrapped = wrap_label(label, width=34, max_lines=2)
        center_y = top_margin + row_index * cell_height + cell_height / 2 + 4
        label_y = center_y if len(wrapped) == 1 else center_y - 6
        for line_index, line in enumerate(wrapped):
            parts.append(
                text_svg(
                    left_margin - 12,
                    label_y + line_index * 13,
                    line,
                    **{"text-anchor": "end", "font-size": 12, "font-family": "Arial", "fill": "#111"},
                )
            )

    for row_index, row in enumerate(ordered_matrix):
        for column_index, value in enumerate(row):
            x = left_margin + column_index * cell_width
            y = top_margin + row_index * cell_height
            parts.append(
                rect_svg(
                    x,
                    y,
                    cell_width,
                    cell_height,
                    fill=heatmap_color(value, scheme, legend_min, legend_max),
                    stroke="#ffffff",
                    **{"stroke-width": 0.6},
                )
            )

    legend_height = 240
    legend_width = 24
    steps = 100
    for step in range(steps):
        if scheme == "diverging":
            value = legend_max - (step / (steps - 1)) * (legend_max - legend_min)
            color = heatmap_color(value, scheme, legend_min, legend_max)
        else:
            value = legend_max - (step / (steps - 1)) * (legend_max - legend_min)
            color = heatmap_color(value, scheme, legend_min, legend_max)
        y = legend_y + step * (legend_height / steps)
        parts.append(rect_svg(legend_x, y, legend_width, legend_height / steps + 1, fill=color, stroke=color))

    parts.append(text_svg(legend_x + legend_width / 2, legend_y - 12, legend_title, **{"text-anchor": "middle", "font-size": 12, "font-weight": "700", "font-family": "Arial"}))
    for tick in legend_ticks:
        if scheme == "diverging":
            ratio = (legend_max - tick) / (legend_max - legend_min)
        else:
            ratio = (legend_max - tick) / (legend_max - legend_min)
        y = legend_y + ratio * legend_height
        parts.append(line_svg(legend_x + legend_width, y, legend_x + legend_width + 8, y, stroke="#333", **{"stroke-width": 1}))
        parts.append(text_svg(legend_x + legend_width + 14, y + 4, str(tick), **{"text-anchor": "start", "font-size": 11, "font-family": "Arial", "fill": "#333"}))

    parts.append("</svg>")
    write_text(output_path, "\n".join(parts))
    return ordered_features, ordered_samples


def render_square_matrix_heatmap(
    output_path: Path,
    title: str,
    subtitle: str,
    labels: list[str],
    matrix: list[list[float]],
    value_min: float,
    value_max: float,
    legend_title: str,
    scheme: str = "vibrant_diverging",
) -> None:
    cell_size = 54
    width = 220 + len(labels) * cell_size + 220
    height = 180 + len(labels) * cell_size + 120
    left_margin = 170
    top_margin = 130
    legend_x = left_margin + len(labels) * cell_size + 45
    legend_y = top_margin
    legend_height = 240
    legend_width = 24

    parts = [svg_header(width, height)]
    parts.append(rect_svg(0, 0, width, height, fill="white"))
    parts.append(text_svg(width / 2, 42, title, **{"text-anchor": "middle", "font-size": 22, "font-weight": "700", "font-family": "Arial"}))
    parts.append(text_svg(width / 2, 70, subtitle, **{"text-anchor": "middle", "font-size": 13, "fill": "#555", "font-family": "Arial"}))

    for index, label in enumerate(labels):
        x = left_margin + index * cell_size + cell_size / 2
        y = top_margin - 16
        parts.append(
            f'<text x="{x:.2f}" y="{y:.2f}" transform="rotate(-45 {x:.2f} {y:.2f})" '
            f'text-anchor="end" font-size="12" font-family="Arial">{html.escape(label)}</text>'
        )
        parts.append(
            text_svg(
                left_margin - 10,
                top_margin + index * cell_size + cell_size / 2 + 4,
                label,
                **{"text-anchor": "end", "font-size": 12, "font-family": "Arial"},
            )
        )

    for row_index, row in enumerate(matrix):
        for col_index, value in enumerate(row):
            x = left_margin + col_index * cell_size
            y = top_margin + row_index * cell_size
            color = heatmap_color(value, scheme, value_min, value_max)
            parts.append(rect_svg(x, y, cell_size, cell_size, fill=color, stroke="#ffffff", **{"stroke-width": 1}))
            parts.append(
                text_svg(
                    x + cell_size / 2,
                    y + cell_size / 2 + 4,
                    f"{value:.2f}",
                    **{"text-anchor": "middle", "font-size": 11, "font-family": "Arial", "fill": "#111"},
                )
            )

    steps = 100
    for step in range(steps):
        value = value_max - (step / (steps - 1)) * (value_max - value_min)
        color = heatmap_color(value, scheme, value_min, value_max)
        y = legend_y + step * (legend_height / steps)
        parts.append(rect_svg(legend_x, y, legend_width, legend_height / steps + 1, fill=color, stroke=color))

    parts.append(text_svg(legend_x + legend_width / 2, legend_y - 12, legend_title, **{"text-anchor": "middle", "font-size": 12, "font-weight": "700", "font-family": "Arial"}))
    for tick in nice_ticks(value_max - value_min, count=5):
        actual = value_min + tick
        if actual > value_max + 1e-9:
            continue
        ratio = (value_max - actual) / (value_max - value_min) if value_max != value_min else 0.0
        y = legend_y + ratio * legend_height
        parts.append(line_svg(legend_x + legend_width, y, legend_x + legend_width + 8, y, stroke="#333", **{"stroke-width": 1}))
        parts.append(text_svg(legend_x + legend_width + 14, y + 4, f"{actual:.2f}", **{"text-anchor": "start", "font-size": 11, "font-family": "Arial"}))

    parts.append("</svg>")
    write_text(output_path, "\n".join(parts))


def build_feature_maps(sample_dirs: list[tuple[str, Path]]) -> tuple[
    dict[str, dict[str, float]],
    dict[str, str],
    dict[str, dict[str, float]],
    dict[str, str],
    dict[str, float],
    dict[str, dict[str, float]],
]:
    pathway_abundance: dict[str, dict[str, float]] = {}
    pathway_names: dict[str, str] = {}
    completeness_pct: dict[str, dict[str, float]] = {}
    completeness_names: dict[str, str] = {}
    total_kos_map: dict[str, float] = {}
    ko_counts: dict[str, dict[str, float]] = {}

    for sample_id, sample_dir in sample_dirs:
        abundance_rows = load_tsv(sample_dir / "annotated_pathways.tsv")
        completeness_rows = load_tsv(sample_dir / "pathway_completeness.tsv")
        ko_rows = load_tsv(sample_dir / "ko_counts.tsv")

        for row in abundance_rows:
            pathway_id = row["pathway_id"].strip()
            pathway_name = row["pathway_name"].strip()
            value = float(row["ko_count"])
            pathway_abundance.setdefault(pathway_id, {})[sample_id] = value
            if pathway_name:
                pathway_names[pathway_id] = pathway_name

        for row in completeness_rows:
            pathway_id = row["pathway_id"].strip()
            pathway_name = row["pathway_name"].strip()
            completeness_pct.setdefault(pathway_id, {})[sample_id] = float(row["completeness_pct"])
            total_kos_map[pathway_id] = float(row["total_kos"])
            if pathway_name:
                completeness_names[pathway_id] = pathway_name

        for row in ko_rows:
            ko_id = row["ko_id"].strip()
            ko_counts.setdefault(ko_id, {})[sample_id] = float(row["count"])

    return pathway_abundance, pathway_names, completeness_pct, completeness_names, total_kos_map, ko_counts


def apply_thresholds(
    sample_ids: list[str],
    pathway_abundance_map: dict[str, dict[str, float]],
    completeness_map: dict[str, dict[str, float]],
    total_kos_map: dict[str, float],
    ko_count_map: dict[str, dict[str, float]],
    min_pathway_abundance_count: float,
    min_ko_count: float,
    min_observed_kos: float,
    min_total_kos: float,
) -> tuple[
    dict[str, dict[str, float]],
    dict[str, dict[str, float]],
    dict[str, dict[str, float]],
    list[dict[str, object]],
]:
    threshold_rows: list[dict[str, object]] = []

    filtered_pathway_abundance: dict[str, dict[str, float]] = {}
    pathway_removed = 0
    pathway_total = 0
    for feature_id, sample_map in pathway_abundance_map.items():
        filtered_samples: dict[str, float] = {}
        for sample_id in sample_ids:
            value = float(sample_map.get(sample_id, 0.0))
            pathway_total += 1
            if value >= min_pathway_abundance_count:
                filtered_samples[sample_id] = value
            elif value > 0:
                pathway_removed += 1
        if filtered_samples:
            filtered_pathway_abundance[feature_id] = filtered_samples
    threshold_rows.append(
        {
            "dataset": "pathway_abundance",
            "threshold_type": "min_pathway_abundance_count",
            "threshold_value": min_pathway_abundance_count,
            "retained_nonzero_values": sum(len(sample_map) for sample_map in filtered_pathway_abundance.values()),
            "removed_nonzero_values": pathway_removed,
            "total_sample_feature_cells": pathway_total,
        }
    )

    filtered_completeness: dict[str, dict[str, float]] = {}
    filtered_observed: dict[str, dict[str, float]] = {}
    completeness_removed = 0
    completeness_total = 0
    for feature_id, sample_map in completeness_map.items():
        filtered_samples: dict[str, float] = {}
        filtered_observed_samples: dict[str, float] = {}
        total_kos = float(total_kos_map.get(feature_id, 0.0))
        for sample_id in sample_ids:
            value = float(sample_map.get(sample_id, 0.0))
            completeness_total += 1
            observed_estimate = (value / 100.0) * total_kos if total_kos > 0 else 0.0
            if total_kos >= min_total_kos and observed_estimate >= min_observed_kos and value > 0:
                filtered_samples[sample_id] = value
                filtered_observed_samples[sample_id] = observed_estimate
            elif value > 0:
                completeness_removed += 1
        if filtered_samples:
            filtered_completeness[feature_id] = filtered_samples
            filtered_observed[feature_id] = filtered_observed_samples
    threshold_rows.append(
        {
            "dataset": "pathway_completeness",
            "threshold_type": "min_total_kos_and_min_observed_kos",
            "threshold_value": f"total>={min_total_kos}; observed>={min_observed_kos}",
            "retained_nonzero_values": sum(len(sample_map) for sample_map in filtered_completeness.values()),
            "removed_nonzero_values": completeness_removed,
            "total_sample_feature_cells": completeness_total,
        }
    )

    filtered_ko_counts: dict[str, dict[str, float]] = {}
    ko_removed = 0
    ko_total = 0
    for feature_id, sample_map in ko_count_map.items():
        filtered_samples: dict[str, float] = {}
        for sample_id in sample_ids:
            value = float(sample_map.get(sample_id, 0.0))
            ko_total += 1
            if value >= min_ko_count:
                filtered_samples[sample_id] = value
            elif value > 0:
                ko_removed += 1
        if filtered_samples:
            filtered_ko_counts[feature_id] = filtered_samples
    threshold_rows.append(
        {
            "dataset": "ko_counts",
            "threshold_type": "min_ko_count",
            "threshold_value": min_ko_count,
            "retained_nonzero_values": sum(len(sample_map) for sample_map in filtered_ko_counts.values()),
            "removed_nonzero_values": ko_removed,
            "total_sample_feature_cells": ko_total,
        }
    )

    return filtered_pathway_abundance, filtered_completeness, filtered_ko_counts, threshold_rows


def build_matrix(
    feature_map: dict[str, dict[str, float]],
    sample_ids: list[str],
    feature_ids: list[str],
) -> list[list[float]]:
    matrix: list[list[float]] = []
    for feature_id in feature_ids:
        row = [float(feature_map.get(feature_id, {}).get(sample_id, 0.0)) for sample_id in sample_ids]
        matrix.append(row)
    return matrix


def write_named_matrix(
    path: Path,
    feature_ids: list[str],
    feature_labels: dict[str, str],
    sample_ids: list[str],
    matrix: list[list[float]],
    extra_column_name: str | None = None,
    extra_values: dict[str, float] | None = None,
) -> None:
    fieldnames = ["feature_id", "feature_name"]
    if extra_column_name:
        fieldnames.append(extra_column_name)
    fieldnames.extend(sample_ids)

    rows: list[dict[str, object]] = []
    for feature_id, values in zip(feature_ids, matrix):
        row: dict[str, object] = {
            "feature_id": feature_id,
            "feature_name": feature_labels.get(feature_id, feature_id),
        }
        if extra_column_name and extra_values is not None:
            row[extra_column_name] = extra_values.get(feature_id, "")
        for sample_id, value in zip(sample_ids, values):
            row[sample_id] = round(value, 6)
        rows.append(row)
    write_tsv(path, fieldnames, rows)


def matrix_to_sample_metrics(
    feature_ids: list[str],
    sample_ids: list[str],
    matrix: list[list[float]],
    log_transform: bool = False,
) -> dict[str, dict[str, list[float]]]:
    sample_vectors: dict[str, list[float]] = {sample_id: [] for sample_id in sample_ids}
    for row in matrix:
        transformed = [log10_1p(value) if log_transform else value for value in row]
        for sample_id, value in zip(sample_ids, transformed):
            sample_vectors[sample_id].append(value)
    return {"vectors": sample_vectors}


def build_pairwise_statistics(
    sample_ids: list[str],
    pathway_matrix: list[list[float]],
    completeness_matrix: list[list[float]],
    ko_matrix: list[list[float]],
) -> list[dict[str, object]]:
    pathway_vectors_raw = {sample_id: [] for sample_id in sample_ids}
    pathway_vectors_log = {sample_id: [] for sample_id in sample_ids}
    completeness_vectors = {sample_id: [] for sample_id in sample_ids}
    ko_vectors_raw = {sample_id: [] for sample_id in sample_ids}
    ko_vectors_log = {sample_id: [] for sample_id in sample_ids}

    for row in pathway_matrix:
        for sample_id, value in zip(sample_ids, row):
            pathway_vectors_raw[sample_id].append(value)
            pathway_vectors_log[sample_id].append(log10_1p(value))

    for row in completeness_matrix:
        for sample_id, value in zip(sample_ids, row):
            completeness_vectors[sample_id].append(value)

    for row in ko_matrix:
        for sample_id, value in zip(sample_ids, row):
            ko_vectors_raw[sample_id].append(value)
            ko_vectors_log[sample_id].append(log10_1p(value))

    rows: list[dict[str, object]] = []
    for left_index in range(len(sample_ids)):
        for right_index in range(left_index + 1, len(sample_ids)):
            sample_a = sample_ids[left_index]
            sample_b = sample_ids[right_index]
            rows.append(
                {
                    "sample_a": sample_a,
                    "sample_b": sample_b,
                    "pathway_abundance_pearson_log": round(
                        pearson_correlation(pathway_vectors_log[sample_a], pathway_vectors_log[sample_b]), 6
                    ),
                    "pathway_abundance_bray_curtis": round(
                        bray_curtis(pathway_vectors_raw[sample_a], pathway_vectors_raw[sample_b]), 6
                    ),
                    "pathway_abundance_jensen_shannon": round(
                        jensen_shannon_divergence(
                            pathway_vectors_raw[sample_a],
                            pathway_vectors_raw[sample_b],
                        ),
                        6,
                    ),
                    "pathway_presence_jaccard": round(
                        jaccard_presence(pathway_vectors_raw[sample_a], pathway_vectors_raw[sample_b]), 6
                    ),
                    "pathway_completeness_pearson": round(
                        pearson_correlation(completeness_vectors[sample_a], completeness_vectors[sample_b]), 6
                    ),
                    "pathway_completeness_mean_abs_diff": round(
                        mean_absolute_difference(completeness_vectors[sample_a], completeness_vectors[sample_b]), 6
                    ),
                    "ko_count_pearson_log": round(
                        pearson_correlation(ko_vectors_log[sample_a], ko_vectors_log[sample_b]), 6
                    ),
                    "ko_count_bray_curtis": round(
                        bray_curtis(ko_vectors_raw[sample_a], ko_vectors_raw[sample_b]), 6
                    ),
                    "ko_count_jensen_shannon": round(
                        jensen_shannon_divergence(
                            ko_vectors_raw[sample_a],
                            ko_vectors_raw[sample_b],
                        ),
                        6,
                    ),
                    "ko_presence_jaccard": round(
                        jaccard_presence(ko_vectors_raw[sample_a], ko_vectors_raw[sample_b]), 6
                    ),
                }
            )
    return rows


def build_square_metric_matrix(
    sample_ids: list[str],
    matrix: list[list[float]],
    metric: str,
    log_transform: bool = False,
) -> list[list[float]]:
    sample_vectors = {sample_id: [] for sample_id in sample_ids}
    for row in matrix:
        transformed = [log10_1p(value) if log_transform else value for value in row]
        for sample_id, value in zip(sample_ids, transformed):
            sample_vectors[sample_id].append(value)

    result: list[list[float]] = []
    for sample_a in sample_ids:
        row: list[float] = []
        for sample_b in sample_ids:
            if metric == "pearson":
                value = pearson_correlation(sample_vectors[sample_a], sample_vectors[sample_b])
            elif metric == "bray_curtis":
                value = bray_curtis(sample_vectors[sample_a], sample_vectors[sample_b])
            elif metric == "jensen_shannon":
                value = jensen_shannon_divergence(sample_vectors[sample_a], sample_vectors[sample_b])
            else:
                raise ValueError(f"Unsupported metric: {metric}")
            row.append(value)
        result.append(row)
    return result


def write_order_file(path: Path, header: str, values: list[str]) -> None:
    rows = [{header: value, "rank": index + 1} for index, value in enumerate(values)]
    write_tsv(path, [header, "rank"], rows)


def render_html_index(output_dir: Path) -> None:
    chart_files = [
        ("cross_sample_pathway_abundance_clustered_heatmap.svg", "Clustered specific-pathway abundance heatmap"),
        ("cross_sample_pathway_completeness_clustered_heatmap.svg", "Clustered pathway completeness heatmap"),
        ("cross_sample_ko_count_clustered_heatmap.svg", "Clustered KO count heatmap"),
        ("cross_sample_pathway_abundance_sample_correlation.svg", "Sample correlation heatmap from pathway abundance"),
        ("cross_sample_pathway_abundance_sample_jsd.svg", "Sample Jensen-Shannon divergence heatmap from pathway abundance"),
        ("cross_sample_ko_count_sample_correlation.svg", "Sample correlation heatmap from KO counts"),
        ("cross_sample_ko_count_sample_jsd.svg", "Sample Jensen-Shannon divergence heatmap from KO counts"),
    ]
    sections = []
    for filename, caption in chart_files:
        sections.append(
            f'<section><h2>{html.escape(caption)}</h2>'
            f'<img src="{html.escape(filename)}" alt="{html.escape(caption)}" style="max-width: 100%; border: 1px solid #ddd;" /></section>'
        )
    html_doc = (
        "<!DOCTYPE html><html><head><meta charset='utf-8'><title>Cross-Sample Metabolic Profiling</title>"
        "<style>body{font-family:Arial,sans-serif;margin:32px;line-height:1.4}section{margin-bottom:40px}</style>"
        "</head><body><h1>Cross-Sample Metabolic Profiling Analysis</h1>"
        + "".join(sections)
        + "</body></html>"
    )
    write_text(output_dir / "cross_sample_chart_index.html", html_doc)


def main() -> int:
    args = parse_args()
    input_root = Path(args.input_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_dirs = [
        (path.name.replace("_charts", ""), path)
        for path in input_root.iterdir()
        if path.is_dir() and path.name.endswith("_charts") and path.name != "1A_charts_test"
    ]
    sample_dirs.sort(key=lambda item: natural_sample_key(item[0]))
    sample_ids = [sample_id for sample_id, _ in sample_dirs]

    (
        pathway_abundance_map,
        pathway_name_map,
        completeness_map,
        completeness_name_map,
        total_kos_map,
        ko_count_map,
    ) = build_feature_maps(sample_dirs)

    (
        pathway_abundance_map,
        completeness_map,
        ko_count_map,
        threshold_rows,
    ) = apply_thresholds(
        sample_ids=sample_ids,
        pathway_abundance_map=pathway_abundance_map,
        completeness_map=completeness_map,
        total_kos_map=total_kos_map,
        ko_count_map=ko_count_map,
        min_pathway_abundance_count=args.min_pathway_abundance_count,
        min_ko_count=args.min_ko_count,
        min_observed_kos=args.min_observed_kos,
        min_total_kos=args.min_total_kos,
    )
    write_tsv(
        output_dir / "cross_sample_threshold_summary.tsv",
        [
            "dataset",
            "threshold_type",
            "threshold_value",
            "retained_nonzero_values",
            "removed_nonzero_values",
            "total_sample_feature_cells",
        ],
        threshold_rows,
    )

    all_specific_pathways = sorted(
        [feature_id for feature_id in pathway_abundance_map if feature_id not in OVERVIEW_PATHWAYS],
        key=lambda feature_id: (
            -sum(pathway_abundance_map[feature_id].values()),
            feature_id,
        ),
    )
    pathway_matrix = build_matrix(pathway_abundance_map, sample_ids, all_specific_pathways)
    write_named_matrix(
        output_dir / "cross_sample_pathway_abundance_matrix.tsv",
        all_specific_pathways,
        {feature_id: pathway_name_map.get(feature_id, feature_id) for feature_id in all_specific_pathways},
        sample_ids,
        pathway_matrix,
    )

    all_completeness_pathways = sorted(
        [feature_id for feature_id in completeness_map if feature_id not in OVERVIEW_PATHWAYS],
        key=lambda feature_id: (
            -sum(completeness_map[feature_id].values()),
            feature_id,
        ),
    )
    completeness_matrix = build_matrix(completeness_map, sample_ids, all_completeness_pathways)
    write_named_matrix(
        output_dir / "cross_sample_pathway_completeness_matrix.tsv",
        all_completeness_pathways,
        {feature_id: completeness_name_map.get(feature_id, feature_id) for feature_id in all_completeness_pathways},
        sample_ids,
        completeness_matrix,
        extra_column_name="total_kos",
        extra_values=total_kos_map,
    )

    all_ko_ids = sorted(
        ko_count_map,
        key=lambda ko_id: (-sum(ko_count_map[ko_id].values()), ko_id),
    )
    ko_matrix = build_matrix(ko_count_map, sample_ids, all_ko_ids)
    write_named_matrix(
        output_dir / "cross_sample_ko_count_matrix.tsv",
        all_ko_ids,
        {ko_id: ko_id for ko_id in all_ko_ids},
        sample_ids,
        ko_matrix,
    )

    pathway_variability_rows: list[dict[str, object]] = []
    for feature_id, row in zip(all_specific_pathways, pathway_matrix):
        transformed = [log10_1p(value) for value in row]
        pathway_variability_rows.append(
            {
                "feature_id": feature_id,
                "feature_name": pathway_name_map.get(feature_id, feature_id),
                "prevalence_count": sum(1 for value in row if value > 0),
                "log10_std": round(std(transformed), 6),
                "raw_sum": round(sum(row), 6),
            }
        )
    pathway_variability_rows.sort(key=lambda row: (-float(row["log10_std"]), -float(row["raw_sum"]), str(row["feature_id"])))
    write_tsv(
        output_dir / "cross_sample_pathway_abundance_variability.tsv",
        ["feature_id", "feature_name", "prevalence_count", "log10_std", "raw_sum"],
        pathway_variability_rows,
    )

    completeness_variability_rows: list[dict[str, object]] = []
    for feature_id, row in zip(all_completeness_pathways, completeness_matrix):
        completeness_variability_rows.append(
            {
                "feature_id": feature_id,
                "feature_name": completeness_name_map.get(feature_id, feature_id),
                "total_kos": total_kos_map.get(feature_id, 0.0),
                "prevalence_count": sum(1 for value in row if value > 0),
                "std": round(std(row), 6),
                "mean": round(mean(row), 6),
            }
        )
    completeness_variability_rows.sort(key=lambda row: (-float(row["std"]), -float(row["mean"]), str(row["feature_id"])))
    write_tsv(
        output_dir / "cross_sample_pathway_completeness_variability.tsv",
        ["feature_id", "feature_name", "total_kos", "prevalence_count", "std", "mean"],
        completeness_variability_rows,
    )

    ko_variability_rows: list[dict[str, object]] = []
    for feature_id, row in zip(all_ko_ids, ko_matrix):
        transformed = [log10_1p(value) for value in row]
        ko_variability_rows.append(
            {
                "feature_id": feature_id,
                "prevalence_count": sum(1 for value in row if value > 0),
                "log10_std": round(std(transformed), 6),
                "raw_sum": round(sum(row), 6),
            }
        )
    ko_variability_rows.sort(key=lambda row: (-float(row["log10_std"]), -float(row["raw_sum"]), str(row["feature_id"])))
    write_tsv(
        output_dir / "cross_sample_ko_variability.tsv",
        ["feature_id", "prevalence_count", "log10_std", "raw_sum"],
        ko_variability_rows,
    )

    top_pathway_ids = [
        row["feature_id"]
        for row in pathway_variability_rows
        if int(row["prevalence_count"]) >= 3
    ][:40]
    top_pathway_matrix_raw = build_matrix(pathway_abundance_map, sample_ids, top_pathway_ids)
    top_pathway_matrix_heat = [zscore_row([log10_1p(value) for value in row]) for row in top_pathway_matrix_raw]

    top_completeness_ids = [
        row["feature_id"]
        for row in completeness_variability_rows
        if float(row["total_kos"]) >= 20 and int(row["prevalence_count"]) >= 3
    ][:40]
    top_completeness_matrix = build_matrix(completeness_map, sample_ids, top_completeness_ids)

    top_ko_ids = [
        row["feature_id"]
        for row in ko_variability_rows
        if int(row["prevalence_count"]) >= 3
    ][:60]
    top_ko_matrix_raw = build_matrix(ko_count_map, sample_ids, top_ko_ids)
    top_ko_matrix_heat = [zscore_row([log10_1p(value) for value in row]) for row in top_ko_matrix_raw]

    kegg_ko_annotations: dict[str, dict[str, str]] = {}
    kegg_ko_cache_path = output_dir / "cross_sample_ko_kegg_annotations.tsv"
    try:
        kegg_ko_annotations = fetch_kegg_ko_labels(top_ko_ids, kegg_ko_cache_path)
    except (error.URLError, OSError, TimeoutError):
        kegg_ko_annotations = load_kegg_ko_cache(kegg_ko_cache_path)

    pathway_labels = {feature_id: pathway_name_map.get(feature_id, feature_id) for feature_id in top_pathway_ids}
    completeness_labels = {feature_id: completeness_name_map.get(feature_id, feature_id) for feature_id in top_completeness_ids}
    ko_labels = {
        feature_id: kegg_ko_annotations.get(feature_id, {}).get("label", feature_id)
        for feature_id in top_ko_ids
    }

    ordered_pathways, ordered_pathway_samples = render_clustered_heatmap(
        output_dir / "cross_sample_pathway_abundance_clustered_heatmap.svg",
        "Cross-Sample Specific Pathway Abundance (Thresholded)",
        (
            "Rows are pathways, columns follow the natural sample order; cells show row-wise z-scores of "
            f"log10(1 + KO count) after removing values < {args.min_pathway_abundance_count:g}"
        ),
        sample_ids,
        top_pathway_ids,
        pathway_labels,
        top_pathway_matrix_heat,
        "vibrant_diverging",
        "Row z-score",
        cluster_columns=False,
    )
    ordered_completeness, ordered_completeness_samples = render_clustered_heatmap(
        output_dir / "cross_sample_pathway_completeness_clustered_heatmap.svg",
        "Cross-Sample Pathway Completeness (Thresholded)",
        (
            "Rows clustered by completeness profiles; columns follow the natural sample order after filtering "
            f"entries with total KOs < {args.min_total_kos:g} or observed KOs < {args.min_observed_kos:g}"
        ),
        sample_ids,
        top_completeness_ids,
        completeness_labels,
        top_completeness_matrix,
        "vibrant_completeness",
        "Completeness %",
        cluster_columns=False,
    )
    ordered_kos, ordered_ko_samples = render_clustered_heatmap(
        output_dir / "cross_sample_ko_count_clustered_heatmap.svg",
        "Cross-Sample KO Count Clustering (Thresholded)",
        (
            "Rows are KO IDs, columns follow the natural sample order; cells show row-wise z-scores of "
            f"log10(1 + KO count) after removing values < {args.min_ko_count:g}"
        ),
        sample_ids,
        top_ko_ids,
        ko_labels,
        top_ko_matrix_heat,
        "vibrant_diverging",
        "Row z-score",
        cluster_columns=False,
    )

    write_tsv(
        output_dir / "cross_sample_ko_heatmap_labels.tsv",
        ["feature_id", "symbol", "name", "label"],
        [
            {
                "feature_id": ko_id,
                "symbol": kegg_ko_annotations.get(ko_id, {}).get("symbol", ""),
                "name": kegg_ko_annotations.get(ko_id, {}).get("name", ""),
                "label": ko_labels.get(ko_id, ko_id),
            }
            for ko_id in ordered_kos
        ],
    )

    write_order_file(output_dir / "cross_sample_pathway_abundance_row_order.tsv", "feature_id", ordered_pathways)
    write_order_file(output_dir / "cross_sample_pathway_abundance_sample_order.tsv", "sample_id", ordered_pathway_samples)
    write_order_file(output_dir / "cross_sample_pathway_completeness_row_order.tsv", "feature_id", ordered_completeness)
    write_order_file(output_dir / "cross_sample_pathway_completeness_sample_order.tsv", "sample_id", ordered_completeness_samples)
    write_order_file(output_dir / "cross_sample_ko_row_order.tsv", "feature_id", ordered_kos)
    write_order_file(output_dir / "cross_sample_ko_sample_order.tsv", "sample_id", ordered_ko_samples)

    pairwise_rows = build_pairwise_statistics(sample_ids, pathway_matrix, completeness_matrix, ko_matrix)
    write_tsv(
        output_dir / "cross_sample_pairwise_statistics.tsv",
        [
            "sample_a",
            "sample_b",
            "pathway_abundance_pearson_log",
            "pathway_abundance_bray_curtis",
            "pathway_abundance_jensen_shannon",
            "pathway_presence_jaccard",
            "pathway_completeness_pearson",
            "pathway_completeness_mean_abs_diff",
            "ko_count_pearson_log",
            "ko_count_bray_curtis",
            "ko_count_jensen_shannon",
            "ko_presence_jaccard",
        ],
        pairwise_rows,
    )

    pathway_sample_corr = build_square_metric_matrix(sample_ids, pathway_matrix, metric="pearson", log_transform=True)
    pathway_sample_jsd = build_square_metric_matrix(sample_ids, pathway_matrix, metric="jensen_shannon")
    ko_sample_corr = build_square_metric_matrix(sample_ids, ko_matrix, metric="pearson", log_transform=True)
    ko_sample_jsd = build_square_metric_matrix(sample_ids, ko_matrix, metric="jensen_shannon")

    pathway_corr_order = list(range(len(sample_ids)))
    ordered_corr_samples = [sample_ids[index] for index in pathway_corr_order]
    pathway_corr_ordered = [
        [pathway_sample_corr[row_index][col_index] for col_index in pathway_corr_order]
        for row_index in pathway_corr_order
    ]
    render_square_matrix_heatmap(
        output_dir / "cross_sample_pathway_abundance_sample_correlation.svg",
        "Sample-to-Sample Correlation from Thresholded Pathway Abundance",
        (
            "Pearson correlation computed on log10(1 + abundance) values across non-overview pathways "
            f"after removing values < {args.min_pathway_abundance_count:g}"
        ),
        ordered_corr_samples,
        pathway_corr_ordered,
        value_min=-1.0,
        value_max=1.0,
        legend_title="Correlation",
        scheme="vibrant_diverging",
    )
    render_square_matrix_heatmap(
        output_dir / "cross_sample_pathway_abundance_sample_jsd.svg",
        "Sample-to-Sample Jensen-Shannon Divergence from Thresholded Pathway Abundance",
        (
            "Jensen-Shannon divergence computed on relative pathway abundance profiles across non-overview pathways "
            f"after removing values < {args.min_pathway_abundance_count:g}; lower values indicate more similar samples"
        ),
        ordered_corr_samples,
        [
            [pathway_sample_jsd[row_index][col_index] for col_index in pathway_corr_order]
            for row_index in pathway_corr_order
        ],
        value_min=0.0,
        value_max=1.0,
        legend_title="JSD",
        scheme="vibrant_sequential",
    )

    ko_corr_order = list(range(len(sample_ids)))
    ordered_ko_corr_samples = [sample_ids[index] for index in ko_corr_order]
    ko_corr_ordered = [
        [ko_sample_corr[row_index][col_index] for col_index in ko_corr_order]
        for row_index in ko_corr_order
    ]
    render_square_matrix_heatmap(
        output_dir / "cross_sample_ko_count_sample_correlation.svg",
        "Sample-to-Sample Correlation from Thresholded KO Profiles",
        (
            "Pearson correlation computed on log10(1 + KO count) values across all detected KO IDs "
            f"after removing values < {args.min_ko_count:g}"
        ),
        ordered_ko_corr_samples,
        ko_corr_ordered,
        value_min=-1.0,
        value_max=1.0,
        legend_title="Correlation",
        scheme="vibrant_diverging",
    )
    render_square_matrix_heatmap(
        output_dir / "cross_sample_ko_count_sample_jsd.svg",
        "Sample-to-Sample Jensen-Shannon Divergence from Thresholded KO Profiles",
        (
            "Jensen-Shannon divergence computed on relative KO count profiles across all detected KO IDs "
            f"after removing values < {args.min_ko_count:g}; lower values indicate more similar samples"
        ),
        ordered_ko_corr_samples,
        [
            [ko_sample_jsd[row_index][col_index] for col_index in ko_corr_order]
            for row_index in ko_corr_order
        ],
        value_min=0.0,
        value_max=1.0,
        legend_title="JSD",
        scheme="vibrant_sequential",
    )

    render_html_index(output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
