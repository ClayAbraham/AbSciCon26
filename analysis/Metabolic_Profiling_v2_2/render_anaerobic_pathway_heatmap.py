from __future__ import annotations

import csv
import html
import math
import xml.etree.ElementTree as et
import zipfile
from pathlib import Path


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS = {"main": MAIN_NS}
WORKBOOK_NAME = "Filters_Samples_Analysis_II.xlsx"
PATHWAY_LIST_NAME = "Filtered Metabolic Pathways.txt"
SAMPLE_ORDER = ["1A", "3A", "5A", "13A", "14A", "42A", "46A", "47A", "49A"]
CATEGORY_ORDER = [
    "Core Metabolism",
    "Anaerobic Metabolism",
    "Xenobiotic Metabolism",
]
PATHWAY_CATEGORY_MAP = {
    "Glycolysis / Gluconeogenesis": "Core Metabolism",
    "Pentose phosphate pathway": "Core Metabolism",
    "Pyruvate metabolism": "Core Metabolism",
    "Citrate cycle (TCA cycle)": "Core Metabolism",
    "Glyoxylate and dicarboxylate metabolism": "Core Metabolism",
    "C5-Branched dibasic acid metabolism": "Core Metabolism",
    "Fatty acid metabolism": "Core Metabolism",
    "Oxidative phosphorylation": "Core Metabolism",
    "Porphyrin metabolism": "Core Metabolism",
    "Purine metabolism": "Core Metabolism",
    "Pyrimidine metabolism": "Core Metabolism",
    "Nucleotide metabolism": "Core Metabolism",
    "Alanine, aspartate and glutamate metabolism": "Core Metabolism",
    "Arginine and proline metabolism": "Core Metabolism",
    "Histidine metabolism": "Core Metabolism",
    "Cysteine and methionine metabolism": "Core Metabolism",
    "Valine, leucine and isoleucine biosynthesis": "Core Metabolism",
    "Valine, leucine and isoleucine degradation": "Core Metabolism",
    "Phenylalanine, tyrosine and tryptophan biosynthesis": "Core Metabolism",
    "Carbon fixation by Calvin cycle": "Anaerobic Metabolism",
    "Other carbon fixation pathways": "Anaerobic Metabolism",
    "Propanoate metabolism": "Anaerobic Metabolism",
    "Butanoate metabolism": "Anaerobic Metabolism",
    "Methane metabolism": "Anaerobic Metabolism",
    "Nitrogen metabolism": "Anaerobic Metabolism",
    "Sulfur metabolism": "Anaerobic Metabolism",
    "Benzoate degradation": "Xenobiotic Metabolism",
    "Bisphenol degradation": "Xenobiotic Metabolism",
    "Nitrotoluene degradation": "Xenobiotic Metabolism",
    "Styrene degradation": "Xenobiotic Metabolism",
    "Toluene degradation": "Xenobiotic Metabolism",
    "Xylene degradation": "Xenobiotic Metabolism",
}
CATEGORY_DISPLAY_LINES = {
    "Core Metabolism": ["Core", "metabolism"],
    "Anaerobic Metabolism": ["Anaerobic", "metabolism"],
    "Xenobiotic Metabolism": ["Xenobiotic", "metabolism"],
}


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def get_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []
    root = et.fromstring(zf.read("xl/sharedStrings.xml"))
    return ["".join(node.text or "" for node in si.iterfind(".//main:t", NS)) for si in root.findall("main:si", NS)]


def get_sheet_targets(zf: zipfile.ZipFile) -> list[tuple[str, str]]:
    workbook_root = et.fromstring(zf.read("xl/workbook.xml"))
    rels_root = et.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    rel_map = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels_root}
    rid_key = f"{{{REL_NS}}}id"
    return [(sheet.attrib["name"], "xl/" + rel_map[sheet.attrib[rid_key]]) for sheet in workbook_root.find("main:sheets", NS)]


def col_index(cell_ref: str) -> int:
    letters = "".join(ch for ch in cell_ref if ch.isalpha()).upper()
    index = 0
    for ch in letters:
        index = index * 26 + (ord(ch) - ord("A") + 1)
    return index


def cell_value(cell: et.Element, shared_strings: list[str]) -> str:
    cell_type = cell.attrib.get("t")
    value = cell.find("main:v", NS)
    if value is not None:
        raw = value.text or ""
        if cell_type == "s":
            return shared_strings[int(raw)]
        return raw
    inline = cell.find("main:is", NS)
    if inline is not None:
        return "".join(node.text or "" for node in inline.iterfind(".//main:t", NS))
    return ""


def get_visible_sheet_rows(path: Path) -> list[dict[int, str]]:
    with zipfile.ZipFile(path) as zf:
        shared_strings = get_shared_strings(zf)
        _sheet_name, sheet_target = get_sheet_targets(zf)[0]
        root = et.fromstring(zf.read(sheet_target))
        rows: list[dict[int, str]] = []
        for row in root.findall(".//main:sheetData/main:row", NS):
            if row.attrib.get("hidden") == "1":
                continue
            data: dict[int, str] = {}
            for cell in row.findall("main:c", NS):
                data[col_index(cell.attrib["r"])] = cell_value(cell, shared_strings)
            rows.append(data)
    return rows


def load_pathway_list(path: Path) -> list[str]:
    pathways: list[str] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("sample_id\t"):
                break
            if line not in seen:
                pathways.append(line)
                seen.add(line)
    return pathways


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


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


def pairwise_distance(vectors: list[list[float]]) -> dict[tuple[int, int], float]:
    distances: dict[tuple[int, int], float] = {}
    for index_a in range(len(vectors)):
        for index_b in range(index_a + 1, len(vectors)):
            corr = pearson_correlation(vectors[index_a], vectors[index_b])
            distances[(index_a, index_b)] = 1.0 - corr
    return distances


def average_cluster_distance(members_a: list[int], members_b: list[int], distances: dict[tuple[int, int], float]) -> float:
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
    clusters = [{"members": [index], "order": [index], "min_member": index} for index in range(len(vectors))]
    while len(clusters) > 1:
        best_pair = None
        best_score = None
        for left_index in range(len(clusters)):
            for right_index in range(left_index + 1, len(clusters)):
                left = clusters[left_index]
                right = clusters[right_index]
                distance = average_cluster_distance(left["members"], right["members"], distances)
                score = (round(distance, 12), min(left["min_member"], right["min_member"]), max(left["min_member"], right["min_member"]))
                if best_score is None or score < best_score:
                    best_score = score
                    best_pair = (left_index, right_index)
        left_index, right_index = best_pair
        left = clusters[left_index]
        right = clusters[right_index]
        merged = {
            "members": left["members"] + right["members"],
            "order": left["order"] + right["order"],
            "min_member": min(left["min_member"], right["min_member"]),
        }
        clusters = [cluster for index, cluster in enumerate(clusters) if index not in {left_index, right_index}] + [merged]
    return clusters[0]["order"]


def pathway_category(pathway_name: str) -> str:
    return PATHWAY_CATEGORY_MAP.get(pathway_name, "Core Metabolism")


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


def multiline_rotated_text_svg(x: float, y: float, lines: list[str], line_spacing: float = 16.0, **attrs: str | float) -> str:
    attr_text = " ".join(f'{key}="{html.escape(str(value), quote=True)}"' for key, value in attrs.items())
    initial_dy = 0.0 if len(lines) == 1 else -line_spacing * (len(lines) - 1) / 2.0
    tspans: list[str] = []
    for index, line in enumerate(lines):
        dy_value = initial_dy if index == 0 else line_spacing
        tspans.append(f'<tspan x="{x:.2f}" dy="{dy_value:.2f}">{html.escape(line)}</tspan>')
    return f'<text x="{x:.2f}" y="{y:.2f}" {attr_text}>{"".join(tspans)}</text>'


def svg_header(width: int, height: int) -> str:
    return f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'


def text_svg(x: float, y: float, text: str, **attrs: str | float) -> str:
    attr_text = " ".join(f'{key}="{html.escape(str(value), quote=True)}"' for key, value in attrs.items())
    return f'<text x="{x:.2f}" y="{y:.2f}" {attr_text}>{html.escape(text)}</text>'


def rect_svg(x: float, y: float, width: float, height: float, **attrs: str | float) -> str:
    attr_text = " ".join(f'{key}="{html.escape(str(value), quote=True)}"' for key, value in attrs.items())
    return f'<rect x="{x:.2f}" y="{y:.2f}" width="{width:.2f}" height="{height:.2f}" {attr_text} />'


def line_svg(x1: float, y1: float, x2: float, y2: float, **attrs: str | float) -> str:
    attr_text = " ".join(f'{key}="{html.escape(str(value), quote=True)}"' for key, value in attrs.items())
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" {attr_text} />'


def vibrant_sequential_color(value: float, minimum: float, maximum: float) -> str:
    if maximum <= minimum:
        ratio = 0.0
    else:
        ratio = (value - minimum) / (maximum - minimum)
    ratio = max(0.0, min(1.0, ratio))
    stops = [
        (0.00, (19, 35, 84)),
        (0.25, (0, 119, 182)),
        (0.50, (0, 191, 99)),
        (0.75, (255, 214, 10)),
        (1.00, (249, 65, 68)),
    ]
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


def main() -> int:
    root = Path(__file__).resolve().parent
    workbook_path = root / WORKBOOK_NAME
    pathway_list_path = root / PATHWAY_LIST_NAME
    output_dir = root / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    visible_rows = get_visible_sheet_rows(workbook_path)
    header = visible_rows[0]
    header_map = {value.strip(): index for index, value in header.items()}
    sample_col = header_map["sample_id"]
    pathway_id_col = header_map["pathway_id"]
    pathway_name_col = header_map["pathway_name"]
    count_col = header_map["annotated_pathway_ko_count"]

    requested_pathways = load_pathway_list(pathway_list_path)
    requested_set = set(requested_pathways)

    pathway_names: dict[str, str] = {}
    matrix_map: dict[str, dict[str, float]] = {}
    used_rows: list[dict[str, object]] = []
    found_names: set[str] = set()

    for row in visible_rows[1:]:
        pathway_id = row.get(pathway_id_col, "").strip()
        pathway_name = row.get(pathway_name_col, "").strip()
        if not pathway_id or pathway_name not in requested_set:
            continue
        sample_id = row.get(sample_col, "").strip()
        annotated_count = float(row.get(count_col, "0") or 0.0)
        pathway_names[pathway_id] = pathway_name
        matrix_map.setdefault(pathway_id, {})[sample_id] = annotated_count
        used_rows.append(
            {
                "sample_id": sample_id,
                "pathway_id": pathway_id,
                "pathway_name": pathway_name,
                "annotated_pathway_ko_count": annotated_count,
            }
        )
        found_names.add(pathway_name)

    missing_names = [name for name in requested_pathways if name not in found_names]
    write_tsv(output_dir / "anaerobic_pathway_list_not_found.tsv", ["pathway_name"], [{"pathway_name": name} for name in missing_names])

    sample_ids = [sample_id for sample_id in SAMPLE_ORDER if any(sample_id in sample_map for sample_map in matrix_map.values())]
    grouped_pathway_ids: dict[str, list[str]] = {category: [] for category in CATEGORY_ORDER}
    for pathway_id in sorted(matrix_map, key=lambda pid: pathway_names[pid]):
        grouped_pathway_ids.setdefault(pathway_category(pathway_names[pathway_id]), []).append(pathway_id)

    ordered_entries: list[tuple[str, str]] = []
    for category in CATEGORY_ORDER:
        category_pathway_ids = grouped_pathway_ids.get(category, [])
        if not category_pathway_ids:
            continue
        category_vectors = [
            [math.log10(1.0 + matrix_map[pathway_id].get(sample_id, 0.0)) for sample_id in sample_ids]
            for pathway_id in category_pathway_ids
        ]
        category_order = cluster_order(category_vectors)
        ordered_entries.extend((category_pathway_ids[index], category) for index in category_order)

    ordered_pathway_ids = [pathway_id for pathway_id, _category in ordered_entries]
    ordered_matrix = [
        [math.log10(1.0 + matrix_map[pathway_id].get(sample_id, 0.0)) for sample_id in sample_ids]
        for pathway_id in ordered_pathway_ids
    ]

    matrix_rows: list[dict[str, object]] = []
    for pathway_id in ordered_pathway_ids:
        row: dict[str, object] = {"pathway_id": pathway_id, "pathway_name": pathway_names[pathway_id]}
        for sample_id in sample_ids:
            row[sample_id] = round(matrix_map[pathway_id].get(sample_id, 0.0), 6)
        matrix_rows.append(row)

    write_tsv(output_dir / "anaerobic_filtered_pathway_abundance_matrix.tsv", ["pathway_id", "pathway_name", *sample_ids], matrix_rows)
    write_tsv(output_dir / "anaerobic_filtered_visible_rows_used.tsv", ["sample_id", "pathway_id", "pathway_name", "annotated_pathway_ko_count"], used_rows)
    write_tsv(
        output_dir / "anaerobic_filtered_row_order.tsv",
        ["pathway_id", "pathway_name", "category", "rank"],
        [
            {
                "pathway_id": pathway_id,
                "pathway_name": pathway_names[pathway_id],
                "category": category,
                "rank": index + 1,
            }
            for index, (pathway_id, category) in enumerate(ordered_entries)
        ],
    )

    max_log = max((max(row) for row in ordered_matrix), default=1.0)
    cell_width = 72
    cell_height = 34
    group_gap = 14
    left_margin = 710
    top_margin = 180
    legend_width = 36
    legend_height = 320
    heatmap_right_x = left_margin + len(sample_ids) * cell_width
    label_x = heatmap_right_x + 18
    unique_categories = []
    for _pathway_id, category in ordered_entries:
        if not unique_categories or unique_categories[-1] != category:
            unique_categories.append(category)
    total_group_gap = max(0, len(unique_categories) - 1) * group_gap
    height = 240 + len(ordered_pathway_ids) * cell_height + total_group_gap + 150
    legend_x = heatmap_right_x + 420
    legend_y = top_margin
    width = int(legend_x + 260)

    row_positions: list[dict[str, object]] = []
    current_y = top_margin
    previous_category = None
    for pathway_id, category in ordered_entries:
        if previous_category is not None and category != previous_category:
            current_y += group_gap
        row_positions.append({"pathway_id": pathway_id, "category": category, "y": current_y})
        current_y += cell_height
        previous_category = category

    parts = [svg_header(width, height)]
    parts.append(rect_svg(0, 0, width, height, fill="white"))
    parts.append(text_svg(width / 2, 52, "Grouped Anaerobic Pathway Heatmap IV", **{"text-anchor": "middle", "font-size": 30, "font-weight": "700", "font-family": "Arial"}))
    parts.append(text_svg(width / 2, 86, "Rows are grouped as core, anaerobic, and xenobiotic metabolism, then ordered by similarity of log10(1 + annotated_pathway_ko_count) within each group", **{"text-anchor": "middle", "font-size": 17, "fill": "#555", "font-family": "Arial"}))

    for sample_index, sample_id in enumerate(sample_ids):
        x = left_margin + sample_index * cell_width + cell_width / 2
        y = top_margin - 22
        parts.append(f'<text x="{x:.2f}" y="{y:.2f}" transform="rotate(-45 {x:.2f} {y:.2f})" text-anchor="end" font-size="16" font-weight="700" font-family="Arial">{html.escape(sample_id)}</text>')

    category_segments: list[dict[str, float | str]] = []
    for row_info in row_positions:
        category = str(row_info["category"])
        y = float(row_info["y"])
        if not category_segments or category_segments[-1]["category"] != category:
            category_segments.append({"category": category, "start_y": y, "end_y": y + cell_height})
        else:
            category_segments[-1]["end_y"] = y + cell_height

    category_x = left_margin - 52
    for segment in category_segments:
        category = str(segment["category"])
        if category == "Core Metabolism":
            category_lines = ["Core", "metabolism"]
        elif category == "Anaerobic Metabolism":
            category_lines = ["Anaerobic", "metabolism"]
        elif category == "Xenobiotic Metabolism":
            category_lines = ["Xenobiotic", "metabolism"]
        else:
            category_lines = wrap_label(category.lower(), width=18, max_lines=2)
        category_center_y = (float(segment["start_y"]) + float(segment["end_y"])) / 2.0
        parts.append(
            multiline_rotated_text_svg(
                category_x,
                category_center_y,
                category_lines,
                line_spacing=18.0,
                transform=f"rotate(-90 {category_x:.2f} {category_center_y:.2f})",
                **{"text-anchor": "middle", "font-size": 16, "font-weight": "700", "font-family": "Arial", "fill": "#333"},
            )
        )

    previous_category = None
    for row_info in row_positions:
        pathway_id = row_info["pathway_id"]
        category = row_info["category"]
        y = float(row_info["y"])
        if previous_category != category:
            if previous_category is not None:
                separator_y = y - (group_gap / 2)
                parts.append(line_svg(left_margin - 18, separator_y, heatmap_right_x, separator_y, stroke="#b8b8b8", **{"stroke-width": 1.5}))
        wrapped = wrap_label(pathway_names[pathway_id], width=40, max_lines=2)
        center_y = y + cell_height / 2 + 5
        label_y = center_y if len(wrapped) == 1 else center_y - 10
        for line_index, line in enumerate(wrapped):
            parts.append(text_svg(label_x, label_y + line_index * 16, line, **{"text-anchor": "start", "font-size": 14, "font-weight": "600", "font-family": "Arial", "fill": "#111"}))
        previous_category = category

    for row_info, row in zip(row_positions, ordered_matrix):
        y = float(row_info["y"])
        for column_index, value in enumerate(row):
            x = left_margin + column_index * cell_width
            parts.append(rect_svg(x, y, cell_width, cell_height, fill=vibrant_sequential_color(value, 0.0, max_log), stroke="#ffffff", **{"stroke-width": 0.5}))

    steps = 100
    for step in range(steps):
        value = max_log - (step / (steps - 1)) * max_log
        color = vibrant_sequential_color(value, 0.0, max_log)
        y = legend_y + step * (legend_height / steps)
        parts.append(rect_svg(legend_x, y, legend_width, legend_height / steps + 1, fill=color, stroke=color))

    parts.append(text_svg(legend_x + legend_width / 2, legend_y - 16, "log10(1+count)", **{"text-anchor": "middle", "font-size": 15, "font-weight": "700", "font-family": "Arial"}))
    raw_ticks = [0, 10, 100, 1000, 10000]
    tick_pairs: list[tuple[int, float]] = []
    for raw_tick in raw_ticks:
        tick_value = math.log10(1.0 + raw_tick)
        if tick_value <= max_log + 1e-9:
            tick_pairs.append((raw_tick, tick_value))
    max_raw = int(round(10**max_log - 1))
    if tick_pairs and tick_pairs[-1][0] != max_raw:
        tick_pairs.append((max_raw, max_log))
    for raw_tick, tick_value in tick_pairs:
        ratio = (max_log - tick_value) / max_log if max_log > 0 else 0.0
        y = legend_y + ratio * legend_height
        parts.append(line_svg(legend_x + legend_width, y, legend_x + legend_width + 10, y, stroke="#333", **{"stroke-width": 1.25}))
        parts.append(text_svg(legend_x + legend_width + 18, y + 5, str(raw_tick), **{"text-anchor": "start", "font-size": 14, "font-family": "Arial", "fill": "#333"}))

    parts.append("</svg>")
    svg_path = output_dir / "anaerobic_clustered_pathway_heatmap_IV.svg"
    write_text(svg_path, "\n".join(parts))

    html_doc = (
        "<!DOCTYPE html><html><head><meta charset='utf-8'><title>Grouped Anaerobic Pathway Heatmap IV</title>"
        "<style>body{font-family:Arial,sans-serif;margin:32px;line-height:1.4}img{max-width:100%;border:1px solid #ddd}</style>"
        "</head><body><h1>Grouped Anaerobic Pathway Heatmap IV</h1>"
        f"<p>Workbook source: {html.escape(workbook_path.name)}. Pathway list source: {html.escape(pathway_list_path.name)}. Hidden rows in the first worksheet were excluded before plotting. The y-axis is organized into core metabolism, anaerobic metabolism, and xenobiotic metabolism, with similarity-based ordering inside each group.</p>"
        f"<img src='{svg_path.name}' alt='Grouped anaerobic pathway heatmap IV' />"
        "</body></html>"
    )
    write_text(output_dir / "anaerobic_clustered_pathway_heatmap_IV_index.html", html_doc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
