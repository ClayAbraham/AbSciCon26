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
SAMPLE_ORDER = ["1A", "3A", "5A", "13A", "14A", "42A", "46A", "47A", "49A"]
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
    strings: list[str] = []
    for si in root.findall("main:si", NS):
        parts = [node.text or "" for node in si.iterfind(".//main:t", NS)]
        strings.append("".join(parts))
    return strings


def get_sheet_targets(zf: zipfile.ZipFile) -> list[tuple[str, str]]:
    workbook_root = et.fromstring(zf.read("xl/workbook.xml"))
    rels_root = et.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    rel_map = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels_root}
    rid_key = f"{{{REL_NS}}}id"
    sheets: list[tuple[str, str]] = []
    for sheet in workbook_root.find("main:sheets", NS):
        sheets.append((sheet.attrib["name"], "xl/" + rel_map[sheet.attrib[rid_key]]))
    return sheets


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
        first_sheet_name, first_target = get_sheet_targets(zf)[0]
        root = et.fromstring(zf.read(first_target))
        rows: list[dict[int, str]] = []
        for row in root.findall(".//main:sheetData/main:row", NS):
            if row.attrib.get("hidden") == "1":
                continue
            data: dict[int, str] = {}
            for cell in row.findall("main:c", NS):
                data[col_index(cell.attrib["r"])] = cell_value(cell, shared_strings)
            rows.append(data)
    return rows


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


def svg_header(width: int, height: int) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
    )


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
    workbook_candidates = sorted(
        [path for path in root.glob("*.xlsx") if not path.name.startswith("~$")],
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not workbook_candidates:
        raise FileNotFoundError("No Excel workbook was found in Metabolic_Profiling_v2_2.")

    workbook_path = workbook_candidates[0]
    output_dir = root / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    visible_rows = get_visible_sheet_rows(workbook_path)
    if not visible_rows:
        raise ValueError("No visible rows were found in the first worksheet.")

    header = visible_rows[0]
    header_map = {value.strip(): index for index, value in header.items()}

    sample_col = header_map["sample_id"]
    pathway_id_col = header_map["pathway_id"]
    pathway_name_col = header_map["pathway_name"]
    count_col = header_map["annotated_pathway_ko_count"]

    pathway_names: dict[str, str] = {}
    matrix_map: dict[str, dict[str, float]] = {}
    source_rows: list[dict[str, object]] = []

    for row in visible_rows[1:]:
        pathway_id = row.get(pathway_id_col, "").strip()
        if not pathway_id or pathway_id in OVERVIEW_PATHWAYS:
            continue
        sample_id = row.get(sample_col, "").strip()
        pathway_name = row.get(pathway_name_col, "").strip()
        annotated_count = float(row.get(count_col, "0") or 0.0)
        pathway_names[pathway_id] = pathway_name
        matrix_map.setdefault(pathway_id, {})[sample_id] = annotated_count
        source_rows.append(
            {
                "sample_id": sample_id,
                "pathway_id": pathway_id,
                "pathway_name": pathway_name,
                "annotated_pathway_ko_count": annotated_count,
            }
        )

    sample_ids = [sample_id for sample_id in SAMPLE_ORDER if any(sample_id in sample_map for sample_map in matrix_map.values())]
    pathway_ids = sorted(
        matrix_map,
        key=lambda pathway_id: (
            -sum(matrix_map[pathway_id].get(sample_id, 0.0) for sample_id in sample_ids),
            pathway_id,
        ),
    )

    matrix_rows: list[dict[str, object]] = []
    values_log: list[list[float]] = []
    for pathway_id in pathway_ids:
        row_values = [matrix_map[pathway_id].get(sample_id, 0.0) for sample_id in sample_ids]
        row_log = [math.log10(1.0 + value) for value in row_values]
        values_log.append(row_log)
        row_payload: dict[str, object] = {
            "pathway_id": pathway_id,
            "pathway_name": pathway_names[pathway_id],
        }
        for sample_id, value in zip(sample_ids, row_values):
            row_payload[sample_id] = round(value, 6)
        matrix_rows.append(row_payload)

    write_tsv(
        output_dir / "excel_filtered_pathway_abundance_matrix.tsv",
        ["pathway_id", "pathway_name", *sample_ids],
        matrix_rows,
    )
    write_tsv(
        output_dir / "excel_filtered_visible_rows_used.tsv",
        ["sample_id", "pathway_id", "pathway_name", "annotated_pathway_ko_count"],
        source_rows,
    )
    write_tsv(
        output_dir / "excluded_overview_pathways.tsv",
        ["pathway_id"],
        [{"pathway_id": pathway_id} for pathway_id in sorted(OVERVIEW_PATHWAYS)],
    )

    max_log = max((max(row) for row in values_log), default=1.0)
    cell_width = 56
    cell_height = 16
    left_margin = 420
    top_margin = 150
    legend_width = 24
    legend_height = 240
    width = 260 + len(sample_ids) * cell_width + 460
    height = 190 + len(pathway_ids) * cell_height + 120
    legend_x = left_margin + len(sample_ids) * cell_width + 60
    legend_y = top_margin

    parts = [svg_header(width, height)]
    parts.append(rect_svg(0, 0, width, height, fill="white"))
    parts.append(
        text_svg(
            width / 2,
            42,
            "Excel-Filtered Non-Overview Pathway Heatmap",
            **{"text-anchor": "middle", "font-size": 22, "font-weight": "700", "font-family": "Arial"},
        )
    )
    parts.append(
        text_svg(
            width / 2,
            70,
            "Built from visible rows in the first worksheet only; KEGG overview pathways were excluded; cells show log10(1 + annotated_pathway_ko_count)",
            **{"text-anchor": "middle", "font-size": 13, "fill": "#555", "font-family": "Arial"},
        )
    )

    for sample_index, sample_id in enumerate(sample_ids):
        x = left_margin + sample_index * cell_width + cell_width / 2
        y = top_margin - 18
        parts.append(
            f'<text x="{x:.2f}" y="{y:.2f}" transform="rotate(-45 {x:.2f} {y:.2f})" '
            f'text-anchor="end" font-size="12" font-family="Arial">{html.escape(sample_id)}</text>'
        )

    for row_index, pathway_id in enumerate(pathway_ids):
        wrapped = wrap_label(pathway_names[pathway_id], width=34, max_lines=2)
        center_y = top_margin + row_index * cell_height + cell_height / 2 + 4
        label_y = center_y if len(wrapped) == 1 else center_y - 6
        for line_index, line in enumerate(wrapped):
            parts.append(
                text_svg(
                    left_margin - 12,
                    label_y + line_index * 12,
                    line,
                    **{"text-anchor": "end", "font-size": 11, "font-family": "Arial", "fill": "#111"},
                )
            )

    for row_index, row in enumerate(values_log):
        for column_index, value in enumerate(row):
            x = left_margin + column_index * cell_width
            y = top_margin + row_index * cell_height
            parts.append(
                rect_svg(
                    x,
                    y,
                    cell_width,
                    cell_height,
                    fill=vibrant_sequential_color(value, 0.0, max_log),
                    stroke="#ffffff",
                    **{"stroke-width": 0.5},
                )
            )

    steps = 100
    for step in range(steps):
        value = max_log - (step / (steps - 1)) * max_log
        color = vibrant_sequential_color(value, 0.0, max_log)
        y = legend_y + step * (legend_height / steps)
        parts.append(rect_svg(legend_x, y, legend_width, legend_height / steps + 1, fill=color, stroke=color))

    parts.append(
        text_svg(
            legend_x + legend_width / 2,
            legend_y - 12,
            "log10(1+count)",
            **{"text-anchor": "middle", "font-size": 12, "font-weight": "700", "font-family": "Arial"},
        )
    )

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
        parts.append(line_svg(legend_x + legend_width, y, legend_x + legend_width + 8, y, stroke="#333", **{"stroke-width": 1}))
        parts.append(
            text_svg(
                legend_x + legend_width + 14,
                y + 4,
                str(raw_tick),
                **{"text-anchor": "start", "font-size": 11, "font-family": "Arial", "fill": "#333"},
            )
        )

    parts.append("</svg>")
    svg_path = output_dir / "excel_filtered_pathway_abundance_heatmap.svg"
    write_text(svg_path, "\n".join(parts))

    html_doc = (
        "<!DOCTYPE html><html><head><meta charset='utf-8'><title>Excel Filtered Pathway Heatmap</title>"
        "<style>body{font-family:Arial,sans-serif;margin:32px;line-height:1.4}img{max-width:100%;border:1px solid #ddd}</style>"
        "</head><body><h1>Excel-Filtered Non-Overview Pathway Heatmap</h1>"
        f"<p>Workbook source: {html.escape(workbook_path.name)}. First worksheet only. Hidden rows from the worksheet filter were excluded before plotting.</p>"
        f"<img src='{svg_path.name}' alt='Excel-filtered non-overview pathway heatmap' />"
        "</body></html>"
    )
    write_text(output_dir / "excel_filtered_pathway_heatmap_index.html", html_doc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
