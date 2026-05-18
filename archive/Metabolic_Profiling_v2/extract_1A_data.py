from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
import zipfile
import xml.etree.ElementTree as et
from pathlib import Path


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS = {"main": MAIN_NS}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract 1A workbook tables into TSV files.")
    parser.add_argument("--workbook", required=True, help="Path to the XLSX workbook.")
    parser.add_argument("--output-dir", required=True, help="Directory for extracted TSV files.")
    return parser.parse_args()


def choose_sheet_name(sheet_targets: dict[str, str], *candidates: str) -> str:
    for candidate in candidates:
        if candidate in sheet_targets:
            return candidate
    available = ", ".join(sorted(sheet_targets))
    raise KeyError(f"None of the candidate sheets {candidates!r} were found. Available sheets: {available}")


def get_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []
    root = et.fromstring(zf.read("xl/sharedStrings.xml"))
    strings: list[str] = []
    for si in root.findall("main:si", NS):
        parts = [node.text or "" for node in si.iterfind(".//main:t", NS)]
        strings.append("".join(parts))
    return strings


def get_sheet_targets(zf: zipfile.ZipFile) -> dict[str, str]:
    workbook_root = et.fromstring(zf.read("xl/workbook.xml"))
    rels_root = et.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    rel_map = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels_root}
    sheet_map: dict[str, str] = {}
    rid_key = f"{{{REL_NS}}}id"
    for sheet in workbook_root.find("main:sheets", NS):
        sheet_map[sheet.attrib["name"]] = "xl/" + rel_map[sheet.attrib[rid_key]]
    return sheet_map


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


def get_sheet_rows(zf: zipfile.ZipFile, target: str, shared_strings: list[str]) -> list[dict[int, str]]:
    root = et.fromstring(zf.read(target))
    rows: list[dict[int, str]] = []
    for row in root.findall(".//main:sheetData/main:row", NS):
        data: dict[int, str] = {}
        for cell in row.findall("main:c", NS):
            data[col_index(cell.attrib["r"])] = cell_value(cell, shared_strings)
        rows.append(data)
    return rows


def to_float(value: str) -> float:
    return float(value.strip())


def normalize_header(value: str) -> str:
    return "".join(ch.lower() for ch in value if ch.isalnum())


def header_index_map(row: dict[int, str]) -> dict[str, int]:
    mapping: dict[str, int] = {}
    for column_index, value in row.items():
        normalized = normalize_header(value)
        if normalized and normalized not in mapping:
            mapping[normalized] = column_index
    return mapping


def is_numeric_string(value: str) -> bool:
    try:
        float(value.strip())
        return True
    except ValueError:
        return False


def infer_annotated_columns(header_row: dict[int, str], sample_row: dict[int, str]) -> tuple[int, int, int]:
    header = header_index_map(header_row)
    pathway_id_pattern = re.compile(r"^(map|ko)\d{5}$", re.IGNORECASE)

    id_col = None
    count_col = None
    name_col = None

    for column_index, value in sorted(sample_row.items()):
        stripped = value.strip()
        if pathway_id_pattern.match(stripped):
            id_col = column_index
        elif is_numeric_string(stripped) and count_col is None:
            count_col = column_index
        elif stripped and name_col is None:
            name_col = column_index

    if id_col is None:
        id_col = header.get("pathwayid", header.get("mapid", 1))
    if count_col is None:
        count_col = header.get("kocount", header.get("count", 2))
    if name_col is None:
        name_col = header.get("name", 3)

    return id_col, count_col, name_col


def infer_ko_columns(header_row: dict[int, str], sample_row: dict[int, str]) -> tuple[int, int]:
    header = header_index_map(header_row)
    ko_pattern = re.compile(r"^K\d{5}$", re.IGNORECASE)

    id_col = None
    count_col = None
    for column_index, value in sorted(sample_row.items()):
        stripped = value.strip()
        if ko_pattern.match(stripped):
            id_col = column_index
        elif is_numeric_string(stripped):
            count_col = column_index

    if id_col is None:
        id_col = header.get("koids", 1)
    if count_col is None:
        count_col = header.get("counts", header.get("count", 2))

    return id_col, count_col


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summary_row(dataset: str, values: list[float]) -> dict[str, object]:
    variance = statistics.variance(values) if len(values) > 1 else 0.0
    return {
        "dataset": dataset,
        "count": len(values),
        "sum": round(sum(values), 4),
        "mean": round(statistics.mean(values), 4),
        "median": round(statistics.median(values), 4),
        "std": round(math.sqrt(variance), 4),
        "min": round(min(values), 4),
        "max": round(max(values), 4),
    }


def main() -> int:
    args = parse_args()
    workbook = Path(args.workbook).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(workbook) as zf:
        shared_strings = get_shared_strings(zf)
        sheet_targets = get_sheet_targets(zf)

        annotated_sheet = choose_sheet_name(sheet_targets, "Annotated_Pathways", "annotated_pathways", "annotated_pathway")
        completeness_sheet = choose_sheet_name(
            sheet_targets,
            "pathway_completeness_annotation",
            "pathway_completeness_annotated",
        )

        annotated_rows = get_sheet_rows(zf, sheet_targets[annotated_sheet], shared_strings)
        completeness_rows = get_sheet_rows(zf, sheet_targets[completeness_sheet], shared_strings)
        ko_rows = get_sheet_rows(zf, sheet_targets["ko_counts"], shared_strings)

    annotated_records: list[dict[str, object]] = []
    annotated_id_col, annotated_count_col, annotated_name_start_col = infer_annotated_columns(
        annotated_rows[0],
        annotated_rows[1],
    )
    for row in annotated_rows[1:]:
        pathway_id = row.get(annotated_id_col, "").strip()
        ko_count = row.get(annotated_count_col, "").strip()
        if not pathway_id or not ko_count:
            continue
        name = " ".join(
            value.strip()
            for col, value in sorted(row.items())
            if col >= annotated_name_start_col and col not in {annotated_id_col, annotated_count_col} and value.strip()
        )
        annotated_records.append(
            {
                "pathway_id": pathway_id,
                "pathway_name": name,
                "ko_count": to_float(ko_count),
            }
        )

    completeness_records: list[dict[str, object]] = []
    completeness_header = header_index_map(completeness_rows[0])
    completeness_id_col = completeness_header.get("pathwayid", 1)
    completeness_observed_col = completeness_header.get("observedkos", 2)
    completeness_total_col = completeness_header.get("totalkos", completeness_header.get("totalkos", 3))
    completeness_pct_col = completeness_header.get("completeness", completeness_header.get("completenesspct", 4))
    completeness_name_col = completeness_header.get("pathwayname", 5)
    for row in completeness_rows[1:]:
        pathway_id = row.get(completeness_id_col, "").strip()
        observed = row.get(completeness_observed_col, "").strip()
        total = row.get(completeness_total_col, "").strip()
        completeness = row.get(completeness_pct_col, "").strip()
        pathway_name = row.get(completeness_name_col, "").strip()
        if not pathway_id or not completeness:
            continue
        completeness_records.append(
            {
                "pathway_id": pathway_id,
                "pathway_name": pathway_name,
                "observed_kos": to_float(observed),
                "total_kos": to_float(total),
                "completeness_pct": to_float(completeness),
            }
        )

    ko_records: list[dict[str, object]] = []
    ko_id_col, ko_count_col = infer_ko_columns(ko_rows[0], ko_rows[1])
    for row in ko_rows[1:]:
        ko_id = row.get(ko_id_col, "").strip()
        count = row.get(ko_count_col, "").strip()
        if not ko_id or not count:
            continue
        ko_records.append({"ko_id": ko_id, "count": to_float(count)})

    annotated_records.sort(key=lambda item: float(item["ko_count"]), reverse=True)
    completeness_records.sort(key=lambda item: float(item["completeness_pct"]), reverse=True)
    ko_records.sort(key=lambda item: float(item["count"]), reverse=True)

    write_tsv(
        output_dir / "annotated_pathways.tsv",
        ["pathway_id", "pathway_name", "ko_count"],
        annotated_records,
    )
    write_tsv(
        output_dir / "pathway_completeness.tsv",
        ["pathway_id", "pathway_name", "observed_kos", "total_kos", "completeness_pct"],
        completeness_records,
    )
    write_tsv(output_dir / "ko_counts.tsv", ["ko_id", "count"], ko_records)

    summary_rows = [
        summary_row("annotated_pathways_ko_count", [float(row["ko_count"]) for row in annotated_records]),
        summary_row("pathway_completeness_pct", [float(row["completeness_pct"]) for row in completeness_records]),
        summary_row("ko_counts", [float(row["count"]) for row in ko_records]),
        {
            "dataset": "pathway_completeness_thresholds",
            "count": len(completeness_records),
            "sum": sum(1 for row in completeness_records if float(row["completeness_pct"]) >= 50.0),
            "mean": sum(1 for row in completeness_records if float(row["completeness_pct"]) >= 75.0),
            "median": sum(1 for row in completeness_records if float(row["completeness_pct"]) >= 90.0),
            "std": "",
            "min": "",
            "max": "",
        },
    ]
    write_tsv(
        output_dir / "summary_statistics.tsv",
        ["dataset", "count", "sum", "mean", "median", "std", "min", "max"],
        summary_rows,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
