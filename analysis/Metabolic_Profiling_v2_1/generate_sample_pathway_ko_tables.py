from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from urllib import error, request


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export per-sample pathway-to-observed-KO tables using KEGG pathway membership."
    )
    parser.add_argument(
        "--input-root",
        default="output",
        help="Directory containing sample *_charts folders.",
    )
    parser.add_argument(
        "--cache-dir",
        default="output/kegg_reference",
        help="Directory for cached KEGG pathway and KO annotations.",
    )
    parser.add_argument(
        "--pathway-id",
        action="append",
        default=[],
        help="Optional pathway MapID filter. Repeat to restrict output to specific pathways.",
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
    try:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    except PermissionError:
        fallback_path = path.with_name(f"{path.stem}_refreshed{path.suffix}")
        with fallback_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)


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


def load_ko_cache(path: Path) -> dict[str, dict[str, str]]:
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


def write_ko_cache(path: Path, cache: dict[str, dict[str, str]]) -> None:
    rows = [cache[ko_id] for ko_id in sorted(cache)]
    write_tsv(path, ["ko_id", "symbol", "name", "label"], rows)


def load_pathway_cache(path: Path) -> dict[str, set[str]]:
    if not path.exists():
        return {}
    rows = load_tsv(path)
    mapping: dict[str, set[str]] = {}
    for row in rows:
        pathway_id = row.get("pathway_id", "").strip()
        ko_id = row.get("ko_id", "").strip()
        if pathway_id and ko_id:
            mapping.setdefault(pathway_id, set()).add(ko_id)
    return mapping


def write_pathway_cache(path: Path, mapping: dict[str, set[str]]) -> None:
    rows: list[dict[str, object]] = []
    for pathway_id in sorted(mapping):
        for ko_id in sorted(mapping[pathway_id]):
            rows.append({"pathway_id": pathway_id, "ko_id": ko_id})
    write_tsv(path, ["pathway_id", "ko_id"], rows)


def fetch_pathway_ko_links(
    pathway_ids: list[str],
    cache_path: Path,
    batch_size: int = 10,
) -> dict[str, set[str]]:
    mapping = load_pathway_cache(cache_path)
    missing = [pathway_id for pathway_id in pathway_ids if pathway_id not in mapping]

    for start in range(0, len(missing), batch_size):
        batch = missing[start : start + batch_size]
        if not batch:
            continue
        url = "https://rest.kegg.jp/link/ko/" + "+".join(batch)
        with request.urlopen(url, timeout=30) as response:
            payload = response.read().decode("utf-8", errors="replace")
        for line in payload.splitlines():
            if not line.strip():
                continue
            left, right = line.split("\t", 1)
            pathway_id = left.split(":", 1)[-1].strip()
            ko_id = right.split(":", 1)[-1].strip()
            if pathway_id and ko_id:
                mapping.setdefault(pathway_id, set()).add(ko_id)

        for pathway_id in batch:
            mapping.setdefault(pathway_id, set())

    write_pathway_cache(cache_path, mapping)
    return {pathway_id: mapping.get(pathway_id, set()) for pathway_id in pathway_ids}


def fetch_ko_annotations(
    ko_ids: list[str],
    cache_path: Path,
    batch_size: int = 50,
) -> dict[str, dict[str, str]]:
    cache = load_ko_cache(cache_path)
    missing = [ko_id for ko_id in ko_ids if ko_id not in cache]

    for start in range(0, len(missing), batch_size):
        batch = missing[start : start + batch_size]
        if not batch:
            continue
        url = "https://rest.kegg.jp/get/" + "+".join(f"ko:{ko_id}" for ko_id in batch)
        with request.urlopen(url, timeout=30) as response:
            payload = response.read().decode("utf-8", errors="replace")
        cache.update(parse_kegg_ko_entries(payload))

    write_ko_cache(cache_path, cache)
    return {ko_id: cache[ko_id] for ko_id in ko_ids if ko_id in cache}


def load_sample_metadata(sample_dir: Path) -> tuple[dict[str, dict[str, object]], dict[str, float], dict[str, float]]:
    completeness_rows = load_tsv(sample_dir / "pathway_completeness.tsv")
    ko_rows = load_tsv(sample_dir / "ko_counts.tsv")
    annotated_candidates = sorted(sample_dir.glob("*_annotated_pathways.tsv"))
    default_annotated_path = sample_dir / "annotated_pathways.tsv"
    if default_annotated_path.exists():
        annotated_candidates.append(default_annotated_path)

    pathway_meta: dict[str, dict[str, object]] = {}
    for row in completeness_rows:
        pathway_id = row["pathway_id"].strip()
        pathway_meta[pathway_id] = {
            "pathway_name": row["pathway_name"].strip(),
            "observed_kos": float(row["observed_kos"]),
            "total_kos": float(row["total_kos"]),
            "completeness_pct": float(row["completeness_pct"]),
            "annotated_ko_count": 0.0,
        }

    if annotated_candidates:
        annotated_path = next((path for path in annotated_candidates if path.exists()), annotated_candidates[0])
        for row in load_tsv(annotated_path):
            pathway_id = row["pathway_id"].strip()
            pathway_meta.setdefault(
                pathway_id,
                {
                    "pathway_name": row["pathway_name"].strip(),
                    "observed_kos": 0.0,
                    "total_kos": 0.0,
                    "completeness_pct": 0.0,
                    "annotated_ko_count": 0.0,
                },
            )
            pathway_meta[pathway_id]["annotated_ko_count"] = float(row["ko_count"])

    ko_counts = {row["ko_id"].strip(): float(row["count"]) for row in ko_rows}
    annotated_counts = {
        pathway_id: float(meta["annotated_ko_count"])
        for pathway_id, meta in pathway_meta.items()
    }
    return pathway_meta, ko_counts, annotated_counts


def main() -> int:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    input_root = (script_dir / args.input_root).resolve()
    cache_dir = (script_dir / args.cache_dir).resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)

    sample_dirs = [
        (path.name.replace("_charts", ""), path)
        for path in input_root.iterdir()
        if path.is_dir() and path.name.endswith("_charts")
    ]
    sample_dirs.sort(key=lambda item: natural_sample_key(item[0]))

    pathway_filter = {pathway_id.strip() for pathway_id in args.pathway_id if pathway_id.strip()}

    all_pathway_ids: set[str] = set()
    sample_data: dict[str, dict[str, object]] = {}
    for sample_id, sample_dir in sample_dirs:
        pathway_meta, ko_counts, annotated_counts = load_sample_metadata(sample_dir)
        selected_pathway_ids = set(pathway_meta)
        if pathway_filter:
            selected_pathway_ids &= pathway_filter
        sample_data[sample_id] = {
            "sample_dir": sample_dir,
            "pathway_meta": pathway_meta,
            "ko_counts": ko_counts,
            "annotated_counts": annotated_counts,
            "selected_pathway_ids": selected_pathway_ids,
        }
        all_pathway_ids.update(selected_pathway_ids)

    pathway_cache_path = cache_dir / "kegg_pathway_ko_links.tsv"
    ko_cache_path = cache_dir / "kegg_ko_annotations.tsv"

    try:
        pathway_links = fetch_pathway_ko_links(sorted(all_pathway_ids), pathway_cache_path)
    except (error.URLError, OSError, TimeoutError):
        pathway_links = load_pathway_cache(pathway_cache_path)

    all_observed_kos: set[str] = set()
    for payload in sample_data.values():
        observed = set(payload["ko_counts"])
        for pathway_id in payload["selected_pathway_ids"]:
            all_observed_kos.update(pathway_links.get(pathway_id, set()) & observed)

    try:
        ko_annotations = fetch_ko_annotations(sorted(all_observed_kos), ko_cache_path)
    except (error.URLError, OSError, TimeoutError):
        ko_annotations = load_ko_cache(ko_cache_path)

    combined_rows: list[dict[str, object]] = []
    combined_summary_rows: list[dict[str, object]] = []

    detail_fields = [
        "sample_id",
        "pathway_id",
        "pathway_name",
        "annotated_pathway_ko_count",
        "observed_kos",
        "total_kos",
        "completeness_pct",
        "ko_id",
        "ko_symbol",
        "ko_name",
        "ko_label",
        "sample_ko_count",
    ]
    summary_fields = [
        "sample_id",
        "pathway_id",
        "pathway_name",
        "annotated_pathway_ko_count",
        "observed_kos",
        "total_kos",
        "completeness_pct",
        "matched_observed_kos",
        "matched_ko_count_sum",
    ]

    for sample_id, _sample_dir in sample_dirs:
        payload = sample_data[sample_id]
        sample_dir = payload["sample_dir"]
        pathway_meta = payload["pathway_meta"]
        ko_counts = payload["ko_counts"]
        selected_pathway_ids = sorted(payload["selected_pathway_ids"])

        detail_rows: list[dict[str, object]] = []
        summary_rows: list[dict[str, object]] = []

        for pathway_id in selected_pathway_ids:
            meta = pathway_meta[pathway_id]
            observed_ko_ids = sorted(pathway_links.get(pathway_id, set()) & set(ko_counts))
            matched_count_sum = sum(ko_counts[ko_id] for ko_id in observed_ko_ids)

            summary_rows.append(
                {
                    "sample_id": sample_id,
                    "pathway_id": pathway_id,
                    "pathway_name": meta["pathway_name"],
                    "annotated_pathway_ko_count": round(float(meta["annotated_ko_count"]), 6),
                    "observed_kos": round(float(meta["observed_kos"]), 6),
                    "total_kos": round(float(meta["total_kos"]), 6),
                    "completeness_pct": round(float(meta["completeness_pct"]), 6),
                    "matched_observed_kos": len(observed_ko_ids),
                    "matched_ko_count_sum": round(float(matched_count_sum), 6),
                }
            )

            for ko_id in observed_ko_ids:
                annotation = ko_annotations.get(
                    ko_id,
                    {"symbol": "", "name": "", "label": ko_id},
                )
                detail_rows.append(
                    {
                        "sample_id": sample_id,
                        "pathway_id": pathway_id,
                        "pathway_name": meta["pathway_name"],
                        "annotated_pathway_ko_count": round(float(meta["annotated_ko_count"]), 6),
                        "observed_kos": round(float(meta["observed_kos"]), 6),
                        "total_kos": round(float(meta["total_kos"]), 6),
                        "completeness_pct": round(float(meta["completeness_pct"]), 6),
                        "ko_id": ko_id,
                        "ko_symbol": annotation.get("symbol", ""),
                        "ko_name": annotation.get("name", ""),
                        "ko_label": annotation.get("label", ko_id),
                        "sample_ko_count": round(float(ko_counts[ko_id]), 6),
                    }
                )

        detail_rows.sort(
            key=lambda row: (
                row["pathway_id"],
                -float(row["sample_ko_count"]),
                row["ko_id"],
            )
        )
        summary_rows.sort(
            key=lambda row: (
                -float(row["completeness_pct"]),
                row["pathway_id"],
            )
        )

        write_tsv(sample_dir / f"{sample_id}_observed_pathway_kos.tsv", detail_fields, detail_rows)
        write_tsv(sample_dir / f"{sample_id}_pathway_ko_summary.tsv", summary_fields, summary_rows)
        combined_rows.extend(detail_rows)
        combined_summary_rows.extend(summary_rows)

    combined_rows.sort(
        key=lambda row: (
            natural_sample_key(str(row["sample_id"])),
            str(row["pathway_id"]),
            -float(row["sample_ko_count"]),
            str(row["ko_id"]),
        )
    )
    combined_summary_rows.sort(
        key=lambda row: (
            natural_sample_key(str(row["sample_id"])),
            -float(row["completeness_pct"]),
            str(row["pathway_id"]),
        )
    )

    write_tsv(input_root / "all_samples_observed_pathway_kos.tsv", detail_fields, combined_rows)
    write_tsv(input_root / "all_samples_pathway_ko_summary.tsv", summary_fields, combined_summary_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
