from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


SEARCH_LIST_PATTERNS = [
    "psychrophilic_gene_search*.txt",
    "Pyschrophilic_gene_search*.txt",
]
OUTPUT_DIRNAME = "output"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Search all-sample pathway KO tables for psychrophilic terms."
    )
    parser.add_argument(
        "--source-table",
        type=Path,
        default=None,
        help=(
            "Path to all_samples_observed_pathway_kos_refreshed.tsv. "
            "Defaults to ../Metabolic_Profiling_v2_1/output/all_samples_observed_pathway_kos_refreshed.tsv."
        ),
    )
    parser.add_argument(
        "--search-file",
        type=Path,
        default=None,
        help="Optional explicit psychrophilic search term file.",
    )
    return parser.parse_args()


def find_search_file(root: Path) -> Path:
    candidates: list[Path] = []
    seen: set[Path] = set()
    for pattern in SEARCH_LIST_PATTERNS:
        for candidate in root.glob(pattern):
            if candidate.is_file() and candidate not in seen:
                candidates.append(candidate)
                seen.add(candidate)
    if candidates:
        return max(candidates, key=lambda path: (path.stat().st_mtime, path.name.casefold()))
    raise FileNotFoundError("Could not find psychrophilic gene search list in the v2_2 folder.")


def search_output_suffix(path: Path) -> str:
    suffix_match = re.search(r"(_[A-Za-z0-9]+)$", path.stem)
    if suffix_match:
        return suffix_match.group(1)
    return ""


def load_search_terms(path: Path) -> list[str]:
    terms: list[str] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            term = raw_line.strip()
            if not term:
                continue
            lowered = term.casefold()
            if lowered in seen:
                continue
            seen.add(lowered)
            terms.append(term)
    return terms


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def normalize_text(value: str) -> str:
    return " ".join(value.casefold().strip().split())


def looks_like_gene_symbol(term: str) -> bool:
    if " " in term:
        return False
    return any(character.isupper() for character in term) or any(character.isdigit() for character in term)


def tokenize_symbols(value: str) -> set[str]:
    tokens = [token for token in re.split(r"[^A-Za-z0-9_-]+", value) if token]
    return {token.casefold() for token in tokens}


def phrase_match(term: str, text: str) -> bool:
    normalized_term = normalize_text(term)
    normalized_text = normalize_text(text)
    if not normalized_term or not normalized_text:
        return False
    pattern = r"(?<![a-z0-9])" + re.escape(normalized_term) + r"(?![a-z0-9])"
    return re.search(pattern, normalized_text) is not None


def detect_ko_match(row: dict[str, str], term: str) -> list[str]:
    matched_fields: list[str] = []
    normalized_term = normalize_text(term)
    if normalized_term == normalize_text(row.get("ko_id", "")):
        matched_fields.append("ko_id")

    symbol_tokens = tokenize_symbols(row.get("ko_symbol", ""))
    if normalized_term in symbol_tokens:
        matched_fields.append("ko_symbol")

    if looks_like_gene_symbol(term):
        return matched_fields

    if phrase_match(term, row.get("ko_name", "")):
        matched_fields.append("ko_name")
    if phrase_match(term, row.get("ko_label", "")):
        matched_fields.append("ko_label")
    return matched_fields


def detect_pathway_match(row: dict[str, str], term: str) -> list[str]:
    stripped_term = term.strip()
    allow_pathway_match = (" " in stripped_term) or ("-" in stripped_term) or stripped_term.islower()
    if not allow_pathway_match:
        return []
    matched_fields: list[str] = []
    normalized_term = normalize_text(term)
    if normalized_term == normalize_text(row.get("pathway_id", "")):
        matched_fields.append("pathway_id")
    if phrase_match(term, row.get("pathway_name", "")):
        matched_fields.append("pathway_name")
    return matched_fields


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    target_path = path
    try:
        handle = target_path.open("w", encoding="utf-8", newline="")
    except PermissionError:
        target_path = path.with_name(f"{path.stem}_refreshed{path.suffix}")
        handle = target_path.open("w", encoding="utf-8", newline="")
    with handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def float_or_zero(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parent
    output_dir = root / OUTPUT_DIRNAME
    output_dir.mkdir(parents=True, exist_ok=True)

    search_file = args.search_file or find_search_file(root)
    output_suffix = search_output_suffix(search_file)
    search_terms = load_search_terms(search_file)
    source_table = args.source_table or (
        root.parent
        / "Metabolic_Profiling_v2_1"
        / "output"
        / "all_samples_observed_pathway_kos_refreshed.tsv"
    )
    rows = load_rows(source_table)

    ko_match_rows: list[dict[str, object]] = []
    pathway_match_map: dict[tuple[str, str, str, str], dict[str, object]] = {}
    per_sample_term: dict[tuple[str, str, str], dict[str, object]] = {}
    per_sample_overall: dict[str, dict[str, object]] = {}
    not_found_terms = {term.casefold(): term for term in search_terms}

    for row in rows:
        sample_id = row.get("sample_id", "")
        pathway_id = row.get("pathway_id", "")
        pathway_name = row.get("pathway_name", "")
        ko_id = row.get("ko_id", "")
        ko_symbol = row.get("ko_symbol", "")
        ko_name = row.get("ko_name", "")
        ko_label = row.get("ko_label", "")
        sample_ko_count = float_or_zero(row.get("sample_ko_count", "0"))
        annotated_pathway_ko_count = float_or_zero(row.get("annotated_pathway_ko_count", "0"))
        observed_kos = float_or_zero(row.get("observed_kos", "0"))
        total_kos = float_or_zero(row.get("total_kos", "0"))
        completeness_pct = float_or_zero(row.get("completeness_pct", "0"))

        for term in search_terms:
            ko_fields = detect_ko_match(row, term)
            pathway_fields = detect_pathway_match(row, term)
            if not ko_fields and not pathway_fields:
                continue
            not_found_terms.pop(term.casefold(), None)

            sample_summary = per_sample_overall.setdefault(
                sample_id,
                {
                    "sample_id": sample_id,
                    "matched_ko_terms": set(),
                    "matched_pathway_terms": set(),
                    "matched_ko_rows": 0,
                    "matched_pathway_rows": 0,
                    "unique_kos": set(),
                    "unique_pathways": set(),
                    "sample_ko_count_sum": 0.0,
                },
            )

            if ko_fields:
                ko_match_rows.append(
                    {
                        "search_term": term,
                        "matched_fields": ",".join(ko_fields),
                        "sample_id": sample_id,
                        "pathway_id": pathway_id,
                        "pathway_name": pathway_name,
                        "annotated_pathway_ko_count": annotated_pathway_ko_count,
                        "observed_kos": observed_kos,
                        "total_kos": total_kos,
                        "completeness_pct": completeness_pct,
                        "ko_id": ko_id,
                        "ko_symbol": ko_symbol,
                        "ko_name": ko_name,
                        "ko_label": ko_label,
                        "sample_ko_count": sample_ko_count,
                    }
                )

                term_key = (sample_id, term, "ko")
                term_summary = per_sample_term.setdefault(
                    term_key,
                    {
                        "sample_id": sample_id,
                        "search_term": term,
                        "match_scope": "ko",
                        "matched_rows": 0,
                        "unique_kos": set(),
                        "unique_ko_labels": {},
                        "unique_pathways": set(),
                        "matched_fields": set(),
                        "sample_ko_count_sum": 0.0,
                        "max_sample_ko_count": 0.0,
                    },
                )
                term_summary["matched_rows"] += 1
                if ko_id:
                    term_summary["unique_kos"].add(ko_id)
                    if ko_label:
                        term_summary["unique_ko_labels"][ko_id] = ko_label
                if pathway_id:
                    term_summary["unique_pathways"].add(f"{pathway_id}|{pathway_name}")
                for field_name in ko_fields:
                    term_summary["matched_fields"].add(field_name)
                term_summary["sample_ko_count_sum"] += sample_ko_count
                term_summary["max_sample_ko_count"] = max(term_summary["max_sample_ko_count"], sample_ko_count)

                sample_summary["matched_ko_terms"].add(term)
                sample_summary["matched_ko_rows"] += 1
                if ko_id:
                    sample_summary["unique_kos"].add(ko_id)
                if pathway_id:
                    sample_summary["unique_pathways"].add(f"{pathway_id}|{pathway_name}")
                sample_summary["sample_ko_count_sum"] += sample_ko_count

            if pathway_fields:
                pathway_key = (sample_id, term, pathway_id, pathway_name)
                if pathway_key not in pathway_match_map:
                    pathway_match_map[pathway_key] = {
                        "search_term": term,
                        "matched_fields": ",".join(pathway_fields),
                        "sample_id": sample_id,
                        "pathway_id": pathway_id,
                        "pathway_name": pathway_name,
                        "annotated_pathway_ko_count": annotated_pathway_ko_count,
                        "observed_kos": observed_kos,
                        "total_kos": total_kos,
                        "completeness_pct": completeness_pct,
                    }

                    term_key = (sample_id, term, "pathway")
                    term_summary = per_sample_term.setdefault(
                        term_key,
                        {
                            "sample_id": sample_id,
                            "search_term": term,
                            "match_scope": "pathway",
                            "matched_rows": 0,
                            "unique_kos": set(),
                            "unique_ko_labels": {},
                            "unique_pathways": set(),
                            "matched_fields": set(),
                            "sample_ko_count_sum": 0.0,
                            "max_sample_ko_count": 0.0,
                        },
                    )
                    term_summary["matched_rows"] += 1
                    term_summary["unique_pathways"].add(f"{pathway_id}|{pathway_name}")
                    for field_name in pathway_fields:
                        term_summary["matched_fields"].add(field_name)
                    term_summary["sample_ko_count_sum"] += annotated_pathway_ko_count
                    term_summary["max_sample_ko_count"] = max(term_summary["max_sample_ko_count"], annotated_pathway_ko_count)

                    sample_summary["matched_pathway_terms"].add(term)
                    sample_summary["matched_pathway_rows"] += 1
                    if pathway_id:
                        sample_summary["unique_pathways"].add(f"{pathway_id}|{pathway_name}")

    ko_detail_fields = [
        "search_term",
        "matched_fields",
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
    ko_match_rows.sort(key=lambda row: (str(row["search_term"]).casefold(), str(row["sample_id"]), str(row["pathway_id"]), str(row["ko_id"])))
    write_tsv(output_dir / f"psychrophilic_summary_detailed{output_suffix}.tsv", ko_detail_fields, ko_match_rows)

    pathway_rows = sorted(
        pathway_match_map.values(),
        key=lambda row: (str(row["search_term"]).casefold(), str(row["sample_id"]), str(row["pathway_id"])),
    )
    write_tsv(
        output_dir / f"psychrophilic_pathway_matches{output_suffix}.tsv",
        [
            "search_term",
            "matched_fields",
            "sample_id",
            "pathway_id",
            "pathway_name",
            "annotated_pathway_ko_count",
            "observed_kos",
            "total_kos",
            "completeness_pct",
        ],
        pathway_rows,
    )

    sample_term_rows: list[dict[str, object]] = []
    for (_sample_id, _term, _scope), summary in sorted(per_sample_term.items(), key=lambda item: (item[0][0], item[0][2], item[0][1].casefold())):
        sample_term_rows.append(
            {
                "sample_id": summary["sample_id"],
                "search_term": summary["search_term"],
                "match_scope": summary["match_scope"],
                "matched_rows": summary["matched_rows"],
                "unique_ko_count": len(summary["unique_kos"]),
                "unique_kos": "; ".join(sorted(summary["unique_kos"])),
                "unique_ko_labels": "; ".join(
                    f"{ko_id} = {summary['unique_ko_labels'][ko_id]}"
                    for ko_id in sorted(summary["unique_ko_labels"])
                ),
                "unique_pathway_count": len(summary["unique_pathways"]),
                "unique_pathways": "; ".join(sorted(summary["unique_pathways"])),
                "matched_fields": ", ".join(sorted(summary["matched_fields"])),
                "sample_ko_count_sum": round(summary["sample_ko_count_sum"], 6),
                "max_sample_ko_count": round(summary["max_sample_ko_count"], 6),
            }
        )
    write_tsv(
        output_dir / f"psychrophilic_summary_by_sample_and_term{output_suffix}.tsv",
        [
            "sample_id",
            "search_term",
            "match_scope",
            "matched_rows",
            "unique_ko_count",
            "unique_kos",
            "unique_ko_labels",
            "unique_pathway_count",
            "unique_pathways",
            "matched_fields",
            "sample_ko_count_sum",
            "max_sample_ko_count",
        ],
        sample_term_rows,
    )

    overall_rows: list[dict[str, object]] = []
    for sample_id, summary in sorted(per_sample_overall.items()):
        overall_rows.append(
            {
                "sample_id": sample_id,
                "matched_ko_term_count": len(summary["matched_ko_terms"]),
                "matched_ko_terms": "; ".join(sorted(summary["matched_ko_terms"], key=str.casefold)),
                "matched_ko_rows": summary["matched_ko_rows"],
                "unique_ko_count": len(summary["unique_kos"]),
                "unique_kos": "; ".join(sorted(summary["unique_kos"])),
                "matched_pathway_term_count": len(summary["matched_pathway_terms"]),
                "matched_pathway_terms": "; ".join(sorted(summary["matched_pathway_terms"], key=str.casefold)),
                "matched_pathway_rows": summary["matched_pathway_rows"],
                "unique_pathway_count": len(summary["unique_pathways"]),
                "unique_pathways": "; ".join(sorted(summary["unique_pathways"])),
                "sample_ko_count_sum": round(summary["sample_ko_count_sum"], 6),
            }
        )
    write_tsv(
        output_dir / f"psychrophilic_summary_overall{output_suffix}.tsv",
        [
            "sample_id",
            "matched_ko_term_count",
            "matched_ko_terms",
            "matched_ko_rows",
            "unique_ko_count",
            "unique_kos",
            "matched_pathway_term_count",
            "matched_pathway_terms",
            "matched_pathway_rows",
            "unique_pathway_count",
            "unique_pathways",
            "sample_ko_count_sum",
        ],
        overall_rows,
    )

    write_tsv(
        output_dir / f"psychrophilic_search_terms_not_found{output_suffix}.tsv",
        ["search_term"],
        [{"search_term": term} for term in sorted(not_found_terms.values(), key=str.casefold)],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
