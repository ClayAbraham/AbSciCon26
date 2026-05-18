from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


TARGET_PATHWAYS = {
    "map00650": "Butanoate metabolism",
    "map00640": "Propanoate metabolism",
    "map00680": "Methane metabolism",
    "map00710": "Carbon fixation by Calvin cycle",
    "map00633": "Nitrotoluene degradation",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export all observed KOs from the selected target pathways across all samples."
    )
    parser.add_argument(
        "--input",
        default="output/all_samples_observed_pathway_kos.tsv",
        help="Combined observed pathway KO table.",
    )
    parser.add_argument(
        "--output-dir",
        default="output",
        help="Directory for TSV outputs.",
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


def main() -> int:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    input_path = (script_dir / args.input).resolve()
    output_dir = (script_dir / args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    all_rows = load_tsv(input_path)
    target_rows = [row for row in all_rows if row["pathway_id"] in TARGET_PATHWAYS]

    sample_ids = sorted({row["sample_id"] for row in target_rows}, key=natural_sample_key)
    pathway_ids = list(TARGET_PATHWAYS)

    target_rows.sort(
        key=lambda row: (
            natural_sample_key(row["sample_id"]),
            row["pathway_id"],
            -float(row["sample_ko_count"]),
            row["ko_id"],
        )
    )

    ko_labels: dict[str, str] = {}
    ko_pathways: dict[str, set[str]] = defaultdict(set)
    ko_samples: dict[str, set[str]] = defaultdict(set)
    ko_counts_by_sample: dict[str, dict[str, float]] = defaultdict(dict)
    ko_pathway_sample_hits: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {pathway_id: set() for pathway_id in pathway_ids}
    )
    pathway_sample_unique_kos: dict[tuple[str, str], set[str]] = defaultdict(set)
    pathway_sample_total_counts: dict[tuple[str, str], float] = defaultdict(float)

    for row in target_rows:
        sample_id = row["sample_id"]
        pathway_id = row["pathway_id"]
        ko_id = row["ko_id"]
        ko_count = float(row["sample_ko_count"])

        ko_labels[ko_id] = row["ko_label"]
        ko_pathways[ko_id].add(pathway_id)
        ko_samples[ko_id].add(sample_id)
        ko_counts_by_sample[ko_id][sample_id] = max(ko_count, ko_counts_by_sample[ko_id].get(sample_id, 0.0))
        ko_pathway_sample_hits[ko_id][pathway_id].add(sample_id)
        pathway_sample_unique_kos[(sample_id, pathway_id)].add(ko_id)
        pathway_sample_total_counts[(sample_id, pathway_id)] += ko_count

    ko_summary_rows: list[dict[str, object]] = []
    for ko_id in sorted(
        ko_pathways,
        key=lambda value: (
            -len(ko_pathways[value]),
            -len(ko_samples[value]),
            value,
        ),
    ):
        row: dict[str, object] = {
            "ko_id": ko_id,
            "ko_label": ko_labels.get(ko_id, ko_id),
            "pathway_count": len(ko_pathways[ko_id]),
            "pathway_ids": "; ".join(sorted(ko_pathways[ko_id])),
            "pathway_names": "; ".join(TARGET_PATHWAYS[pathway_id] for pathway_id in sorted(ko_pathways[ko_id])),
            "sample_count": len(ko_samples[ko_id]),
            "sample_ids": "; ".join(sorted(ko_samples[ko_id], key=natural_sample_key)),
            "present_in_all_9_samples": "yes" if len(ko_samples[ko_id]) == len(sample_ids) else "no",
        }
        for pathway_id in pathway_ids:
            row[f"{pathway_id}_samples_detected"] = len(ko_pathway_sample_hits[ko_id][pathway_id])
        for sample_id in sample_ids:
            row[sample_id] = round(ko_counts_by_sample[ko_id].get(sample_id, 0.0), 6)
        ko_summary_rows.append(row)

    pathway_sample_summary_rows: list[dict[str, object]] = []
    for sample_id in sample_ids:
        for pathway_id in pathway_ids:
            unique_ko_ids = sorted(pathway_sample_unique_kos[(sample_id, pathway_id)])
            pathway_sample_summary_rows.append(
                {
                    "sample_id": sample_id,
                    "pathway_id": pathway_id,
                    "pathway_name": TARGET_PATHWAYS[pathway_id],
                    "unique_ko_count": len(unique_ko_ids),
                    "ko_ids": "; ".join(unique_ko_ids),
                    "total_ko_count_sum": round(pathway_sample_total_counts[(sample_id, pathway_id)], 6),
                }
            )

    pathway_pair_rows: list[dict[str, object]] = []
    pathway_ko_sets = {
        pathway_id: {row["ko_id"] for row in target_rows if row["pathway_id"] == pathway_id}
        for pathway_id in pathway_ids
    }
    for index, pathway_a in enumerate(pathway_ids):
        for pathway_b in pathway_ids[index + 1 :]:
            overlap = pathway_ko_sets[pathway_a] & pathway_ko_sets[pathway_b]
            pathway_pair_rows.append(
                {
                    "pathway_a_id": pathway_a,
                    "pathway_a_name": TARGET_PATHWAYS[pathway_a],
                    "pathway_b_id": pathway_b,
                    "pathway_b_name": TARGET_PATHWAYS[pathway_b],
                    "shared_ko_count": len(overlap),
                    "shared_ko_ids": "; ".join(sorted(overlap)),
                }
            )

    target_detail_fields = [
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
    ko_summary_fields = [
        "ko_id",
        "ko_label",
        "pathway_count",
        "pathway_ids",
        "pathway_names",
        "sample_count",
        "sample_ids",
        "present_in_all_9_samples",
    ]
    ko_summary_fields.extend(f"{pathway_id}_samples_detected" for pathway_id in pathway_ids)
    ko_summary_fields.extend(sample_ids)

    write_tsv(output_dir / "target_pathway_all_kos.tsv", target_detail_fields, target_rows)
    write_tsv(output_dir / "target_pathway_ko_summary.tsv", ko_summary_fields, ko_summary_rows)
    write_tsv(
        output_dir / "target_pathway_pathway_sample_summary.tsv",
        ["sample_id", "pathway_id", "pathway_name", "unique_ko_count", "ko_ids", "total_ko_count_sum"],
        pathway_sample_summary_rows,
    )
    write_tsv(
        output_dir / "target_pathway_pairwise_overlap_summary.tsv",
        [
            "pathway_a_id",
            "pathway_a_name",
            "pathway_b_id",
            "pathway_b_name",
            "shared_ko_count",
            "shared_ko_ids",
        ],
        pathway_pair_rows,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
