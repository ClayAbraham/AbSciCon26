from __future__ import annotations

import csv
import html
import math
from pathlib import Path


SAMPLE_ORDER = ["1A", "3A", "5A", "13A", "14A", "42A", "46A", "47A", "49A"]
MIN_MAX_SAMPLE_KO_COUNT = 10.0
TERM_ORDER_CANDIDATES = [
    "psychrophilic_heatmap_vIII_row_order.tsv",
    "psychrophilic_summary_by_sample_and_term_vIII_cutoff10.tsv",
]
DETAILED_INPUT_CANDIDATES = [
    "psychrophilic_summary_detailed_II.tsv",
    "psychrophilic_summary_detailed.tsv",
]
CATEGORY_ORDER = [
    "Membrane and Lipid Remodeling",
    "Osmoprotectants and Compatible Solutes",
    "Biofilm and Extracellular Matrix",
    "Carbon Storage",
    "Other Stress-Linked Genes",
]
SUBGROUP_ORDER = {
    "Membrane and Lipid Remodeling": [
        "Desaturases",
        "Fatty Acid Biosynthesis",
        "Other Membrane Remodeling",
    ],
    "Osmoprotectants and Compatible Solutes": [
        "Betaine and Osmoprotection",
        "Ectoine System",
        "Trehalose System",
        "Other Compatible Solutes",
    ],
    "Biofilm and Extracellular Matrix": [
        "Extracellular Polysaccharides",
        "Other Biofilm Functions",
    ],
    "Carbon Storage": [
        "Storage Polymers",
    ],
    "Other Stress-Linked Genes": [
        "Other Stress-Linked Genes",
    ],
}
TERM_CATEGORY_MAP = {
    "desaturase": "Membrane and Lipid Remodeling",
    "omega-6 desaturase": "Membrane and Lipid Remodeling",
    "fabA": "Membrane and Lipid Remodeling",
    "fabB": "Membrane and Lipid Remodeling",
    "fabD": "Membrane and Lipid Remodeling",
    "fabF": "Membrane and Lipid Remodeling",
    "fabG": "Membrane and Lipid Remodeling",
    "fabH": "Membrane and Lipid Remodeling",
    "osmoprotectant": "Osmoprotectants and Compatible Solutes",
    "betA": "Osmoprotectants and Compatible Solutes",
    "betB": "Osmoprotectants and Compatible Solutes",
    "betaine": "Osmoprotectants and Compatible Solutes",
    "glycine betaine": "Osmoprotectants and Compatible Solutes",
    "opuA": "Osmoprotectants and Compatible Solutes",
    "ectA": "Osmoprotectants and Compatible Solutes",
    "ectB": "Osmoprotectants and Compatible Solutes",
    "ectC": "Osmoprotectants and Compatible Solutes",
    "ectD": "Osmoprotectants and Compatible Solutes",
    "ectoine": "Osmoprotectants and Compatible Solutes",
    "otsA": "Osmoprotectants and Compatible Solutes",
    "otsB": "Osmoprotectants and Compatible Solutes",
    "trehalose": "Osmoprotectants and Compatible Solutes",
    "trehalose synthase": "Osmoprotectants and Compatible Solutes",
    "treS": "Osmoprotectants and Compatible Solutes",
    "treY": "Osmoprotectants and Compatible Solutes",
    "treZ": "Osmoprotectants and Compatible Solutes",
    "alginate": "Biofilm and Extracellular Matrix",
    "biofilm": "Biofilm and Extracellular Matrix",
    "capsular polysaccharide": "Biofilm and Extracellular Matrix",
    "cellulose synthase": "Biofilm and Extracellular Matrix",
    "exopolysaccharide": "Biofilm and Extracellular Matrix",
    "pel": "Biofilm and Extracellular Matrix",
    "carbon storage": "Carbon Storage",
    "glycogen": "Carbon Storage",
}
TERM_SUBGROUP_MAP = {
    "desaturase": "Desaturases",
    "omega-6 desaturase": "Desaturases",
    "fabA": "Fatty Acid Biosynthesis",
    "fabB": "Fatty Acid Biosynthesis",
    "fabD": "Fatty Acid Biosynthesis",
    "fabF": "Fatty Acid Biosynthesis",
    "fabG": "Fatty Acid Biosynthesis",
    "fabH": "Fatty Acid Biosynthesis",
    "osmoprotectant": "Betaine and Osmoprotection",
    "betA": "Betaine and Osmoprotection",
    "betB": "Betaine and Osmoprotection",
    "betaine": "Betaine and Osmoprotection",
    "glycine betaine": "Betaine and Osmoprotection",
    "opuA": "Betaine and Osmoprotection",
    "ectA": "Ectoine System",
    "ectB": "Ectoine System",
    "ectC": "Ectoine System",
    "ectD": "Ectoine System",
    "ectoine": "Ectoine System",
    "otsA": "Trehalose System",
    "otsB": "Trehalose System",
    "trehalose": "Trehalose System",
    "trehalose synthase": "Trehalose System",
    "treS": "Trehalose System",
    "treY": "Trehalose System",
    "treZ": "Trehalose System",
    "alginate": "Extracellular Polysaccharides",
    "capsular polysaccharide": "Extracellular Polysaccharides",
    "cellulose synthase": "Extracellular Polysaccharides",
    "exopolysaccharide": "Extracellular Polysaccharides",
    "pel": "Extracellular Polysaccharides",
    "biofilm": "Other Biofilm Functions",
    "carbon storage": "Storage Polymers",
    "glycogen": "Storage Polymers",
}
SUBGROUP_TERM_ORDER = {
    "Desaturases": ["desaturase", "omega-6 desaturase"],
    "Fatty Acid Biosynthesis": ["fabA", "fabB", "fabD", "fabF", "fabG", "fabH"],
    "Betaine and Osmoprotection": ["osmoprotectant", "betA", "betB", "betaine", "glycine betaine", "opuA"],
    "Ectoine System": ["ectA", "ectB", "ectC", "ectD", "ectoine"],
    "Trehalose System": ["otsA", "otsB", "trehalose", "trehalose synthase", "treS", "treY", "treZ"],
    "Extracellular Polysaccharides": ["alginate", "capsular polysaccharide", "cellulose synthase", "exopolysaccharide", "pel"],
    "Other Biofilm Functions": ["biofilm"],
    "Storage Polymers": ["carbon storage", "glycogen"],
}
CATEGORY_DISPLAY_LINES = {
    "Membrane and Lipid Remodeling": ["Membrane and Lipid", "Remodeling"],
    "Osmoprotectants and Compatible Solutes": ["Osmoprotectants", "and Compatible Solutes"],
    "Biofilm and Extracellular Matrix": ["Biofilm and", "Extracellular Matrix"],
    "Carbon Storage": ["Carbon", "Storage"],
    "Other Stress-Linked Genes": ["Other Stress-Linked", "Genes"],
}


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def first_existing_path(root: Path, candidates: list[str]) -> Path:
    for candidate_name in candidates:
        candidate = root / candidate_name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not find any of these files in {root}: {', '.join(candidates)}")


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


def term_category(term: str) -> str:
    return TERM_CATEGORY_MAP.get(term, "Other Stress-Linked Genes")


def term_subgroup(term: str, category: str) -> str:
    fallback_map = {
        "Membrane and Lipid Remodeling": "Other Membrane Remodeling",
        "Osmoprotectants and Compatible Solutes": "Other Compatible Solutes",
        "Biofilm and Extracellular Matrix": "Other Biofilm Functions",
        "Carbon Storage": "Storage Polymers",
        "Other Stress-Linked Genes": "Other Stress-Linked Genes",
    }
    return TERM_SUBGROUP_MAP.get(term, fallback_map.get(category, "Other Stress-Linked Genes"))


def subgroup_sorted_terms(category: str, terms: list[str], matrix_map: dict[str, dict[str, float]], sample_ids: list[str]) -> list[str]:
    grouped_terms: dict[str, list[str]] = {}
    for term in terms:
        grouped_terms.setdefault(term_subgroup(term, category), []).append(term)

    ordered_terms: list[str] = []
    seen_terms: set[str] = set()
    for subgroup in SUBGROUP_ORDER.get(category, []):
        subgroup_terms = grouped_terms.get(subgroup, [])
        if not subgroup_terms:
            continue
        manual_order = SUBGROUP_TERM_ORDER.get(subgroup, [])
        manual_rank = {term: index for index, term in enumerate(manual_order)}
        known_terms = [term for term in subgroup_terms if term in manual_rank]
        unknown_terms = [term for term in subgroup_terms if term not in manual_rank]
        ordered_terms.extend(sorted(known_terms, key=lambda term: manual_rank[term]))
        if len(unknown_terms) > 1:
            vectors = [
                [math.log10(1.0 + matrix_map[term].get(sample_id, 0.0)) for sample_id in sample_ids]
                for term in unknown_terms
            ]
            order = cluster_order(vectors)
            ordered_terms.extend(unknown_terms[index] for index in order)
        else:
            ordered_terms.extend(unknown_terms)
        seen_terms.update(subgroup_terms)

    leftover_terms = [term for term in terms if term not in seen_terms]
    if len(leftover_terms) > 1:
        vectors = [
            [math.log10(1.0 + matrix_map[term].get(sample_id, 0.0)) for sample_id in leftover_terms]
            for term in leftover_terms
        ]
        order = cluster_order(vectors)
        ordered_terms.extend(leftover_terms[index] for index in order)
    else:
        ordered_terms.extend(leftover_terms)
    return ordered_terms


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


def looks_like_gene_symbol(term: str) -> bool:
    if " " in term:
        return False
    return any(character.isupper() for character in term) or any(character.isdigit() for character in term)


def first_ko_annotation(unique_ko_labels: str) -> str:
    if not unique_ko_labels:
        return ""
    first_label = unique_ko_labels.split("; K", 1)[0].strip()
    if first_label.startswith("K") and " = " in first_label:
        return first_label.split(" = ", 1)[1].strip()
    return first_label


def ko_annotation_entries(unique_ko_labels: str) -> list[str]:
    if not unique_ko_labels:
        return []
    entries: list[str] = []
    for index, chunk in enumerate(unique_ko_labels.split("; K")):
        entry = chunk.strip()
        if not entry:
            continue
        if index > 0:
            entry = "K" + entry
        entries.append(entry)
    return entries


def matching_ko_annotation(term: str, unique_ko_labels: str) -> str:
    term_lower = term.lower()
    for entry in ko_annotation_entries(unique_ko_labels):
        annotation = entry.split(" = ", 1)[1].strip() if entry.startswith("K") and " = " in entry else entry
        symbol_text = annotation.split(";", 1)[0]
        symbol_tokens = {token.strip().lower() for token in symbol_text.split(",")}
        if term_lower in symbol_tokens:
            return annotation
    return ""


def annotation_display_name(annotation: str) -> str:
    if not annotation:
        return ""
    annotation = annotation.strip().strip('"')
    if ";" in annotation:
        _symbol_text, enzyme_text = annotation.split(";", 1)
        annotation = enzyme_text.strip()
    if " [EC:" in annotation:
        annotation = annotation.split(" [EC:", 1)[0].rstrip()
    return " ".join(annotation.split())


def clean_enzyme_name(term: str, annotation: str) -> str:
    if not annotation:
        return ""
    annotation = annotation.strip().strip('"')
    if ";" in annotation:
        symbol_text, enzyme_text = annotation.split(";", 1)
        symbol_tokens = {token.strip().lower() for token in symbol_text.split(",")}
        if term.lower() in symbol_tokens:
            annotation = enzyme_text.strip()
    return annotation_display_name(annotation)


def compact_chart_label_name(annotation: str) -> str:
    if not annotation:
        return ""
    compact = annotation
    replacements = [
        ("3-hydroxyacyl-[acyl-carrier protein]", "3-hydroxyacyl-ACP"),
        ("3-hydroxyacyl-[acyl-carrier-protein]", "3-hydroxyacyl-ACP"),
        ("3-oxoacyl-[acyl-carrier protein]", "3-oxoacyl-ACP"),
        ("3-oxoacyl-[acyl-carrier-protein]", "3-oxoacyl-ACP"),
        ("acyl-[acyl-carrier protein]", "acyl-ACP"),
        ("acyl-[acyl-carrier-protein]", "acyl-ACP"),
        ("[acyl-carrier protein]", "ACP"),
        ("[acyl-carrier-protein]", "ACP"),
    ]
    for old, new in replacements:
        compact = compact.replace(old, new)
    return compact


def summarize_ko_annotations(unique_ko_labels: str, max_names: int = 2) -> str:
    names: list[str] = []
    seen: set[str] = set()
    for entry in ko_annotation_entries(unique_ko_labels):
        annotation = entry.split(" = ", 1)[1].strip() if entry.startswith("K") and " = " in entry else entry
        name = annotation_display_name(annotation)
        key = name.lower()
        if not name or key in seen:
            continue
        seen.add(key)
        names.append(name)
    if not names:
        return ""
    if len(names) <= max_names:
        return "; ".join(names)
    return f"{'; '.join(names[:max_names])} (+{len(names) - max_names} more)"


def build_display_label(term: str, metadata: dict[str, str]) -> str:
    unique_ko_labels = metadata.get("unique_ko_labels", "")
    if looks_like_gene_symbol(term):
        annotation = matching_ko_annotation(term, unique_ko_labels)
        if not annotation:
            unique_ko_count = int(metadata.get("unique_ko_count", "0") or 0)
            if unique_ko_count == 1:
                annotation = first_ko_annotation(unique_ko_labels)
        annotation = clean_enzyme_name(term, annotation)
        if annotation:
            return f"{term} | {annotation}"
    summary = summarize_ko_annotations(unique_ko_labels)
    if summary:
        return f"{term} | {summary}"
    return term


def build_ko_display_label(term: str, metadata: dict[str, str]) -> str:
    ko_id = metadata.get("ko_id", "").strip()
    ko_symbol = metadata.get("ko_symbol", "").strip()
    ko_name = metadata.get("ko_name", "").strip()
    symbol_parts = [part.strip() for part in ko_symbol.split(",") if part.strip()]
    primary_symbol = symbol_parts[0] if symbol_parts else ""
    enzyme_name = annotation_display_name(f"{ko_symbol}; {ko_name}" if ko_symbol and ko_name else (ko_name or metadata.get("ko_label", "")))
    if looks_like_gene_symbol(term):
        if ko_id and enzyme_name:
            return f"{term} | {ko_id} | {enzyme_name}"
        if enzyme_name:
            return f"{term} | {enzyme_name}"
        if ko_id:
            return f"{term} | {ko_id}"
        return term
    ko_part = " ".join(part for part in [ko_id, primary_symbol] if part).strip()
    if ko_part and enzyme_name:
        return f"{term} | {ko_part} | {enzyme_name}"
    if enzyme_name:
        return f"{term} | {enzyme_name}"
    if ko_part:
        return f"{term} | {ko_part}"
    return term


def subgroup_category(subgroup: str) -> str:
    for category, subgroups in SUBGROUP_ORDER.items():
        if subgroup in subgroups:
            return category
    return "Other Stress-Linked Genes"


def primary_symbol(ko_symbol: str) -> str:
    parts = [part.strip() for part in ko_symbol.split(",") if part.strip()]
    return parts[0] if parts else ""


def resolve_ko_subgroup(metadata: dict[str, object]) -> str:
    symbol = primary_symbol(str(metadata.get("ko_symbol", ""))).casefold()
    name = str(metadata.get("ko_name", "")).casefold()
    source_terms = {str(term).casefold() for term in metadata.get("source_terms", set())}

    if "desaturase" in name or symbol.startswith("des") or symbol in {"scd", "fab2", "fad6"}:
        return "Desaturases"
    if symbol.startswith("fab") or symbol in {"mcat", "mct1", "oar1", "oxsm", "cem1"}:
        return "Fatty Acid Biosynthesis"
    if symbol in {"opua", "opubd", "opuc", "prov", "prow", "prox"} or "osmoprotectant transport" in name or "glycine betaine/proline transport" in name:
        return "Betaine and Osmoprotection"
    if symbol in {"ecta", "ectb", "ectc", "ectd", "doea"} or "ectoine" in name or "diaminobutyrate" in name:
        return "Ectoine System"
    if symbol in {"otsa", "otsb", "tres", "trey", "trez"} or "trehalose" in name:
        return "Trehalose System"
    if symbol.startswith("alg") or symbol.startswith("exo") or symbol.startswith("pga") or symbol.startswith("kps") or symbol == "bcsa" or symbol.startswith("pel"):
        return "Extracellular Polysaccharides"
    if "biofilm" in name:
        return "Other Biofilm Functions"
    if symbol in {"csra", "glgp"} or "carbon storage" in name or "glycogen" in name:
        return "Storage Polymers"

    for term in source_terms:
        if term in {"osmoprotectant", "betaine", "glycine betaine", "opua"}:
            return "Betaine and Osmoprotection"
        if term in {"ecta", "ectb", "ectc", "ectd", "ectoine"}:
            return "Ectoine System"
        if term in {"otsa", "otsb", "trehalose", "trehalose synthase", "tres", "trey", "trez"}:
            return "Trehalose System"
        if term in {"alginate", "exopolysaccharide", "capsular polysaccharide", "cellulose synthase", "pel"}:
            return "Extracellular Polysaccharides"
        if term == "biofilm":
            return "Other Biofilm Functions"
        if term in {"carbon storage", "glycogen"}:
            return "Storage Polymers"
        if term in {"desaturase", "omega-6 desaturase"}:
            return "Desaturases"
        if term in {"faba", "fabb", "fabd", "fabf", "fabg", "fabh"}:
            return "Fatty Acid Biosynthesis"
    return "Other Stress-Linked Genes"


def build_ko_resolved_label(metadata: dict[str, object]) -> str:
    ko_id = str(metadata.get("ko_id", "")).strip()
    ko_symbol = str(metadata.get("ko_symbol", "")).strip()
    ko_name = annotation_display_name(str(metadata.get("ko_name", "")).strip() or str(metadata.get("ko_label", "")).strip())
    ko_name = compact_chart_label_name(ko_name)
    symbol = primary_symbol(ko_symbol)
    left = " ".join(part for part in [ko_id, symbol] if part).strip()
    if left and ko_name:
        return f"{left} | {ko_name}"
    return left or ko_name or ko_id or "Unlabeled KO"


def multiline_rotated_text_svg(x: float, y: float, lines: list[str], line_spacing: float = 16.0, **attrs: str | float) -> str:
    attr_text = " ".join(f'{key}="{html.escape(str(value), quote=True)}"' for key, value in attrs.items())
    initial_dy = 0.0 if len(lines) == 1 else -line_spacing * (len(lines) - 1) / 2.0
    tspans: list[str] = []
    for index, line in enumerate(lines):
        dy_value = initial_dy if index == 0 else line_spacing
        tspans.append(f'<tspan x="{x:.2f}" dy="{dy_value:.2f}">{html.escape(line)}</tspan>')
    tspans_text = "".join(tspans)
    return f'<text x="{x:.2f}" y="{y:.2f}" {attr_text}>{tspans_text}</text>'


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
    output_dir = root / "output"

    term_order_path = first_existing_path(output_dir, TERM_ORDER_CANDIDATES)
    detailed_input_path = first_existing_path(output_dir, DETAILED_INPUT_CANDIDATES)

    term_order_rows = load_rows(term_order_path)
    detailed_rows = load_rows(detailed_input_path)

    ordered_terms_with_category: list[tuple[str, str]] = []
    seen_terms: set[str] = set()
    for row in term_order_rows:
        term = row.get("search_term", "").strip()
        if not term or term in seen_terms:
            continue
        seen_terms.add(term)
        category = row.get("category", "").strip() or term_category(term)
        ordered_terms_with_category.append((term, category))

    ordered_terms = [term for term, _category in ordered_terms_with_category]
    term_set = set(ordered_terms)
    term_category_lookup = {term: category for term, category in ordered_terms_with_category}
    term_rank_lookup = {term: index + 1 for index, (term, _category) in enumerate(ordered_terms_with_category)}

    matrix_map: dict[str, dict[str, float]] = {}
    row_metadata: dict[str, dict[str, object]] = {}
    for row in detailed_rows:
        sample_id = row.get("sample_id", "").strip()
        search_term = row.get("search_term", "").strip()
        ko_id = row.get("ko_id", "").strip()
        if search_term not in term_set or not sample_id or not ko_id:
            continue
        value = float(row.get("sample_ko_count", "0") or 0.0)
        sample_map = matrix_map.setdefault(ko_id, {})
        sample_map[sample_id] = max(sample_map.get(sample_id, 0.0), value)
        metadata = row_metadata.setdefault(
            ko_id,
            {
                "ko_id": ko_id,
                "ko_symbol": row.get("ko_symbol", "").strip(),
                "ko_name": row.get("ko_name", "").strip(),
                "ko_label": row.get("ko_label", "").strip(),
                "matched_fields": set(),
                "source_terms": set(),
                "source_term_ranks": set(),
            },
        )
        if not metadata["ko_symbol"] and row.get("ko_symbol", "").strip():
            metadata["ko_symbol"] = row.get("ko_symbol", "").strip()
        if not metadata["ko_name"] and row.get("ko_name", "").strip():
            metadata["ko_name"] = row.get("ko_name", "").strip()
        if not metadata["ko_label"] and row.get("ko_label", "").strip():
            metadata["ko_label"] = row.get("ko_label", "").strip()
        metadata["source_terms"].add(search_term)
        metadata["source_term_ranks"].add(term_rank_lookup.get(search_term, 9999))
        if row.get("matched_fields", "").strip():
            metadata["matched_fields"].add(row.get("matched_fields", "").strip())

    filtered_ko_ids = [
        ko_id
        for ko_id, sample_map in matrix_map.items()
        if max(sample_map.values(), default=0.0) >= MIN_MAX_SAMPLE_KO_COUNT
    ]
    sample_ids = [
        sample_id
        for sample_id in SAMPLE_ORDER
        if any(sample_id in matrix_map[ko_id] for ko_id in filtered_ko_ids)
    ]

    grouped_ko_ids: dict[str, dict[str, list[str]]] = {
        category: {subgroup: [] for subgroup in SUBGROUP_ORDER.get(category, [])}
        for category in CATEGORY_ORDER
    }
    for ko_id in filtered_ko_ids:
        subgroup = resolve_ko_subgroup(row_metadata[ko_id])
        category = subgroup_category(subgroup)
        row_metadata[ko_id]["subgroup"] = subgroup
        row_metadata[ko_id]["category"] = category
        grouped_ko_ids.setdefault(category, {})
        grouped_ko_ids[category].setdefault(subgroup, []).append(ko_id)

    ordered_entries: list[tuple[str, str, str]] = []
    for category in CATEGORY_ORDER:
        category_groups = grouped_ko_ids.get(category, {})
        for subgroup in SUBGROUP_ORDER.get(category, []):
            ko_ids = category_groups.get(subgroup, [])
            if not ko_ids:
                continue
            ko_ids = sorted(
                ko_ids,
                key=lambda ko_id: (
                    min(row_metadata[ko_id]["source_term_ranks"]) if row_metadata[ko_id]["source_term_ranks"] else 9999,
                    primary_symbol(str(row_metadata[ko_id].get("ko_symbol", ""))).casefold(),
                    ko_id,
                ),
            )
            vectors = [
                [math.log10(1.0 + matrix_map[ko_id].get(sample_id, 0.0)) for sample_id in sample_ids]
                for ko_id in ko_ids
            ]
            ko_order = cluster_order(vectors) if len(ko_ids) > 1 else list(range(len(ko_ids)))
            ordered_entries.extend((ko_ids[index], category, subgroup) for index in ko_order)

    ordered_keys = [ko_id for ko_id, _category, _subgroup in ordered_entries]
    ordered_matrix = [
        [math.log10(1.0 + matrix_map[ko_id].get(sample_id, 0.0)) for sample_id in sample_ids]
        for ko_id in ordered_keys
    ]

    matrix_rows: list[dict[str, object]] = []
    for ko_id in ordered_keys:
        metadata = row_metadata[ko_id]
        row: dict[str, object] = {
            "search_terms": "; ".join(sorted(metadata["source_terms"])),
            "ko_id": metadata["ko_id"],
            "ko_symbol": metadata["ko_symbol"],
            "ko_name": metadata["ko_name"],
            "display_label": build_ko_resolved_label(metadata),
            "category": metadata["category"],
            "subgroup": metadata["subgroup"],
            "max_sample_ko_count": round(max(matrix_map[ko_id].values(), default=0.0), 6),
        }
        for sample_id in sample_ids:
            row[sample_id] = round(matrix_map[ko_id].get(sample_id, 0.0), 6)
        matrix_rows.append(row)

    write_tsv(
        output_dir / "psychrophilic_heatmap_vIV_matrix.tsv",
        ["search_terms", "ko_id", "ko_symbol", "ko_name", "display_label", "category", "subgroup", "max_sample_ko_count", *sample_ids],
        matrix_rows,
    )

    write_tsv(
        output_dir / "psychrophilic_heatmap_vIV_row_order.tsv",
        ["ko_id", "ko_symbol", "ko_name", "display_label", "category", "subgroup", "source_terms", "source_term_ranks", "rank", "max_sample_ko_count"],
        [
            {
                "ko_id": row_metadata[ko_id]["ko_id"],
                "ko_symbol": row_metadata[ko_id]["ko_symbol"],
                "ko_name": row_metadata[ko_id]["ko_name"],
                "display_label": build_ko_resolved_label(row_metadata[ko_id]),
                "category": category,
                "subgroup": subgroup,
                "source_terms": "; ".join(sorted(row_metadata[ko_id]["source_terms"])),
                "source_term_ranks": "; ".join(str(rank) for rank in sorted(row_metadata[ko_id]["source_term_ranks"])),
                "rank": index + 1,
                "max_sample_ko_count": round(max(matrix_map[ko_id].values(), default=0.0), 6),
            }
            for index, (ko_id, category, subgroup) in enumerate(ordered_entries)
        ],
    )

    max_log = max((max(row) for row in ordered_matrix), default=1.0)
    cell_width = 72
    cell_height = 40
    group_gap = 22
    left_margin = 980
    top_margin = 180
    legend_width = 36
    legend_height = 320
    heatmap_right_x = left_margin + len(sample_ids) * cell_width
    label_x = heatmap_right_x + 18
    unique_categories = []
    for _ko_id, category, _subgroup in ordered_entries:
        if not unique_categories or unique_categories[-1] != category:
            unique_categories.append(category)
    total_group_gap = max(0, len(unique_categories) - 1) * group_gap
    height = 240 + len(ordered_keys) * cell_height + total_group_gap + 150
    legend_x = heatmap_right_x + 430
    legend_y = top_margin
    width = int(legend_x + 280)

    row_positions: list[dict[str, object]] = []
    current_y = top_margin
    previous_category = None
    for ko_id, category, subgroup in ordered_entries:
        if previous_category is not None and category != previous_category:
            current_y += group_gap
        row_positions.append({"ko_id": ko_id, "category": category, "subgroup": subgroup, "y": current_y})
        current_y += cell_height
        previous_category = category

    category_segments: list[dict[str, float | str]] = []
    for row_info in row_positions:
        category = str(row_info["category"])
        y = float(row_info["y"])
        if not category_segments or category_segments[-1]["category"] != category:
            category_segments.append({"category": category, "start_y": y, "end_y": y + cell_height})
        else:
            category_segments[-1]["end_y"] = y + cell_height

    parts = [svg_header(width, height)]
    parts.append(rect_svg(0, 0, width, height, fill="white"))
    parts.append(text_svg(width / 2, 52, "Psychrophilic Gene Summary Heatmap vIV", **{"text-anchor": "middle", "font-size": 30, "font-weight": "700", "font-family": "Arial"}))
    parts.append(
        text_svg(
            width / 2,
            86,
            "KO-resolved rows from the vIII psychrophilic term set; samples fixed in genomic order; KOs grouped by biological response family and annotated by function",
            **{"text-anchor": "middle", "font-size": 17, "fill": "#555", "font-family": "Arial"},
        )
    )

    for sample_index, sample_id in enumerate(sample_ids):
        x = left_margin + sample_index * cell_width + cell_width / 2
        y = top_margin - 22
        parts.append(f'<text x="{x:.2f}" y="{y:.2f}" transform="rotate(-45 {x:.2f} {y:.2f})" text-anchor="end" font-size="16" font-weight="700" font-family="Arial">{html.escape(sample_id)}</text>')

    category_base_x = left_margin - 52
    for segment_index, segment in enumerate(category_segments):
        category = str(segment["category"])
        category_lines = CATEGORY_DISPLAY_LINES.get(category, wrap_label(category, width=20, max_lines=2))
        category_center_y = (float(segment["start_y"]) + float(segment["end_y"])) / 2.0
        category_x = category_base_x
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
    for index, row_info in enumerate(row_positions):
        ko_id = str(row_info["ko_id"])
        category = str(row_info["category"])
        y = float(row_info["y"])
        if previous_category != category:
            if previous_category is not None:
                separator_y = y - (group_gap / 2)
                parts.append(line_svg(left_margin - 18, separator_y, heatmap_right_x, separator_y, stroke="#b8b8b8", **{"stroke-width": 1.5}))

        display_label = build_ko_resolved_label(row_metadata[ko_id])
        wrapped = wrap_label(display_label, width=44, max_lines=4)
        center_y = y + cell_height / 2 + 5
        stagger_offset = -3 if (index % 2 == 0) else 3
        label_y = center_y + stagger_offset if len(wrapped) == 1 else center_y - 8.5 * (len(wrapped) - 1) + stagger_offset
        for line_index, line in enumerate(wrapped):
            parts.append(
                text_svg(
                    label_x,
                    label_y + line_index * 17,
                    line,
                    **{"text-anchor": "start", "font-size": 13, "font-weight": "600", "font-family": "Arial", "fill": "#111"},
                )
            )
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

    parts.append(text_svg(legend_x + legend_width / 2, legend_y - 42, "KO Count", **{"text-anchor": "middle", "font-size": 15, "font-weight": "700", "font-family": "Arial"}))
    parts.append(text_svg(legend_x + legend_width / 2, legend_y - 24, "log10(1 + value)", **{"text-anchor": "middle", "font-size": 13, "font-family": "Arial", "fill": "#444"}))
    parts.append(text_svg(legend_x + legend_width / 2, legend_y - 8, "per-KO max >= 10", **{"text-anchor": "middle", "font-size": 12, "font-family": "Arial", "fill": "#666"}))

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
    parts.append(rect_svg(legend_x, legend_y, legend_width, legend_height, fill="none", stroke="#666", **{"stroke-width": 0.8}))

    parts.append("</svg>")
    svg_path = output_dir / "psychrophilic_heatmap_vIV.svg"
    write_text(svg_path, "\n".join(parts))

    html_doc = (
        "<!DOCTYPE html><html><head><meta charset='utf-8'><title>Psychrophilic Gene Summary Heatmap vIV</title>"
        "<style>body{font-family:Arial,sans-serif;margin:32px;line-height:1.4}img{max-width:100%;border:1px solid #ddd}</style>"
        "</head><body><h1>Psychrophilic Gene Summary Heatmap vIV</h1>"
        f"<p>Source tables: {html.escape(term_order_path.name)} and {html.escape(detailed_input_path.name)}. Rows are unique KO-resolved entries within the psychrophilic term set, deduplicated across overlapping search terms and grouped by biological function. Each KO is required to have max_sample_ko_count &gt;= 10 before plotting. Values shown are log10(1 + sample_ko_count).</p>"
        f"<img src='{svg_path.name}' alt='Psychrophilic gene summary heatmap vIV' />"
        "</body></html>"
    )
    write_text(output_dir / "psychrophilic_heatmap_vIV_index.html", html_doc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
