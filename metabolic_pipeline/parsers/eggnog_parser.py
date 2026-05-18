"""eggNOG-mapper parser for precomputed annotations."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def _find_first_column(columns: list[str], candidates: list[str]) -> str | None:
    lowered = {c.lower(): c for c in columns}
    for name in candidates:
        if name.lower() in lowered:
            return lowered[name.lower()]
    return None


def parse_eggnog_file(path: Path, sample_id: str) -> pd.DataFrame:
    """Parse one eggNOG annotation file."""
    if not path.exists():
        raise FileNotFoundError(f"eggNOG file not found: {path}")
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    if df.empty:
        return pd.DataFrame(
            columns=[
                "sample_id",
                "protein_id",
                "eggnog_og",
                "ec_number",
                "preferred_name",
                "product",
                "ko_id",
                "evalue",
                "score",
                "source_tool",
            ]
        )

    columns = list(df.columns)
    protein_col = _find_first_column(columns, ["query", "#query", "query_name", "protein_id"])
    if protein_col is None:
        protein_col = columns[0]
    ko_col = _find_first_column(columns, ["KEGG_ko", "KEGG_KO", "ko_id", "KEGG_ko(s)"])
    og_col = _find_first_column(columns, ["eggNOG_OGs", "OGs"])
    ec_col = _find_first_column(columns, ["EC", "ec", "ec_number"])
    pref_col = _find_first_column(columns, ["Preferred_name", "preferred_name"])
    prod_col = _find_first_column(columns, ["Description", "product", "description"])
    eval_col = _find_first_column(columns, ["evalue", "e-value", "seed_ortholog_evalue"])
    score_col = _find_first_column(columns, ["score", "bitscore", "seed_ortholog_score"])

    out = pd.DataFrame(
        {
            "sample_id": sample_id,
            "protein_id": df[protein_col].astype(str),
            "eggnog_og": df[og_col].astype(str) if og_col else None,
            "ec_number": df[ec_col].astype(str) if ec_col else None,
            "preferred_name": df[pref_col].astype(str) if pref_col else None,
            "product": df[prod_col].astype(str) if prod_col else None,
            "ko_id": df[ko_col].astype(str) if ko_col else None,
            "evalue": df[eval_col].astype(str) if eval_col else None,
            "score": df[score_col].astype(str) if score_col else None,
            "source_tool": "eggnog",
        }
    )
    return out

