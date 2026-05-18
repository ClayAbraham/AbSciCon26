"""Parsers for Prodigal outputs."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from metabolic_pipeline.parsers.fasta_parser import parse_faa_file

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")


def parse_prodigal_gff(path: Path, sample_id: str) -> pd.DataFrame:
    """Parse a Prodigal GFF file into protein-to-contig mappings."""
    if not path.exists():
        raise FileNotFoundError(f"Prodigal GFF not found: {path}")

    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                continue
            contig_id = fields[0]
            attributes = fields[8]
            match = ID_RE.search(attributes)
            if match is None:
                continue
            rows.append(
                {
                    "sample_id": sample_id,
                    "protein_id": match.group(1),
                    "contig_id": contig_id,
                }
            )

    if not rows:
        return pd.DataFrame(columns=["sample_id", "protein_id", "contig_id"])
    return pd.DataFrame(rows).drop_duplicates(subset=["sample_id", "protein_id"], keep="first")


def parse_prodigal_outputs(protein_faa: Path, gff_path: Path, sample_id: str) -> pd.DataFrame:
    """Parse Prodigal amino acid FASTA and overlay contig IDs from the GFF."""
    proteins = parse_faa_file(protein_faa, sample_id=sample_id)
    gff_map = parse_prodigal_gff(gff_path, sample_id=sample_id)
    if gff_map.empty:
        return proteins

    merged = proteins.drop(columns=["contig_id"]).merge(
        gff_map,
        on=["sample_id", "protein_id"],
        how="left",
    )
    return merged[["sample_id", "protein_id", "contig_id", "raw_header", "sequence_length"]]
