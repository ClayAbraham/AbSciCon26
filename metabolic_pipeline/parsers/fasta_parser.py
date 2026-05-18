"""FASTA parsing utilities."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


def _extract_protein_id(header_line: str) -> str:
    header = header_line.lstrip(">").strip()
    if not header:
        raise ValueError("Encountered malformed FASTA header with no identifier.")
    return header.split()[0]


def _extract_contig_id(header_line: str, protein_id: str) -> str | None:
    header = header_line.lstrip(">").strip()
    contig_match = re.search(r"(?:contig|contig_id)[=:]([^\s;|]+)", header, flags=re.IGNORECASE)
    if contig_match:
        return contig_match.group(1)
    if "_" in protein_id:
        return protein_id.rsplit("_", 1)[0]
    return None


def parse_faa_file(faa_path: Path, sample_id: str) -> pd.DataFrame:
    """Parse one protein FASTA file into metadata rows."""
    if not faa_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {faa_path}")

    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    current_header: str | None = None
    seq_chunks: list[str] = []

    def finalize_record() -> None:
        nonlocal current_header, seq_chunks
        if current_header is None:
            return
        protein_id = _extract_protein_id(current_header)
        if protein_id in seen:
            raise ValueError(
                f"Duplicate protein ID detected within sample '{sample_id}': {protein_id}"
            )
        seen.add(protein_id)
        sequence = "".join(seq_chunks).strip()
        rows.append(
            {
                "sample_id": sample_id,
                "protein_id": protein_id,
                "contig_id": _extract_contig_id(current_header, protein_id),
                "raw_header": current_header.lstrip(">").strip(),
                "sequence_length": len(sequence),
            }
        )
        current_header = None
        seq_chunks = []

    with faa_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                finalize_record()
                current_header = line
            else:
                if current_header is None:
                    raise ValueError(
                        f"Malformed FASTA in {faa_path}: sequence encountered before header."
                    )
                seq_chunks.append(line.strip())
    finalize_record()

    if not rows:
        raise ValueError(f"Empty FASTA file: {faa_path}")
    return pd.DataFrame(rows)


def parse_faa_files(faa_paths: list[Path]) -> pd.DataFrame:
    """Parse many FASTA files and concatenate rows."""
    frames: list[pd.DataFrame] = []
    for path in faa_paths:
        frames.append(parse_faa_file(path, sample_id=path.stem))
    return pd.concat(frames, ignore_index=True)

