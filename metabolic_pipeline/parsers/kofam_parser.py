"""KOfamScan output parser."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

KO_RE = re.compile(r"K\d{5}")


def _parse_kofam_line(line: str) -> dict[str, object] | None:
    text = line.strip()
    if not text or text.startswith("#"):
        return None

    passed = text.startswith("*")
    normalized = text.lstrip("*").strip()
    fields = re.split(r"\s+", normalized, maxsplit=6)
    if len(fields) < 2:
        return None

    protein_id = fields[0]
    ko_candidate = next((f for f in fields if KO_RE.fullmatch(f)), None)
    if not ko_candidate:
        return None

    threshold = float(fields[2]) if len(fields) > 2 and _is_float(fields[2]) else None
    score = float(fields[3]) if len(fields) > 3 and _is_float(fields[3]) else None
    evalue = fields[4] if len(fields) > 4 else None
    ko_description = fields[6] if len(fields) > 6 else (fields[5] if len(fields) > 5 else None)

    return {
        "protein_id": protein_id,
        "ko_id": ko_candidate,
        "ko_description": ko_description,
        "score": score,
        "threshold": threshold,
        "evalue": evalue,
        "passed_threshold": bool(passed),
        "source_tool": "kofamscan",
        "raw_annotation": line.rstrip("\n"),
    }


def _is_float(value: str) -> bool:
    try:
        float(value)
        return True
    except Exception:
        return False


def parse_kofam_output(path: Path, sample_id: str) -> pd.DataFrame:
    """Parse one KOfam output file."""
    if not path.exists():
        raise FileNotFoundError(f"Missing annotation output: {path}")
    rows: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parsed = _parse_kofam_line(line)
            if parsed is None:
                continue
            parsed["sample_id"] = sample_id
            rows.append(parsed)
    if not rows:
        return pd.DataFrame(
            columns=[
                "sample_id",
                "protein_id",
                "ko_id",
                "ko_description",
                "score",
                "threshold",
                "evalue",
                "passed_threshold",
                "source_tool",
                "raw_annotation",
            ]
        )
    return pd.DataFrame(rows)
