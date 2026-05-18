"""Pathway definition loading and validation."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

KO_RE = re.compile(r"^K\d{5}$")


@dataclass(frozen=True)
class PathwayDefinition:
    """A single user-defined pathway."""

    pathway_id: str
    pathway_name: str
    description: str
    required_kos: list[str]
    optional_kos: list[str]
    marker_kos: list[str]
    notes: str


def _validate_ko_list(values: object, field_name: str, pathway_id: str) -> list[str]:
    if values is None:
        return []
    if not isinstance(values, list):
        raise ValueError(f"{field_name} must be a list in pathway {pathway_id}.")
    out: list[str] = []
    for val in values:
        text = str(val)
        if not KO_RE.match(text):
            raise ValueError(f"Invalid KO format in {pathway_id}/{field_name}: {text}")
        out.append(text)
    return out


def load_pathway_definitions(path: Path) -> list[PathwayDefinition]:
    """Load and validate pathway definitions from JSON."""
    if not path.exists():
        raise FileNotFoundError(f"Pathway definition file not found: {path}")
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Malformed pathway JSON: {exc}") from exc

    if not isinstance(raw, list):
        raise ValueError("Malformed pathway JSON: top-level value must be a list.")

    definitions: list[PathwayDefinition] = []
    for item in raw:
        if not isinstance(item, dict):
            raise ValueError("Malformed pathway JSON: each pathway entry must be an object.")
        pathway_id = str(item.get("pathway_id", "")).strip()
        pathway_name = str(item.get("pathway_name", "")).strip()
        if not pathway_id or not pathway_name:
            raise ValueError("Malformed pathway JSON: pathway_id and pathway_name are required.")
        definitions.append(
            PathwayDefinition(
                pathway_id=pathway_id,
                pathway_name=pathway_name,
                description=str(item.get("description", "")),
                required_kos=_validate_ko_list(item.get("required_kos", []), "required_kos", pathway_id),
                optional_kos=_validate_ko_list(item.get("optional_kos", []), "optional_kos", pathway_id),
                marker_kos=_validate_ko_list(item.get("marker_kos", []), "marker_kos", pathway_id),
                notes=str(item.get("notes", "")),
            )
        )
    return definitions

