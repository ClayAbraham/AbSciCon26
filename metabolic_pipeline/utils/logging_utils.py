"""Logging utilities."""

from __future__ import annotations

import logging
from pathlib import Path


def setup_logger(output_dir: Path) -> logging.Logger:
    """Create a file and console logger."""
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "run.log"
    logger = logging.getLogger("metabolic_pipeline")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(fmt)
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)
    return logger

