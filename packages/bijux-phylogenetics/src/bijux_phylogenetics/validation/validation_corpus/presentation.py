from __future__ import annotations

from dataclasses import asdict, is_dataclass
import json
from pathlib import Path


def normalize_jsonable(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value):
        return {key: normalize_jsonable(item) for key, item in asdict(value).items()}
    if isinstance(value, dict):
        return {str(key): normalize_jsonable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [normalize_jsonable(item) for item in value]
    return value


def write_validation_corpus_json(path: Path, report: object) -> Path:
    """Write a validation-corpus or dashboard report as deterministic JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = normalize_jsonable(report)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path
