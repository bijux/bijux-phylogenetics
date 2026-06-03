from __future__ import annotations

import csv
import json
from pathlib import Path


def section(name: str, payload: object) -> tuple[str, str]:
    return name, json.dumps(payload, default=str, indent=2, sort_keys=True)


def truncate_report_rows(
    rows: list[object],
    *,
    limit: int | None,
    section_name: str,
    truncated_sections: list[str],
) -> tuple[list[object], int]:
    if limit is None or len(rows) <= limit:
        return rows, 0
    truncated_sections.append(section_name)
    return rows[:limit], len(rows) - limit


def preview_report_rows(rows: list[object], *, limit: int = 5) -> list[object]:
    if limit <= 0:
        return []
    return rows[:limit]


def write_json_artifact(path: Path, payload: object) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_tabular_artifact(path: Path, rows: list[dict[str, object]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value, sort_keys=True)
                    if isinstance(value, (dict, list))
                    else value
                    for key, value in row.items()
                }
            )
    return path


def write_machine_manifest(path: Path, payload: dict[str, object]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def report_sidecar_path(out_path: Path) -> Path:
    return out_path.with_suffix(".json")
