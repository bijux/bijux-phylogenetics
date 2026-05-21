from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any

from ....artifacts import preview_report_rows, truncate_report_rows


def artifact_link(
    *,
    artifact_paths: dict[str, Path],
    artifact_key: str,
    out_path: Path,
) -> str:
    """Build a stable relative artifact path for reviewer-facing sections."""
    return artifact_paths[artifact_key].relative_to(out_path.parent).as_posix()


def truncate_dataclass_rows(
    *,
    rows: list[Any],
    limit: int,
    section_name: str,
    truncated_sections: list[str],
) -> tuple[list[dict[str, Any]], int]:
    """Serialize dataclass rows and apply the report-row budget."""
    return truncate_report_rows(
        [asdict(row) for row in rows],
        limit=limit,
        section_name=section_name,
        truncated_sections=truncated_sections,
    )


def preview_payload(
    *,
    rows: list[dict[str, Any]],
    row_count: int,
    truncated_row_count: int,
    preview_limit: int,
) -> dict[str, Any]:
    """Build the standard preview payload for linked uncertainty tables."""
    return {
        "row_count": row_count,
        "truncated_row_count": truncated_row_count,
        "preview_row_count": min(len(rows), preview_limit),
        "preview_rows": preview_report_rows(rows, limit=preview_limit),
    }
