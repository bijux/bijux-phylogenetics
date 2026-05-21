from __future__ import annotations

from pathlib import Path


def _write_tsv(path: Path, rows: list[list[str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )
    return path


def _format_float(value: float) -> str:
    return format(value, ".12g")


def _format_optional_float(value: float | None) -> str:
    return "" if value is None else format(value, ".12g")


def _format_optional_bool(value: bool | None) -> str:
    if value is None:
        return ""
    return str(value).lower()
