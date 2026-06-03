from __future__ import annotations

from hashlib import sha256
from html import escape
from pathlib import Path


def _html_list(items: list[str]) -> str:
    return "<ul>" + "".join(f"<li>{item}</li>" for item in items) + "</ul>"


def _table(
    headers: list[str],
    rows: list[list[str]],
    *,
    max_rows: int | None = None,
) -> str:
    rendered_rows = rows if max_rows is None else rows[:max_rows]
    head = "".join(f"<th>{escape(header)}</th>" for header in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{escape(cell)}</td>" for cell in row) + "</tr>"
        for row in rendered_rows
    )
    truncation_note = ""
    if max_rows is not None and len(rows) > max_rows:
        truncation_note = (
            f"<p><em>Showing the first {max_rows} of {len(rows)} rows. "
            "Use the linked TSV artifacts for the full table.</em></p>"
        )
    return (
        truncation_note
        + f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
    )


def _support_range_text(
    minimum_support: float | None, maximum_support: float | None
) -> str:
    if minimum_support is None or maximum_support is None:
        return "not available"
    return f"{_format_number(minimum_support)}-{_format_number(maximum_support)}"


def _format_number(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".12g")


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()
