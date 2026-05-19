from __future__ import annotations


def _format_number(value: float) -> str:
    return format(value, ".15g")


def _format_optional_number(value: float | None) -> str:
    if value is None:
        return ""
    return _format_number(value)


def _format_optional_bool(value: bool | None) -> str:
    if value is None:
        return ""
    return str(value).lower()


def _optional_string(value: object) -> str | None:
    return value if isinstance(value, str) and value else None
