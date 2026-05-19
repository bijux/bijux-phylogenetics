from __future__ import annotations


def _format_number(value: float) -> str:
    return format(value, ".15g")


def _format_optional_bool(value: bool | None) -> str:
    if value is None:
        return ""
    return str(value).lower()


def _optional_string(value: object) -> str | None:
    return value if isinstance(value, str) and value else None


def _optional_float(value: object) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _optional_bool(value: object) -> bool | None:
    if isinstance(value, bool):
        return value
    return None
