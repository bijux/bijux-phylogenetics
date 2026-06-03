from __future__ import annotations


def stable_float(value: float | None) -> float:
    """Normalize reviewer-facing floats to one deterministic string-backed value."""
    if value is None:
        raise ValueError("expected a float value")
    return float(format(round(float(value), 15), ".15g"))
