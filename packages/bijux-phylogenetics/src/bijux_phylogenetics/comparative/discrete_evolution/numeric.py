from __future__ import annotations


def _quantile(sorted_values: list[float], fraction: float) -> float:
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return float(format(sorted_values[0], ".15g"))
    index = max(
        0, min(len(sorted_values) - 1, int(round(fraction * (len(sorted_values) - 1))))
    )
    return float(format(sorted_values[index], ".15g"))
