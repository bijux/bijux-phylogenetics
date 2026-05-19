from __future__ import annotations


def transition_recovery_match(
    *,
    true_changed: bool,
    estimated_changed: bool,
    true_source: str,
    true_target: str,
    estimated_source: str,
    estimated_target: str,
    true_event_count: int,
    estimated_event_count: int,
) -> bool:
    return (
        true_changed == estimated_changed
        and true_source == estimated_source
        and true_target == estimated_target
        and true_event_count == estimated_event_count
    )


def format_transition(source_state: str, target_state: str) -> str:
    return f"{source_state}->{target_state}"


def evaluate_threshold(
    *,
    comparator: str,
    threshold: str,
    observed_value: bool | float,
) -> bool:
    if comparator == "==":
        if isinstance(observed_value, bool):
            return observed_value is (threshold.strip().lower() == "true")
        return float(observed_value) == float(threshold)
    threshold_value = float(threshold)
    numeric_observed = float(observed_value)
    if comparator == "<=":
        return numeric_observed <= threshold_value
    if comparator == ">=":
        return numeric_observed >= threshold_value
    raise ValueError(f"unsupported threshold comparator '{comparator}'")


def mean(values) -> float:
    materialized = list(values)
    if not materialized:
        return 0.0
    return sum(materialized) / len(materialized)


def format_number(value: float) -> str:
    return format(value, ".15g")


def format_observed_value(value: bool | float) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    return format_number(value)
