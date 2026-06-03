from __future__ import annotations

import math


def validate_increasing_parameter_bounds(
    *,
    parameter_name: str,
    lower_bound: float,
    upper_bound: float,
    owner_name: str,
) -> tuple[float, float]:
    validated_lower_bound = _validate_finite_bound(
        parameter_name=parameter_name,
        bound_name="lower bound",
        bound_value=lower_bound,
        owner_name=owner_name,
    )
    validated_upper_bound = _validate_finite_bound(
        parameter_name=parameter_name,
        bound_name="upper bound",
        bound_value=upper_bound,
        owner_name=owner_name,
    )
    if validated_upper_bound <= validated_lower_bound:
        raise ValueError(
            f"{owner_name} requires strictly increasing bounds for '{parameter_name}'"
        )
    return validated_lower_bound, validated_upper_bound


def validate_parameter_within_bounds(
    *,
    parameter_name: str,
    value: float,
    lower_bound: float,
    upper_bound: float,
    owner_name: str,
) -> float:
    validated_value = _validate_finite_bound(
        parameter_name=parameter_name,
        bound_name="value",
        bound_value=value,
        owner_name=owner_name,
    )
    validated_lower_bound, validated_upper_bound = validate_increasing_parameter_bounds(
        parameter_name=parameter_name,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        owner_name=owner_name,
    )
    if not (validated_lower_bound <= validated_value <= validated_upper_bound):
        raise ValueError(
            f"{owner_name} requires '{parameter_name}' to lie within [{validated_lower_bound}, {validated_upper_bound}]"
        )
    return validated_value


def validate_positive_parameter_bounds(
    *,
    parameter_name: str,
    lower_bound: float,
    upper_bound: float,
    owner_name: str,
) -> tuple[float, float]:
    validated_lower_bound, validated_upper_bound = validate_increasing_parameter_bounds(
        parameter_name=parameter_name,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        owner_name=owner_name,
    )
    if validated_lower_bound <= 0.0:
        raise ValueError(
            f"{owner_name} requires a positive lower bound for '{parameter_name}'"
        )
    return validated_lower_bound, validated_upper_bound


def _validate_finite_bound(
    *,
    parameter_name: str,
    bound_name: str,
    bound_value: float,
    owner_name: str,
) -> float:
    if not math.isfinite(bound_value):
        raise ValueError(
            f"{owner_name} requires a finite {bound_name} for '{parameter_name}'"
        )
    return float(bound_value)
