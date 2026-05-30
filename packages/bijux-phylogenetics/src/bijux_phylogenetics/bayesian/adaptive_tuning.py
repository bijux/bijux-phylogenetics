from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.runtime.errors import PhylogeneticsError


@dataclass(frozen=True, slots=True)
class AdaptiveTuningController:
    """One validated burn-in tuning policy for a scalar proposal scale."""

    proposal_name: str
    scale_parameter_name: str
    initial_scale: float
    target_acceptance_rate: float
    burnin_iteration_count: int
    adaptation_window_size: int
    decrease_factor: float
    increase_factor: float
    minimum_scale: float
    maximum_scale: float


def build_adaptive_tuning_controller(
    *,
    proposal_name: str,
    scale_parameter_name: str,
    initial_scale: float,
    target_acceptance_rate: float,
    burnin_iteration_count: int,
    adaptation_window_size: int,
    decrease_factor: float = 0.5,
    increase_factor: float = 2.0,
    minimum_scale: float = 1e-6,
    maximum_scale: float = 1e6,
) -> AdaptiveTuningController:
    """Build one validated controller for burn-in-only proposal adaptation."""
    validated_proposal_name = _validate_nonblank_name(
        value=proposal_name,
        field_name="proposal_name",
        owner_name="adaptive tuning controller",
    )
    validated_scale_parameter_name = _validate_nonblank_name(
        value=scale_parameter_name,
        field_name="scale_parameter_name",
        owner_name="adaptive tuning controller",
    )
    validated_initial_scale = _validate_positive_finite_float(
        value=initial_scale,
        field_name="initial_scale",
        owner_name="adaptive tuning controller",
    )
    validated_target_acceptance_rate = _validate_probability_rate(
        value=target_acceptance_rate,
        field_name="target_acceptance_rate",
        owner_name="adaptive tuning controller",
    )
    validated_burnin_iteration_count = _validate_nonnegative_integer(
        value=burnin_iteration_count,
        field_name="burnin_iteration_count",
        owner_name="adaptive tuning controller",
    )
    validated_adaptation_window_size = _validate_positive_integer(
        value=adaptation_window_size,
        field_name="adaptation_window_size",
        owner_name="adaptive tuning controller",
    )
    validated_decrease_factor = _validate_open_interval_float(
        value=decrease_factor,
        field_name="decrease_factor",
        owner_name="adaptive tuning controller",
        lower_bound=0.0,
        upper_bound=1.0,
    )
    validated_increase_factor = _validate_greater_than_float(
        value=increase_factor,
        field_name="increase_factor",
        owner_name="adaptive tuning controller",
        lower_bound=1.0,
    )
    validated_minimum_scale = _validate_positive_finite_float(
        value=minimum_scale,
        field_name="minimum_scale",
        owner_name="adaptive tuning controller",
    )
    validated_maximum_scale = _validate_positive_finite_float(
        value=maximum_scale,
        field_name="maximum_scale",
        owner_name="adaptive tuning controller",
    )
    if validated_minimum_scale > validated_maximum_scale:
        raise PhylogeneticsError(
            "adaptive tuning controller requires 'minimum_scale' to be less than or equal to 'maximum_scale'",
            code="adaptive_tuning_scale_bounds_invalid",
        )
    if not validated_minimum_scale <= validated_initial_scale <= validated_maximum_scale:
        raise PhylogeneticsError(
            "adaptive tuning controller requires 'initial_scale' to lie within the configured scale bounds",
            code="adaptive_tuning_initial_scale_out_of_bounds",
            details={
                "initial_scale": validated_initial_scale,
                "minimum_scale": validated_minimum_scale,
                "maximum_scale": validated_maximum_scale,
            },
        )
    return AdaptiveTuningController(
        proposal_name=validated_proposal_name,
        scale_parameter_name=validated_scale_parameter_name,
        initial_scale=validated_initial_scale,
        target_acceptance_rate=validated_target_acceptance_rate,
        burnin_iteration_count=validated_burnin_iteration_count,
        adaptation_window_size=validated_adaptation_window_size,
        decrease_factor=validated_decrease_factor,
        increase_factor=validated_increase_factor,
        minimum_scale=validated_minimum_scale,
        maximum_scale=validated_maximum_scale,
    )


def _validate_nonblank_name(
    *,
    value: str,
    field_name: str,
    owner_name: str,
) -> str:
    if not isinstance(value, str) or not value.strip():
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one nonblank name",
            code="adaptive_tuning_name_invalid",
            details={"field_name": field_name},
        )
    return value.strip()


def _validate_positive_integer(
    *,
    value: int,
    field_name: str,
    owner_name: str,
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one integer",
            code="adaptive_tuning_integer_required",
            details={"field_name": field_name},
        )
    if value <= 0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be positive",
            code="adaptive_tuning_positive_integer_required",
            details={"field_name": field_name},
        )
    return value


def _validate_nonnegative_integer(
    *,
    value: int,
    field_name: str,
    owner_name: str,
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one integer",
            code="adaptive_tuning_integer_required",
            details={"field_name": field_name},
        )
    if value < 0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be nonnegative",
            code="adaptive_tuning_nonnegative_integer_required",
            details={"field_name": field_name},
        )
    return value


def _validate_positive_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be numeric",
            code="adaptive_tuning_float_required",
            details={"field_name": field_name},
        )
    normalized_value = float(value)
    if not math.isfinite(normalized_value):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be finite",
            code="adaptive_tuning_finite_float_required",
            details={"field_name": field_name},
        )
    if normalized_value <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be strictly positive",
            code="adaptive_tuning_positive_float_required",
            details={"field_name": field_name},
        )
    return normalized_value


def _validate_probability_rate(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be numeric",
            code="adaptive_tuning_float_required",
            details={"field_name": field_name},
        )
    normalized_value = float(value)
    if not math.isfinite(normalized_value):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be finite",
            code="adaptive_tuning_finite_float_required",
            details={"field_name": field_name},
        )
    if not 0.0 < normalized_value < 1.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to lie strictly between 0 and 1",
            code="adaptive_tuning_probability_rate_invalid",
            details={"field_name": field_name},
        )
    return normalized_value


def _validate_open_interval_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
    lower_bound: float,
    upper_bound: float,
) -> float:
    normalized_value = _validate_positive_finite_float(
        value=value,
        field_name=field_name,
        owner_name=owner_name,
    )
    if not lower_bound < normalized_value < upper_bound:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to lie strictly between {lower_bound} and {upper_bound}",
            code="adaptive_tuning_interval_float_invalid",
            details={
                "field_name": field_name,
                "lower_bound": lower_bound,
                "upper_bound": upper_bound,
            },
        )
    return normalized_value


def _validate_greater_than_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
    lower_bound: float,
) -> float:
    normalized_value = _validate_positive_finite_float(
        value=value,
        field_name=field_name,
        owner_name=owner_name,
    )
    if normalized_value <= lower_bound:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be greater than {lower_bound}",
            code="adaptive_tuning_lower_bound_invalid",
            details={"field_name": field_name, "lower_bound": lower_bound},
        )
    return normalized_value
