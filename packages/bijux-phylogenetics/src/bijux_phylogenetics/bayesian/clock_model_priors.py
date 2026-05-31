from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.bayesian.required_values import require_present
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

CLOCK_MODEL_SCALAR_PRIOR_FAMILIES = (
    "exponential",
    "fixed",
    "gamma",
    "lognormal",
)

_FIXED_PRIOR_TOLERANCE = 1e-12


@dataclass(frozen=True, slots=True)
class ClockModelScalarPriorModel:
    """One validated scalar prior over one positive clock-model parameter."""

    family: str
    rate: float | None = None
    shape: float | None = None
    scale: float | None = None
    log_mean: float | None = None
    log_standard_deviation: float | None = None
    fixed_value: float | None = None
    fixed_tolerance: float | None = None

    def parameter_values(self) -> dict[str, float]:
        parameter_values: dict[str, float] = {}
        if self.rate is not None:
            parameter_values["rate"] = self.rate
        if self.shape is not None:
            parameter_values["shape"] = self.shape
        if self.scale is not None:
            parameter_values["scale"] = self.scale
        if self.log_mean is not None:
            parameter_values["log_mean"] = self.log_mean
        if self.log_standard_deviation is not None:
            parameter_values["log_standard_deviation"] = self.log_standard_deviation
        if self.fixed_value is not None:
            parameter_values["fixed_value"] = self.fixed_value
        if self.fixed_tolerance is not None:
            parameter_values["fixed_tolerance"] = self.fixed_tolerance
        return parameter_values


def build_exponential_clock_model_scalar_prior(
    *,
    rate: float,
) -> ClockModelScalarPriorModel:
    """Build one exponential prior over one positive clock-model parameter."""
    return ClockModelScalarPriorModel(
        family="exponential",
        rate=_validate_positive_finite_value(
            parameter_name="rate",
            value=rate,
            owner_name="exponential clock-model scalar prior",
        ),
    )


def build_gamma_clock_model_scalar_prior(
    *,
    shape: float,
    scale: float,
) -> ClockModelScalarPriorModel:
    """Build one gamma prior over one positive clock-model parameter."""
    return ClockModelScalarPriorModel(
        family="gamma",
        shape=_validate_positive_finite_value(
            parameter_name="shape",
            value=shape,
            owner_name="gamma clock-model scalar prior",
        ),
        scale=_validate_positive_finite_value(
            parameter_name="scale",
            value=scale,
            owner_name="gamma clock-model scalar prior",
        ),
    )


def build_lognormal_clock_model_scalar_prior(
    *,
    log_mean: float,
    log_standard_deviation: float,
) -> ClockModelScalarPriorModel:
    """Build one lognormal prior over one positive clock-model parameter."""
    return ClockModelScalarPriorModel(
        family="lognormal",
        log_mean=_validate_finite_value(
            parameter_name="log_mean",
            value=log_mean,
            owner_name="lognormal clock-model scalar prior",
        ),
        log_standard_deviation=_validate_positive_finite_value(
            parameter_name="log_standard_deviation",
            value=log_standard_deviation,
            owner_name="lognormal clock-model scalar prior",
        ),
    )


def build_fixed_clock_model_scalar_prior(
    *,
    fixed_value: float,
    fixed_tolerance: float = _FIXED_PRIOR_TOLERANCE,
) -> ClockModelScalarPriorModel:
    """Build one fixed-value prior over one positive clock-model parameter."""
    return ClockModelScalarPriorModel(
        family="fixed",
        fixed_value=_validate_positive_finite_value(
            parameter_name="fixed_value",
            value=fixed_value,
            owner_name="fixed clock-model scalar prior",
        ),
        fixed_tolerance=_validate_nonnegative_finite_value(
            parameter_name="fixed_tolerance",
            value=fixed_tolerance,
            owner_name="fixed clock-model scalar prior",
        ),
    )


def evaluate_clock_model_scalar_log_prior(
    *,
    parameter_value: float,
    prior_model: ClockModelScalarPriorModel,
    parameter_name: str,
) -> float:
    """Evaluate one positive clock-model parameter log prior."""
    validated_value = _validate_positive_finite_value(
        parameter_name=parameter_name,
        value=parameter_value,
        owner_name="clock-model scalar prior evaluation",
    )
    if prior_model.family == "exponential":
        rate = require_present(
            prior_model.rate,
            owner_name="clock-model scalar prior evaluation",
            field_name="rate",
        )
        return math.log(rate) - (rate * validated_value)
    if prior_model.family == "gamma":
        shape = require_present(
            prior_model.shape,
            owner_name="clock-model scalar prior evaluation",
            field_name="shape",
        )
        scale = require_present(
            prior_model.scale,
            owner_name="clock-model scalar prior evaluation",
            field_name="scale",
        )
        return (
            ((shape - 1.0) * math.log(validated_value))
            - (validated_value / scale)
            - math.lgamma(shape)
            - (shape * math.log(scale))
        )
    if prior_model.family == "lognormal":
        log_mean = require_present(
            prior_model.log_mean,
            owner_name="clock-model scalar prior evaluation",
            field_name="log_mean",
        )
        log_standard_deviation = require_present(
            prior_model.log_standard_deviation,
            owner_name="clock-model scalar prior evaluation",
            field_name="log_standard_deviation",
        )
        return (
            -math.log(validated_value)
            - math.log(log_standard_deviation)
            - (0.5 * math.log(2.0 * math.pi))
            - (
                ((math.log(validated_value) - log_mean) ** 2)
                / (2.0 * (log_standard_deviation**2))
            )
        )
    if prior_model.family == "fixed":
        fixed_value = require_present(
            prior_model.fixed_value,
            owner_name="clock-model scalar prior evaluation",
            field_name="fixed_value",
        )
        fixed_tolerance = require_present(
            prior_model.fixed_tolerance,
            owner_name="clock-model scalar prior evaluation",
            field_name="fixed_tolerance",
        )
        return (
            0.0
            if math.isclose(
                validated_value,
                fixed_value,
                rel_tol=0.0,
                abs_tol=fixed_tolerance,
            )
            else -math.inf
        )
    raise PhylogeneticsError(
        "clock-model scalar prior family is unsupported",
        code="clock_model_scalar_prior_family_invalid",
        details={
            "family": prior_model.family,
            "allowed_families": list(CLOCK_MODEL_SCALAR_PRIOR_FAMILIES),
        },
    )


def _validate_positive_finite_value(
    *,
    parameter_name: str,
    value: float,
    owner_name: str,
) -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{parameter_name}' to be positive and finite",
            code="clock_model_scalar_prior_positive_value_invalid",
            details={parameter_name: value},
        )
    return float(format(value, ".15g"))


def _validate_nonnegative_finite_value(
    *,
    parameter_name: str,
    value: float,
    owner_name: str,
) -> float:
    if not math.isfinite(value) or value < 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{parameter_name}' to be non-negative and finite",
            code="clock_model_scalar_prior_nonnegative_value_invalid",
            details={parameter_name: value},
        )
    return float(format(value, ".15g"))


def _validate_finite_value(
    *,
    parameter_name: str,
    value: float,
    owner_name: str,
) -> float:
    if not math.isfinite(value):
        raise PhylogeneticsError(
            f"{owner_name} requires '{parameter_name}' to be finite",
            code="clock_model_scalar_prior_finite_value_invalid",
            details={parameter_name: value},
        )
    return float(format(value, ".15g"))
