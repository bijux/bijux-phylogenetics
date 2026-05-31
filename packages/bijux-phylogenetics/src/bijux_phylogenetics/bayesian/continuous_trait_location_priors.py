from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.bayesian.required_values import require_present
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

CONTINUOUS_TRAIT_LOCATION_PRIOR_FAMILIES = ("fixed", "normal")
_FIXED_PRIOR_TOLERANCE = 1e-12


@dataclass(frozen=True, slots=True)
class ContinuousTraitLocationPriorModel:
    """One validated prior over one real-valued continuous-trait location parameter."""

    family: str
    mean: float | None = None
    standard_deviation: float | None = None
    fixed_value: float | None = None
    fixed_tolerance: float | None = None

    def parameter_values(self) -> dict[str, float]:
        parameter_values: dict[str, float] = {}
        if self.mean is not None:
            parameter_values["mean"] = self.mean
        if self.standard_deviation is not None:
            parameter_values["standard_deviation"] = self.standard_deviation
        if self.fixed_value is not None:
            parameter_values["fixed_value"] = self.fixed_value
        if self.fixed_tolerance is not None:
            parameter_values["fixed_tolerance"] = self.fixed_tolerance
        return parameter_values


def build_normal_continuous_trait_location_prior(
    *,
    mean: float,
    standard_deviation: float,
) -> ContinuousTraitLocationPriorModel:
    """Build one normal prior over one continuous-trait location parameter."""
    return ContinuousTraitLocationPriorModel(
        family="normal",
        mean=_validate_finite_value(
            parameter_name="mean",
            value=mean,
            owner_name="normal continuous-trait location prior",
        ),
        standard_deviation=_validate_positive_finite_value(
            parameter_name="standard_deviation",
            value=standard_deviation,
            owner_name="normal continuous-trait location prior",
        ),
    )


def build_fixed_continuous_trait_location_prior(
    *,
    fixed_value: float,
    fixed_tolerance: float = _FIXED_PRIOR_TOLERANCE,
) -> ContinuousTraitLocationPriorModel:
    """Build one fixed-value prior over one continuous-trait location parameter."""
    return ContinuousTraitLocationPriorModel(
        family="fixed",
        fixed_value=_validate_finite_value(
            parameter_name="fixed_value",
            value=fixed_value,
            owner_name="fixed continuous-trait location prior",
        ),
        fixed_tolerance=_validate_nonnegative_finite_value(
            parameter_name="fixed_tolerance",
            value=fixed_tolerance,
            owner_name="fixed continuous-trait location prior",
        ),
    )


def evaluate_continuous_trait_location_log_prior(
    *,
    parameter_value: float,
    prior_model: ContinuousTraitLocationPriorModel,
    parameter_name: str,
) -> float:
    """Evaluate one real-valued continuous-trait location log prior."""
    validated_value = _validate_finite_value(
        parameter_name=parameter_name,
        value=parameter_value,
        owner_name="continuous-trait location prior evaluation",
    )
    if prior_model.family == "normal":
        mean = require_present(
            prior_model.mean,
            owner_name="continuous-trait location prior evaluation",
            field_name="mean",
        )
        standard_deviation = require_present(
            prior_model.standard_deviation,
            owner_name="continuous-trait location prior evaluation",
            field_name="standard_deviation",
        )
        return _normal_log_density(
            validated_value,
            mean=mean,
            standard_deviation=standard_deviation,
        )
    if prior_model.family == "fixed":
        fixed_value = require_present(
            prior_model.fixed_value,
            owner_name="continuous-trait location prior evaluation",
            field_name="fixed_value",
        )
        fixed_tolerance = require_present(
            prior_model.fixed_tolerance,
            owner_name="continuous-trait location prior evaluation",
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
        "continuous-trait location prior family is unsupported",
        code="continuous_trait_location_prior_family_invalid",
        details={
            "family": prior_model.family,
            "allowed_families": list(CONTINUOUS_TRAIT_LOCATION_PRIOR_FAMILIES),
        },
    )


def _normal_log_density(
    parameter_value: float,
    *,
    mean: float,
    standard_deviation: float,
) -> float:
    centered_value = (parameter_value - mean) / standard_deviation
    return (
        -math.log(standard_deviation)
        - (0.5 * math.log(2.0 * math.pi))
        - (0.5 * (centered_value**2))
    )


def _validate_positive_finite_value(
    *,
    parameter_name: str,
    value: float,
    owner_name: str,
) -> float:
    validated_value = _validate_finite_value(
        parameter_name=parameter_name,
        value=value,
        owner_name=owner_name,
    )
    if validated_value <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires {parameter_name} > 0",
            code="continuous_trait_location_prior_parameter_nonpositive",
            details={
                "parameter_name": parameter_name,
                "parameter_value": value,
                "owner_name": owner_name,
            },
        )
    return validated_value


def _validate_nonnegative_finite_value(
    *,
    parameter_name: str,
    value: float,
    owner_name: str,
) -> float:
    validated_value = _validate_finite_value(
        parameter_name=parameter_name,
        value=value,
        owner_name=owner_name,
    )
    if validated_value < 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires {parameter_name} >= 0",
            code="continuous_trait_location_prior_parameter_negative",
            details={
                "parameter_name": parameter_name,
                "parameter_value": value,
                "owner_name": owner_name,
            },
        )
    return validated_value


def _validate_finite_value(
    *,
    parameter_name: str,
    value: float,
    owner_name: str,
) -> float:
    validated_value = float(value)
    if not math.isfinite(validated_value):
        raise PhylogeneticsError(
            f"{owner_name} requires finite {parameter_name}",
            code="continuous_trait_location_prior_parameter_nonfinite",
            details={
                "parameter_name": parameter_name,
                "parameter_value": value,
                "owner_name": owner_name,
            },
        )
    return validated_value


__all__ = [
    "CONTINUOUS_TRAIT_LOCATION_PRIOR_FAMILIES",
    "ContinuousTraitLocationPriorModel",
    "build_fixed_continuous_trait_location_prior",
    "build_normal_continuous_trait_location_prior",
    "evaluate_continuous_trait_location_log_prior",
]
