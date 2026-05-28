from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.runtime.errors import PhylogeneticsError

DISCRETE_TRAIT_RATE_PRIOR_FAMILIES = (
    "exponential",
    "gamma",
    "lognormal",
)
DISCRETE_TRAIT_RATE_PRIOR_MODELS = (
    "equal-rates",
    "symmetric",
    "all-rates-different",
)


@dataclass(frozen=True, slots=True)
class DiscreteTraitRatePriorModel:
    """One validated prior over one positive discrete-trait transition rate."""

    family: str
    rate: float | None = None
    shape: float | None = None
    scale: float | None = None
    log_mean: float | None = None
    log_standard_deviation: float | None = None

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
        return parameter_values


def build_exponential_discrete_trait_rate_prior(
    *,
    rate: float,
) -> DiscreteTraitRatePriorModel:
    """Build one exponential prior over one discrete-trait transition rate."""
    return DiscreteTraitRatePriorModel(
        family="exponential",
        rate=_validate_positive_finite_value(
            parameter_name="rate",
            value=rate,
            owner_name="exponential discrete-trait rate prior",
        ),
    )


def build_gamma_discrete_trait_rate_prior(
    *,
    shape: float,
    scale: float,
) -> DiscreteTraitRatePriorModel:
    """Build one gamma prior over one discrete-trait transition rate."""
    return DiscreteTraitRatePriorModel(
        family="gamma",
        shape=_validate_positive_finite_value(
            parameter_name="shape",
            value=shape,
            owner_name="gamma discrete-trait rate prior",
        ),
        scale=_validate_positive_finite_value(
            parameter_name="scale",
            value=scale,
            owner_name="gamma discrete-trait rate prior",
        ),
    )


def build_lognormal_discrete_trait_rate_prior(
    *,
    log_mean: float,
    log_standard_deviation: float,
) -> DiscreteTraitRatePriorModel:
    """Build one lognormal prior over one discrete-trait transition rate."""
    return DiscreteTraitRatePriorModel(
        family="lognormal",
        log_mean=_validate_finite_value(
            parameter_name="log_mean",
            value=log_mean,
            owner_name="lognormal discrete-trait rate prior",
        ),
        log_standard_deviation=_validate_positive_finite_value(
            parameter_name="log_standard_deviation",
            value=log_standard_deviation,
            owner_name="lognormal discrete-trait rate prior",
        ),
    )


def evaluate_discrete_trait_rate_value_log_prior(
    *,
    rate_value: float,
    prior_model: DiscreteTraitRatePriorModel,
) -> float:
    """Evaluate one discrete-trait transition-rate prior density."""
    validated_rate_value = _validate_positive_finite_value(
        parameter_name="rate_value",
        value=rate_value,
        owner_name="discrete-trait rate prior evaluation",
    )
    if prior_model.family == "exponential":
        assert prior_model.rate is not None
        return math.log(prior_model.rate) - (prior_model.rate * validated_rate_value)
    if prior_model.family == "gamma":
        assert prior_model.shape is not None
        assert prior_model.scale is not None
        return _gamma_log_density(
            validated_rate_value,
            shape=prior_model.shape,
            scale=prior_model.scale,
        )
    if prior_model.family == "lognormal":
        assert prior_model.log_mean is not None
        assert prior_model.log_standard_deviation is not None
        return _lognormal_log_density(
            validated_rate_value,
            log_mean=prior_model.log_mean,
            log_standard_deviation=prior_model.log_standard_deviation,
        )
    raise PhylogeneticsError(
        "discrete-trait rate prior family is unsupported",
        code="discrete_trait_rate_prior_family_invalid",
        details={
            "family": prior_model.family,
            "allowed_families": list(DISCRETE_TRAIT_RATE_PRIOR_FAMILIES),
        },
    )


def _gamma_log_density(
    rate_value: float,
    *,
    shape: float,
    scale: float,
) -> float:
    return (
        ((shape - 1.0) * math.log(rate_value))
        - (rate_value / scale)
        - math.lgamma(shape)
        - (shape * math.log(scale))
    )


def _lognormal_log_density(
    rate_value: float,
    *,
    log_mean: float,
    log_standard_deviation: float,
) -> float:
    return (
        -math.log(rate_value)
        - math.log(log_standard_deviation)
        - (0.5 * math.log(2.0 * math.pi))
        - (
            ((math.log(rate_value) - log_mean) ** 2)
            / (2.0 * (log_standard_deviation**2))
        )
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
            code="discrete_trait_rate_prior_parameter_nonpositive",
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
            code="discrete_trait_rate_prior_parameter_nonfinite",
            details={
                "parameter_name": parameter_name,
                "parameter_value": value,
                "owner_name": owner_name,
            },
        )
    return validated_value
