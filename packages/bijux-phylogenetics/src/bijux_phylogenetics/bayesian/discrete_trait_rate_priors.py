from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import math

from bijux_phylogenetics.ancestral.discrete import DiscreteTransitionRateRow
from bijux_phylogenetics.bayesian.discrete_trait_rate_parameters import (
    DISCRETE_TRAIT_RATE_PARAMETER_MODELS,
    parameterize_discrete_trait_rate_rows,
)
from bijux_phylogenetics.bayesian.required_values import require_present
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

DISCRETE_TRAIT_RATE_PRIOR_FAMILIES = (
    "exponential",
    "gamma",
    "lognormal",
)
DISCRETE_TRAIT_RATE_PRIOR_MODELS = DISCRETE_TRAIT_RATE_PARAMETER_MODELS


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


@dataclass(frozen=True, slots=True)
class DiscreteTraitRatePriorRow:
    """One discrete-trait rate-parameter prior contribution."""

    parameter_name: str
    transition_pairs: list[tuple[str, str]]
    rate_value: float
    log_prior_contribution: float


@dataclass(frozen=True, slots=True)
class DiscreteTraitRatePriorEvaluationReport:
    """One discrete-trait transition-rate prior evaluation report."""

    model: str
    family: str
    parameter_count: int
    parameter_values: dict[str, float]
    total_log_prior: float
    rows: list[DiscreteTraitRatePriorRow]


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
        rate = require_present(
            prior_model.rate,
            owner_name="discrete-trait rate prior evaluation",
            field_name="rate",
        )
        return math.log(rate) - (rate * validated_rate_value)
    if prior_model.family == "gamma":
        shape = require_present(
            prior_model.shape,
            owner_name="discrete-trait rate prior evaluation",
            field_name="shape",
        )
        scale = require_present(
            prior_model.scale,
            owner_name="discrete-trait rate prior evaluation",
            field_name="scale",
        )
        return _gamma_log_density(
            validated_rate_value,
            shape=shape,
            scale=scale,
        )
    if prior_model.family == "lognormal":
        log_mean = require_present(
            prior_model.log_mean,
            owner_name="discrete-trait rate prior evaluation",
            field_name="log_mean",
        )
        log_standard_deviation = require_present(
            prior_model.log_standard_deviation,
            owner_name="discrete-trait rate prior evaluation",
            field_name="log_standard_deviation",
        )
        return _lognormal_log_density(
            validated_rate_value,
            log_mean=log_mean,
            log_standard_deviation=log_standard_deviation,
        )
    raise PhylogeneticsError(
        "discrete-trait rate prior family is unsupported",
        code="discrete_trait_rate_prior_family_invalid",
        details={
            "family": prior_model.family,
            "allowed_families": list(DISCRETE_TRAIT_RATE_PRIOR_FAMILIES),
        },
    )


def evaluate_discrete_trait_rate_log_prior(
    *,
    model: str,
    transition_rate_rows: Sequence[DiscreteTransitionRateRow],
    prior_model: DiscreteTraitRatePriorModel,
) -> DiscreteTraitRatePriorEvaluationReport:
    """Evaluate one prior over the governed rate parameters of one Mk fit surface."""
    try:
        parameterization = parameterize_discrete_trait_rate_rows(
            model=model,
            transition_rate_rows=transition_rate_rows,
        )
    except PhylogeneticsError as error:
        raise _translate_parameterization_error_for_prior(error) from error
    rows = [
        DiscreteTraitRatePriorRow(
            parameter_name=group.parameter_name,
            transition_pairs=group.transition_pairs,
            rate_value=group.rate_value,
            log_prior_contribution=float(
                format(
                    evaluate_discrete_trait_rate_value_log_prior(
                        rate_value=group.rate_value,
                        prior_model=prior_model,
                    ),
                    ".15g",
                )
            ),
        )
        for group in parameterization.groups
    ]
    total_log_prior = math.fsum(row.log_prior_contribution for row in rows)
    return DiscreteTraitRatePriorEvaluationReport(
        model=parameterization.model,
        family=prior_model.family,
        parameter_count=len(rows),
        parameter_values=prior_model.parameter_values(),
        total_log_prior=float(format(total_log_prior, ".15g")),
        rows=rows,
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


def _translate_parameterization_error_for_prior(
    error: PhylogeneticsError,
) -> PhylogeneticsError:
    code = error.code
    if code == "discrete_trait_rate_parameter_model_unsupported":
        return PhylogeneticsError(
            "discrete-trait rate priors support only ER, SYM, and ARD transition-rate surfaces",
            code="discrete_trait_rate_prior_model_unsupported",
            details=error.details,
        )
    if code == "discrete_trait_rate_parameter_no_allowed_rates":
        return PhylogeneticsError(
            "discrete-trait rate prior evaluation requires at least one allowed transition rate",
            code="discrete_trait_rate_prior_no_allowed_rates",
            details=error.details,
        )
    if code == "discrete_trait_rate_parameter_equal_rates_inconsistent":
        return PhylogeneticsError(
            "equal-rates prior evaluation requires one shared transition rate",
            code="discrete_trait_rate_prior_equal_rates_inconsistent",
            details=error.details,
        )
    if code == "discrete_trait_rate_parameter_symmetric_pair_incomplete":
        return PhylogeneticsError(
            "symmetric prior evaluation requires one forward and one reverse transition per rate parameter",
            code="discrete_trait_rate_prior_symmetric_pair_incomplete",
            details=error.details,
        )
    if code == "discrete_trait_rate_parameter_symmetric_pair_mismatched":
        return PhylogeneticsError(
            "symmetric prior evaluation requires matched forward and reverse transition rates",
            code="discrete_trait_rate_prior_symmetric_pair_mismatched",
            details=error.details,
        )
    return error


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
