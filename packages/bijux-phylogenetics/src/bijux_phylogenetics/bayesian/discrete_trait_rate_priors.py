from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import math

from bijux_phylogenetics.ancestral.discrete import DiscreteTransitionRateRow
from bijux_phylogenetics.ancestral.discrete.policy import resolve_discrete_model_name
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
_RATE_GROUP_TOLERANCE = 1e-12


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


def evaluate_discrete_trait_rate_log_prior(
    *,
    model: str,
    transition_rate_rows: Sequence[DiscreteTransitionRateRow],
    prior_model: DiscreteTraitRatePriorModel,
) -> DiscreteTraitRatePriorEvaluationReport:
    """Evaluate one prior over the governed rate parameters of one Mk fit surface."""
    resolved_model = _resolve_supported_discrete_trait_rate_model(model)
    rate_groups = _resolve_rate_parameter_groups(
        model=resolved_model,
        transition_rate_rows=transition_rate_rows,
    )
    rows = [
        DiscreteTraitRatePriorRow(
            parameter_name=parameter_name,
            transition_pairs=transition_pairs,
            rate_value=float(format(rate_value, ".15g")),
            log_prior_contribution=float(
                format(
                    evaluate_discrete_trait_rate_value_log_prior(
                        rate_value=rate_value,
                        prior_model=prior_model,
                    ),
                    ".15g",
                )
            ),
        )
        for parameter_name, transition_pairs, rate_value in rate_groups
    ]
    total_log_prior = math.fsum(row.log_prior_contribution for row in rows)
    return DiscreteTraitRatePriorEvaluationReport(
        model=resolved_model,
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


def _resolve_supported_discrete_trait_rate_model(model: str) -> str:
    resolved_model = resolve_discrete_model_name(model)
    if resolved_model not in DISCRETE_TRAIT_RATE_PRIOR_MODELS:
        raise PhylogeneticsError(
            "discrete-trait rate priors support only ER, SYM, and ARD transition-rate surfaces",
            code="discrete_trait_rate_prior_model_unsupported",
            details={
                "model": model,
                "resolved_model": resolved_model,
                "allowed_models": list(DISCRETE_TRAIT_RATE_PRIOR_MODELS),
            },
        )
    return resolved_model


def _resolve_rate_parameter_groups(
    *,
    model: str,
    transition_rate_rows: Sequence[DiscreteTransitionRateRow],
) -> list[tuple[str, list[tuple[str, str]], float]]:
    allowed_rows = [
        row for row in transition_rate_rows if row.transition_allowed
    ]
    if not allowed_rows:
        raise PhylogeneticsError(
            "discrete-trait rate prior evaluation requires at least one allowed transition rate",
            code="discrete_trait_rate_prior_no_allowed_rates",
        )
    validated_rows = [
        (
            row.source_state,
            row.target_state,
            _validate_positive_finite_value(
                parameter_name=f"{row.source_state}->{row.target_state} rate",
                value=row.rate,
                owner_name="discrete-trait rate prior evaluation",
            ),
        )
        for row in allowed_rows
    ]
    if model == "equal-rates":
        shared_rate = validated_rows[0][2]
        if any(
            not _rates_close(rate_value, shared_rate)
            for _source_state, _target_state, rate_value in validated_rows[1:]
        ):
            raise PhylogeneticsError(
                "equal-rates prior evaluation requires one shared transition rate",
                code="discrete_trait_rate_prior_equal_rates_inconsistent",
            )
        return [
            (
                "shared-rate",
                sorted(
                    (source_state, target_state)
                    for source_state, target_state, _rate_value in validated_rows
                ),
                shared_rate,
            )
        ]
    if model == "symmetric":
        pair_groups: dict[tuple[str, str], list[tuple[str, str, float]]] = {}
        for source_state, target_state, rate_value in validated_rows:
            pair_groups.setdefault(
                tuple(sorted((source_state, target_state))),
                [],
            ).append((source_state, target_state, rate_value))
        resolved_groups: list[tuple[str, list[tuple[str, str]], float]] = []
        for pair_key in sorted(pair_groups):
            grouped_rows = pair_groups[pair_key]
            transition_pairs = sorted(
                (source_state, target_state)
                for source_state, target_state, _rate_value in grouped_rows
            )
            if len(grouped_rows) != 2:
                raise PhylogeneticsError(
                    "symmetric prior evaluation requires one forward and one reverse transition per rate parameter",
                    code="discrete_trait_rate_prior_symmetric_pair_incomplete",
                    details={"pair": list(pair_key)},
                )
            first_rate = grouped_rows[0][2]
            if any(
                not _rates_close(rate_value, first_rate)
                for _source_state, _target_state, rate_value in grouped_rows[1:]
            ):
                raise PhylogeneticsError(
                    "symmetric prior evaluation requires matched forward and reverse transition rates",
                    code="discrete_trait_rate_prior_symmetric_pair_mismatched",
                    details={"pair": list(pair_key)},
                )
            resolved_groups.append(
                (
                    f"{pair_key[0]}<->{pair_key[1]}",
                    transition_pairs,
                    first_rate,
                )
            )
        return resolved_groups
    return [
        (
            f"{source_state}->{target_state}",
            [(source_state, target_state)],
            rate_value,
        )
        for source_state, target_state, rate_value in sorted(validated_rows)
    ]


def _rates_close(left: float, right: float) -> bool:
    return math.isclose(
        left,
        right,
        rel_tol=_RATE_GROUP_TOLERANCE,
        abs_tol=_RATE_GROUP_TOLERANCE,
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
