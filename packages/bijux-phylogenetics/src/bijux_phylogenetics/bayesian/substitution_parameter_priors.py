from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math

import numpy

from bijux_phylogenetics.bayesian.required_values import require_present
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    DNA_EXCHANGEABILITY_LABELS,
    parameterize_dna_base_frequency_simplex,
    parameterize_dna_exchangeability_simplex,
)
from bijux_phylogenetics.phylo.likelihood.gamma import validate_discrete_gamma_alpha
from bijux_phylogenetics.phylo.likelihood.invariant import (
    validate_invariant_proportion,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

POSITIVE_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES = (
    "exponential",
    "fixed",
    "gamma",
    "lognormal",
)
PROBABILITY_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES = ("beta", "fixed")
SIMPLEX_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES = ("dirichlet", "fixed")
SUBSTITUTION_PARAMETER_PRIOR_TARGETS = (
    "base-frequencies",
    "exchangeabilities",
    "gamma-alpha",
    "invariant-proportion",
    "kappa",
)
_FIXED_PRIOR_TOLERANCE = 1e-12
_DNA_BASE_FREQUENCY_COMPONENT_NAMES = ("A", "C", "G", "T")


@dataclass(frozen=True, slots=True)
class PositiveSubstitutionParameterPriorModel:
    """One validated prior over one strictly positive substitution parameter."""

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


@dataclass(frozen=True, slots=True)
class ProbabilitySubstitutionParameterPriorModel:
    """One validated prior over one probability-valued substitution parameter."""

    family: str
    alpha: float | None = None
    beta: float | None = None
    fixed_value: float | None = None
    fixed_tolerance: float | None = None

    def parameter_values(self) -> dict[str, float]:
        parameter_values: dict[str, float] = {}
        if self.alpha is not None:
            parameter_values["alpha"] = self.alpha
        if self.beta is not None:
            parameter_values["beta"] = self.beta
        if self.fixed_value is not None:
            parameter_values["fixed_value"] = self.fixed_value
        if self.fixed_tolerance is not None:
            parameter_values["fixed_tolerance"] = self.fixed_tolerance
        return parameter_values


@dataclass(frozen=True, slots=True)
class SimplexSubstitutionParameterPriorModel:
    """One validated prior over one normalized simplex-valued substitution parameter."""

    family: str
    component_names: tuple[str, ...]
    concentration_parameters: tuple[float, ...] | None = None
    fixed_values: tuple[float, ...] | None = None
    fixed_tolerance: float | None = None

    def parameter_values(self) -> dict[str, float]:
        parameter_values: dict[str, float] = {}
        if self.concentration_parameters is not None:
            parameter_values.update(
                {
                    f"alpha_{component_name}": concentration_parameter
                    for component_name, concentration_parameter in zip(
                        self.component_names,
                        self.concentration_parameters,
                        strict=True,
                    )
                }
            )
        if self.fixed_values is not None:
            parameter_values.update(
                {
                    f"fixed_{component_name}": fixed_value
                    for component_name, fixed_value in zip(
                        self.component_names,
                        self.fixed_values,
                        strict=True,
                    )
                }
            )
        if self.fixed_tolerance is not None:
            parameter_values["fixed_tolerance"] = self.fixed_tolerance
        return parameter_values


@dataclass(frozen=True, slots=True)
class SubstitutionParameterPriorBundle:
    """One explicit prior bundle over supported substitution-parameter targets."""

    kappa_prior: PositiveSubstitutionParameterPriorModel | None = None
    exchangeability_prior: SimplexSubstitutionParameterPriorModel | None = None
    base_frequency_prior: SimplexSubstitutionParameterPriorModel | None = None
    gamma_alpha_prior: PositiveSubstitutionParameterPriorModel | None = None
    invariant_proportion_prior: ProbabilitySubstitutionParameterPriorModel | None = None

    def prior_count(self) -> int:
        return sum(
            prior is not None
            for prior in (
                self.kappa_prior,
                self.exchangeability_prior,
                self.base_frequency_prior,
                self.gamma_alpha_prior,
                self.invariant_proportion_prior,
            )
        )


@dataclass(frozen=True, slots=True)
class SubstitutionParameterPriorRow:
    """One target-level substitution-parameter prior contribution."""

    target_name: str
    family: str
    component_values: dict[str, float]
    hyperparameter_values: dict[str, float]
    log_prior_contribution: float


@dataclass(frozen=True, slots=True)
class SubstitutionParameterPriorEvaluationReport:
    """One substitution-parameter prior evaluation report."""

    prior_count: int
    total_log_prior: float
    rows: list[SubstitutionParameterPriorRow]


def build_exponential_positive_substitution_parameter_prior(
    *,
    rate: float,
) -> PositiveSubstitutionParameterPriorModel:
    """Build one exponential prior on one positive substitution parameter."""
    return PositiveSubstitutionParameterPriorModel(
        family="exponential",
        rate=_validate_positive_finite_value(
            parameter_name="rate",
            value=rate,
            owner_name="exponential substitution-parameter prior",
        ),
    )


def build_gamma_positive_substitution_parameter_prior(
    *,
    shape: float,
    scale: float,
) -> PositiveSubstitutionParameterPriorModel:
    """Build one gamma prior on one positive substitution parameter."""
    return PositiveSubstitutionParameterPriorModel(
        family="gamma",
        shape=_validate_positive_finite_value(
            parameter_name="shape",
            value=shape,
            owner_name="gamma substitution-parameter prior",
        ),
        scale=_validate_positive_finite_value(
            parameter_name="scale",
            value=scale,
            owner_name="gamma substitution-parameter prior",
        ),
    )


def build_lognormal_positive_substitution_parameter_prior(
    *,
    log_mean: float,
    log_standard_deviation: float,
) -> PositiveSubstitutionParameterPriorModel:
    """Build one lognormal prior on one positive substitution parameter."""
    return PositiveSubstitutionParameterPriorModel(
        family="lognormal",
        log_mean=_validate_finite_value(
            parameter_name="log_mean",
            value=log_mean,
            owner_name="lognormal substitution-parameter prior",
        ),
        log_standard_deviation=_validate_positive_finite_value(
            parameter_name="log_standard_deviation",
            value=log_standard_deviation,
            owner_name="lognormal substitution-parameter prior",
        ),
    )


def build_fixed_positive_substitution_parameter_prior(
    *,
    fixed_value: float,
    fixed_tolerance: float = _FIXED_PRIOR_TOLERANCE,
) -> PositiveSubstitutionParameterPriorModel:
    """Build one fixed-value prior on one positive substitution parameter."""
    return PositiveSubstitutionParameterPriorModel(
        family="fixed",
        fixed_value=_validate_positive_finite_value(
            parameter_name="fixed_value",
            value=fixed_value,
            owner_name="fixed positive substitution-parameter prior",
        ),
        fixed_tolerance=_validate_nonnegative_finite_value(
            parameter_name="fixed_tolerance",
            value=fixed_tolerance,
            owner_name="fixed positive substitution-parameter prior",
        ),
    )


def build_beta_probability_substitution_parameter_prior(
    *,
    alpha: float,
    beta: float,
) -> ProbabilitySubstitutionParameterPriorModel:
    """Build one beta prior on one probability-valued substitution parameter."""
    return ProbabilitySubstitutionParameterPriorModel(
        family="beta",
        alpha=_validate_positive_finite_value(
            parameter_name="alpha",
            value=alpha,
            owner_name="beta substitution-parameter prior",
        ),
        beta=_validate_positive_finite_value(
            parameter_name="beta",
            value=beta,
            owner_name="beta substitution-parameter prior",
        ),
    )


def build_fixed_probability_substitution_parameter_prior(
    *,
    fixed_value: float,
    fixed_tolerance: float = _FIXED_PRIOR_TOLERANCE,
) -> ProbabilitySubstitutionParameterPriorModel:
    """Build one fixed-value prior on one probability-valued substitution parameter."""
    validated_fixed_value = validate_invariant_proportion(
        fixed_value,
        model_name="fixed probability substitution-parameter prior",
    )
    return ProbabilitySubstitutionParameterPriorModel(
        family="fixed",
        fixed_value=validated_fixed_value,
        fixed_tolerance=_validate_nonnegative_finite_value(
            parameter_name="fixed_tolerance",
            value=fixed_tolerance,
            owner_name="fixed probability substitution-parameter prior",
        ),
    )


def build_dirichlet_simplex_substitution_parameter_prior(
    *,
    expected_component_names: Sequence[str],
    concentration_parameters: Mapping[str, float] | Sequence[float],
) -> SimplexSubstitutionParameterPriorModel:
    """Build one Dirichlet prior on one normalized simplex-valued substitution parameter."""
    component_names = tuple(expected_component_names)
    concentration_values = _resolve_named_component_values(
        component_names=component_names,
        raw_values=concentration_parameters,
        owner_name="Dirichlet simplex substitution-parameter prior",
    )
    return SimplexSubstitutionParameterPriorModel(
        family="dirichlet",
        component_names=component_names,
        concentration_parameters=tuple(
            _validate_positive_finite_value(
                parameter_name=f"concentration_parameters[{component_name!r}]",
                value=value,
                owner_name="Dirichlet simplex substitution-parameter prior",
            )
            for component_name, value in zip(
                component_names,
                concentration_values,
                strict=True,
            )
        ),
    )


def build_fixed_simplex_substitution_parameter_prior(
    *,
    expected_component_names: Sequence[str],
    fixed_values: Mapping[str, float] | Sequence[float],
    fixed_tolerance: float = _FIXED_PRIOR_TOLERANCE,
) -> SimplexSubstitutionParameterPriorModel:
    """Build one fixed prior on one normalized simplex-valued substitution parameter."""
    component_names = tuple(expected_component_names)
    validated_fixed_values = _validate_simplex_values(
        component_names=component_names,
        raw_values=fixed_values,
        owner_name="fixed simplex substitution-parameter prior",
    )
    return SimplexSubstitutionParameterPriorModel(
        family="fixed",
        component_names=component_names,
        fixed_values=validated_fixed_values,
        fixed_tolerance=_validate_nonnegative_finite_value(
            parameter_name="fixed_tolerance",
            value=fixed_tolerance,
            owner_name="fixed simplex substitution-parameter prior",
        ),
    )


def build_substitution_parameter_prior_bundle(
    *,
    kappa_prior: PositiveSubstitutionParameterPriorModel | None = None,
    exchangeability_prior: SimplexSubstitutionParameterPriorModel | None = None,
    base_frequency_prior: SimplexSubstitutionParameterPriorModel | None = None,
    gamma_alpha_prior: PositiveSubstitutionParameterPriorModel | None = None,
    invariant_proportion_prior: ProbabilitySubstitutionParameterPriorModel
    | None = None,
) -> SubstitutionParameterPriorBundle:
    """Build one explicit prior bundle over supported substitution parameters."""
    if exchangeability_prior is not None:
        _require_expected_component_names(
            prior_model=exchangeability_prior,
            expected_component_names=DNA_EXCHANGEABILITY_LABELS,
            target_name="exchangeabilities",
        )
    if base_frequency_prior is not None:
        _require_expected_component_names(
            prior_model=base_frequency_prior,
            expected_component_names=_DNA_BASE_FREQUENCY_COMPONENT_NAMES,
            target_name="base frequencies",
        )
    bundle = SubstitutionParameterPriorBundle(
        kappa_prior=kappa_prior,
        exchangeability_prior=exchangeability_prior,
        base_frequency_prior=base_frequency_prior,
        gamma_alpha_prior=gamma_alpha_prior,
        invariant_proportion_prior=invariant_proportion_prior,
    )
    if bundle.prior_count() == 0:
        raise PhylogeneticsError(
            "substitution-parameter prior bundle requires at least one explicit prior",
            code="substitution_parameter_prior_bundle_empty",
        )
    return bundle


def evaluate_substitution_parameter_log_prior(
    *,
    prior_bundle: SubstitutionParameterPriorBundle,
    kappa: float | None = None,
    exchangeabilities: (
        Mapping[tuple[str, str], float]
        | Mapping[str, float]
        | numpy.ndarray
        | Sequence[float]
        | None
    ) = None,
    base_frequencies: (
        Mapping[str, float] | numpy.ndarray | Sequence[float] | None
    ) = None,
    gamma_alpha: float | None = None,
    invariant_proportion: float | None = None,
) -> SubstitutionParameterPriorEvaluationReport:
    """Evaluate one substitution-parameter prior bundle on one explicit parameter set."""
    rows: list[SubstitutionParameterPriorRow] = []
    if prior_bundle.kappa_prior is not None:
        validated_kappa = _validate_positive_parameter_value(
            parameter_name="kappa",
            value=kappa,
        )
        rows.append(
            _evaluate_positive_prior_row(
                target_name="kappa",
                realized_value=validated_kappa,
                prior_model=prior_bundle.kappa_prior,
            )
        )
    if prior_bundle.exchangeability_prior is not None:
        if exchangeabilities is None:
            raise PhylogeneticsError(
                "substitution-parameter prior evaluation requires exchangeabilities when one prior is configured",
                code="substitution_parameter_prior_missing_exchangeabilities",
            )
        normalized_exchangeability_simplex = parameterize_dna_exchangeability_simplex(
            exchangeabilities
        )
        rows.append(
            _evaluate_simplex_prior_row(
                target_name="exchangeabilities",
                component_names=DNA_EXCHANGEABILITY_LABELS,
                realized_values=normalized_exchangeability_simplex.constrained_mapping(),
                prior_model=prior_bundle.exchangeability_prior,
            )
        )
    if prior_bundle.base_frequency_prior is not None:
        if base_frequencies is None:
            raise PhylogeneticsError(
                "substitution-parameter prior evaluation requires base frequencies when one prior is configured",
                code="substitution_parameter_prior_missing_base_frequencies",
            )
        normalized_base_frequency_simplex = parameterize_dna_base_frequency_simplex(
            base_frequencies
        )
        rows.append(
            _evaluate_simplex_prior_row(
                target_name="base-frequencies",
                component_names=_DNA_BASE_FREQUENCY_COMPONENT_NAMES,
                realized_values=normalized_base_frequency_simplex.constrained_mapping(),
                prior_model=prior_bundle.base_frequency_prior,
            )
        )
    if prior_bundle.gamma_alpha_prior is not None:
        if gamma_alpha is None:
            raise PhylogeneticsError(
                "substitution-parameter prior evaluation requires gamma_alpha when one prior is configured",
                code="substitution_parameter_prior_missing_gamma_alpha",
            )
        validated_gamma_alpha = validate_discrete_gamma_alpha(gamma_alpha)
        rows.append(
            _evaluate_positive_prior_row(
                target_name="gamma-alpha",
                realized_value=validated_gamma_alpha,
                prior_model=prior_bundle.gamma_alpha_prior,
            )
        )
    if prior_bundle.invariant_proportion_prior is not None:
        if invariant_proportion is None:
            raise PhylogeneticsError(
                "substitution-parameter prior evaluation requires invariant_proportion when one prior is configured",
                code="substitution_parameter_prior_missing_invariant_proportion",
            )
        validated_invariant_proportion = validate_invariant_proportion(
            invariant_proportion,
            model_name="substitution-parameter prior evaluation",
        )
        rows.append(
            _evaluate_probability_prior_row(
                target_name="invariant-proportion",
                realized_value=validated_invariant_proportion,
                prior_model=prior_bundle.invariant_proportion_prior,
            )
        )
    total_log_prior = math.fsum(row.log_prior_contribution for row in rows)
    return SubstitutionParameterPriorEvaluationReport(
        prior_count=prior_bundle.prior_count(),
        total_log_prior=float(format(total_log_prior, ".15g")),
        rows=rows,
    )


def _evaluate_positive_prior_row(
    *,
    target_name: str,
    realized_value: float,
    prior_model: PositiveSubstitutionParameterPriorModel,
) -> SubstitutionParameterPriorRow:
    log_prior_contribution = _evaluate_positive_log_prior(
        value=realized_value,
        prior_model=prior_model,
    )
    return SubstitutionParameterPriorRow(
        target_name=target_name,
        family=prior_model.family,
        component_values={"value": float(format(realized_value, ".15g"))},
        hyperparameter_values=prior_model.parameter_values(),
        log_prior_contribution=float(format(log_prior_contribution, ".15g")),
    )


def _evaluate_probability_prior_row(
    *,
    target_name: str,
    realized_value: float,
    prior_model: ProbabilitySubstitutionParameterPriorModel,
) -> SubstitutionParameterPriorRow:
    log_prior_contribution = _evaluate_probability_log_prior(
        value=realized_value,
        prior_model=prior_model,
    )
    return SubstitutionParameterPriorRow(
        target_name=target_name,
        family=prior_model.family,
        component_values={"value": float(format(realized_value, ".15g"))},
        hyperparameter_values=prior_model.parameter_values(),
        log_prior_contribution=float(format(log_prior_contribution, ".15g")),
    )


def _evaluate_simplex_prior_row(
    *,
    target_name: str,
    component_names: Sequence[str],
    realized_values: Mapping[str, float],
    prior_model: SimplexSubstitutionParameterPriorModel,
) -> SubstitutionParameterPriorRow:
    log_prior_contribution = _evaluate_simplex_log_prior(
        component_names=component_names,
        realized_values=realized_values,
        prior_model=prior_model,
    )
    return SubstitutionParameterPriorRow(
        target_name=target_name,
        family=prior_model.family,
        component_values={
            component_name: float(format(realized_values[component_name], ".15g"))
            for component_name in component_names
        },
        hyperparameter_values=prior_model.parameter_values(),
        log_prior_contribution=float(format(log_prior_contribution, ".15g")),
    )


def _evaluate_positive_log_prior(
    *,
    value: float,
    prior_model: PositiveSubstitutionParameterPriorModel,
) -> float:
    if prior_model.family == "exponential":
        rate = require_present(
            prior_model.rate,
            owner_name="positive substitution-parameter prior evaluation",
            field_name="rate",
        )
        return math.log(rate) - (rate * value)
    if prior_model.family == "gamma":
        shape = require_present(
            prior_model.shape,
            owner_name="positive substitution-parameter prior evaluation",
            field_name="shape",
        )
        scale = require_present(
            prior_model.scale,
            owner_name="positive substitution-parameter prior evaluation",
            field_name="scale",
        )
        return _gamma_log_density(value, shape=shape, scale=scale)
    if prior_model.family == "lognormal":
        log_mean = require_present(
            prior_model.log_mean,
            owner_name="positive substitution-parameter prior evaluation",
            field_name="log_mean",
        )
        log_standard_deviation = require_present(
            prior_model.log_standard_deviation,
            owner_name="positive substitution-parameter prior evaluation",
            field_name="log_standard_deviation",
        )
        return _lognormal_log_density(
            value,
            log_mean=log_mean,
            log_standard_deviation=log_standard_deviation,
        )
    if prior_model.family == "fixed":
        fixed_value = require_present(
            prior_model.fixed_value,
            owner_name="positive substitution-parameter prior evaluation",
            field_name="fixed_value",
        )
        fixed_tolerance = require_present(
            prior_model.fixed_tolerance,
            owner_name="positive substitution-parameter prior evaluation",
            field_name="fixed_tolerance",
        )
        return (
            0.0
            if math.isclose(
                value,
                fixed_value,
                rel_tol=0.0,
                abs_tol=fixed_tolerance,
            )
            else -math.inf
        )
    raise PhylogeneticsError(
        "positive substitution-parameter prior family is unsupported",
        code="positive_substitution_parameter_prior_family_invalid",
        details={
            "family": prior_model.family,
            "allowed_families": list(POSITIVE_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES),
        },
    )


def _evaluate_probability_log_prior(
    *,
    value: float,
    prior_model: ProbabilitySubstitutionParameterPriorModel,
) -> float:
    if prior_model.family == "beta":
        alpha = require_present(
            prior_model.alpha,
            owner_name="probability substitution-parameter prior evaluation",
            field_name="alpha",
        )
        beta = require_present(
            prior_model.beta,
            owner_name="probability substitution-parameter prior evaluation",
            field_name="beta",
        )
        return _beta_log_density(value, alpha=alpha, beta=beta)
    if prior_model.family == "fixed":
        fixed_value = require_present(
            prior_model.fixed_value,
            owner_name="probability substitution-parameter prior evaluation",
            field_name="fixed_value",
        )
        fixed_tolerance = require_present(
            prior_model.fixed_tolerance,
            owner_name="probability substitution-parameter prior evaluation",
            field_name="fixed_tolerance",
        )
        return (
            0.0
            if math.isclose(
                value,
                fixed_value,
                rel_tol=0.0,
                abs_tol=fixed_tolerance,
            )
            else -math.inf
        )
    raise PhylogeneticsError(
        "probability substitution-parameter prior family is unsupported",
        code="probability_substitution_parameter_prior_family_invalid",
        details={
            "family": prior_model.family,
            "allowed_families": list(PROBABILITY_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES),
        },
    )


def _evaluate_simplex_log_prior(
    *,
    component_names: Sequence[str],
    realized_values: Mapping[str, float],
    prior_model: SimplexSubstitutionParameterPriorModel,
) -> float:
    ordered_realized_values = _validate_simplex_values(
        component_names=tuple(component_names),
        raw_values=realized_values,
        owner_name="substitution-parameter prior evaluation",
    )
    if prior_model.family == "dirichlet":
        concentration_parameters = require_present(
            prior_model.concentration_parameters,
            owner_name="simplex substitution-parameter prior evaluation",
            field_name="concentration_parameters",
        )
        return _dirichlet_log_density(
            values=ordered_realized_values,
            concentration_parameters=concentration_parameters,
        )
    if prior_model.family == "fixed":
        fixed_values = require_present(
            prior_model.fixed_values,
            owner_name="simplex substitution-parameter prior evaluation",
            field_name="fixed_values",
        )
        fixed_tolerance = require_present(
            prior_model.fixed_tolerance,
            owner_name="simplex substitution-parameter prior evaluation",
            field_name="fixed_tolerance",
        )
        return (
            0.0
            if all(
                math.isclose(
                    observed_value,
                    expected_value,
                    rel_tol=0.0,
                    abs_tol=fixed_tolerance,
                )
                for observed_value, expected_value in zip(
                    ordered_realized_values,
                    fixed_values,
                    strict=True,
                )
            )
            else -math.inf
        )
    raise PhylogeneticsError(
        "simplex substitution-parameter prior family is unsupported",
        code="simplex_substitution_parameter_prior_family_invalid",
        details={
            "family": prior_model.family,
            "allowed_families": list(SIMPLEX_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES),
        },
    )


def _validate_positive_parameter_value(
    *,
    parameter_name: str,
    value: float | None,
) -> float:
    if value is None:
        raise PhylogeneticsError(
            f"substitution-parameter prior evaluation requires {parameter_name}",
            code="substitution_parameter_prior_missing_positive_parameter",
            details={"parameter_name": parameter_name},
        )
    return _validate_positive_finite_value(
        parameter_name=parameter_name,
        value=value,
        owner_name="substitution-parameter prior evaluation",
    )


def _validate_simplex_values(
    *,
    component_names: tuple[str, ...],
    raw_values: Mapping[str, float] | Sequence[float],
    owner_name: str,
) -> tuple[float, ...]:
    if component_names == _DNA_BASE_FREQUENCY_COMPONENT_NAMES:
        return parameterize_dna_base_frequency_simplex(raw_values).constrained_values
    if component_names == DNA_EXCHANGEABILITY_LABELS:
        return parameterize_dna_exchangeability_simplex(raw_values).constrained_values
    raise PhylogeneticsError(
        f"{owner_name} does not support simplex components {component_names!r}",
        code="substitution_parameter_prior_simplex_components_unsupported",
        details={"component_names": list(component_names)},
    )


def _resolve_named_component_values(
    *,
    component_names: tuple[str, ...],
    raw_values: Mapping[str, float] | Sequence[float],
    owner_name: str,
) -> tuple[float, ...]:
    if isinstance(raw_values, Mapping):
        if set(raw_values) != set(component_names):
            raise PhylogeneticsError(
                f"{owner_name} requires values for exactly {list(component_names)}",
                code="substitution_parameter_prior_component_names_invalid",
                details={
                    "expected_component_names": list(component_names),
                    "observed_component_names": sorted(
                        str(name) for name in raw_values
                    ),
                },
            )
        return tuple(
            float(raw_values[component_name]) for component_name in component_names
        )
    resolved_values = tuple(float(value) for value in raw_values)
    if len(resolved_values) != len(component_names):
        raise PhylogeneticsError(
            f"{owner_name} requires exactly {len(component_names)} component values",
            code="substitution_parameter_prior_component_dimension_invalid",
            details={
                "expected_dimension": len(component_names),
                "observed_dimension": len(resolved_values),
            },
        )
    return resolved_values


def _require_expected_component_names(
    *,
    prior_model: SimplexSubstitutionParameterPriorModel,
    expected_component_names: Sequence[str],
    target_name: str,
) -> None:
    if prior_model.component_names != tuple(expected_component_names):
        raise PhylogeneticsError(
            f"substitution-parameter prior bundle requires one supported component order for {target_name}",
            code="substitution_parameter_prior_component_order_invalid",
            details={
                "target_name": target_name,
                "expected_component_names": list(expected_component_names),
                "observed_component_names": list(prior_model.component_names),
            },
        )


def _gamma_log_density(value: float, *, shape: float, scale: float) -> float:
    if value == 0.0:
        if shape < 1.0:
            return math.inf
        if shape == 1.0:
            return -math.log(scale)
        return -math.inf
    return (
        (shape - 1.0) * math.log(value)
        - (value / scale)
        - math.lgamma(shape)
        - (shape * math.log(scale))
    )


def _lognormal_log_density(
    value: float,
    *,
    log_mean: float,
    log_standard_deviation: float,
) -> float:
    return (
        -math.log(value)
        - math.log(log_standard_deviation)
        - 0.5 * math.log(2.0 * math.pi)
        - ((math.log(value) - log_mean) ** 2) / (2.0 * (log_standard_deviation**2))
    )


def _beta_log_density(value: float, *, alpha: float, beta: float) -> float:
    if value == 0.0:
        if alpha < 1.0:
            return math.inf
        if alpha == 1.0:
            return math.log(beta)
        return -math.inf
    return (
        math.lgamma(alpha + beta)
        - math.lgamma(alpha)
        - math.lgamma(beta)
        + ((alpha - 1.0) * math.log(value))
        + ((beta - 1.0) * math.log1p(-value))
    )


def _dirichlet_log_density(
    *,
    values: tuple[float, ...],
    concentration_parameters: tuple[float, ...],
) -> float:
    return (
        math.lgamma(math.fsum(concentration_parameters))
        - math.fsum(math.lgamma(alpha) for alpha in concentration_parameters)
        + math.fsum(
            (alpha - 1.0) * math.log(value)
            for alpha, value in zip(concentration_parameters, values, strict=True)
        )
    )


def _validate_positive_finite_value(
    *,
    parameter_name: str,
    value: float,
    owner_name: str,
) -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires a strictly positive finite {parameter_name}",
            code="substitution_parameter_prior_value_invalid",
            details={
                "parameter_name": parameter_name,
                "value": value,
                "owner_name": owner_name,
            },
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
            f"{owner_name} requires a non-negative finite {parameter_name}",
            code="substitution_parameter_prior_value_invalid",
            details={
                "parameter_name": parameter_name,
                "value": value,
                "owner_name": owner_name,
            },
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
            f"{owner_name} requires a finite {parameter_name}",
            code="substitution_parameter_prior_value_invalid",
            details={
                "parameter_name": parameter_name,
                "value": value,
                "owner_name": owner_name,
            },
        )
    return float(format(value, ".15g"))


__all__ = [
    "POSITIVE_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES",
    "PROBABILITY_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES",
    "SIMPLEX_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES",
    "SUBSTITUTION_PARAMETER_PRIOR_TARGETS",
    "PositiveSubstitutionParameterPriorModel",
    "ProbabilitySubstitutionParameterPriorModel",
    "SimplexSubstitutionParameterPriorModel",
    "SubstitutionParameterPriorBundle",
    "SubstitutionParameterPriorEvaluationReport",
    "SubstitutionParameterPriorRow",
    "build_beta_probability_substitution_parameter_prior",
    "build_dirichlet_simplex_substitution_parameter_prior",
    "build_exponential_positive_substitution_parameter_prior",
    "build_fixed_positive_substitution_parameter_prior",
    "build_fixed_probability_substitution_parameter_prior",
    "build_fixed_simplex_substitution_parameter_prior",
    "build_gamma_positive_substitution_parameter_prior",
    "build_lognormal_positive_substitution_parameter_prior",
    "build_substitution_parameter_prior_bundle",
    "evaluate_substitution_parameter_log_prior",
]
