from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.bayesian.required_values import require_present
from bijux_phylogenetics.comparative.evolutionary_modes import (
    ContinuousEvolutionaryModeFitReport,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

CONTINUOUS_TRAIT_SCALAR_PRIOR_FAMILIES = (
    "exponential",
    "fixed",
    "gamma",
    "lognormal",
)
CONTINUOUS_TRAIT_PROBABILITY_PRIOR_FAMILIES = ("beta", "fixed")
CONTINUOUS_TRAIT_PRIOR_MODES = (
    "brownian",
    "white-noise",
    "pagel-lambda",
    "pagel-kappa",
    "pagel-delta",
    "ornstein-uhlenbeck",
    "early-burst",
)
CONTINUOUS_TRAIT_PRIOR_TARGETS = (
    "alpha",
    "delta",
    "kappa",
    "lambda",
    "rate_change",
    "sigma",
)
_FIXED_PRIOR_TOLERANCE = 1e-12
_MODE_PARAMETER_TARGETS: dict[str, str | None] = {
    "brownian": None,
    "white-noise": None,
    "pagel-lambda": "lambda",
    "pagel-kappa": "kappa",
    "pagel-delta": "delta",
    "ornstein-uhlenbeck": "alpha",
    "early-burst": "rate_change",
}


@dataclass(frozen=True, slots=True)
class ContinuousTraitScalarPriorModel:
    """One validated scalar prior over one continuous-trait model parameter."""

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
class ContinuousTraitProbabilityPriorModel:
    """One validated probability prior over one bounded continuous-trait parameter."""

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
class ContinuousTraitModelPriorBundle:
    """One explicit prior bundle over the owned continuous-trait mode parameters."""

    sigma_prior: ContinuousTraitScalarPriorModel | None = None
    lambda_prior: ContinuousTraitProbabilityPriorModel | None = None
    kappa_prior: ContinuousTraitScalarPriorModel | None = None
    delta_prior: ContinuousTraitScalarPriorModel | None = None
    alpha_prior: ContinuousTraitScalarPriorModel | None = None
    rate_change_prior: ContinuousTraitScalarPriorModel | None = None

    def prior_count(self) -> int:
        return sum(
            prior is not None
            for prior in (
                self.sigma_prior,
                self.lambda_prior,
                self.kappa_prior,
                self.delta_prior,
                self.alpha_prior,
                self.rate_change_prior,
            )
        )


@dataclass(frozen=True, slots=True)
class ContinuousTraitModelPriorRow:
    """One continuous-trait model prior contribution."""

    target_name: str
    fit_parameter_source: str
    family: str
    parameter_value: float
    hyperparameter_values: dict[str, float]
    log_prior_contribution: float


@dataclass(frozen=True, slots=True)
class ContinuousTraitModelPriorEvaluationReport:
    """One evaluated prior bundle over one owned continuous-trait fit report."""

    tree_path: str
    traits_path: str
    trait: str
    mode: str
    log_likelihood: float
    total_log_prior: float
    posterior_score: float
    prior_count: int
    rows: list[ContinuousTraitModelPriorRow]


def build_exponential_continuous_trait_scalar_prior(
    *,
    rate: float,
) -> ContinuousTraitScalarPriorModel:
    """Build one exponential prior over one scalar continuous-trait parameter."""
    return ContinuousTraitScalarPriorModel(
        family="exponential",
        rate=_validate_positive_finite_value(
            parameter_name="rate",
            value=rate,
            owner_name="exponential continuous-trait scalar prior",
        ),
    )


def build_gamma_continuous_trait_scalar_prior(
    *,
    shape: float,
    scale: float,
) -> ContinuousTraitScalarPriorModel:
    """Build one gamma prior over one scalar continuous-trait parameter."""
    return ContinuousTraitScalarPriorModel(
        family="gamma",
        shape=_validate_positive_finite_value(
            parameter_name="shape",
            value=shape,
            owner_name="gamma continuous-trait scalar prior",
        ),
        scale=_validate_positive_finite_value(
            parameter_name="scale",
            value=scale,
            owner_name="gamma continuous-trait scalar prior",
        ),
    )


def build_lognormal_continuous_trait_scalar_prior(
    *,
    log_mean: float,
    log_standard_deviation: float,
) -> ContinuousTraitScalarPriorModel:
    """Build one lognormal prior over one scalar continuous-trait parameter."""
    return ContinuousTraitScalarPriorModel(
        family="lognormal",
        log_mean=_validate_finite_value(
            parameter_name="log_mean",
            value=log_mean,
            owner_name="lognormal continuous-trait scalar prior",
        ),
        log_standard_deviation=_validate_positive_finite_value(
            parameter_name="log_standard_deviation",
            value=log_standard_deviation,
            owner_name="lognormal continuous-trait scalar prior",
        ),
    )


def build_fixed_continuous_trait_scalar_prior(
    *,
    fixed_value: float,
    fixed_tolerance: float = _FIXED_PRIOR_TOLERANCE,
) -> ContinuousTraitScalarPriorModel:
    """Build one fixed-value prior over one scalar continuous-trait parameter."""
    return ContinuousTraitScalarPriorModel(
        family="fixed",
        fixed_value=_validate_nonnegative_finite_value(
            parameter_name="fixed_value",
            value=fixed_value,
            owner_name="fixed continuous-trait scalar prior",
        ),
        fixed_tolerance=_validate_nonnegative_finite_value(
            parameter_name="fixed_tolerance",
            value=fixed_tolerance,
            owner_name="fixed continuous-trait scalar prior",
        ),
    )


def build_beta_continuous_trait_probability_prior(
    *,
    alpha: float,
    beta: float,
) -> ContinuousTraitProbabilityPriorModel:
    """Build one beta prior over one probability-bounded continuous-trait parameter."""
    return ContinuousTraitProbabilityPriorModel(
        family="beta",
        alpha=_validate_positive_finite_value(
            parameter_name="alpha",
            value=alpha,
            owner_name="beta continuous-trait probability prior",
        ),
        beta=_validate_positive_finite_value(
            parameter_name="beta",
            value=beta,
            owner_name="beta continuous-trait probability prior",
        ),
    )


def build_fixed_continuous_trait_probability_prior(
    *,
    fixed_value: float,
    fixed_tolerance: float = _FIXED_PRIOR_TOLERANCE,
) -> ContinuousTraitProbabilityPriorModel:
    """Build one fixed-value prior over one probability-bounded parameter."""
    return ContinuousTraitProbabilityPriorModel(
        family="fixed",
        fixed_value=_validate_probability_value(
            parameter_name="fixed_value",
            value=fixed_value,
            owner_name="fixed continuous-trait probability prior",
        ),
        fixed_tolerance=_validate_nonnegative_finite_value(
            parameter_name="fixed_tolerance",
            value=fixed_tolerance,
            owner_name="fixed continuous-trait probability prior",
        ),
    )


def build_continuous_trait_model_prior_bundle(
    *,
    sigma_prior: ContinuousTraitScalarPriorModel | None = None,
    lambda_prior: ContinuousTraitProbabilityPriorModel | None = None,
    kappa_prior: ContinuousTraitScalarPriorModel | None = None,
    delta_prior: ContinuousTraitScalarPriorModel | None = None,
    alpha_prior: ContinuousTraitScalarPriorModel | None = None,
    rate_change_prior: ContinuousTraitScalarPriorModel | None = None,
) -> ContinuousTraitModelPriorBundle:
    """Build one explicit prior bundle over the owned continuous-trait mode surface."""
    return ContinuousTraitModelPriorBundle(
        sigma_prior=sigma_prior,
        lambda_prior=lambda_prior,
        kappa_prior=kappa_prior,
        delta_prior=delta_prior,
        alpha_prior=alpha_prior,
        rate_change_prior=rate_change_prior,
    )


def evaluate_continuous_trait_scalar_log_prior(
    *,
    parameter_value: float,
    prior_model: ContinuousTraitScalarPriorModel,
    parameter_name: str,
    allow_zero: bool,
) -> float:
    """Evaluate one scalar continuous-trait parameter log prior."""
    validated_value = _validate_scalar_parameter_value(
        parameter_name=parameter_name,
        value=parameter_value,
        allow_zero=allow_zero,
    )
    if prior_model.family == "exponential":
        rate = require_present(
            prior_model.rate,
            owner_name="continuous-trait scalar prior evaluation",
            field_name="rate",
        )
        return math.log(rate) - (rate * validated_value)
    if prior_model.family == "gamma":
        shape = require_present(
            prior_model.shape,
            owner_name="continuous-trait scalar prior evaluation",
            field_name="shape",
        )
        scale = require_present(
            prior_model.scale,
            owner_name="continuous-trait scalar prior evaluation",
            field_name="scale",
        )
        return _gamma_log_density(
            validated_value,
            shape=shape,
            scale=scale,
            parameter_name=parameter_name,
            allow_zero=allow_zero,
        )
    if prior_model.family == "lognormal":
        if validated_value == 0.0:
            raise PhylogeneticsError(
                "lognormal continuous-trait scalar priors cannot score a boundary value of zero",
                code="continuous_trait_scalar_prior_zero_lognormal",
                details={"parameter_name": parameter_name},
            )
        log_mean = require_present(
            prior_model.log_mean,
            owner_name="continuous-trait scalar prior evaluation",
            field_name="log_mean",
        )
        log_standard_deviation = require_present(
            prior_model.log_standard_deviation,
            owner_name="continuous-trait scalar prior evaluation",
            field_name="log_standard_deviation",
        )
        return _lognormal_log_density(
            validated_value,
            log_mean=log_mean,
            log_standard_deviation=log_standard_deviation,
        )
    if prior_model.family == "fixed":
        fixed_value = require_present(
            prior_model.fixed_value,
            owner_name="continuous-trait scalar prior evaluation",
            field_name="fixed_value",
        )
        fixed_tolerance = require_present(
            prior_model.fixed_tolerance,
            owner_name="continuous-trait scalar prior evaluation",
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
        "continuous-trait scalar prior family is unsupported",
        code="continuous_trait_scalar_prior_family_invalid",
        details={
            "family": prior_model.family,
            "allowed_families": list(CONTINUOUS_TRAIT_SCALAR_PRIOR_FAMILIES),
        },
    )


def evaluate_continuous_trait_probability_log_prior(
    *,
    parameter_value: float,
    prior_model: ContinuousTraitProbabilityPriorModel,
    parameter_name: str,
) -> float:
    """Evaluate one bounded probability continuous-trait parameter log prior."""
    validated_value = _validate_probability_value(
        parameter_name=parameter_name,
        value=parameter_value,
        owner_name="continuous-trait probability prior evaluation",
    )
    if prior_model.family == "beta":
        alpha = require_present(
            prior_model.alpha,
            owner_name="continuous-trait probability prior evaluation",
            field_name="alpha",
        )
        beta = require_present(
            prior_model.beta,
            owner_name="continuous-trait probability prior evaluation",
            field_name="beta",
        )
        return _beta_log_density(
            validated_value,
            alpha=alpha,
            beta=beta,
            parameter_name=parameter_name,
        )
    if prior_model.family == "fixed":
        fixed_value = require_present(
            prior_model.fixed_value,
            owner_name="continuous-trait probability prior evaluation",
            field_name="fixed_value",
        )
        fixed_tolerance = require_present(
            prior_model.fixed_tolerance,
            owner_name="continuous-trait probability prior evaluation",
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
        "continuous-trait probability prior family is unsupported",
        code="continuous_trait_probability_prior_family_invalid",
        details={
            "family": prior_model.family,
            "allowed_families": list(CONTINUOUS_TRAIT_PROBABILITY_PRIOR_FAMILIES),
        },
    )


def evaluate_continuous_trait_model_log_prior(
    *,
    fit_report: ContinuousEvolutionaryModeFitReport,
    prior_bundle: ContinuousTraitModelPriorBundle,
) -> ContinuousTraitModelPriorEvaluationReport:
    """Evaluate one prior bundle on one owned continuous-trait fit report."""
    resolved_mode = _resolve_supported_continuous_trait_prior_mode(fit_report.mode)
    mode_parameter_rows = _resolve_mode_parameter_rows(
        fit_report=fit_report,
        prior_bundle=prior_bundle,
        mode=resolved_mode,
    )
    sigma_rows = _resolve_sigma_rows(fit_report=fit_report, prior_bundle=prior_bundle)
    rows = [*mode_parameter_rows, *sigma_rows]
    total_log_prior = math.fsum(row.log_prior_contribution for row in rows)
    posterior_score = fit_report.log_likelihood + total_log_prior
    return ContinuousTraitModelPriorEvaluationReport(
        tree_path=str(fit_report.tree_path),
        traits_path=str(fit_report.traits_path),
        trait=fit_report.trait,
        mode=resolved_mode,
        log_likelihood=float(format(fit_report.log_likelihood, ".15g")),
        total_log_prior=float(format(total_log_prior, ".15g")),
        posterior_score=float(format(posterior_score, ".15g")),
        prior_count=len(rows),
        rows=rows,
    )


def _resolve_mode_parameter_rows(
    *,
    fit_report: ContinuousEvolutionaryModeFitReport,
    prior_bundle: ContinuousTraitModelPriorBundle,
    mode: str,
) -> list[ContinuousTraitModelPriorRow]:
    target_name = _MODE_PARAMETER_TARGETS[mode]
    if target_name is None:
        if (
            fit_report.parameter_name is not None
            or fit_report.parameter_value is not None
        ):
            raise PhylogeneticsError(
                "non-parameterized continuous-trait modes must not carry an auxiliary parameter value",
                code="continuous_trait_prior_unexpected_mode_parameter",
                details={
                    "mode": mode,
                    "parameter_name": fit_report.parameter_name,
                    "parameter_value": fit_report.parameter_value,
                },
            )
        return []
    if fit_report.parameter_name != target_name or fit_report.parameter_value is None:
        raise PhylogeneticsError(
            "continuous-trait fit report parameter surface is inconsistent with its mode",
            code="continuous_trait_prior_mode_parameter_mismatch",
            details={
                "mode": mode,
                "expected_parameter_name": target_name,
                "actual_parameter_name": fit_report.parameter_name,
                "actual_parameter_value": fit_report.parameter_value,
            },
        )
    prior_model = getattr(prior_bundle, f"{target_name}_prior")
    if prior_model is None:
        return []
    if target_name == "lambda":
        log_prior = evaluate_continuous_trait_probability_log_prior(
            parameter_value=fit_report.parameter_value,
            prior_model=prior_model,
            parameter_name=target_name,
        )
    else:
        if not isinstance(prior_model, ContinuousTraitScalarPriorModel):
            raise PhylogeneticsError(
                "continuous-trait prior bundle target requires a scalar prior model",
                code="continuous_trait_scalar_prior_model_invalid",
                details={"target_name": target_name},
            )
        log_prior = evaluate_continuous_trait_scalar_log_prior(
            parameter_value=fit_report.parameter_value,
            prior_model=prior_model,
            parameter_name=target_name,
            allow_zero=True,
        )
    return [
        ContinuousTraitModelPriorRow(
            target_name=target_name,
            fit_parameter_source="parameter_value",
            family=prior_model.family,
            parameter_value=float(format(fit_report.parameter_value, ".15g")),
            hyperparameter_values=prior_model.parameter_values(),
            log_prior_contribution=float(format(log_prior, ".15g")),
        )
    ]


def _resolve_sigma_rows(
    *,
    fit_report: ContinuousEvolutionaryModeFitReport,
    prior_bundle: ContinuousTraitModelPriorBundle,
) -> list[ContinuousTraitModelPriorRow]:
    if prior_bundle.sigma_prior is None:
        return []
    log_prior = evaluate_continuous_trait_scalar_log_prior(
        parameter_value=fit_report.rate,
        prior_model=prior_bundle.sigma_prior,
        parameter_name="sigma",
        allow_zero=False,
    )
    return [
        ContinuousTraitModelPriorRow(
            target_name="sigma",
            fit_parameter_source="rate",
            family=prior_bundle.sigma_prior.family,
            parameter_value=float(format(fit_report.rate, ".15g")),
            hyperparameter_values=prior_bundle.sigma_prior.parameter_values(),
            log_prior_contribution=float(format(log_prior, ".15g")),
        )
    ]


def _gamma_log_density(
    parameter_value: float,
    *,
    shape: float,
    scale: float,
    parameter_name: str,
    allow_zero: bool,
) -> float:
    if parameter_value == 0.0:
        if not allow_zero:
            raise PhylogeneticsError(
                "gamma continuous-trait scalar prior cannot score a zero value on positive support",
                code="continuous_trait_scalar_prior_zero_gamma_positive_support",
                details={"parameter_name": parameter_name},
            )
        if shape < 1.0:
            raise PhylogeneticsError(
                "gamma continuous-trait scalar prior with shape below one is singular at zero and is unsupported for boundary-valued fits",
                code="continuous_trait_scalar_prior_zero_gamma_singular",
                details={"parameter_name": parameter_name, "shape": shape},
            )
        if math.isclose(shape, 1.0, rel_tol=0.0, abs_tol=1e-12):
            return -math.log(scale)
        return -math.inf
    return (
        ((shape - 1.0) * math.log(parameter_value))
        - (parameter_value / scale)
        - math.lgamma(shape)
        - (shape * math.log(scale))
    )


def _lognormal_log_density(
    parameter_value: float,
    *,
    log_mean: float,
    log_standard_deviation: float,
) -> float:
    return (
        -math.log(parameter_value)
        - math.log(log_standard_deviation)
        - (0.5 * math.log(2.0 * math.pi))
        - (
            ((math.log(parameter_value) - log_mean) ** 2)
            / (2.0 * (log_standard_deviation**2))
        )
    )


def _beta_log_density(
    parameter_value: float,
    *,
    alpha: float,
    beta: float,
    parameter_name: str,
) -> float:
    normalization = math.lgamma(alpha + beta) - math.lgamma(alpha) - math.lgamma(beta)
    if parameter_value == 0.0:
        if alpha < 1.0:
            raise PhylogeneticsError(
                "beta continuous-trait probability prior with alpha below one is singular at zero and is unsupported for boundary-valued fits",
                code="continuous_trait_probability_prior_zero_beta_singular",
                details={"parameter_name": parameter_name, "alpha": alpha},
            )
        return (
            normalization
            if math.isclose(alpha, 1.0, rel_tol=0.0, abs_tol=1e-12)
            else -math.inf
        )
    if parameter_value == 1.0:
        if beta < 1.0:
            raise PhylogeneticsError(
                "beta continuous-trait probability prior with beta below one is singular at one and is unsupported for boundary-valued fits",
                code="continuous_trait_probability_prior_one_beta_singular",
                details={"parameter_name": parameter_name, "beta": beta},
            )
        return (
            normalization
            if math.isclose(beta, 1.0, rel_tol=0.0, abs_tol=1e-12)
            else -math.inf
        )
    return (
        normalization
        + ((alpha - 1.0) * math.log(parameter_value))
        + ((beta - 1.0) * math.log1p(-parameter_value))
    )


def _resolve_supported_continuous_trait_prior_mode(mode: str) -> str:
    if mode not in CONTINUOUS_TRAIT_PRIOR_MODES:
        raise PhylogeneticsError(
            "continuous-trait model priors support only the owned evolutionary mode fit surface",
            code="continuous_trait_prior_mode_unsupported",
            details={
                "mode": mode,
                "allowed_modes": list(CONTINUOUS_TRAIT_PRIOR_MODES),
            },
        )
    return mode


def _validate_scalar_parameter_value(
    *,
    parameter_name: str,
    value: float,
    allow_zero: bool,
) -> float:
    if allow_zero:
        return _validate_nonnegative_finite_value(
            parameter_name=parameter_name,
            value=value,
            owner_name="continuous-trait scalar prior evaluation",
        )
    return _validate_positive_finite_value(
        parameter_name=parameter_name,
        value=value,
        owner_name="continuous-trait scalar prior evaluation",
    )


def _validate_finite_value(
    *,
    parameter_name: str,
    value: float,
    owner_name: str,
) -> float:
    if not math.isfinite(value):
        raise PhylogeneticsError(
            f"{owner_name} requires {parameter_name} to be finite",
            code="continuous_trait_prior_value_not_finite",
            details={"parameter_name": parameter_name, "value": value},
        )
    return value


def _validate_positive_finite_value(
    *,
    parameter_name: str,
    value: float,
    owner_name: str,
) -> float:
    validated = _validate_finite_value(
        parameter_name=parameter_name,
        value=value,
        owner_name=owner_name,
    )
    if validated <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires {parameter_name} to be strictly positive",
            code="continuous_trait_prior_value_not_positive",
            details={"parameter_name": parameter_name, "value": value},
        )
    return validated


def _validate_nonnegative_finite_value(
    *,
    parameter_name: str,
    value: float,
    owner_name: str,
) -> float:
    validated = _validate_finite_value(
        parameter_name=parameter_name,
        value=value,
        owner_name=owner_name,
    )
    if validated < 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires {parameter_name} to be nonnegative",
            code="continuous_trait_prior_value_negative",
            details={"parameter_name": parameter_name, "value": value},
        )
    return validated


def _validate_probability_value(
    *,
    parameter_name: str,
    value: float,
    owner_name: str,
) -> float:
    validated = _validate_finite_value(
        parameter_name=parameter_name,
        value=value,
        owner_name=owner_name,
    )
    if not 0.0 <= validated <= 1.0:
        raise PhylogeneticsError(
            f"{owner_name} requires {parameter_name} to lie within [0, 1]",
            code="continuous_trait_prior_value_not_probability",
            details={"parameter_name": parameter_name, "value": value},
        )
    return validated
