from __future__ import annotations

from dataclasses import dataclass
from statistics import mean

from bijux_phylogenetics.phylo.likelihood.dna import validate_positive_kappa
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .metropolis_hastings import (
    MetropolisHastingsRunReport,
    list_reversible_jump_model_switch_families,
    validate_reversible_jump_model_switch_family,
)
from .posterior_sets.diagnostics import highest_posterior_density_interval
from .state import BayesianPhylogeneticState

_DEFAULT_MODEL_FAMILY = "nucleotide-substitution-model"
_NUCLEOTIDE_MODELS = ("JC69", "K80")
_TRANSITION_TRANSVERSION_RATIO = "transition-transversion-ratio"


@dataclass(frozen=True, slots=True)
class PosteriorModelSupportRow:
    """One posterior model-support summary row."""

    model_name: str
    posterior_probability: float
    supporting_sample_count: int
    total_sample_count: int


@dataclass(frozen=True, slots=True)
class PosteriorModelEstimateRow:
    """One posterior estimate summary row for one supported model."""

    estimate_name: str
    model_name: str
    posterior_probability: float
    posterior_mean: float
    hpd_95_lower: float
    hpd_95_upper: float
    supporting_sample_count: int
    total_sample_count: int


@dataclass(frozen=True, slots=True)
class PosteriorModelAveragedEstimateRow:
    """One model-averaged posterior estimate summary row."""

    estimate_name: str
    posterior_mean: float
    hpd_95_lower: float
    hpd_95_upper: float
    sample_count: int
    contributing_models: list[str]


@dataclass(frozen=True, slots=True)
class PosteriorModelAveragingReport:
    """One posterior model-averaging summary across sampled model states."""

    sample_count: int
    model_family: str
    sampled_models: list[str]
    support_rows: list[PosteriorModelSupportRow]
    per_model_estimate_rows: list[PosteriorModelEstimateRow]
    model_averaged_estimate_rows: list[PosteriorModelAveragedEstimateRow]
    warnings: list[str]


def summarize_metropolis_hastings_model_averaged_estimates(
    *,
    run_report: MetropolisHastingsRunReport,
    model_family: str = _DEFAULT_MODEL_FAMILY,
) -> PosteriorModelAveragingReport:
    """Summarize posterior model support and model-averaged estimates from one chain."""
    if not isinstance(run_report, MetropolisHastingsRunReport):
        raise PhylogeneticsError(
            "posterior model-averaging summary requires one MetropolisHastingsRunReport",
            code="posterior_model_averaging_run_report_type_invalid",
        )
    return summarize_posterior_model_averaged_estimates(
        sampled_states=run_report.sampled_states,
        model_family=model_family,
    )


def summarize_posterior_model_averaged_estimates(
    *,
    sampled_states: list[BayesianPhylogeneticState]
    | tuple[BayesianPhylogeneticState, ...],
    model_family: str = _DEFAULT_MODEL_FAMILY,
) -> PosteriorModelAveragingReport:
    """Summarize posterior model support and model-averaged estimates."""
    validated_model_family = validate_reversible_jump_model_switch_family(model_family)
    validated_sampled_states = _validate_sampled_states(sampled_states)
    if validated_model_family != _DEFAULT_MODEL_FAMILY:
        raise AssertionError(
            "posterior model-averaging summary reached one unsupported family handler"
        )
    estimate_values_by_model: dict[str, list[float]] = {
        model_name: [] for model_name in _NUCLEOTIDE_MODELS
    }
    for sampled_state in validated_sampled_states:
        model_name, estimate_value = _resolve_nucleotide_model_averaging_sample(
            sampled_state
        )
        estimate_values_by_model[model_name].append(estimate_value)
    sample_count = len(validated_sampled_states)
    support_rows = _build_support_rows(
        estimate_values_by_model=estimate_values_by_model,
        sample_count=sample_count,
    )
    per_model_estimate_rows = _build_per_model_estimate_rows(
        estimate_values_by_model=estimate_values_by_model,
        support_rows=support_rows,
        sample_count=sample_count,
    )
    model_averaged_estimate_rows = _build_model_averaged_estimate_rows(
        estimate_values_by_model=estimate_values_by_model,
        sample_count=sample_count,
    )
    return PosteriorModelAveragingReport(
        sample_count=sample_count,
        model_family=validated_model_family,
        sampled_models=[support_row.model_name for support_row in support_rows],
        support_rows=support_rows,
        per_model_estimate_rows=per_model_estimate_rows,
        model_averaged_estimate_rows=model_averaged_estimate_rows,
        warnings=_build_posterior_model_averaging_warnings(
            support_rows=support_rows,
            model_family=validated_model_family,
        ),
    )


def _resolve_nucleotide_model_averaging_sample(
    sampled_state: BayesianPhylogeneticState,
) -> tuple[str, float]:
    model_name = sampled_state.model_parameters.categorical_parameters.get(
        "substitution-model"
    )
    if model_name is None:
        raise PhylogeneticsError(
            "posterior model-averaging summary requires one 'substitution-model' categorical parameter in every sampled state",
            code="posterior_model_averaging_model_name_missing",
        )
    normalized_model_name = model_name.strip().upper()
    if normalized_model_name not in _NUCLEOTIDE_MODELS:
        raise PhylogeneticsError(
            "posterior model-averaging summary currently supports only JC69 and K80 within the nucleotide-substitution-model family",
            code="posterior_model_averaging_model_name_unsupported",
            details={
                "model_name": model_name,
                "allowed_model_names": list(_NUCLEOTIDE_MODELS),
                "declared_model_families": list(
                    list_reversible_jump_model_switch_families()
                ),
            },
        )
    if normalized_model_name == "JC69":
        if "kappa" in sampled_state.model_parameters.scalar_parameters:
            raise PhylogeneticsError(
                "posterior model-averaging summary requires JC69 sampled states to omit the standalone 'kappa' scalar parameter",
                code="posterior_model_averaging_jc69_kappa_unexpected",
            )
        return normalized_model_name, 1.0
    raw_kappa = sampled_state.model_parameters.scalar_parameters.get("kappa")
    if raw_kappa is None:
        raise PhylogeneticsError(
            "posterior model-averaging summary requires K80 sampled states to include one positive 'kappa' scalar parameter",
            code="posterior_model_averaging_k80_kappa_missing",
        )
    try:
        validated_kappa = float(
            format(
                validate_positive_kappa(
                    raw_kappa,
                    model_name="posterior model-averaging summary",
                ),
                ".15g",
            )
        )
    except ValueError as error:
        raise PhylogeneticsError(
            str(error),
            code="posterior_model_averaging_k80_kappa_invalid",
        ) from error
    return normalized_model_name, validated_kappa


def _validate_sampled_states(
    sampled_states: list[BayesianPhylogeneticState]
    | tuple[BayesianPhylogeneticState, ...],
) -> tuple[BayesianPhylogeneticState, ...]:
    validated_sampled_states = tuple(sampled_states)
    if not validated_sampled_states:
        raise PhylogeneticsError(
            "posterior model-averaging summary requires at least one sampled Bayesian phylogenetic state",
            code="posterior_model_averaging_sampled_states_empty",
        )
    if any(
        not isinstance(sampled_state, BayesianPhylogeneticState)
        for sampled_state in validated_sampled_states
    ):
        raise PhylogeneticsError(
            "posterior model-averaging summary requires every sampled state to be one BayesianPhylogeneticState",
            code="posterior_model_averaging_sampled_state_type_invalid",
        )
    return validated_sampled_states


def _build_support_rows(
    *,
    estimate_values_by_model: dict[str, list[float]],
    sample_count: int,
) -> list[PosteriorModelSupportRow]:
    rows: list[PosteriorModelSupportRow] = []
    for model_name, values in sorted(
        estimate_values_by_model.items(),
        key=lambda item: (-len(item[1]), item[0]),
    ):
        if not values:
            continue
        rows.append(
            PosteriorModelSupportRow(
                model_name=model_name,
                posterior_probability=float(format(len(values) / sample_count, ".15g")),
                supporting_sample_count=len(values),
                total_sample_count=sample_count,
            )
        )
    return rows


def _build_per_model_estimate_rows(
    *,
    estimate_values_by_model: dict[str, list[float]],
    support_rows: list[PosteriorModelSupportRow],
    sample_count: int,
) -> list[PosteriorModelEstimateRow]:
    support_by_model = {row.model_name: row for row in support_rows}
    rows: list[PosteriorModelEstimateRow] = []
    for model_name in [support_row.model_name for support_row in support_rows]:
        values = estimate_values_by_model[model_name]
        hpd_95_lower, hpd_95_upper = highest_posterior_density_interval(values)
        support_row = support_by_model[model_name]
        rows.append(
            PosteriorModelEstimateRow(
                estimate_name=_TRANSITION_TRANSVERSION_RATIO,
                model_name=model_name,
                posterior_probability=support_row.posterior_probability,
                posterior_mean=float(format(mean(values), ".15g")),
                hpd_95_lower=float(format(hpd_95_lower, ".15g")),
                hpd_95_upper=float(format(hpd_95_upper, ".15g")),
                supporting_sample_count=len(values),
                total_sample_count=sample_count,
            )
        )
    return rows


def _build_model_averaged_estimate_rows(
    *,
    estimate_values_by_model: dict[str, list[float]],
    sample_count: int,
) -> list[PosteriorModelAveragedEstimateRow]:
    pooled_values = [
        estimate_value
        for model_name in _NUCLEOTIDE_MODELS
        for estimate_value in estimate_values_by_model[model_name]
    ]
    if not pooled_values:
        raise PhylogeneticsError(
            "posterior model-averaging summary requires at least one supported sampled value",
            code="posterior_model_averaging_estimate_values_empty",
        )
    hpd_95_lower, hpd_95_upper = highest_posterior_density_interval(pooled_values)
    return [
        PosteriorModelAveragedEstimateRow(
            estimate_name=_TRANSITION_TRANSVERSION_RATIO,
            posterior_mean=float(format(mean(pooled_values), ".15g")),
            hpd_95_lower=float(format(hpd_95_lower, ".15g")),
            hpd_95_upper=float(format(hpd_95_upper, ".15g")),
            sample_count=sample_count,
            contributing_models=[
                model_name
                for model_name in _NUCLEOTIDE_MODELS
                if estimate_values_by_model[model_name]
            ],
        )
    ]


def _build_posterior_model_averaging_warnings(
    *,
    support_rows: list[PosteriorModelSupportRow],
    model_family: str,
) -> list[str]:
    if len(support_rows) > 1:
        return []
    if not support_rows:
        return []
    return [
        f"posterior model-averaging summary for {model_family} observed only one sampled model state, so model-averaged estimates equal the retained model-specific summary"
    ]


__all__ = [
    "PosteriorModelAveragedEstimateRow",
    "PosteriorModelAveragingReport",
    "PosteriorModelEstimateRow",
    "PosteriorModelSupportRow",
    "summarize_metropolis_hastings_model_averaged_estimates",
    "summarize_posterior_model_averaged_estimates",
]
