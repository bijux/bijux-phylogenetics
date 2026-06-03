from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math
from statistics import mean

from bijux_phylogenetics.bayesian.continuous_trait_location_priors import (
    ContinuousTraitLocationPriorModel,
    evaluate_continuous_trait_location_log_prior,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    ContinuousTraitScalarPriorModel,
    evaluate_continuous_trait_scalar_log_prior,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsRunReport,
    propose_clock_rate_move,
    propose_continuous_trait_location_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.posterior_sets.diagnostics import (
    highest_posterior_density_interval,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
    build_bayesian_prior_component_state,
)
from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    log_determinant,
    quadratic_form,
    stable_covariance,
)
from bijux_phylogenetics.comparative.common import build_brownian_covariance_matrix
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

BROWNIAN_CONTINUOUS_TRAIT_MODELS = ("brownian",)

_BROWNIAN_MODEL_NAME = "brownian"
_MODEL_PARAMETER_NAME = "continuous-trait-model"
_ROOT_STATE_PARAMETER_NAME = "root-state"
_SIGMA_SQUARED_PARAMETER_NAME = "sigma-squared"


@dataclass(frozen=True, slots=True)
class BrownianContinuousTraitModelDefinition:
    """One validated fixed-topology Bayesian Brownian continuous-trait model."""

    root_state_prior: ContinuousTraitLocationPriorModel
    sigma_squared_prior: ContinuousTraitScalarPriorModel
    initial_root_state: float | None = None
    initial_sigma_squared: float | None = None


@dataclass(frozen=True, slots=True)
class BrownianContinuousTraitProposalSchedule:
    """One validated proposal schedule for Bayesian Brownian trait sampling."""

    root_state_move_weight: float
    root_state_standard_deviation: float
    sigma_squared_move_weight: float
    sigma_squared_log_scale_standard_deviation: float


@dataclass(frozen=True, slots=True)
class BrownianContinuousTraitPosteriorRow:
    """One sampled posterior row from a Bayesian Brownian trait chain."""

    sample_index: int
    iteration_index: int
    topology_id: str
    model_name: str
    total_log_prior: float
    log_likelihood: float
    posterior_log_score: float
    prior_component_log_priors: dict[str, float]
    root_state: float
    sigma_squared: float


@dataclass(frozen=True, slots=True)
class BrownianContinuousTraitParameterSummary:
    """One posterior summary for a Brownian continuous-trait parameter."""

    parameter_name: str
    sample_count: int
    posterior_mean: float
    hpd_95_lower: float
    hpd_95_upper: float


@dataclass(frozen=True, slots=True)
class BrownianContinuousTraitRunReport:
    """One completed fixed-topology Bayesian Brownian trait posterior run."""

    model_definition: BrownianContinuousTraitModelDefinition
    proposal_schedule: BrownianContinuousTraitProposalSchedule
    taxa: list[str]
    tip_values: dict[str, float]
    chain_report: MetropolisHastingsRunReport
    posterior_rows: list[BrownianContinuousTraitPosteriorRow]
    parameter_summaries: list[BrownianContinuousTraitParameterSummary]


def build_brownian_continuous_trait_model_definition(
    *,
    root_state_prior: ContinuousTraitLocationPriorModel,
    sigma_squared_prior: ContinuousTraitScalarPriorModel,
    initial_root_state: float | None = None,
    initial_sigma_squared: float | None = None,
) -> BrownianContinuousTraitModelDefinition:
    """Build one validated Bayesian Brownian continuous-trait model definition."""
    if not isinstance(root_state_prior, ContinuousTraitLocationPriorModel):
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait model requires one ContinuousTraitLocationPriorModel for root_state_prior",
            code="brownian_continuous_trait_root_state_prior_type_invalid",
        )
    if root_state_prior.family == "fixed":
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait model requires one non-fixed root-state prior so root state remains sampled",
            code="brownian_continuous_trait_root_state_prior_fixed",
        )
    if not isinstance(sigma_squared_prior, ContinuousTraitScalarPriorModel):
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait model requires one ContinuousTraitScalarPriorModel for sigma_squared_prior",
            code="brownian_continuous_trait_sigma_squared_prior_type_invalid",
        )
    if sigma_squared_prior.family == "fixed":
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait model requires one non-fixed sigma-squared prior so diffusion variance remains sampled",
            code="brownian_continuous_trait_sigma_squared_prior_fixed",
        )
    return BrownianContinuousTraitModelDefinition(
        root_state_prior=root_state_prior,
        sigma_squared_prior=sigma_squared_prior,
        initial_root_state=_validate_optional_finite_float(
            value=initial_root_state,
            field_name="initial_root_state",
            owner_name="bayesian Brownian continuous-trait model",
        ),
        initial_sigma_squared=_validate_optional_positive_finite_float(
            value=initial_sigma_squared,
            field_name="initial_sigma_squared",
            owner_name="bayesian Brownian continuous-trait model",
        ),
    )


def build_brownian_continuous_trait_proposal_schedule(
    *,
    model_definition: BrownianContinuousTraitModelDefinition,
    root_state_move_weight: float,
    root_state_standard_deviation: float,
    sigma_squared_move_weight: float,
    sigma_squared_log_scale_standard_deviation: float,
) -> BrownianContinuousTraitProposalSchedule:
    """Build one validated proposal schedule for Bayesian Brownian trait sampling."""
    if not isinstance(model_definition, BrownianContinuousTraitModelDefinition):
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait proposal schedule requires one BrownianContinuousTraitModelDefinition",
            code="brownian_continuous_trait_proposal_schedule_model_definition_type_invalid",
        )
    validated_root_state_move_weight = _validate_positive_finite_float(
        value=root_state_move_weight,
        field_name="root_state_move_weight",
        owner_name="bayesian Brownian continuous-trait proposal schedule",
    )
    validated_sigma_squared_move_weight = _validate_positive_finite_float(
        value=sigma_squared_move_weight,
        field_name="sigma_squared_move_weight",
        owner_name="bayesian Brownian continuous-trait proposal schedule",
    )
    return BrownianContinuousTraitProposalSchedule(
        root_state_move_weight=validated_root_state_move_weight,
        root_state_standard_deviation=_validate_positive_finite_float(
            value=root_state_standard_deviation,
            field_name="root_state_standard_deviation",
            owner_name="bayesian Brownian continuous-trait proposal schedule",
        ),
        sigma_squared_move_weight=validated_sigma_squared_move_weight,
        sigma_squared_log_scale_standard_deviation=_validate_positive_finite_float(
            value=sigma_squared_log_scale_standard_deviation,
            field_name="sigma_squared_log_scale_standard_deviation",
            owner_name="bayesian Brownian continuous-trait proposal schedule",
        ),
    )


def run_brownian_continuous_trait_metropolis_hastings(
    *,
    tree: PhyloTree,
    tip_values: Mapping[str, float],
    model_definition: BrownianContinuousTraitModelDefinition,
    proposal_schedule: BrownianContinuousTraitProposalSchedule,
    iteration_count: int,
    sample_every: int = 1,
    seed: int = 0,
) -> BrownianContinuousTraitRunReport:
    """Run one fixed-topology Bayesian Brownian trait posterior sampler."""
    if not isinstance(tree, PhyloTree):
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait runner requires one PhyloTree",
            code="brownian_continuous_trait_tree_type_invalid",
        )
    if not isinstance(model_definition, BrownianContinuousTraitModelDefinition):
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait runner requires one BrownianContinuousTraitModelDefinition",
            code="brownian_continuous_trait_model_definition_type_invalid",
        )
    if not isinstance(proposal_schedule, BrownianContinuousTraitProposalSchedule):
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait runner requires one BrownianContinuousTraitProposalSchedule",
            code="brownian_continuous_trait_proposal_schedule_type_invalid",
        )
    normalized_tree = tree.copy()
    normalized_tree.rooted = tree.rooted
    normalized_tree.refresh()
    if normalized_tree.rooted is not True:
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait runner requires one rooted tree",
            code="brownian_continuous_trait_tree_rooting_invalid",
        )
    taxa, normalized_tip_values = _normalize_tip_values(
        tree=normalized_tree,
        tip_values=tip_values,
    )
    covariance = stable_covariance(
        build_brownian_covariance_matrix(normalized_tree, taxa)
    )
    inverse_covariance = invert_matrix(covariance)
    covariance_log_determinant = log_determinant(covariance)
    ordered_values = [normalized_tip_values[taxon] for taxon in taxa]
    initial_root_state = (
        model_definition.initial_root_state
        if model_definition.initial_root_state is not None
        else _estimate_initial_root_state(
            ordered_values=ordered_values,
            inverse_covariance=inverse_covariance,
        )
    )
    initial_sigma_squared = (
        model_definition.initial_sigma_squared
        if model_definition.initial_sigma_squared is not None
        else _estimate_initial_sigma_squared(
            ordered_values=ordered_values,
            inverse_covariance=inverse_covariance,
            root_state=initial_root_state,
        )
    )
    initial_model_parameters = build_bayesian_model_parameter_state(
        categorical_parameters={_MODEL_PARAMETER_NAME: _BROWNIAN_MODEL_NAME},
        scalar_parameters={
            _ROOT_STATE_PARAMETER_NAME: initial_root_state,
            _SIGMA_SQUARED_PARAMETER_NAME: initial_sigma_squared,
        },
    )
    initial_state = score_bayesian_phylogenetic_state(
        tree=normalized_tree,
        model_parameters=initial_model_parameters,
        update_prior_components=lambda state: (
            _build_brownian_continuous_trait_prior_components(
                state=state,
                model_definition=model_definition,
                fixed_topology_id=None,
            )
        ),
        update_log_likelihood=lambda state: (
            _evaluate_brownian_continuous_trait_log_likelihood(
                state=state,
                fixed_topology_id=None,
                ordered_values=ordered_values,
                inverse_covariance=inverse_covariance,
                covariance_log_determinant=covariance_log_determinant,
            )
        ),
    )
    fixed_topology_id = initial_state.tree.topology_id
    chain_report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: (
            _propose_brownian_continuous_trait_state(
                current_state=current_state,
                rng=rng,
                proposal_schedule=proposal_schedule,
            )
        ),
        update_prior_components=lambda state: (
            _build_brownian_continuous_trait_prior_components(
                state=state,
                model_definition=model_definition,
                fixed_topology_id=fixed_topology_id,
            )
        ),
        update_log_likelihood=lambda state: (
            _evaluate_brownian_continuous_trait_log_likelihood(
                state=state,
                fixed_topology_id=fixed_topology_id,
                ordered_values=ordered_values,
                inverse_covariance=inverse_covariance,
                covariance_log_determinant=covariance_log_determinant,
            )
        ),
        iteration_count=iteration_count,
        sample_every=sample_every,
        seed=seed,
    )
    posterior_rows = _build_brownian_continuous_trait_posterior_rows(
        chain_report=chain_report,
        fixed_topology_id=fixed_topology_id,
    )
    return BrownianContinuousTraitRunReport(
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        taxa=taxa,
        tip_values=normalized_tip_values,
        chain_report=chain_report,
        posterior_rows=posterior_rows,
        parameter_summaries=_build_parameter_summaries(posterior_rows),
    )


def _normalize_tip_values(
    *,
    tree: PhyloTree,
    tip_values: Mapping[str, float],
) -> tuple[list[str], dict[str, float]]:
    if not isinstance(tip_values, Mapping):
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait runner requires one mapping of tip values",
            code="brownian_continuous_trait_tip_values_type_invalid",
        )
    taxa = list(tree.tip_names)
    if len(taxa) < 3:
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait runner requires at least three taxa",
            code="brownian_continuous_trait_taxon_count_invalid",
            details={"taxon_count": len(taxa)},
        )
    tree_taxa = set(taxa)
    tip_value_taxa = set(tip_values)
    missing_taxa = sorted(tree_taxa - tip_value_taxa)
    extra_taxa = sorted(tip_value_taxa - tree_taxa)
    if missing_taxa or extra_taxa:
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait runner requires tip values to match the tree tip set exactly",
            code="brownian_continuous_trait_tip_value_taxa_mismatch",
            details={
                "missing_taxa": missing_taxa,
                "extra_taxa": extra_taxa,
            },
        )
    normalized_tip_values = {
        taxon: _validate_finite_float(
            value=tip_values[taxon],
            field_name=taxon,
            owner_name="bayesian Brownian continuous-trait tip values",
        )
        for taxon in taxa
    }
    if len({format(value, ".12g") for value in normalized_tip_values.values()}) == 1:
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait runner requires non-constant tip values so sigma-squared remains identifiable under positive-support sampling",
            code="brownian_continuous_trait_tip_values_constant",
        )
    return taxa, normalized_tip_values


def _estimate_initial_root_state(
    *,
    ordered_values: list[float],
    inverse_covariance: list[list[float]],
) -> float:
    ones = [1.0] * len(ordered_values)
    denominator = quadratic_form(ones, inverse_covariance)
    if not math.isfinite(denominator) or denominator <= 0.0:
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait runner could not estimate one finite positive generalized-least-squares denominator for the initial root state",
            code="brownian_continuous_trait_initial_root_state_denominator_invalid",
            details={"denominator": denominator},
        )
    numerator = sum(
        inverse_covariance[row_index][column_index] * ordered_values[column_index]
        for row_index in range(len(ordered_values))
        for column_index in range(len(ordered_values))
    )
    return _validate_finite_float(
        value=numerator / denominator,
        field_name="initial_root_state",
        owner_name="bayesian Brownian continuous-trait runner",
    )


def _estimate_initial_sigma_squared(
    *,
    ordered_values: list[float],
    inverse_covariance: list[list[float]],
    root_state: float,
) -> float:
    residuals = [value - root_state for value in ordered_values]
    sigma_squared = quadratic_form(residuals, inverse_covariance) / len(ordered_values)
    if not math.isfinite(sigma_squared) or sigma_squared <= 0.0:
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait runner requires one strictly positive initial sigma-squared value",
            code="brownian_continuous_trait_initial_sigma_squared_invalid",
            details={"sigma_squared": sigma_squared},
        )
    return sigma_squared


def _build_brownian_continuous_trait_prior_components(
    *,
    state: BayesianPhylogeneticState,
    model_definition: BrownianContinuousTraitModelDefinition,
    fixed_topology_id: str | None,
) -> list[BayesianPriorComponentState]:
    _require_brownian_continuous_trait_state_consistency(
        state=state,
        fixed_topology_id=fixed_topology_id,
    )
    root_state = state.model_parameters.scalar_parameters[_ROOT_STATE_PARAMETER_NAME]
    sigma_squared = state.model_parameters.scalar_parameters[
        _SIGMA_SQUARED_PARAMETER_NAME
    ]
    root_state_log_prior = evaluate_continuous_trait_location_log_prior(
        parameter_value=root_state,
        prior_model=model_definition.root_state_prior,
        parameter_name=_ROOT_STATE_PARAMETER_NAME,
    )
    sigma_squared_log_prior = evaluate_continuous_trait_scalar_log_prior(
        parameter_value=sigma_squared,
        prior_model=model_definition.sigma_squared_prior,
        parameter_name=_SIGMA_SQUARED_PARAMETER_NAME,
        allow_zero=False,
    )
    return [
        build_bayesian_prior_component_state(
            component_name="continuous-trait:root-state",
            family=model_definition.root_state_prior.family,
            log_prior=root_state_log_prior,
            parameter_values=model_definition.root_state_prior.parameter_values(),
        ),
        build_bayesian_prior_component_state(
            component_name="continuous-trait:sigma-squared",
            family=model_definition.sigma_squared_prior.family,
            log_prior=sigma_squared_log_prior,
            parameter_values=model_definition.sigma_squared_prior.parameter_values(),
        ),
    ]


def _evaluate_brownian_continuous_trait_log_likelihood(
    *,
    state: BayesianPhylogeneticState,
    fixed_topology_id: str | None,
    ordered_values: list[float],
    inverse_covariance: list[list[float]],
    covariance_log_determinant: float,
) -> float:
    _require_brownian_continuous_trait_state_consistency(
        state=state,
        fixed_topology_id=fixed_topology_id,
    )
    root_state = _validate_finite_float(
        value=state.model_parameters.scalar_parameters[_ROOT_STATE_PARAMETER_NAME],
        field_name=_ROOT_STATE_PARAMETER_NAME,
        owner_name="bayesian Brownian continuous-trait log likelihood",
    )
    sigma_squared = _validate_positive_finite_float(
        value=state.model_parameters.scalar_parameters[_SIGMA_SQUARED_PARAMETER_NAME],
        field_name=_SIGMA_SQUARED_PARAMETER_NAME,
        owner_name="bayesian Brownian continuous-trait log likelihood",
    )
    residuals = [value - root_state for value in ordered_values]
    quadratic_term = quadratic_form(residuals, inverse_covariance)
    return -0.5 * (
        len(ordered_values) * math.log(2.0 * math.pi)
        + len(ordered_values) * math.log(sigma_squared)
        + covariance_log_determinant
        + (quadratic_term / sigma_squared)
    )


def _propose_brownian_continuous_trait_state(
    *,
    current_state: BayesianPhylogeneticState,
    rng,
    proposal_schedule: BrownianContinuousTraitProposalSchedule,
):
    weighted_moves = [
        (
            proposal_schedule.root_state_move_weight,
            lambda: propose_continuous_trait_location_move(
                current_state,
                rng,
                standard_deviation=proposal_schedule.root_state_standard_deviation,
                parameter_name=_ROOT_STATE_PARAMETER_NAME,
            ),
        ),
        (
            proposal_schedule.sigma_squared_move_weight,
            lambda: propose_clock_rate_move(
                current_state,
                rng,
                log_scale_standard_deviation=(
                    proposal_schedule.sigma_squared_log_scale_standard_deviation
                ),
                parameter_name=_SIGMA_SQUARED_PARAMETER_NAME,
            ),
        ),
    ]
    total_weight = math.fsum(weight for weight, _move in weighted_moves)
    move_threshold = rng.random() * total_weight
    cumulative_weight = 0.0
    for weight, move in weighted_moves:
        cumulative_weight += weight
        if move_threshold <= cumulative_weight:
            return move()
    return weighted_moves[-1][1]()


def _build_brownian_continuous_trait_posterior_rows(
    *,
    chain_report: MetropolisHastingsRunReport,
    fixed_topology_id: str,
) -> list[BrownianContinuousTraitPosteriorRow]:
    posterior_rows: list[BrownianContinuousTraitPosteriorRow] = []
    for sample_index, state in enumerate(chain_report.sampled_states):
        if state.tree.topology_id != fixed_topology_id:
            raise PhylogeneticsError(
                "bayesian Brownian continuous-trait posterior trace detected one topology change in sampled states",
                code="brownian_continuous_trait_trace_topology_changed",
                details={
                    "expected_topology_id": fixed_topology_id,
                    "observed_topology_id": state.tree.topology_id,
                    "sample_index": sample_index,
                },
            )
        posterior_rows.append(
            BrownianContinuousTraitPosteriorRow(
                sample_index=sample_index,
                iteration_index=sample_index * chain_report.sample_every,
                topology_id=state.tree.topology_id,
                model_name=_BROWNIAN_MODEL_NAME,
                total_log_prior=state.total_log_prior,
                log_likelihood=state.log_likelihood,
                posterior_log_score=state.posterior_log_score,
                prior_component_log_priors={
                    component.component_name: component.log_prior
                    for component in state.prior_components
                },
                root_state=state.model_parameters.scalar_parameters[
                    _ROOT_STATE_PARAMETER_NAME
                ],
                sigma_squared=state.model_parameters.scalar_parameters[
                    _SIGMA_SQUARED_PARAMETER_NAME
                ],
            )
        )
    return posterior_rows


def _build_parameter_summaries(
    posterior_rows: list[BrownianContinuousTraitPosteriorRow],
) -> list[BrownianContinuousTraitParameterSummary]:
    root_state_values = [row.root_state for row in posterior_rows]
    sigma_squared_values = [row.sigma_squared for row in posterior_rows]
    return [
        _build_parameter_summary(
            parameter_name=_ROOT_STATE_PARAMETER_NAME,
            values=root_state_values,
        ),
        _build_parameter_summary(
            parameter_name=_SIGMA_SQUARED_PARAMETER_NAME,
            values=sigma_squared_values,
        ),
    ]


def _build_parameter_summary(
    *,
    parameter_name: str,
    values: list[float],
) -> BrownianContinuousTraitParameterSummary:
    hpd_95_lower, hpd_95_upper = highest_posterior_density_interval(values)
    return BrownianContinuousTraitParameterSummary(
        parameter_name=parameter_name,
        sample_count=len(values),
        posterior_mean=round(mean(values), 6),
        hpd_95_lower=round(hpd_95_lower, 6),
        hpd_95_upper=round(hpd_95_upper, 6),
    )


def _require_brownian_continuous_trait_state_consistency(
    *,
    state: BayesianPhylogeneticState,
    fixed_topology_id: str | None,
) -> None:
    model_name = state.model_parameters.categorical_parameters.get(
        _MODEL_PARAMETER_NAME
    )
    if model_name != _BROWNIAN_MODEL_NAME:
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait model requires every sampled state to preserve the configured model label",
            code="brownian_continuous_trait_state_model_label_invalid",
            details={
                "expected_model_name": _BROWNIAN_MODEL_NAME,
                "observed_model_name": model_name,
            },
        )
    if fixed_topology_id is not None and state.tree.topology_id != fixed_topology_id:
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait model requires topology to remain unchanged across sampled states",
            code="brownian_continuous_trait_state_topology_changed",
            details={
                "expected_topology_id": fixed_topology_id,
                "observed_topology_id": state.tree.topology_id,
            },
        )
    if _ROOT_STATE_PARAMETER_NAME not in state.model_parameters.scalar_parameters:
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait model requires one root-state scalar parameter in every sampled state",
            code="brownian_continuous_trait_state_root_state_missing",
        )
    if _SIGMA_SQUARED_PARAMETER_NAME not in state.model_parameters.scalar_parameters:
        raise PhylogeneticsError(
            "bayesian Brownian continuous-trait model requires one sigma-squared scalar parameter in every sampled state",
            code="brownian_continuous_trait_state_sigma_squared_missing",
        )


def _validate_positive_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    validated_value = _validate_finite_float(
        value=value,
        field_name=field_name,
        owner_name=owner_name,
    )
    if validated_value <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires {field_name} > 0",
            code="brownian_continuous_trait_parameter_nonpositive",
            details={
                "field_name": field_name,
                "field_value": value,
                "owner_name": owner_name,
            },
        )
    return validated_value


def _validate_optional_positive_finite_float(
    *,
    value: float | None,
    field_name: str,
    owner_name: str,
) -> float | None:
    return (
        _validate_positive_finite_float(
            value=value,
            field_name=field_name,
            owner_name=owner_name,
        )
        if value is not None
        else None
    )


def _validate_optional_finite_float(
    *,
    value: float | None,
    field_name: str,
    owner_name: str,
) -> float | None:
    return (
        _validate_finite_float(
            value=value,
            field_name=field_name,
            owner_name=owner_name,
        )
        if value is not None
        else None
    )


def _validate_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    validated_value = float(value)
    if not math.isfinite(validated_value):
        raise PhylogeneticsError(
            f"{owner_name} requires finite {field_name}",
            code="brownian_continuous_trait_parameter_nonfinite",
            details={
                "field_name": field_name,
                "field_value": value,
                "owner_name": owner_name,
            },
        )
    return validated_value


__all__ = [
    "BROWNIAN_CONTINUOUS_TRAIT_MODELS",
    "BrownianContinuousTraitModelDefinition",
    "BrownianContinuousTraitParameterSummary",
    "BrownianContinuousTraitPosteriorRow",
    "BrownianContinuousTraitProposalSchedule",
    "BrownianContinuousTraitRunReport",
    "build_brownian_continuous_trait_model_definition",
    "build_brownian_continuous_trait_proposal_schedule",
    "run_brownian_continuous_trait_metropolis_hastings",
]
