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
from bijux_phylogenetics.comparative.common import (
    build_brownian_covariance_matrix,
    build_ou_covariance_matrix,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

ORNSTEIN_UHLENBECK_CONTINUOUS_TRAIT_MODELS = ("ornstein-uhlenbeck",)

_MODEL_NAME = "ornstein-uhlenbeck"
_MODEL_PARAMETER_NAME = "continuous-trait-model"
_ALPHA_PARAMETER_NAME = "alpha"
_OPTIMUM_PARAMETER_NAME = "optimum"
_SIGMA_SQUARED_PARAMETER_NAME = "sigma-squared"
_MINIMUM_IDENTIFIABILITY_SAMPLE_SIZE = 5
_BOUNDARY_ALPHA_SCALE_FRACTION = 0.25
_BROAD_ALPHA_POSTERIOR_SCALE_MULTIPLIER = 1.0


@dataclass(frozen=True, slots=True)
class OrnsteinUhlenbeckContinuousTraitModelDefinition:
    """One validated fixed-topology Bayesian OU continuous-trait model."""

    alpha_prior: ContinuousTraitScalarPriorModel
    optimum_prior: ContinuousTraitLocationPriorModel
    sigma_squared_prior: ContinuousTraitScalarPriorModel
    initial_alpha: float | None = None
    initial_optimum: float | None = None
    initial_sigma_squared: float | None = None


@dataclass(frozen=True, slots=True)
class OrnsteinUhlenbeckContinuousTraitProposalSchedule:
    """One validated proposal schedule for Bayesian OU trait sampling."""

    alpha_move_weight: float
    alpha_log_scale_standard_deviation: float
    optimum_move_weight: float
    optimum_standard_deviation: float
    sigma_squared_move_weight: float
    sigma_squared_log_scale_standard_deviation: float


@dataclass(frozen=True, slots=True)
class OrnsteinUhlenbeckContinuousTraitPosteriorRow:
    """One sampled posterior row from a Bayesian OU trait chain."""

    sample_index: int
    iteration_index: int
    topology_id: str
    model_name: str
    total_log_prior: float
    log_likelihood: float
    posterior_log_score: float
    prior_component_log_priors: dict[str, float]
    alpha: float
    optimum: float
    sigma_squared: float


@dataclass(frozen=True, slots=True)
class OrnsteinUhlenbeckContinuousTraitParameterSummary:
    """One posterior summary for a Bayesian OU continuous-trait parameter."""

    parameter_name: str
    sample_count: int
    posterior_mean: float
    hpd_95_lower: float
    hpd_95_upper: float


@dataclass(frozen=True, slots=True)
class OrnsteinUhlenbeckContinuousTraitIdentifiabilityWarning:
    """One posterior-based warning that the OU parameter surface is weakly identified."""

    kind: str
    message: str
    observed_value: float
    threshold: float


@dataclass(frozen=True, slots=True)
class OrnsteinUhlenbeckContinuousTraitRunReport:
    """One completed fixed-topology Bayesian OU trait posterior run."""

    model_definition: OrnsteinUhlenbeckContinuousTraitModelDefinition
    proposal_schedule: OrnsteinUhlenbeckContinuousTraitProposalSchedule
    taxa: list[str]
    tip_values: dict[str, float]
    chain_report: MetropolisHastingsRunReport
    posterior_rows: list[OrnsteinUhlenbeckContinuousTraitPosteriorRow]
    parameter_summaries: list[OrnsteinUhlenbeckContinuousTraitParameterSummary]
    identifiability_warnings: list[
        OrnsteinUhlenbeckContinuousTraitIdentifiabilityWarning
    ]


def build_ornstein_uhlenbeck_continuous_trait_model_definition(
    *,
    alpha_prior: ContinuousTraitScalarPriorModel,
    optimum_prior: ContinuousTraitLocationPriorModel,
    sigma_squared_prior: ContinuousTraitScalarPriorModel,
    initial_alpha: float | None = None,
    initial_optimum: float | None = None,
    initial_sigma_squared: float | None = None,
) -> OrnsteinUhlenbeckContinuousTraitModelDefinition:
    """Build one validated Bayesian OU continuous-trait model definition."""
    if not isinstance(alpha_prior, ContinuousTraitScalarPriorModel):
        raise PhylogeneticsError(
            "bayesian OU continuous-trait model requires one ContinuousTraitScalarPriorModel for alpha_prior",
            code="ou_continuous_trait_alpha_prior_type_invalid",
        )
    if alpha_prior.family == "fixed":
        raise PhylogeneticsError(
            "bayesian OU continuous-trait model requires one non-fixed alpha prior so pull-to-optimum strength remains sampled",
            code="ou_continuous_trait_alpha_prior_fixed",
        )
    if not isinstance(optimum_prior, ContinuousTraitLocationPriorModel):
        raise PhylogeneticsError(
            "bayesian OU continuous-trait model requires one ContinuousTraitLocationPriorModel for optimum_prior",
            code="ou_continuous_trait_optimum_prior_type_invalid",
        )
    if optimum_prior.family == "fixed":
        raise PhylogeneticsError(
            "bayesian OU continuous-trait model requires one non-fixed optimum prior so the optimum remains sampled",
            code="ou_continuous_trait_optimum_prior_fixed",
        )
    if not isinstance(sigma_squared_prior, ContinuousTraitScalarPriorModel):
        raise PhylogeneticsError(
            "bayesian OU continuous-trait model requires one ContinuousTraitScalarPriorModel for sigma_squared_prior",
            code="ou_continuous_trait_sigma_squared_prior_type_invalid",
        )
    if sigma_squared_prior.family == "fixed":
        raise PhylogeneticsError(
            "bayesian OU continuous-trait model requires one non-fixed sigma-squared prior so diffusion variance remains sampled",
            code="ou_continuous_trait_sigma_squared_prior_fixed",
        )
    return OrnsteinUhlenbeckContinuousTraitModelDefinition(
        alpha_prior=alpha_prior,
        optimum_prior=optimum_prior,
        sigma_squared_prior=sigma_squared_prior,
        initial_alpha=_validate_optional_positive_finite_float(
            value=initial_alpha,
            field_name="initial_alpha",
            owner_name="bayesian OU continuous-trait model",
        ),
        initial_optimum=_validate_optional_finite_float(
            value=initial_optimum,
            field_name="initial_optimum",
            owner_name="bayesian OU continuous-trait model",
        ),
        initial_sigma_squared=_validate_optional_positive_finite_float(
            value=initial_sigma_squared,
            field_name="initial_sigma_squared",
            owner_name="bayesian OU continuous-trait model",
        ),
    )


def build_ornstein_uhlenbeck_continuous_trait_proposal_schedule(
    *,
    model_definition: OrnsteinUhlenbeckContinuousTraitModelDefinition,
    alpha_move_weight: float,
    alpha_log_scale_standard_deviation: float,
    optimum_move_weight: float,
    optimum_standard_deviation: float,
    sigma_squared_move_weight: float,
    sigma_squared_log_scale_standard_deviation: float,
) -> OrnsteinUhlenbeckContinuousTraitProposalSchedule:
    """Build one validated proposal schedule for Bayesian OU trait sampling."""
    if not isinstance(
        model_definition,
        OrnsteinUhlenbeckContinuousTraitModelDefinition,
    ):
        raise PhylogeneticsError(
            "bayesian OU continuous-trait proposal schedule requires one OrnsteinUhlenbeckContinuousTraitModelDefinition",
            code="ou_continuous_trait_proposal_schedule_model_definition_type_invalid",
        )
    return OrnsteinUhlenbeckContinuousTraitProposalSchedule(
        alpha_move_weight=_validate_positive_finite_float(
            value=alpha_move_weight,
            field_name="alpha_move_weight",
            owner_name="bayesian OU continuous-trait proposal schedule",
        ),
        alpha_log_scale_standard_deviation=_validate_positive_finite_float(
            value=alpha_log_scale_standard_deviation,
            field_name="alpha_log_scale_standard_deviation",
            owner_name="bayesian OU continuous-trait proposal schedule",
        ),
        optimum_move_weight=_validate_positive_finite_float(
            value=optimum_move_weight,
            field_name="optimum_move_weight",
            owner_name="bayesian OU continuous-trait proposal schedule",
        ),
        optimum_standard_deviation=_validate_positive_finite_float(
            value=optimum_standard_deviation,
            field_name="optimum_standard_deviation",
            owner_name="bayesian OU continuous-trait proposal schedule",
        ),
        sigma_squared_move_weight=_validate_positive_finite_float(
            value=sigma_squared_move_weight,
            field_name="sigma_squared_move_weight",
            owner_name="bayesian OU continuous-trait proposal schedule",
        ),
        sigma_squared_log_scale_standard_deviation=_validate_positive_finite_float(
            value=sigma_squared_log_scale_standard_deviation,
            field_name="sigma_squared_log_scale_standard_deviation",
            owner_name="bayesian OU continuous-trait proposal schedule",
        ),
    )


def run_ornstein_uhlenbeck_continuous_trait_metropolis_hastings(
    *,
    tree: PhyloTree,
    tip_values: Mapping[str, float],
    model_definition: OrnsteinUhlenbeckContinuousTraitModelDefinition,
    proposal_schedule: OrnsteinUhlenbeckContinuousTraitProposalSchedule,
    iteration_count: int,
    sample_every: int = 1,
    seed: int = 0,
) -> OrnsteinUhlenbeckContinuousTraitRunReport:
    """Run one fixed-topology Bayesian OU trait posterior sampler."""
    if not isinstance(tree, PhyloTree):
        raise PhylogeneticsError(
            "bayesian OU continuous-trait runner requires one PhyloTree",
            code="ou_continuous_trait_tree_type_invalid",
        )
    if not isinstance(
        model_definition,
        OrnsteinUhlenbeckContinuousTraitModelDefinition,
    ):
        raise PhylogeneticsError(
            "bayesian OU continuous-trait runner requires one OrnsteinUhlenbeckContinuousTraitModelDefinition",
            code="ou_continuous_trait_model_definition_type_invalid",
        )
    if not isinstance(
        proposal_schedule,
        OrnsteinUhlenbeckContinuousTraitProposalSchedule,
    ):
        raise PhylogeneticsError(
            "bayesian OU continuous-trait runner requires one OrnsteinUhlenbeckContinuousTraitProposalSchedule",
            code="ou_continuous_trait_proposal_schedule_type_invalid",
        )
    normalized_tree = tree.copy()
    normalized_tree.rooted = tree.rooted
    normalized_tree.refresh()
    if normalized_tree.rooted is not True:
        raise PhylogeneticsError(
            "bayesian OU continuous-trait runner requires one rooted tree",
            code="ou_continuous_trait_tree_rooting_invalid",
        )
    taxa, normalized_tip_values = _normalize_tip_values(
        tree=normalized_tree,
        tip_values=tip_values,
    )
    tip_value_vector = [normalized_tip_values[taxon] for taxon in taxa]
    reference_alpha_scale = _reference_alpha_scale(
        tree=normalized_tree,
        taxa=taxa,
    )
    initial_alpha = (
        model_definition.initial_alpha
        if model_definition.initial_alpha is not None
        else reference_alpha_scale
    )
    initial_covariance = _build_ou_covariance(
        tree=normalized_tree,
        taxa=taxa,
        alpha=initial_alpha,
    )
    initial_optimum = (
        model_definition.initial_optimum
        if model_definition.initial_optimum is not None
        else _estimate_gls_optimum(
            tip_values=tip_value_vector,
            inverse_covariance=invert_matrix(initial_covariance),
        )
    )
    initial_sigma_squared = (
        model_definition.initial_sigma_squared
        if model_definition.initial_sigma_squared is not None
        else _estimate_sigma_squared(
            tip_values=tip_value_vector,
            inverse_covariance=invert_matrix(initial_covariance),
            optimum=initial_optimum,
        )
    )
    initial_state = score_bayesian_phylogenetic_state(
        tree=normalized_tree,
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={_MODEL_PARAMETER_NAME: _MODEL_NAME},
            scalar_parameters={
                _ALPHA_PARAMETER_NAME: initial_alpha,
                _OPTIMUM_PARAMETER_NAME: initial_optimum,
                _SIGMA_SQUARED_PARAMETER_NAME: initial_sigma_squared,
            },
        ),
        update_prior_components=lambda state: _build_ou_prior_components(
            state=state,
            model_definition=model_definition,
            fixed_topology_id=None,
        ),
        update_log_likelihood=lambda state: _evaluate_ou_log_likelihood(
            state=state,
            fixed_topology_id=None,
            tree=normalized_tree,
            taxa=taxa,
            tip_values=tip_value_vector,
        ),
    )
    fixed_topology_id = initial_state.tree.topology_id
    chain_report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: _propose_ou_state(
            current_state=current_state,
            rng=rng,
            proposal_schedule=proposal_schedule,
        ),
        update_prior_components=lambda state: _build_ou_prior_components(
            state=state,
            model_definition=model_definition,
            fixed_topology_id=fixed_topology_id,
        ),
        update_log_likelihood=lambda state: _evaluate_ou_log_likelihood(
            state=state,
            fixed_topology_id=fixed_topology_id,
            tree=normalized_tree,
            taxa=taxa,
            tip_values=tip_value_vector,
        ),
        iteration_count=iteration_count,
        sample_every=sample_every,
        seed=seed,
    )
    posterior_rows = _build_ou_posterior_rows(
        chain_report=chain_report,
        fixed_topology_id=fixed_topology_id,
    )
    parameter_summaries = _build_parameter_summaries(posterior_rows)
    return OrnsteinUhlenbeckContinuousTraitRunReport(
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        taxa=taxa,
        tip_values=normalized_tip_values,
        chain_report=chain_report,
        posterior_rows=posterior_rows,
        parameter_summaries=parameter_summaries,
        identifiability_warnings=_build_identifiability_warnings(
            taxa=taxa,
            parameter_summaries=parameter_summaries,
            reference_alpha_scale=reference_alpha_scale,
        ),
    )


def _normalize_tip_values(
    *,
    tree: PhyloTree,
    tip_values: Mapping[str, float],
) -> tuple[list[str], dict[str, float]]:
    if not isinstance(tip_values, Mapping):
        raise PhylogeneticsError(
            "bayesian OU continuous-trait runner requires one mapping of tip values",
            code="ou_continuous_trait_tip_values_type_invalid",
        )
    taxa = list(tree.tip_names)
    if len(taxa) < 3:
        raise PhylogeneticsError(
            "bayesian OU continuous-trait runner requires at least three taxa",
            code="ou_continuous_trait_taxon_count_invalid",
            details={"taxon_count": len(taxa)},
        )
    tree_taxa = set(taxa)
    tip_value_taxa = set(tip_values)
    missing_taxa = sorted(tree_taxa - tip_value_taxa)
    extra_taxa = sorted(tip_value_taxa - tree_taxa)
    if missing_taxa or extra_taxa:
        raise PhylogeneticsError(
            "bayesian OU continuous-trait runner requires tip values to match the tree tip set exactly",
            code="ou_continuous_trait_tip_value_taxa_mismatch",
            details={
                "missing_taxa": missing_taxa,
                "extra_taxa": extra_taxa,
            },
        )
    normalized_tip_values = {
        taxon: _validate_finite_float(
            value=tip_values[taxon],
            field_name=taxon,
            owner_name="bayesian OU continuous-trait tip values",
        )
        for taxon in taxa
    }
    if len({format(value, ".12g") for value in normalized_tip_values.values()}) == 1:
        raise PhylogeneticsError(
            "bayesian OU continuous-trait runner requires non-constant tip values so alpha and sigma-squared remain identifiable under positive-support sampling",
            code="ou_continuous_trait_tip_values_constant",
        )
    return taxa, normalized_tip_values


def _reference_alpha_scale(*, tree: PhyloTree, taxa: list[str]) -> float:
    max_depth = max(max(row) for row in build_brownian_covariance_matrix(tree, taxa))
    return _validate_positive_finite_float(
        value=1.0 / max(max_depth, 1e-6),
        field_name="reference_alpha_scale",
        owner_name="bayesian OU continuous-trait runner",
    )


def _build_ou_covariance(
    *,
    tree: PhyloTree,
    taxa: list[str],
    alpha: float,
) -> list[list[float]]:
    return stable_covariance(build_ou_covariance_matrix(tree, taxa, alpha=alpha))


def _estimate_gls_optimum(
    *,
    tip_values: list[float],
    inverse_covariance: list[list[float]],
) -> float:
    ones = [1.0] * len(tip_values)
    denominator = quadratic_form(ones, inverse_covariance)
    if not math.isfinite(denominator) or denominator <= 0.0:
        raise PhylogeneticsError(
            "bayesian OU continuous-trait runner could not estimate one finite positive generalized-least-squares denominator for the initial optimum",
            code="ou_continuous_trait_initial_optimum_denominator_invalid",
            details={"denominator": denominator},
        )
    numerator = sum(
        ones[row_index]
        * sum(
            inverse_covariance[row_index][column_index] * tip_values[column_index]
            for column_index in range(len(tip_values))
        )
        for row_index in range(len(tip_values))
    )
    return _validate_finite_float(
        value=numerator / denominator,
        field_name="initial_optimum",
        owner_name="bayesian OU continuous-trait runner",
    )


def _estimate_sigma_squared(
    *,
    tip_values: list[float],
    inverse_covariance: list[list[float]],
    optimum: float,
) -> float:
    residuals = [value - optimum for value in tip_values]
    sigma_squared = quadratic_form(residuals, inverse_covariance) / len(tip_values)
    if not math.isfinite(sigma_squared) or sigma_squared <= 0.0:
        raise PhylogeneticsError(
            "bayesian OU continuous-trait runner requires one strictly positive initial sigma-squared value",
            code="ou_continuous_trait_initial_sigma_squared_invalid",
            details={"sigma_squared": sigma_squared},
        )
    return sigma_squared


def _build_ou_prior_components(
    *,
    state: BayesianPhylogeneticState,
    model_definition: OrnsteinUhlenbeckContinuousTraitModelDefinition,
    fixed_topology_id: str | None,
) -> list[BayesianPriorComponentState]:
    _require_ou_state_consistency(
        state=state,
        fixed_topology_id=fixed_topology_id,
    )
    alpha = state.model_parameters.scalar_parameters[_ALPHA_PARAMETER_NAME]
    optimum = state.model_parameters.scalar_parameters[_OPTIMUM_PARAMETER_NAME]
    sigma_squared = state.model_parameters.scalar_parameters[
        _SIGMA_SQUARED_PARAMETER_NAME
    ]
    alpha_log_prior = evaluate_continuous_trait_scalar_log_prior(
        parameter_value=alpha,
        prior_model=model_definition.alpha_prior,
        parameter_name=_ALPHA_PARAMETER_NAME,
        allow_zero=False,
    )
    optimum_log_prior = evaluate_continuous_trait_location_log_prior(
        parameter_value=optimum,
        prior_model=model_definition.optimum_prior,
        parameter_name=_OPTIMUM_PARAMETER_NAME,
    )
    sigma_squared_log_prior = evaluate_continuous_trait_scalar_log_prior(
        parameter_value=sigma_squared,
        prior_model=model_definition.sigma_squared_prior,
        parameter_name=_SIGMA_SQUARED_PARAMETER_NAME,
        allow_zero=False,
    )
    return [
        build_bayesian_prior_component_state(
            component_name="continuous-trait:alpha",
            family=model_definition.alpha_prior.family,
            log_prior=alpha_log_prior,
            parameter_values=model_definition.alpha_prior.parameter_values(),
        ),
        build_bayesian_prior_component_state(
            component_name="continuous-trait:optimum",
            family=model_definition.optimum_prior.family,
            log_prior=optimum_log_prior,
            parameter_values=model_definition.optimum_prior.parameter_values(),
        ),
        build_bayesian_prior_component_state(
            component_name="continuous-trait:sigma-squared",
            family=model_definition.sigma_squared_prior.family,
            log_prior=sigma_squared_log_prior,
            parameter_values=model_definition.sigma_squared_prior.parameter_values(),
        ),
    ]


def _evaluate_ou_log_likelihood(
    *,
    state: BayesianPhylogeneticState,
    fixed_topology_id: str | None,
    tree: PhyloTree,
    taxa: list[str],
    tip_values: list[float],
) -> float:
    _require_ou_state_consistency(
        state=state,
        fixed_topology_id=fixed_topology_id,
    )
    alpha = _validate_positive_finite_float(
        value=state.model_parameters.scalar_parameters[_ALPHA_PARAMETER_NAME],
        field_name=_ALPHA_PARAMETER_NAME,
        owner_name="bayesian OU continuous-trait log likelihood",
    )
    optimum = _validate_finite_float(
        value=state.model_parameters.scalar_parameters[_OPTIMUM_PARAMETER_NAME],
        field_name=_OPTIMUM_PARAMETER_NAME,
        owner_name="bayesian OU continuous-trait log likelihood",
    )
    sigma_squared = _validate_positive_finite_float(
        value=state.model_parameters.scalar_parameters[_SIGMA_SQUARED_PARAMETER_NAME],
        field_name=_SIGMA_SQUARED_PARAMETER_NAME,
        owner_name="bayesian OU continuous-trait log likelihood",
    )
    covariance = _build_ou_covariance(tree=tree, taxa=taxa, alpha=alpha)
    inverse_covariance = invert_matrix(covariance)
    residuals = [value - optimum for value in tip_values]
    quadratic_term = quadratic_form(residuals, inverse_covariance)
    return -0.5 * (
        len(tip_values) * math.log(2.0 * math.pi)
        + len(tip_values) * math.log(sigma_squared)
        + log_determinant(covariance)
        + (quadratic_term / sigma_squared)
    )


def _propose_ou_state(
    *,
    current_state: BayesianPhylogeneticState,
    rng,
    proposal_schedule: OrnsteinUhlenbeckContinuousTraitProposalSchedule,
):
    weighted_moves = [
        (
            proposal_schedule.alpha_move_weight,
            lambda: propose_clock_rate_move(
                current_state,
                rng,
                log_scale_standard_deviation=(
                    proposal_schedule.alpha_log_scale_standard_deviation
                ),
                parameter_name=_ALPHA_PARAMETER_NAME,
            ),
        ),
        (
            proposal_schedule.optimum_move_weight,
            lambda: propose_continuous_trait_location_move(
                current_state,
                rng,
                standard_deviation=proposal_schedule.optimum_standard_deviation,
                parameter_name=_OPTIMUM_PARAMETER_NAME,
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


def _build_ou_posterior_rows(
    *,
    chain_report: MetropolisHastingsRunReport,
    fixed_topology_id: str,
) -> list[OrnsteinUhlenbeckContinuousTraitPosteriorRow]:
    posterior_rows: list[OrnsteinUhlenbeckContinuousTraitPosteriorRow] = []
    for sample_index, state in enumerate(chain_report.sampled_states):
        if state.tree.topology_id != fixed_topology_id:
            raise PhylogeneticsError(
                "bayesian OU continuous-trait posterior trace detected one topology change in sampled states",
                code="ou_continuous_trait_trace_topology_changed",
                details={
                    "expected_topology_id": fixed_topology_id,
                    "observed_topology_id": state.tree.topology_id,
                    "sample_index": sample_index,
                },
            )
        posterior_rows.append(
            OrnsteinUhlenbeckContinuousTraitPosteriorRow(
                sample_index=sample_index,
                iteration_index=sample_index * chain_report.sample_every,
                topology_id=state.tree.topology_id,
                model_name=_MODEL_NAME,
                total_log_prior=state.total_log_prior,
                log_likelihood=state.log_likelihood,
                posterior_log_score=state.posterior_log_score,
                prior_component_log_priors={
                    component.component_name: component.log_prior
                    for component in state.prior_components
                },
                alpha=state.model_parameters.scalar_parameters[_ALPHA_PARAMETER_NAME],
                optimum=state.model_parameters.scalar_parameters[
                    _OPTIMUM_PARAMETER_NAME
                ],
                sigma_squared=state.model_parameters.scalar_parameters[
                    _SIGMA_SQUARED_PARAMETER_NAME
                ],
            )
        )
    return posterior_rows


def _build_parameter_summaries(
    posterior_rows: list[OrnsteinUhlenbeckContinuousTraitPosteriorRow],
) -> list[OrnsteinUhlenbeckContinuousTraitParameterSummary]:
    return [
        _build_parameter_summary(
            parameter_name=_ALPHA_PARAMETER_NAME,
            values=[row.alpha for row in posterior_rows],
        ),
        _build_parameter_summary(
            parameter_name=_OPTIMUM_PARAMETER_NAME,
            values=[row.optimum for row in posterior_rows],
        ),
        _build_parameter_summary(
            parameter_name=_SIGMA_SQUARED_PARAMETER_NAME,
            values=[row.sigma_squared for row in posterior_rows],
        ),
    ]


def _build_parameter_summary(
    *,
    parameter_name: str,
    values: list[float],
) -> OrnsteinUhlenbeckContinuousTraitParameterSummary:
    hpd_95_lower, hpd_95_upper = highest_posterior_density_interval(values)
    return OrnsteinUhlenbeckContinuousTraitParameterSummary(
        parameter_name=parameter_name,
        sample_count=len(values),
        posterior_mean=round(mean(values), 6),
        hpd_95_lower=round(hpd_95_lower, 6),
        hpd_95_upper=round(hpd_95_upper, 6),
    )


def _build_identifiability_warnings(
    *,
    taxa: list[str],
    parameter_summaries: list[OrnsteinUhlenbeckContinuousTraitParameterSummary],
    reference_alpha_scale: float,
) -> list[OrnsteinUhlenbeckContinuousTraitIdentifiabilityWarning]:
    alpha_summary = next(
        summary
        for summary in parameter_summaries
        if summary.parameter_name == _ALPHA_PARAMETER_NAME
    )
    warnings: list[OrnsteinUhlenbeckContinuousTraitIdentifiabilityWarning] = []
    if len(taxa) < _MINIMUM_IDENTIFIABILITY_SAMPLE_SIZE:
        warnings.append(
            OrnsteinUhlenbeckContinuousTraitIdentifiabilityWarning(
                kind="small_sample_size",
                message="OU alpha is hard to identify with fewer than five taxa",
                observed_value=float(len(taxa)),
                threshold=float(_MINIMUM_IDENTIFIABILITY_SAMPLE_SIZE),
            )
        )
    boundary_threshold = reference_alpha_scale * _BOUNDARY_ALPHA_SCALE_FRACTION
    if alpha_summary.hpd_95_lower <= boundary_threshold:
        warnings.append(
            OrnsteinUhlenbeckContinuousTraitIdentifiabilityWarning(
                kind="boundary_alpha",
                message="OU alpha posterior reaches the Brownian boundary region near zero and may not be well identified",
                observed_value=alpha_summary.hpd_95_lower,
                threshold=boundary_threshold,
            )
        )
    posterior_width = alpha_summary.hpd_95_upper - alpha_summary.hpd_95_lower
    broad_posterior_threshold = (
        reference_alpha_scale * _BROAD_ALPHA_POSTERIOR_SCALE_MULTIPLIER
    )
    if posterior_width >= broad_posterior_threshold:
        warnings.append(
            OrnsteinUhlenbeckContinuousTraitIdentifiabilityWarning(
                kind="broad_alpha_posterior",
                message="OU alpha posterior remains broad relative to the tree depth scale, so pull-to-optimum strength is weakly identified",
                observed_value=posterior_width,
                threshold=broad_posterior_threshold,
            )
        )
    if alpha_summary.posterior_mean <= reference_alpha_scale:
        warnings.append(
            OrnsteinUhlenbeckContinuousTraitIdentifiabilityWarning(
                kind="weak_pull_to_optimum",
                message="OU alpha posterior is concentrated on weak pull-to-optimum values that are difficult to distinguish from Brownian motion",
                observed_value=alpha_summary.posterior_mean,
                threshold=reference_alpha_scale,
            )
        )
    return warnings


def _require_ou_state_consistency(
    *,
    state: BayesianPhylogeneticState,
    fixed_topology_id: str | None,
) -> None:
    model_name = state.model_parameters.categorical_parameters.get(
        _MODEL_PARAMETER_NAME
    )
    if model_name != _MODEL_NAME:
        raise PhylogeneticsError(
            "bayesian OU continuous-trait model requires every sampled state to preserve the configured model label",
            code="ou_continuous_trait_state_model_label_invalid",
            details={
                "expected_model_name": _MODEL_NAME,
                "observed_model_name": model_name,
            },
        )
    if fixed_topology_id is not None and state.tree.topology_id != fixed_topology_id:
        raise PhylogeneticsError(
            "bayesian OU continuous-trait model requires topology to remain unchanged across sampled states",
            code="ou_continuous_trait_state_topology_changed",
            details={
                "expected_topology_id": fixed_topology_id,
                "observed_topology_id": state.tree.topology_id,
            },
        )
    for parameter_name in (
        _ALPHA_PARAMETER_NAME,
        _OPTIMUM_PARAMETER_NAME,
        _SIGMA_SQUARED_PARAMETER_NAME,
    ):
        if parameter_name not in state.model_parameters.scalar_parameters:
            raise PhylogeneticsError(
                "bayesian OU continuous-trait model requires one complete scalar parameter state in every sampled state",
                code="ou_continuous_trait_state_scalar_parameter_missing",
                details={"parameter_name": parameter_name},
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
            code="ou_continuous_trait_parameter_nonpositive",
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
            code="ou_continuous_trait_parameter_nonfinite",
            details={
                "field_name": field_name,
                "field_value": value,
                "owner_name": owner_name,
            },
        )
    return validated_value


__all__ = [
    "ORNSTEIN_UHLENBECK_CONTINUOUS_TRAIT_MODELS",
    "OrnsteinUhlenbeckContinuousTraitIdentifiabilityWarning",
    "OrnsteinUhlenbeckContinuousTraitModelDefinition",
    "OrnsteinUhlenbeckContinuousTraitParameterSummary",
    "OrnsteinUhlenbeckContinuousTraitPosteriorRow",
    "OrnsteinUhlenbeckContinuousTraitProposalSchedule",
    "OrnsteinUhlenbeckContinuousTraitRunReport",
    "build_ornstein_uhlenbeck_continuous_trait_model_definition",
    "build_ornstein_uhlenbeck_continuous_trait_proposal_schedule",
    "run_ornstein_uhlenbeck_continuous_trait_metropolis_hastings",
]
