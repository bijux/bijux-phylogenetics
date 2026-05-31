from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math

import numpy

from bijux_phylogenetics.ancestral.common import node_descendant_taxa, node_signature
from bijux_phylogenetics.ancestral.discrete.likelihood.likelihood_math import (
    tree_log_likelihood,
)
from bijux_phylogenetics.ancestral.discrete.likelihood.posterior_probabilities import (
    estimate_marginal_state_probabilities,
)
from bijux_phylogenetics.ancestral.discrete.likelihood.rate_matrix import (
    DiscreteTransitionRateRow,
    build_transition_rate_rows,
    rate_matrix_from_log_parameters,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_allowed_transition_pairs,
    resolve_discrete_model_name,
    resolve_root_prior,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_parameters import (
    DISCRETE_TRAIT_RATE_PARAMETER_MODELS,
    parameterize_discrete_trait_rate_rows,
    resolve_discrete_trait_rate_rows,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    DiscreteTraitRatePriorModel,
    evaluate_discrete_trait_rate_log_prior,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsRunReport,
    propose_discrete_trait_rate_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
    build_bayesian_prior_component_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

DISCRETE_TRAIT_MK_MODELS = DISCRETE_TRAIT_RATE_PARAMETER_MODELS
DISCRETE_TRAIT_MK_ROOT_PRIOR_MODES = ("equal", "empirical", "fixed")

_DISCRETE_TRAIT_RATE_PARAMETER_NAME = "discrete-trait-rates"
_DISCRETE_TRAIT_STATE_ORDERING = "unordered"


@dataclass(frozen=True, slots=True)
class DiscreteTraitMkModelDefinition:
    """One validated fixed-topology Bayesian Mk model definition."""

    transition_model_name: str
    rate_prior: DiscreteTraitRatePriorModel
    root_prior_mode: str
    fixed_root_state: str | None = None
    initial_rate: float = 1.0


@dataclass(frozen=True, slots=True)
class DiscreteTraitMkProposalSchedule:
    """One validated proposal schedule for Bayesian Mk transition-rate sampling."""

    transition_model_name: str
    rate_log_scale_standard_deviation: float


@dataclass(frozen=True, slots=True)
class DiscreteTraitMkNodeStateSummary:
    """One internal-node posterior state summary for a sampled Mk state."""

    node_id: str
    node_name: str | None
    descendant_taxa: list[str]
    most_likely_state: str
    state_probabilities: dict[str, float]


@dataclass(frozen=True, slots=True)
class DiscreteTraitMkPosteriorRow:
    """One sampled posterior row from a Bayesian Mk chain."""

    sample_index: int
    iteration_index: int
    topology_id: str
    transition_model_name: str
    total_log_prior: float
    log_likelihood: float
    posterior_log_score: float
    prior_component_log_priors: dict[str, float]
    rate_parameters: dict[str, float]
    node_state_summaries: list[DiscreteTraitMkNodeStateSummary]


@dataclass(frozen=True, slots=True)
class DiscreteTraitMkRunReport:
    """One completed fixed-topology Bayesian Mk posterior run."""

    model_definition: DiscreteTraitMkModelDefinition
    proposal_schedule: DiscreteTraitMkProposalSchedule
    state_order: list[str]
    taxa: list[str]
    chain_report: MetropolisHastingsRunReport
    posterior_rows: list[DiscreteTraitMkPosteriorRow]


def build_discrete_trait_mk_model_definition(
    *,
    transition_model_name: str,
    rate_prior: DiscreteTraitRatePriorModel,
    root_prior_mode: str = "equal",
    fixed_root_state: str | None = None,
    initial_rate: float = 1.0,
) -> DiscreteTraitMkModelDefinition:
    """Build one validated Bayesian discrete-trait Mk model definition."""
    validated_transition_model_name = _validate_transition_model_name(
        transition_model_name
    )
    if not isinstance(rate_prior, DiscreteTraitRatePriorModel):
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk model requires one DiscreteTraitRatePriorModel",
            code="discrete_trait_mk_rate_prior_type_invalid",
        )
    validated_root_prior_mode = _validate_root_prior_mode(root_prior_mode)
    validated_fixed_root_state = _validate_optional_nonblank_state_name(
        fixed_root_state,
        field_name="fixed_root_state",
        owner_name="bayesian discrete-trait Mk model",
    )
    if validated_root_prior_mode != "fixed" and validated_fixed_root_state is not None:
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk model accepts 'fixed_root_state' only when root_prior_mode is 'fixed'",
            code="discrete_trait_mk_fixed_root_state_unused",
            details={"root_prior_mode": validated_root_prior_mode},
        )
    if validated_root_prior_mode == "fixed" and validated_fixed_root_state is None:
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk model requires 'fixed_root_state' when root_prior_mode is 'fixed'",
            code="discrete_trait_mk_fixed_root_state_missing",
        )
    return DiscreteTraitMkModelDefinition(
        transition_model_name=validated_transition_model_name,
        rate_prior=rate_prior,
        root_prior_mode=validated_root_prior_mode,
        fixed_root_state=validated_fixed_root_state,
        initial_rate=_validate_positive_finite_float(
            value=initial_rate,
            field_name="initial_rate",
            owner_name="bayesian discrete-trait Mk model",
        ),
    )


def build_discrete_trait_mk_proposal_schedule(
    *,
    model_definition: DiscreteTraitMkModelDefinition,
    rate_log_scale_standard_deviation: float,
) -> DiscreteTraitMkProposalSchedule:
    """Build one validated proposal schedule for Bayesian Mk sampling."""
    if not isinstance(model_definition, DiscreteTraitMkModelDefinition):
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk proposal schedule requires one DiscreteTraitMkModelDefinition",
            code="discrete_trait_mk_proposal_schedule_model_definition_type_invalid",
        )
    return DiscreteTraitMkProposalSchedule(
        transition_model_name=model_definition.transition_model_name,
        rate_log_scale_standard_deviation=_validate_positive_finite_float(
            value=rate_log_scale_standard_deviation,
            field_name="rate_log_scale_standard_deviation",
            owner_name="bayesian discrete-trait Mk proposal schedule",
        ),
    )


def run_discrete_trait_mk_metropolis_hastings(
    *,
    tree: PhyloTree,
    tip_states: Mapping[str, str],
    model_definition: DiscreteTraitMkModelDefinition,
    proposal_schedule: DiscreteTraitMkProposalSchedule,
    iteration_count: int,
    sample_every: int = 1,
    seed: int = 0,
) -> DiscreteTraitMkRunReport:
    """Run one Bayesian Mk Metropolis-Hastings chain on a fixed rooted tree."""
    if not isinstance(tree, PhyloTree):
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk runner requires one PhyloTree",
            code="discrete_trait_mk_tree_type_invalid",
        )
    if not isinstance(model_definition, DiscreteTraitMkModelDefinition):
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk runner requires one DiscreteTraitMkModelDefinition",
            code="discrete_trait_mk_model_definition_type_invalid",
        )
    if not isinstance(proposal_schedule, DiscreteTraitMkProposalSchedule):
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk runner requires one DiscreteTraitMkProposalSchedule",
            code="discrete_trait_mk_proposal_schedule_type_invalid",
        )
    if (
        proposal_schedule.transition_model_name
        != model_definition.transition_model_name
    ):
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk runner requires the proposal schedule and model definition to use the same transition model",
            code="discrete_trait_mk_model_schedule_mismatch",
            details={
                "model_definition": model_definition.transition_model_name,
                "proposal_schedule": proposal_schedule.transition_model_name,
            },
        )
    normalized_tree = tree.copy()
    normalized_tree.rooted = tree.rooted
    normalized_tree.refresh()
    if len(normalized_tree.root.children) != 2:
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk runner requires one rooted bifurcating tree",
            code="discrete_trait_mk_tree_rooting_invalid",
        )
    normalized_tip_states, state_order, state_counts = _normalize_tip_states(
        tree=normalized_tree,
        tip_states=tip_states,
    )
    allowed_transition_pairs = resolve_allowed_transition_pairs(
        state_order,
        model=model_definition.transition_model_name,
        state_ordering=_DISCRETE_TRAIT_STATE_ORDERING,
        allowed_transition_pairs=None,
    )
    initial_transition_rate_rows = _build_initial_transition_rate_rows(
        state_order=state_order,
        model_definition=model_definition,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    initial_parameterization = parameterize_discrete_trait_rate_rows(
        model=model_definition.transition_model_name,
        transition_rate_rows=initial_transition_rate_rows,
    )
    initial_model_parameters = build_bayesian_model_parameter_state(
        categorical_parameters={
            "discrete-trait-model": model_definition.transition_model_name,
            "root-prior-mode": model_definition.root_prior_mode,
        },
        vector_parameters={
            _DISCRETE_TRAIT_RATE_PARAMETER_NAME: initial_parameterization.parameter_values
        },
    )
    chain_report = run_metropolis_hastings_sampler(
        initial_state=score_bayesian_phylogenetic_state(
            tree=normalized_tree,
            model_parameters=initial_model_parameters,
            update_prior_components=lambda state: (
                _build_discrete_trait_mk_prior_components(
                    state=state,
                    model_definition=model_definition,
                    state_order=state_order,
                    allowed_transition_pairs=allowed_transition_pairs,
                )
            ),
            update_log_likelihood=lambda state: (
                _evaluate_discrete_trait_mk_log_likelihood(
                    state=state,
                    tip_states=normalized_tip_states,
                    model_definition=model_definition,
                    state_order=state_order,
                    state_counts=state_counts,
                    allowed_transition_pairs=allowed_transition_pairs,
                )
            ),
        ),
        propose_state=lambda current_state, rng: propose_discrete_trait_rate_move(
            current_state,
            rng,
            log_scale_standard_deviation=(
                proposal_schedule.rate_log_scale_standard_deviation
            ),
            parameter_name=_DISCRETE_TRAIT_RATE_PARAMETER_NAME,
        ),
        update_prior_components=lambda state: _build_discrete_trait_mk_prior_components(
            state=state,
            model_definition=model_definition,
            state_order=state_order,
            allowed_transition_pairs=allowed_transition_pairs,
        ),
        update_log_likelihood=lambda state: _evaluate_discrete_trait_mk_log_likelihood(
            state=state,
            tip_states=normalized_tip_states,
            model_definition=model_definition,
            state_order=state_order,
            state_counts=state_counts,
            allowed_transition_pairs=allowed_transition_pairs,
        ),
        iteration_count=iteration_count,
        sample_every=sample_every,
        seed=seed,
    )
    posterior_rows = _build_discrete_trait_mk_posterior_rows(
        chain_report=chain_report,
        transition_model_name=model_definition.transition_model_name,
        tip_states=normalized_tip_states,
        state_order=state_order,
        state_counts=state_counts,
        model_definition=model_definition,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    return DiscreteTraitMkRunReport(
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        state_order=state_order,
        taxa=normalized_tree.tip_names,
        chain_report=chain_report,
        posterior_rows=posterior_rows,
    )


def _build_initial_transition_rate_rows(
    *,
    state_order: list[str],
    model_definition: DiscreteTraitMkModelDefinition,
    allowed_transition_pairs: set[tuple[int, int]],
) -> list[DiscreteTransitionRateRow]:
    parameter_count = _resolve_parameter_count(
        model_name=model_definition.transition_model_name,
        state_order=state_order,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    initial_log_parameters = numpy.full(
        parameter_count,
        math.log(model_definition.initial_rate),
        dtype=float,
    )
    initial_rate_matrix = rate_matrix_from_log_parameters(
        initial_log_parameters,
        state_order=state_order,
        model=model_definition.transition_model_name,
        state_ordering=_DISCRETE_TRAIT_STATE_ORDERING,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    return build_transition_rate_rows(
        state_order=state_order,
        state_ordering=_DISCRETE_TRAIT_STATE_ORDERING,
        rate_matrix=initial_rate_matrix,
        allowed_transition_pairs=allowed_transition_pairs,
    )


def _build_discrete_trait_mk_prior_components(
    *,
    state: BayesianPhylogeneticState,
    model_definition: DiscreteTraitMkModelDefinition,
    state_order: list[str],
    allowed_transition_pairs: set[tuple[int, int]],
) -> list[BayesianPriorComponentState]:
    transition_rate_rows = _resolve_transition_rate_rows_from_state(
        state=state,
        model_definition=model_definition,
        state_order=state_order,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    prior_report = evaluate_discrete_trait_rate_log_prior(
        model=model_definition.transition_model_name,
        transition_rate_rows=transition_rate_rows,
        prior_model=model_definition.rate_prior,
    )
    return [
        build_bayesian_prior_component_state(
            component_name=f"discrete-trait-rate:{row.parameter_name}",
            family=prior_report.family,
            log_prior=row.log_prior_contribution,
            parameter_values={
                "rate_value": row.rate_value,
                "transition_pair_count": float(len(row.transition_pairs)),
            },
        )
        for row in prior_report.rows
    ]


def _evaluate_discrete_trait_mk_log_likelihood(
    *,
    state: BayesianPhylogeneticState,
    tip_states: dict[str, str],
    model_definition: DiscreteTraitMkModelDefinition,
    state_order: list[str],
    state_counts: dict[str, int],
    allowed_transition_pairs: set[tuple[int, int]],
) -> float:
    transition_rate_rows = _resolve_transition_rate_rows_from_state(
        state=state,
        model_definition=model_definition,
        state_order=state_order,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    rate_matrix = _build_rate_matrix_from_transition_rows(
        state_order=state_order,
        transition_rate_rows=transition_rate_rows,
    )
    root_prior = _resolve_root_prior_for_state(
        state_order=state_order,
        state_counts=state_counts,
        model_definition=model_definition,
    )
    return float(
        format(
            tree_log_likelihood(
                state.tree.to_tree(),
                tip_states,
                state_order=state_order,
                rate_matrix=rate_matrix,
                root_prior=root_prior,
            ),
            ".15g",
        )
    )


def _build_discrete_trait_mk_posterior_rows(
    *,
    chain_report: MetropolisHastingsRunReport,
    transition_model_name: str,
    tip_states: dict[str, str],
    state_order: list[str],
    state_counts: dict[str, int],
    model_definition: DiscreteTraitMkModelDefinition,
    allowed_transition_pairs: set[tuple[int, int]],
) -> list[DiscreteTraitMkPosteriorRow]:
    posterior_rows: list[DiscreteTraitMkPosteriorRow] = []
    for sample_index, state in enumerate(chain_report.sampled_states):
        node_state_summaries = _build_node_state_summaries(
            state=state,
            tip_states=tip_states,
            state_order=state_order,
            state_counts=state_counts,
            model_definition=model_definition,
            allowed_transition_pairs=allowed_transition_pairs,
        )
        posterior_rows.append(
            DiscreteTraitMkPosteriorRow(
                sample_index=sample_index,
                iteration_index=sample_index * chain_report.sample_every,
                topology_id=state.tree.topology_id,
                transition_model_name=transition_model_name,
                total_log_prior=state.total_log_prior,
                log_likelihood=state.log_likelihood,
                posterior_log_score=state.posterior_log_score,
                prior_component_log_priors={
                    component.component_name: component.log_prior
                    for component in state.prior_components
                },
                rate_parameters=dict(
                    state.model_parameters.vector_parameters[
                        _DISCRETE_TRAIT_RATE_PARAMETER_NAME
                    ]
                ),
                node_state_summaries=node_state_summaries,
            )
        )
    return posterior_rows


def _build_node_state_summaries(
    *,
    state: BayesianPhylogeneticState,
    tip_states: dict[str, str],
    state_order: list[str],
    state_counts: dict[str, int],
    model_definition: DiscreteTraitMkModelDefinition,
    allowed_transition_pairs: set[tuple[int, int]],
) -> list[DiscreteTraitMkNodeStateSummary]:
    transition_rate_rows = _resolve_transition_rate_rows_from_state(
        state=state,
        model_definition=model_definition,
        state_order=state_order,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    rate_matrix = _build_rate_matrix_from_transition_rows(
        state_order=state_order,
        transition_rate_rows=transition_rate_rows,
    )
    root_prior = _resolve_root_prior_for_state(
        state_order=state_order,
        state_counts=state_counts,
        model_definition=model_definition,
    )
    tree = state.tree.to_tree()
    posterior_by_node = estimate_marginal_state_probabilities(
        tree,
        tip_states,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )
    summaries: list[DiscreteTraitMkNodeStateSummary] = []
    for node in tree.iter_internal_nodes(order="preorder"):
        summary_node_id = node_signature(node)
        state_probabilities = posterior_by_node[summary_node_id]
        summaries.append(
            DiscreteTraitMkNodeStateSummary(
                node_id=summary_node_id,
                node_name=node.name,
                descendant_taxa=node_descendant_taxa(node),
                most_likely_state=_select_most_likely_state(
                    state_probabilities,
                    state_order=state_order,
                ),
                state_probabilities=dict(state_probabilities),
            )
        )
    return summaries


def _resolve_transition_rate_rows_from_state(
    *,
    state: BayesianPhylogeneticState,
    model_definition: DiscreteTraitMkModelDefinition,
    state_order: list[str],
    allowed_transition_pairs: set[tuple[int, int]],
) -> list[DiscreteTransitionRateRow]:
    _validate_state_model_consistency(
        state=state,
        model_definition=model_definition,
    )
    transition_rate_rows = _build_initial_transition_rate_rows(
        state_order=state_order,
        model_definition=model_definition,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    try:
        return resolve_discrete_trait_rate_rows(
            model=model_definition.transition_model_name,
            transition_rate_rows=transition_rate_rows,
            parameter_values=state.model_parameters.vector_parameters[
                _DISCRETE_TRAIT_RATE_PARAMETER_NAME
            ],
        )
    except KeyError as error:
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk state requires one 'discrete-trait-rates' vector parameter",
            code="discrete_trait_mk_state_rate_parameters_missing",
        ) from error


def _validate_state_model_consistency(
    *,
    state: BayesianPhylogeneticState,
    model_definition: DiscreteTraitMkModelDefinition,
) -> None:
    state_model_name = state.model_parameters.categorical_parameters.get(
        "discrete-trait-model"
    )
    if state_model_name != model_definition.transition_model_name:
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk state carries one transition model that does not match the model definition",
            code="discrete_trait_mk_state_model_mismatch",
            details={
                "state_model_name": state_model_name,
                "model_definition": model_definition.transition_model_name,
            },
        )


def _build_rate_matrix_from_transition_rows(
    *,
    state_order: list[str],
    transition_rate_rows: list[DiscreteTransitionRateRow],
) -> numpy.ndarray:
    state_index = {state: index for index, state in enumerate(state_order)}
    rate_matrix = numpy.zeros((len(state_order), len(state_order)), dtype=float)
    for row in transition_rate_rows:
        if not row.transition_allowed:
            continue
        left_index = state_index[row.source_state]
        right_index = state_index[row.target_state]
        rate_matrix[left_index, right_index] = row.rate
    for state_position in range(len(state_order)):
        rate_matrix[state_position, state_position] = -float(
            numpy.sum(rate_matrix[state_position, :])
        )
    return rate_matrix


def _resolve_root_prior_for_state(
    *,
    state_order: list[str],
    state_counts: dict[str, int],
    model_definition: DiscreteTraitMkModelDefinition,
) -> numpy.ndarray:
    try:
        return resolve_root_prior(
            state_order,
            state_counts=state_counts,
            mode=model_definition.root_prior_mode,
            fixed_root_state=model_definition.fixed_root_state,
            default_root_prior=None,
        )
    except ValueError as error:
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk root-prior configuration is invalid for the analyzed state vocabulary",
            code="discrete_trait_mk_root_prior_invalid",
            details={
                "root_prior_mode": model_definition.root_prior_mode,
                "fixed_root_state": model_definition.fixed_root_state,
            },
        ) from error


def _resolve_parameter_count(
    *,
    model_name: str,
    state_order: list[str],
    allowed_transition_pairs: set[tuple[int, int]],
) -> int:
    if model_name == "equal-rates":
        return 1
    if model_name == "symmetric":
        return sum(
            1
            for left_index in range(len(state_order))
            for right_index in range(left_index + 1, len(state_order))
            if (left_index, right_index) in allowed_transition_pairs
            and (right_index, left_index) in allowed_transition_pairs
        )
    return len(allowed_transition_pairs)


def _normalize_tip_states(
    *,
    tree: PhyloTree,
    tip_states: Mapping[str, str],
) -> tuple[dict[str, str], list[str], dict[str, int]]:
    if not isinstance(tip_states, Mapping):
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk runner requires one mapping of tip states",
            code="discrete_trait_mk_tip_states_type_invalid",
        )
    normalized_tip_states = {
        _validate_nonblank_state_name(
            taxon_name,
            field_name=f"tip_states[{taxon_name!r}]",
            owner_name="bayesian discrete-trait Mk runner",
        ): _validate_nonblank_state_name(
            state_name,
            field_name=f"tip_states[{taxon_name!r}]",
            owner_name="bayesian discrete-trait Mk runner",
        )
        for taxon_name, state_name in tip_states.items()
    }
    expected_taxa = set(tree.tip_names)
    provided_taxa = set(normalized_tip_states)
    missing_taxa = sorted(expected_taxa - provided_taxa)
    unexpected_taxa = sorted(provided_taxa - expected_taxa)
    if missing_taxa or unexpected_taxa:
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk runner requires tip states to match the tree tip set exactly",
            code="discrete_trait_mk_tip_states_taxa_mismatch",
            details={
                "missing_taxa": missing_taxa,
                "unexpected_taxa": unexpected_taxa,
            },
        )
    observed_states = sorted(set(normalized_tip_states.values()))
    if len(tree.tip_names) < 2:
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk runner requires at least two named tips",
            code="discrete_trait_mk_tip_count_invalid",
        )
    if len(observed_states) < 2:
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk runner requires at least two observed states",
            code="discrete_trait_mk_observed_state_count_invalid",
        )
    state_counts = {
        state_name: sum(
            1
            for candidate_state in normalized_tip_states.values()
            if candidate_state == state_name
        )
        for state_name in observed_states
    }
    return normalized_tip_states, observed_states, state_counts


def _select_most_likely_state(
    state_probabilities: Mapping[str, float],
    *,
    state_order: list[str],
) -> str:
    return max(
        state_order,
        key=lambda state: (
            state_probabilities.get(state, float("-inf")),
            -state_order.index(state),
        ),
    )


def _validate_transition_model_name(model_name: str) -> str:
    try:
        resolved_model_name = resolve_discrete_model_name(model_name.strip().lower())
    except ValueError as error:
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk model supports only ER, SYM, and ARD transition models",
            code="discrete_trait_mk_model_invalid",
            details={"model_name": model_name},
        ) from error
    if resolved_model_name not in DISCRETE_TRAIT_MK_MODELS:
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk model supports only ER, SYM, and ARD transition models",
            code="discrete_trait_mk_model_unsupported",
            details={
                "model_name": model_name,
                "resolved_model_name": resolved_model_name,
                "allowed_model_names": list(DISCRETE_TRAIT_MK_MODELS),
            },
        )
    return resolved_model_name


def _validate_root_prior_mode(root_prior_mode: str) -> str:
    normalized_root_prior_mode = root_prior_mode.strip().lower()
    if normalized_root_prior_mode not in DISCRETE_TRAIT_MK_ROOT_PRIOR_MODES:
        raise PhylogeneticsError(
            "bayesian discrete-trait Mk model root_prior_mode must be one of equal, empirical, or fixed",
            code="discrete_trait_mk_root_prior_mode_invalid",
            details={
                "root_prior_mode": root_prior_mode,
                "allowed_root_prior_modes": list(DISCRETE_TRAIT_MK_ROOT_PRIOR_MODES),
            },
        )
    return normalized_root_prior_mode


def _validate_optional_nonblank_state_name(
    value: str | None,
    *,
    field_name: str,
    owner_name: str,
) -> str | None:
    if value is None:
        return None
    return _validate_nonblank_state_name(
        value,
        field_name=field_name,
        owner_name=owner_name,
    )


def _validate_nonblank_state_name(
    value: str,
    *,
    field_name: str,
    owner_name: str,
) -> str:
    if not isinstance(value, str):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one non-blank string",
            code="discrete_trait_mk_state_name_type_invalid",
            details={"field_name": field_name},
        )
    normalized_value = value.strip()
    if not normalized_value:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one non-blank string",
            code="discrete_trait_mk_state_name_blank",
            details={"field_name": field_name},
        )
    return normalized_value


def _validate_positive_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    try:
        validated_value = float(value)
    except (TypeError, ValueError) as error:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one finite float",
            code="discrete_trait_mk_float_type_invalid",
            details={"field_name": field_name},
        ) from error
    if not math.isfinite(validated_value) or validated_value <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be positive and finite",
            code="discrete_trait_mk_positive_float_invalid",
            details={"field_name": field_name, "value": value},
        )
    return float(format(validated_value, ".15g"))


__all__ = [
    "DISCRETE_TRAIT_MK_MODELS",
    "DISCRETE_TRAIT_MK_ROOT_PRIOR_MODES",
    "DiscreteTraitMkModelDefinition",
    "DiscreteTraitMkNodeStateSummary",
    "DiscreteTraitMkPosteriorRow",
    "DiscreteTraitMkProposalSchedule",
    "DiscreteTraitMkRunReport",
    "build_discrete_trait_mk_model_definition",
    "build_discrete_trait_mk_proposal_schedule",
    "run_discrete_trait_mk_metropolis_hastings",
]
