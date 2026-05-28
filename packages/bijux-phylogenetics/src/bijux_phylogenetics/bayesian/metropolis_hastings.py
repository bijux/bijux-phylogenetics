from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import math
import random

from bijux_phylogenetics.bayesian.state import (
    BayesianModelParameterState,
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

BayesianPriorUpdate = Callable[
    [BayesianPhylogeneticState],
    list[BayesianPriorComponentState],
]
BayesianLikelihoodUpdate = Callable[[BayesianPhylogeneticState], float]
BayesianStateProposal = Callable[
    [BayesianPhylogeneticState, random.Random],
    "MetropolisHastingsProposal",
]


@dataclass(frozen=True, slots=True)
class MetropolisHastingsProposal:
    """One validated proposal result for a Metropolis-Hastings iteration."""

    changed_fields: tuple[str, ...]
    log_forward_density: float
    log_reverse_density: float
    is_valid: bool
    invalid_reason: str | None
    proposed_tree: object | None = None
    proposed_model_parameters: BayesianModelParameterState | None = None

    @property
    def log_hastings_ratio(self) -> float:
        return self.log_reverse_density - self.log_forward_density


@dataclass(frozen=True, slots=True)
class MetropolisHastingsStepRow:
    """One iteration-level Metropolis-Hastings trace row."""

    iteration_index: int
    proposal_changed_fields: tuple[str, ...]
    proposal_valid: bool
    proposal_invalid_reason: str | None
    log_forward_density: float
    log_reverse_density: float
    accepted: bool
    log_hastings_ratio: float
    current_posterior_log_score: float
    proposed_posterior_log_score: float | None
    log_acceptance_ratio: float | None
    recorded_posterior_log_score: float


@dataclass(frozen=True, slots=True)
class MetropolisHastingsRunReport:
    """One completed Metropolis-Hastings chain report."""

    iteration_count: int
    sample_every: int
    seed: int
    accepted_count: int
    rejected_count: int
    acceptance_rate: float
    initial_state: BayesianPhylogeneticState
    final_state: BayesianPhylogeneticState
    sampled_states: list[BayesianPhylogeneticState]
    step_rows: list[MetropolisHastingsStepRow]


def build_metropolis_hastings_proposal(
    *,
    changed_fields: list[str] | tuple[str, ...],
    log_forward_density: float,
    log_reverse_density: float,
    is_valid: bool,
    invalid_reason: str | None = None,
    proposed_tree: object | None = None,
    proposed_model_parameters: BayesianModelParameterState | None = None,
) -> MetropolisHastingsProposal:
    """Build one validated Metropolis-Hastings proposal result."""
    return _validate_metropolis_hastings_proposal(
        MetropolisHastingsProposal(
            changed_fields=tuple(changed_fields),
            log_forward_density=log_forward_density,
            log_reverse_density=log_reverse_density,
            is_valid=is_valid,
            invalid_reason=invalid_reason,
            proposed_tree=proposed_tree,
            proposed_model_parameters=proposed_model_parameters,
        )
    )


def propose_branch_length_scaling_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    log_scale_standard_deviation: float,
) -> MetropolisHastingsProposal:
    """Propose one multiplicative change to one positive branch length."""
    validated_log_scale_standard_deviation = _validate_positive_finite_float(
        value=log_scale_standard_deviation,
        field_name="log_scale_standard_deviation",
        owner_name="branch-length scaling proposal",
    )
    current_tree = current_state.tree.to_tree()
    candidate_branch_rows = [
        (child.node_id, float(child.branch_length))
        for _parent, child in current_tree.iter_edges()
        if child.node_id is not None
        and child.branch_length is not None
        and math.isfinite(child.branch_length)
        and child.branch_length > 0.0
    ]
    if not candidate_branch_rows:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.branch_lengths",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="tree has no strictly positive finite branch lengths",
        )
    selected_branch_id, current_branch_length = candidate_branch_rows[
        rng.randrange(len(candidate_branch_rows))
    ]
    try:
        scale_factor = math.exp(
            rng.gauss(0.0, validated_log_scale_standard_deviation)
        )
    except OverflowError:
        return build_metropolis_hastings_proposal(
            changed_fields=(f"tree.branch_length:{selected_branch_id}",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="branch-length scaling factor overflowed",
        )
    proposed_branch_length = current_branch_length * scale_factor
    if not math.isfinite(proposed_branch_length) or proposed_branch_length <= 0.0:
        return build_metropolis_hastings_proposal(
            changed_fields=(f"tree.branch_length:{selected_branch_id}",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="branch-length scaling produced a non-positive finite branch length",
        )
    proposed_tree = _copy_tree_with_scaled_branch_length(
        current_tree=current_tree,
        branch_id=selected_branch_id,
        proposed_branch_length=proposed_branch_length,
    )
    branch_selection_log_probability = -math.log(len(candidate_branch_rows))
    log_forward_density = branch_selection_log_probability + _lognormal_scaling_density(
        current_branch_length=current_branch_length,
        proposed_branch_length=proposed_branch_length,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    log_reverse_density = branch_selection_log_probability + _lognormal_scaling_density(
        current_branch_length=proposed_branch_length,
        proposed_branch_length=current_branch_length,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    return build_metropolis_hastings_proposal(
        changed_fields=(f"tree.branch_length:{selected_branch_id}",),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=proposed_tree,
        proposed_model_parameters=current_state.model_parameters,
    )


def propose_global_tree_height_scaling_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    log_scale_standard_deviation: float,
) -> MetropolisHastingsProposal:
    """Propose one coherent multiplicative scale change across every branch."""
    validated_log_scale_standard_deviation = _validate_positive_finite_float(
        value=log_scale_standard_deviation,
        field_name="log_scale_standard_deviation",
        owner_name="global tree-height scaling proposal",
    )
    current_tree = current_state.tree.to_tree()
    current_branch_lengths = [
        float(child.branch_length)
        for _parent, child in current_tree.iter_edges()
        if child.branch_length is not None
    ]
    if not current_branch_lengths:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.branch_lengths",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="global tree-height scaling requires explicit branch lengths on every edge",
        )
    if len(current_branch_lengths) != sum(1 for _ in current_tree.iter_edges()):
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.branch_lengths",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="global tree-height scaling requires explicit branch lengths on every edge",
        )
    if any(
        not math.isfinite(branch_length) or branch_length <= 0.0
        for branch_length in current_branch_lengths
    ):
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.branch_lengths",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="global tree-height scaling requires strictly positive finite branch lengths on every edge",
        )
    try:
        scale_factor = math.exp(
            rng.gauss(0.0, validated_log_scale_standard_deviation)
        )
    except OverflowError:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.branch_lengths",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="global tree-height scaling factor overflowed",
        )
    if not math.isfinite(scale_factor) or scale_factor <= 0.0:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.branch_lengths",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="global tree-height scaling factor was not strictly positive and finite",
        )
    proposed_tree = _copy_tree_with_globally_scaled_branch_lengths(
        current_tree=current_tree,
        scale_factor=scale_factor,
    )
    proposed_branch_lengths = [
        float(child.branch_length)
        for _parent, child in proposed_tree.iter_edges()
        if child.branch_length is not None
    ]
    log_forward_density = _global_tree_scaling_density(
        current_branch_lengths=current_branch_lengths,
        proposed_branch_lengths=proposed_branch_lengths,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    log_reverse_density = _global_tree_scaling_density(
        current_branch_lengths=proposed_branch_lengths,
        proposed_branch_lengths=current_branch_lengths,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    return build_metropolis_hastings_proposal(
        changed_fields=("tree.branch_lengths",),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=proposed_tree,
        proposed_model_parameters=current_state.model_parameters,
    )


def score_bayesian_phylogenetic_state(
    *,
    tree,
    model_parameters: BayesianModelParameterState,
    update_prior_components: BayesianPriorUpdate,
    update_log_likelihood: BayesianLikelihoodUpdate,
) -> BayesianPhylogeneticState:
    """Score one Bayesian state from owned tree and model-parameter surfaces."""
    provisional_state = build_bayesian_phylogenetic_state(
        tree=tree,
        model_parameters=model_parameters,
        prior_components=[
            BayesianPriorComponentState(
                component_name="provisional",
                family=None,
                log_prior=0.0,
            )
        ],
        log_likelihood=0.0,
    )
    updated_prior_components = update_prior_components(provisional_state)
    if not updated_prior_components:
        raise PhylogeneticsError(
            "metropolis-hastings state scoring requires at least one prior component",
            code="metropolis_hastings_prior_components_empty",
        )
    return build_bayesian_phylogenetic_state(
        tree=provisional_state.tree.to_tree(),
        model_parameters=provisional_state.model_parameters,
        prior_components=updated_prior_components,
        log_likelihood=update_log_likelihood(provisional_state),
    )


def run_metropolis_hastings_sampler(
    *,
    initial_state: BayesianPhylogeneticState,
    propose_state: BayesianStateProposal,
    update_prior_components: BayesianPriorUpdate,
    update_log_likelihood: BayesianLikelihoodUpdate,
    iteration_count: int,
    sample_every: int = 1,
    seed: int = 0,
) -> MetropolisHastingsRunReport:
    """Run one native Metropolis-Hastings chain over Bayesian phylogenetic states."""
    validated_iteration_count = _validate_positive_integer(
        value=iteration_count,
        field_name="iteration_count",
        owner_name="metropolis-hastings sampler",
    )
    validated_sample_every = _validate_positive_integer(
        value=sample_every,
        field_name="sample_every",
        owner_name="metropolis-hastings sampler",
    )
    validated_seed = _validate_integer_seed(seed)
    rng = random.Random(validated_seed)  # nosec B311
    current_state = initial_state
    step_rows: list[MetropolisHastingsStepRow] = []
    sampled_states: list[BayesianPhylogeneticState] = [current_state]
    accepted_count = 0

    for iteration_index in range(1, validated_iteration_count + 1):
        previous_state = current_state
        proposal = _validate_metropolis_hastings_proposal(
            propose_state(current_state, rng)
        )
        proposed_state: BayesianPhylogeneticState | None = None
        log_acceptance_ratio: float | None = None
        accepted = False
        if proposal.is_valid:
            proposed_state = score_bayesian_phylogenetic_state(
                tree=proposal.proposed_tree,
                model_parameters=proposal.proposed_model_parameters,
                update_prior_components=update_prior_components,
                update_log_likelihood=update_log_likelihood,
            )
            log_acceptance_ratio = (
                proposed_state.posterior_log_score
                - current_state.posterior_log_score
                + proposal.log_hastings_ratio
            )
            accepted = _accept_metropolis_hastings_proposal(
                log_acceptance_ratio=log_acceptance_ratio,
                rng=rng,
            )
            if accepted:
                current_state = proposed_state
                accepted_count += 1
        if iteration_index % validated_sample_every == 0:
            sampled_states.append(current_state)
        step_rows.append(
            MetropolisHastingsStepRow(
                iteration_index=iteration_index,
                proposal_changed_fields=proposal.changed_fields,
                proposal_valid=proposal.is_valid,
                proposal_invalid_reason=proposal.invalid_reason,
                log_forward_density=proposal.log_forward_density,
                log_reverse_density=proposal.log_reverse_density,
                accepted=accepted,
                log_hastings_ratio=proposal.log_hastings_ratio,
                current_posterior_log_score=previous_state.posterior_log_score,
                proposed_posterior_log_score=(
                    proposed_state.posterior_log_score
                    if proposed_state is not None
                    else None
                ),
                log_acceptance_ratio=log_acceptance_ratio,
                recorded_posterior_log_score=current_state.posterior_log_score,
            )
        )

    rejected_count = validated_iteration_count - accepted_count
    return MetropolisHastingsRunReport(
        iteration_count=validated_iteration_count,
        sample_every=validated_sample_every,
        seed=validated_seed,
        accepted_count=accepted_count,
        rejected_count=rejected_count,
        acceptance_rate=accepted_count / validated_iteration_count,
        initial_state=initial_state,
        final_state=current_state,
        sampled_states=sampled_states,
        step_rows=step_rows,
    )


def _accept_metropolis_hastings_proposal(
    *,
    log_acceptance_ratio: float,
    rng: random.Random,
) -> bool:
    if log_acceptance_ratio >= 0.0:
        return True
    return math.log(rng.random()) <= log_acceptance_ratio


def _validate_positive_integer(
    *,
    value: int,
    field_name: str,
    owner_name: str,
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one integer",
            code="metropolis_hastings_integer_required",
            details={"field_name": field_name},
        )
    if value <= 0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be positive",
            code="metropolis_hastings_positive_integer_required",
            details={"field_name": field_name},
        )
    return value


def _validate_metropolis_hastings_proposal(
    proposal: MetropolisHastingsProposal,
) -> MetropolisHastingsProposal:
    if not isinstance(proposal, MetropolisHastingsProposal):
        raise PhylogeneticsError(
            "metropolis-hastings sampler requires proposals to return one MetropolisHastingsProposal",
            code="metropolis_hastings_proposal_type_invalid",
        )
    validated_changed_fields = _validate_changed_fields(proposal.changed_fields)
    validated_log_forward_density = _validate_finite_float(
        value=proposal.log_forward_density,
        field_name="log_forward_density",
        owner_name="metropolis-hastings proposal",
    )
    validated_log_reverse_density = _validate_finite_float(
        value=proposal.log_reverse_density,
        field_name="log_reverse_density",
        owner_name="metropolis-hastings proposal",
    )
    if not isinstance(proposal.is_valid, bool):
        raise PhylogeneticsError(
            "metropolis-hastings proposal requires 'is_valid' to be boolean",
            code="metropolis_hastings_proposal_validity_type_invalid",
        )
    validated_invalid_reason = _validate_invalid_reason(
        invalid_reason=proposal.invalid_reason,
        is_valid=proposal.is_valid,
    )
    if proposal.is_valid and not isinstance(
        proposal.proposed_model_parameters,
        BayesianModelParameterState,
    ):
        raise PhylogeneticsError(
            "metropolis-hastings proposal requires one BayesianModelParameterState",
            code="metropolis_hastings_model_parameter_state_type_invalid",
        )
    return MetropolisHastingsProposal(
        changed_fields=validated_changed_fields,
        log_forward_density=validated_log_forward_density,
        log_reverse_density=validated_log_reverse_density,
        is_valid=proposal.is_valid,
        invalid_reason=validated_invalid_reason,
        proposed_tree=proposal.proposed_tree,
        proposed_model_parameters=proposal.proposed_model_parameters,
    )


def _validate_integer_seed(seed: int) -> int:
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise PhylogeneticsError(
            "metropolis-hastings sampler requires 'seed' to be one integer",
            code="metropolis_hastings_seed_type_invalid",
        )
    return seed


def _validate_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be numeric",
            code="metropolis_hastings_float_required",
            details={"field_name": field_name},
        )
    normalized_value = float(value)
    if not math.isfinite(normalized_value):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be finite",
            code="metropolis_hastings_finite_float_required",
            details={"field_name": field_name},
        )
    return normalized_value


def _validate_positive_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    normalized_value = _validate_finite_float(
        value=value,
        field_name=field_name,
        owner_name=owner_name,
    )
    if normalized_value <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be strictly positive",
            code="metropolis_hastings_positive_float_required",
            details={"field_name": field_name},
        )
    return normalized_value


def _validate_changed_fields(
    changed_fields: tuple[str, ...],
) -> tuple[str, ...]:
    if not isinstance(changed_fields, tuple):
        raise PhylogeneticsError(
            "metropolis-hastings proposal requires 'changed_fields' to be one tuple",
            code="metropolis_hastings_changed_fields_type_invalid",
        )
    if not changed_fields:
        raise PhylogeneticsError(
            "metropolis-hastings proposal requires at least one changed field",
            code="metropolis_hastings_changed_fields_empty",
        )
    normalized_changed_fields: list[str] = []
    seen_changed_fields: set[str] = set()
    for changed_field in changed_fields:
        if not isinstance(changed_field, str) or not changed_field.strip():
            raise PhylogeneticsError(
                "metropolis-hastings proposal requires every changed field to be nonblank text",
                code="metropolis_hastings_changed_field_invalid",
            )
        normalized_changed_field = changed_field.strip()
        if normalized_changed_field in seen_changed_fields:
            raise PhylogeneticsError(
                "metropolis-hastings proposal requires changed fields to be unique",
                code="metropolis_hastings_changed_fields_duplicate",
                details={"changed_field": normalized_changed_field},
            )
        seen_changed_fields.add(normalized_changed_field)
        normalized_changed_fields.append(normalized_changed_field)
    return tuple(normalized_changed_fields)


def _validate_invalid_reason(
    *,
    invalid_reason: str | None,
    is_valid: bool,
) -> str | None:
    if is_valid:
        if invalid_reason is not None:
            raise PhylogeneticsError(
                "metropolis-hastings proposal does not allow 'invalid_reason' on valid proposals",
                code="metropolis_hastings_invalid_reason_not_allowed",
            )
        return None
    if not isinstance(invalid_reason, str) or not invalid_reason.strip():
        raise PhylogeneticsError(
            "metropolis-hastings proposal requires one nonblank invalid reason when 'is_valid' is false",
            code="metropolis_hastings_invalid_reason_required",
        )
    return invalid_reason.strip()


def _copy_tree_with_scaled_branch_length(
    *,
    current_tree: PhyloTree,
    branch_id: str,
    proposed_branch_length: float,
) -> PhyloTree:
    proposed_tree = current_tree.copy()
    proposed_tree.node_by_id(branch_id).branch_length = proposed_branch_length
    return proposed_tree


def _copy_tree_with_globally_scaled_branch_lengths(
    *,
    current_tree: PhyloTree,
    scale_factor: float,
) -> PhyloTree:
    proposed_tree = current_tree.copy()
    for _parent, child in proposed_tree.iter_edges():
        if child.branch_length is None:
            raise PhylogeneticsError(
                "global tree-height scaling requires explicit branch lengths on every edge",
                code="metropolis_hastings_global_tree_scale_branch_length_missing",
            )
        child.branch_length = child.branch_length * scale_factor
    return proposed_tree


def _lognormal_scaling_density(
    *,
    current_branch_length: float,
    proposed_branch_length: float,
    log_scale_standard_deviation: float,
) -> float:
    log_scale_change = math.log(proposed_branch_length / current_branch_length)
    z_score = log_scale_change / log_scale_standard_deviation
    return (
        -math.log(proposed_branch_length)
        - math.log(log_scale_standard_deviation)
        - (math.log(2.0 * math.pi) / 2.0)
        - ((z_score * z_score) / 2.0)
    )


def _global_tree_scaling_density(
    *,
    current_branch_lengths: list[float],
    proposed_branch_lengths: list[float],
    log_scale_standard_deviation: float,
) -> float:
    if len(current_branch_lengths) != len(proposed_branch_lengths):
        raise PhylogeneticsError(
            "global tree-height scaling requires matching branch dimensions",
            code="metropolis_hastings_global_tree_scale_dimension_mismatch",
        )
    log_scale_change = math.log(proposed_branch_lengths[0] / current_branch_lengths[0])
    z_score = log_scale_change / log_scale_standard_deviation
    return (
        -math.fsum(math.log(branch_length) for branch_length in proposed_branch_lengths)
        - math.log(log_scale_standard_deviation)
        - (math.log(2.0 * math.pi) / 2.0)
        - ((z_score * z_score) / 2.0)
    )
