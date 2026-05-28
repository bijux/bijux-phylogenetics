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
    """One proposed Bayesian state plus its Hastings correction."""

    proposed_tree: object
    proposed_model_parameters: BayesianModelParameterState
    log_hastings_ratio: float


@dataclass(frozen=True, slots=True)
class MetropolisHastingsStepRow:
    """One iteration-level Metropolis-Hastings trace row."""

    iteration_index: int
    accepted: bool
    log_hastings_ratio: float
    current_posterior_log_score: float
    proposed_posterior_log_score: float
    log_acceptance_ratio: float
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
                accepted=accepted,
                log_hastings_ratio=proposal.log_hastings_ratio,
                current_posterior_log_score=previous_state.posterior_log_score,
                proposed_posterior_log_score=proposed_state.posterior_log_score,
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
    if not isinstance(proposal.proposed_model_parameters, BayesianModelParameterState):
        raise PhylogeneticsError(
            "metropolis-hastings proposal requires one BayesianModelParameterState",
            code="metropolis_hastings_model_parameter_state_type_invalid",
        )
    return MetropolisHastingsProposal(
        proposed_tree=proposal.proposed_tree,
        proposed_model_parameters=proposal.proposed_model_parameters,
        log_hastings_ratio=_validate_finite_float(
            value=proposal.log_hastings_ratio,
            field_name="log_hastings_ratio",
            owner_name="metropolis-hastings proposal",
        ),
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
