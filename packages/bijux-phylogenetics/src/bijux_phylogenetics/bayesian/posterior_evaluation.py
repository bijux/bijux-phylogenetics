from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import math

from bijux_phylogenetics.bayesian.state import (
    BayesianLikelihoodContextState,
    BayesianModelParameterState,
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    MetropolisHastingsRandomState,
    build_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

_POSTERIOR_DECOMPOSITION_TOLERANCE = 1e-12

BayesianPriorUpdate = Callable[
    [BayesianPhylogeneticState],
    list[BayesianPriorComponentState],
]
BayesianLikelihoodUpdate = Callable[[BayesianPhylogeneticState], float]


@dataclass(slots=True)
class BayesianPosteriorDecompositionRow:
    """One explicit term inside one native Bayesian posterior decomposition."""

    term_kind: str
    term_name: str
    family: str | None
    log_value: float
    parameter_values: dict[str, float]


@dataclass(slots=True)
class BayesianPosteriorEvaluationReport:
    """One explicit native Bayesian posterior evaluation with decomposition."""

    analysis_name: str
    topology_id: str
    posterior_term_source: str
    likelihood_term_source: str
    prior_term_source: str
    identity_tolerance: float
    log_likelihood: float
    total_log_prior: float
    hastings_independent_log_posterior: float
    decomposition_delta: float
    decomposition_valid: bool
    state: BayesianPhylogeneticState
    rows: list[BayesianPosteriorDecompositionRow]


def evaluate_bayesian_posterior(
    *,
    tree: PhyloTree,
    model_parameters: BayesianModelParameterState,
    likelihood_context: BayesianLikelihoodContextState | None = None,
    update_prior_components: BayesianPriorUpdate,
    update_log_likelihood: BayesianLikelihoodUpdate,
    random_state: MetropolisHastingsRandomState | None = None,
    identity_tolerance: float = _POSTERIOR_DECOMPOSITION_TOLERANCE,
) -> BayesianPosteriorEvaluationReport:
    """Evaluate one native Bayesian posterior from owned tree and parameter inputs."""
    provisional_state = build_bayesian_phylogenetic_state(
        tree=tree,
        model_parameters=model_parameters,
        likelihood_context=likelihood_context,
        prior_components=[
            BayesianPriorComponentState(
                component_name="evaluation-baseline",
                family=None,
                log_prior=0.0,
            )
        ],
        log_likelihood=0.0,
    )
    updated_prior_components = update_prior_components(provisional_state)
    if not updated_prior_components:
        raise PhylogeneticsError(
            "bayesian posterior evaluation requires at least one prior component",
            code="bayesian_posterior_evaluation_prior_components_empty",
        )
    scored_state = build_bayesian_phylogenetic_state(
        tree=provisional_state.tree.to_tree(),
        model_parameters=provisional_state.model_parameters,
        likelihood_context=provisional_state.likelihood_context,
        prior_components=updated_prior_components,
        log_likelihood=update_log_likelihood(provisional_state),
        random_state=random_state,
    )
    return summarize_bayesian_posterior_decomposition(
        scored_state,
        identity_tolerance=identity_tolerance,
    )


def summarize_bayesian_posterior_decomposition(
    state: BayesianPhylogeneticState,
    *,
    identity_tolerance: float = _POSTERIOR_DECOMPOSITION_TOLERANCE,
) -> BayesianPosteriorEvaluationReport:
    """Summarize likelihood, prior, and posterior terms for one scored state."""
    if not isinstance(state, BayesianPhylogeneticState):
        raise PhylogeneticsError(
            "bayesian posterior decomposition requires one BayesianPhylogeneticState",
            code="bayesian_posterior_decomposition_state_type_invalid",
        )
    validated_identity_tolerance = _validate_identity_tolerance(identity_tolerance)
    decomposition_delta = state.posterior_log_score - (
        state.log_likelihood + state.total_log_prior
    )
    return BayesianPosteriorEvaluationReport(
        analysis_name=state.likelihood_context.analysis_name,
        topology_id=state.tree.topology_id,
        posterior_term_source="derived_from_likelihood_and_prior",
        likelihood_term_source="evaluated",
        prior_term_source="summed_prior_components",
        identity_tolerance=validated_identity_tolerance,
        log_likelihood=state.log_likelihood,
        total_log_prior=state.total_log_prior,
        hastings_independent_log_posterior=state.posterior_log_score,
        decomposition_delta=decomposition_delta,
        decomposition_valid=abs(decomposition_delta) <= validated_identity_tolerance,
        state=state,
        rows=[
            BayesianPosteriorDecompositionRow(
                term_kind="likelihood",
                term_name="likelihood",
                family=None,
                log_value=state.log_likelihood,
                parameter_values={},
            ),
            *[
                BayesianPosteriorDecompositionRow(
                    term_kind="prior",
                    term_name=component.component_name,
                    family=component.family,
                    log_value=component.log_prior,
                    parameter_values=dict(component.parameter_values),
                )
                for component in state.prior_components
            ],
        ],
    )


def _validate_identity_tolerance(value: float) -> float:
    if not isinstance(value, int | float) or isinstance(value, bool):
        raise PhylogeneticsError(
            "bayesian posterior decomposition identity_tolerance requires one numeric value",
            code="bayesian_posterior_decomposition_identity_tolerance_type_invalid",
        )
    validated_value = float(value)
    if not math.isfinite(validated_value) or validated_value < 0.0:
        raise PhylogeneticsError(
            "bayesian posterior decomposition identity_tolerance requires one nonnegative finite value",
            code="bayesian_posterior_decomposition_identity_tolerance_invalid",
            details={"identity_tolerance": value},
        )
    return validated_value


__all__ = [
    "BayesianPosteriorDecompositionRow",
    "BayesianPosteriorEvaluationReport",
    "evaluate_bayesian_posterior",
    "summarize_bayesian_posterior_decomposition",
]
