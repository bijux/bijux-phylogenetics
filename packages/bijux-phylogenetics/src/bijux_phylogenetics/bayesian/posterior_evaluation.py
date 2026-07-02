from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.bayesian.state import BayesianPhylogeneticState
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

_POSTERIOR_DECOMPOSITION_TOLERANCE = 1e-12


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
    "summarize_bayesian_posterior_decomposition",
]
