from __future__ import annotations

from pathlib import Path

from .analysis import (
    assess_geographic_state_analysis_readiness as assess_geographic_state_analysis_readiness,
    estimate_ancestral_geographic_states as estimate_ancestral_geographic_states,
    run_discrete_state_transition_model as run_discrete_state_transition_model,
)
from .models import (
    BiogeographicComputedResult as BiogeographicComputedResult,
    BiogeographicInterpretationReport as BiogeographicInterpretationReport,
    DiscreteEvolutionNarrative as DiscreteEvolutionNarrative,
    DiscreteEvolutionReportBuildResult as DiscreteEvolutionReportBuildResult,
    DiscreteModelComparisonReport as DiscreteModelComparisonReport,
    DiscreteModelComparisonRow as DiscreteModelComparisonRow,
    DiscreteStateEvolutionReport as DiscreteStateEvolutionReport,
    DiscreteTransitionReferenceObservation as DiscreteTransitionReferenceObservation,
    DiscreteTransitionReferenceRate as DiscreteTransitionReferenceRate,
    DiscreteTransitionReferenceValidationReport as DiscreteTransitionReferenceValidationReport,
    DominantStateBiasReport as DominantStateBiasReport,
    ModelSensitiveRegionRow as ModelSensitiveRegionRow,
    NodeStateDifference as NodeStateDifference,
    NodeStateEstimate as NodeStateEstimate,
    SparseStateInstabilityReport as SparseStateInstabilityReport,
    StateCodingAuditReport as StateCodingAuditReport,
    StateCodingAuditRow as StateCodingAuditRow,
    StateCodingIssue as StateCodingIssue,
    StateCodingValidationReport as StateCodingValidationReport,
    StateImbalanceReport as StateImbalanceReport,
    StateImbalanceWarning as StateImbalanceWarning,
    TransitionEvent as TransitionEvent,
    TransitionRateRow as TransitionRateRow,
    TransitionRateUncertaintyReport as TransitionRateUncertaintyReport,
    TransitionRateUncertaintyRow as TransitionRateUncertaintyRow,
    TransitionSupportRow as TransitionSupportRow,
)
from .state_coding import (
    audit_discrete_state_coding as audit_discrete_state_coding,
    detect_state_imbalance_problems as detect_state_imbalance_problems,
    validate_discrete_state_coding as validate_discrete_state_coding,
)
from .transition_engine import (
    _estimate_node_states as _estimate_node_states,
    _estimate_transition_support_rows as _estimate_transition_support_rows,
    _fit_transition_matrix as _fit_transition_matrix,
    _fitch_candidate_sets as _fitch_candidate_sets,
    _resolve_state_order as _resolve_state_order,
    _root_prior as _root_prior,
    _stationary_frequencies as _stationary_frequencies,
    _transition_events as _transition_events,
)

_DEFAULT_STATE_COLORS = (
    "#0f766e",
    "#1d4ed8",
    "#c2410c",
    "#7c3aed",
    "#b91c1c",
    "#047857",
    "#a16207",
    "#0f172a",
)


def _quantile(sorted_values: list[float], fraction: float) -> float:
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return float(format(sorted_values[0], ".15g"))
    index = max(
        0, min(len(sorted_values) - 1, int(round(fraction * (len(sorted_values) - 1))))
    )
    return float(format(sorted_values[index], ".15g"))


def _build_model_sensitive_regions(
    differences: list[NodeStateDifference],
) -> list[ModelSensitiveRegionRow]:
    rows: list[ModelSensitiveRegionRow] = []
    for difference in differences:
        if not difference.differs:
            continue
        left_probability = difference.left_probabilities.get(difference.left_state, 0.0)
        right_probability = difference.right_probabilities.get(
            difference.right_state, 0.0
        )
        rows.append(
            ModelSensitiveRegionRow(
                node=difference.node,
                descendant_taxa=difference.descendant_taxa,
                left_state=difference.left_state,
                right_state=difference.right_state,
                sensitivity_score=float(
                    format(abs(left_probability - right_probability), ".15g")
                ),
            )
        )
    return sorted(rows, key=lambda row: (-row.sensitivity_score, row.node))


def compare_discrete_state_models(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    left_model: str = "equal-rates",
    right_model: str = "all-rates-different",
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> DiscreteModelComparisonReport:
    """Compare discrete-state reconstructions across two supported models."""
    left = run_discrete_state_transition_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=left_model,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    right = run_discrete_state_transition_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=right_model,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    right_by_node = {estimate.node: estimate for estimate in right.estimates}
    differences: list[NodeStateDifference] = []
    for left_estimate in left.estimates:
        right_estimate = right_by_node[left_estimate.node]
        differences.append(
            NodeStateDifference(
                node=left_estimate.node,
                descendant_taxa=left_estimate.descendant_taxa,
                left_state=left_estimate.most_likely_state,
                right_state=right_estimate.most_likely_state,
                differs=left_estimate.most_likely_state
                != right_estimate.most_likely_state,
                left_probabilities=left_estimate.state_probabilities,
                right_probabilities=right_estimate.state_probabilities,
            )
        )
    rows = [
        DiscreteModelComparisonRow(
            model=left_model,
            parameter_count=left.transition_model.parameter_count,
            pseudo_log_likelihood=left.transition_model.pseudo_log_likelihood,
            aic=left.transition_model.aic,
            transition_count=left.transition_summary.transition_count,
        ),
        DiscreteModelComparisonRow(
            model=right_model,
            parameter_count=right.transition_model.parameter_count,
            pseudo_log_likelihood=right.transition_model.pseudo_log_likelihood,
            aic=right.transition_model.aic,
            transition_count=right.transition_summary.transition_count,
        ),
    ]
    better_model = min(rows, key=lambda row: row.aic).model
    sensitive_regions = _build_model_sensitive_regions(differences)
    return DiscreteModelComparisonReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=left.taxon_column,
        trait=trait,
        left_model=left_model,
        right_model=right_model,
        better_model=better_model,
        rows=rows,
        node_differences=differences,
        sensitive_region_count=len(sensitive_regions),
        sensitive_regions=sensitive_regions,
    )
