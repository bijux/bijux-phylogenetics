from __future__ import annotations

from bijux_phylogenetics.ancestral.continuous import ContinuousAncestralReport
from bijux_phylogenetics.ancestral.discrete import DiscreteAncestralReport
from bijux_phylogenetics.comparative.continuous import (
    BrownianTraitEvolutionSummaryReport,
    OUTraitEvolutionSummaryReport,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    DiscreteStateEvolutionReport,
)
from bijux_phylogenetics.compare.topology import TreeComparisonReport
from bijux_phylogenetics.ecology import HostSwitchingReport

from .models import (
    KnownAnswerContinuousNodeRecoveryRow,
    KnownAnswerContinuousNodeTruth,
    KnownAnswerDiscreteNodeRecoveryRow,
    KnownAnswerDiscreteNodeTruth,
    KnownAnswerParameterRecoveryRow,
    KnownAnswerRecoveryThreshold,
    KnownAnswerThresholdEvaluationRow,
    KnownAnswerTransitionRecoveryRow,
    KnownAnswerTransitionTruth,
)
from .policy import (
    evaluate_threshold,
    format_observed_value,
    format_transition,
    mean,
    transition_recovery_match,
)


def build_parameter_recovery_rows(
    *,
    true_parameters: dict[str, str],
    brownian_fit: BrownianTraitEvolutionSummaryReport,
    ou_fit: OUTraitEvolutionSummaryReport,
) -> list[KnownAnswerParameterRecoveryRow]:
    true_root_state = float(true_parameters["continuous_root_state"])
    true_sigma_squared = float(true_parameters["continuous_sigma_squared"])
    true_ou_alpha = float(true_parameters["ou_alpha"])
    true_ou_theta = float(true_parameters["ou_theta"])
    true_ou_sigma_squared = float(true_parameters["ou_sigma_squared"])
    return [
        _parameter_recovery_row(
            parameter="continuous_root_state",
            true_value=true_root_state,
            estimated_value=brownian_fit.root_state,
            interpretation="Brownian root-state estimate recovered from observed tip values on the true tree.",
        ),
        _parameter_recovery_row(
            parameter="continuous_sigma_squared",
            true_value=true_sigma_squared,
            estimated_value=brownian_fit.sigma_squared,
            interpretation="Brownian evolutionary rate recovered from observed tip values on the true tree.",
        ),
        _parameter_recovery_row(
            parameter="ou_alpha",
            true_value=true_ou_alpha,
            estimated_value=ou_fit.alpha,
            interpretation="OU alpha recovered from the governed OU tip trait on the true tree.",
        ),
        _parameter_recovery_row(
            parameter="ou_theta",
            true_value=true_ou_theta,
            estimated_value=ou_fit.theta,
            interpretation="OU optimum recovered from the governed OU tip trait on the true tree.",
        ),
        _parameter_recovery_row(
            parameter="ou_sigma_squared",
            true_value=true_ou_sigma_squared,
            estimated_value=ou_fit.sigma_squared,
            interpretation="OU sigma-squared recovered from the governed OU tip trait on the true tree.",
        ),
    ]


def build_continuous_node_recovery_rows(
    *,
    true_nodes: list[KnownAnswerContinuousNodeTruth],
    report: ContinuousAncestralReport,
) -> list[KnownAnswerContinuousNodeRecoveryRow]:
    truth_by_node = {row.node: row for row in true_nodes if not row.is_tip}
    return [
        KnownAnswerContinuousNodeRecoveryRow(
            node=estimate.node,
            descendant_taxa=estimate.descendant_taxa,
            true_value=truth_by_node[estimate.node].true_value,
            estimated_value=estimate.estimate,
            absolute_error=abs(
                estimate.estimate - truth_by_node[estimate.node].true_value
            ),
            standard_error=estimate.standard_error,
            lower_95_interval=estimate.lower_95_interval,
            upper_95_interval=estimate.upper_95_interval,
            confidence=estimate.confidence,
        )
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def build_discrete_node_recovery_rows(
    *,
    true_nodes: list[KnownAnswerDiscreteNodeTruth],
    report: DiscreteAncestralReport,
) -> list[KnownAnswerDiscreteNodeRecoveryRow]:
    truth_by_node = {row.node: row for row in true_nodes if not row.is_tip}
    return [
        KnownAnswerDiscreteNodeRecoveryRow(
            node=estimate.node,
            descendant_taxa=estimate.descendant_taxa,
            true_state=truth_by_node[estimate.node].true_state,
            estimated_state=estimate.most_likely_state,
            true_state_probability=estimate.state_probabilities.get(
                truth_by_node[estimate.node].true_state,
                0.0,
            ),
            confidence=estimate.confidence,
            correct=estimate.most_likely_state
            == truth_by_node[estimate.node].true_state,
            ambiguous=estimate.ambiguous,
        )
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def build_geographic_node_recovery_rows(
    *,
    true_nodes: list[KnownAnswerDiscreteNodeTruth],
    report: DiscreteStateEvolutionReport,
) -> list[KnownAnswerDiscreteNodeRecoveryRow]:
    truth_by_node = {row.node: row for row in true_nodes if not row.is_tip}
    return [
        KnownAnswerDiscreteNodeRecoveryRow(
            node=estimate.node,
            descendant_taxa=estimate.descendant_taxa,
            true_state=truth_by_node[estimate.node].true_state,
            estimated_state=estimate.most_likely_state,
            true_state_probability=estimate.state_probabilities.get(
                truth_by_node[estimate.node].true_state,
                0.0,
            ),
            confidence=max(estimate.state_probabilities.values()),
            correct=estimate.most_likely_state
            == truth_by_node[estimate.node].true_state,
            ambiguous=estimate.ambiguous,
        )
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def build_host_node_recovery_rows(
    *,
    true_nodes: list[KnownAnswerDiscreteNodeTruth],
    report: HostSwitchingReport,
) -> list[KnownAnswerDiscreteNodeRecoveryRow]:
    truth_by_node = {row.node: row for row in true_nodes if not row.is_tip}
    return [
        KnownAnswerDiscreteNodeRecoveryRow(
            node=estimate.node,
            descendant_taxa=estimate.descendant_taxa,
            true_state=truth_by_node[estimate.node].true_state,
            estimated_state=estimate.most_likely_host,
            true_state_probability=estimate.host_probabilities.get(
                truth_by_node[estimate.node].true_state,
                0.0,
            ),
            confidence=estimate.confidence,
            correct=estimate.most_likely_host
            == truth_by_node[estimate.node].true_state,
            ambiguous=estimate.ambiguous,
        )
        for estimate in report.node_rows
        if not estimate.node_name or estimate.node != estimate.node_name
    ]


def build_host_event_recovery_rows(
    *,
    true_rows: list[KnownAnswerTransitionTruth],
    report: HostSwitchingReport,
) -> list[KnownAnswerTransitionRecoveryRow]:
    inferred_by_branch = {
        (row.parent_node, row.child_node): row for row in report.branch_rows
    }
    return [
        KnownAnswerTransitionRecoveryRow(
            parent_node=row.parent_node,
            child_node=row.child_node,
            true_transition=format_transition(row.source_state, row.target_state),
            estimated_transition=format_transition(
                inferred_by_branch[
                    (row.parent_node, row.child_node)
                ].parent_most_likely_host,
                inferred_by_branch[
                    (row.parent_node, row.child_node)
                ].child_most_likely_host,
            ),
            true_changed=row.changed,
            estimated_changed=inferred_by_branch[
                (row.parent_node, row.child_node)
            ].changed,
            true_event_count=row.event_count,
            estimated_event_count=1
            if inferred_by_branch[(row.parent_node, row.child_node)].changed
            else 0,
            correct=transition_recovery_match(
                true_changed=row.changed,
                estimated_changed=inferred_by_branch[
                    (row.parent_node, row.child_node)
                ].changed,
                true_source=row.source_state,
                true_target=row.target_state,
                estimated_source=inferred_by_branch[
                    (row.parent_node, row.child_node)
                ].parent_most_likely_host,
                estimated_target=inferred_by_branch[
                    (row.parent_node, row.child_node)
                ].child_most_likely_host,
                true_event_count=row.event_count,
                estimated_event_count=1
                if inferred_by_branch[(row.parent_node, row.child_node)].changed
                else 0,
            ),
        )
        for row in true_rows
    ]


def build_transition_recovery_rows(
    *,
    true_rows: list[KnownAnswerTransitionTruth],
    report: DiscreteStateEvolutionReport,
) -> list[KnownAnswerTransitionRecoveryRow]:
    inferred_by_branch = {
        (row.parent_node, row.child_node): row
        for row in report.transition_summary.events
    }
    return [
        KnownAnswerTransitionRecoveryRow(
            parent_node=row.parent_node,
            child_node=row.child_node,
            true_transition=format_transition(row.source_state, row.target_state),
            estimated_transition=format_transition(
                inferred_by_branch[(row.parent_node, row.child_node)].source_state,
                inferred_by_branch[(row.parent_node, row.child_node)].target_state,
            ),
            true_changed=row.changed,
            estimated_changed=inferred_by_branch[
                (row.parent_node, row.child_node)
            ].changed,
            true_event_count=row.event_count,
            estimated_event_count=1
            if inferred_by_branch[(row.parent_node, row.child_node)].changed
            else 0,
            correct=transition_recovery_match(
                true_changed=row.changed,
                estimated_changed=inferred_by_branch[
                    (row.parent_node, row.child_node)
                ].changed,
                true_source=row.source_state,
                true_target=row.target_state,
                estimated_source=inferred_by_branch[
                    (row.parent_node, row.child_node)
                ].source_state,
                estimated_target=inferred_by_branch[
                    (row.parent_node, row.child_node)
                ].target_state,
                true_event_count=row.event_count,
                estimated_event_count=1
                if inferred_by_branch[(row.parent_node, row.child_node)].changed
                else 0,
            ),
        )
        for row in true_rows
    ]


def build_threshold_evaluation_rows(
    *,
    thresholds: list[KnownAnswerRecoveryThreshold],
    tree_recovery: TreeComparisonReport,
    parameter_recovery_rows: list[KnownAnswerParameterRecoveryRow],
    discrete_node_recovery_rows: list[KnownAnswerDiscreteNodeRecoveryRow],
    host_node_recovery_rows: list[KnownAnswerDiscreteNodeRecoveryRow],
    host_event_recovery_rows: list[KnownAnswerTransitionRecoveryRow],
    geographic_node_recovery_rows: list[KnownAnswerDiscreteNodeRecoveryRow],
    geographic_event_recovery_rows: list[KnownAnswerTransitionRecoveryRow],
) -> list[KnownAnswerThresholdEvaluationRow]:
    parameter_by_name = {row.parameter: row for row in parameter_recovery_rows}
    metric_values: dict[str, bool | float] = {
        "same_unrooted_topology": tree_recovery.same_unrooted_topology,
        "brownian_root_absolute_error": parameter_by_name[
            "continuous_root_state"
        ].absolute_error,
        "brownian_sigma_squared_absolute_error": parameter_by_name[
            "continuous_sigma_squared"
        ].absolute_error,
        "ou_alpha_absolute_error": parameter_by_name["ou_alpha"].absolute_error,
        "ou_theta_absolute_error": parameter_by_name["ou_theta"].absolute_error,
        "ou_sigma_squared_absolute_error": parameter_by_name[
            "ou_sigma_squared"
        ].absolute_error,
        "discrete_internal_node_accuracy": mean(
            1.0 if row.correct else 0.0 for row in discrete_node_recovery_rows
        ),
        "host_internal_node_accuracy": mean(
            1.0 if row.correct else 0.0 for row in host_node_recovery_rows
        ),
        "host_event_accuracy": mean(
            1.0 if row.correct else 0.0 for row in host_event_recovery_rows
        ),
        "geographic_internal_node_accuracy": mean(
            1.0 if row.correct else 0.0 for row in geographic_node_recovery_rows
        ),
        "geographic_event_accuracy": mean(
            1.0 if row.correct else 0.0 for row in geographic_event_recovery_rows
        ),
    }
    rows: list[KnownAnswerThresholdEvaluationRow] = []
    for threshold in thresholds:
        observed_value = metric_values[threshold.metric]
        rows.append(
            KnownAnswerThresholdEvaluationRow(
                metric=threshold.metric,
                comparator=threshold.comparator,
                threshold=threshold.threshold,
                observed_value=format_observed_value(observed_value),
                passed=evaluate_threshold(
                    comparator=threshold.comparator,
                    threshold=threshold.threshold,
                    observed_value=observed_value,
                ),
                rationale=threshold.rationale,
            )
        )
    return rows


def _parameter_recovery_row(
    *,
    parameter: str,
    true_value: float,
    estimated_value: float,
    interpretation: str,
) -> KnownAnswerParameterRecoveryRow:
    absolute_error = abs(estimated_value - true_value)
    denominator = abs(true_value)
    relative_error = 0.0 if denominator == 0.0 else absolute_error / denominator
    return KnownAnswerParameterRecoveryRow(
        parameter=parameter,
        true_value=true_value,
        estimated_value=estimated_value,
        absolute_error=absolute_error,
        relative_error=relative_error,
        interpretation=interpretation,
    )
