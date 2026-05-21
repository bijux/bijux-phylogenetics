from __future__ import annotations

from pathlib import Path
import shutil

from bijux_phylogenetics.ancestral.continuous import (
    write_continuous_ancestral_summary_table,
    write_continuous_ancestral_uncertainty_table,
)
from bijux_phylogenetics.ancestral.discrete import (
    write_discrete_ancestral_probability_table,
    write_discrete_ancestral_summary_table,
)
from bijux_phylogenetics.comparative.continuous import (
    write_brownian_trait_evolution_summary_table,
    write_ou_trait_evolution_summary_table,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    DiscreteStateEvolutionReport,
    write_node_state_probability_table,
    write_transition_summary_table,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.ecology import (
    write_host_state_node_table,
    write_host_switch_branch_table,
    write_host_switch_summary_table,
)

from .models import (
    KnownAnswerContinuousNodeRecoveryRow,
    KnownAnswerDiscreteNodeRecoveryRow,
    KnownAnswerParameterRecoveryRow,
    KnownAnswerReferenceWorkflowBundle,
    KnownAnswerReferenceWorkflowReport,
    KnownAnswerThresholdEvaluationRow,
    KnownAnswerTransitionRecoveryRow,
)
from .policy import format_number, mean


def write_known_answer_reference_workflow_bundle(
    output_root: Path,
    report: KnownAnswerReferenceWorkflowReport,
) -> KnownAnswerReferenceWorkflowBundle:
    """Write the governed recovery outputs for the packaged simulation panel."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    continuous_mae = mean(
        row.absolute_error for row in report.continuous_node_recovery_rows
    )
    discrete_accuracy = mean(
        1.0 if row.correct else 0.0 for row in report.discrete_node_recovery_rows
    )
    discrete_mean_true_probability = mean(
        row.true_state_probability for row in report.discrete_node_recovery_rows
    )
    host_node_accuracy = mean(
        1.0 if row.correct else 0.0 for row in report.host_node_recovery_rows
    )
    host_event_accuracy = mean(
        1.0 if row.correct else 0.0 for row in report.host_event_recovery_rows
    )
    geographic_node_accuracy = mean(
        1.0 if row.correct else 0.0 for row in report.geographic_node_recovery_rows
    )
    geographic_event_accuracy = mean(
        1.0 if row.correct else 0.0 for row in report.geographic_event_recovery_rows
    )
    workflow_summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report=report,
        continuous_mae=continuous_mae,
        discrete_accuracy=discrete_accuracy,
        discrete_mean_true_probability=discrete_mean_true_probability,
        host_node_accuracy=host_node_accuracy,
        host_event_accuracy=host_event_accuracy,
        geographic_node_accuracy=geographic_node_accuracy,
        geographic_event_accuracy=geographic_event_accuracy,
    )
    distance_tree_path = _write_distance_tree(
        output_root / "recovered-distance-tree.nwk",
        report.distance_tree_newick,
    )
    tree_recovery_path = _write_tree_recovery_table(
        output_root / "tree-recovery.tsv",
        report,
    )
    parameter_recovery_path = _write_parameter_recovery_table(
        output_root / "parameter-recovery.tsv",
        report.parameter_recovery_rows,
    )
    brownian_fit_summary_path = write_brownian_trait_evolution_summary_table(
        output_root / "brownian-fit-summary.tsv",
        report.brownian_fit,
    )
    ou_fit_summary_path = write_ou_trait_evolution_summary_table(
        output_root / "ou-fit-summary.tsv",
        report.ou_fit,
    )
    continuous_ancestral_summary_path = write_continuous_ancestral_summary_table(
        output_root / "continuous-ancestral-summary.tsv",
        report.continuous_ancestral,
    )
    continuous_ancestral_uncertainty_path = (
        write_continuous_ancestral_uncertainty_table(
            output_root / "continuous-ancestral-uncertainty.tsv",
            report.continuous_ancestral,
        )
    )
    continuous_node_recovery_path = _write_continuous_node_recovery_table(
        output_root / "continuous-node-recovery.tsv",
        report.continuous_node_recovery_rows,
    )
    discrete_ancestral_summary_path = write_discrete_ancestral_summary_table(
        output_root / "discrete-ancestral-summary.tsv",
        report.discrete_ancestral,
    )
    discrete_ancestral_probability_path = write_discrete_ancestral_probability_table(
        output_root / "discrete-ancestral-probabilities.tsv",
        report.discrete_ancestral,
    )
    discrete_node_recovery_path = _write_discrete_node_recovery_table(
        output_root / "discrete-node-recovery.tsv",
        report.discrete_node_recovery_rows,
    )
    host_switch_summary_path = write_host_switch_summary_table(
        output_root / "host-switch-summary.tsv",
        report.host_switching,
    )
    host_state_nodes_path = write_host_state_node_table(
        output_root / "host-state-nodes.tsv",
        report.host_switching,
    )
    host_switch_branches_path = write_host_switch_branch_table(
        output_root / "host-switch-branches.tsv",
        report.host_switching,
    )
    host_node_recovery_path = _write_discrete_node_recovery_table(
        output_root / "host-node-recovery.tsv",
        report.host_node_recovery_rows,
    )
    host_event_recovery_path = _write_transition_recovery_table(
        output_root / "host-event-recovery.tsv",
        report.host_event_recovery_rows,
    )
    geographic_ancestral_summary_path = _write_geographic_summary_table(
        output_root / "geographic-ancestral-summary.tsv",
        report.geographic_states,
    )
    geographic_state_probability_path = write_node_state_probability_table(
        output_root / "geographic-state-probabilities.tsv",
        report.geographic_states,
    )
    geographic_transition_summary_path = write_transition_summary_table(
        output_root / "geographic-transition-summary.tsv",
        report.geographic_states,
    )
    geographic_node_recovery_path = _write_discrete_node_recovery_table(
        output_root / "geographic-node-recovery.tsv",
        report.geographic_node_recovery_rows,
    )
    geographic_event_recovery_path = _write_transition_recovery_table(
        output_root / "geographic-event-recovery.tsv",
        report.geographic_event_recovery_rows,
    )
    threshold_evaluation_path = _write_threshold_evaluation_table(
        output_root / "recovery-threshold-evaluation.tsv",
        report.threshold_evaluation_rows,
    )
    return KnownAnswerReferenceWorkflowBundle(
        output_root=output_root,
        rooted_topology_equal=report.tree_recovery.topology_equal,
        same_unrooted_topology=report.tree_recovery.same_unrooted_topology,
        same_taxa_different_rooting=report.tree_recovery.same_taxa_different_rooting,
        robinson_foulds_distance=report.tree_recovery.robinson_foulds_distance,
        continuous_internal_node_mean_absolute_error=continuous_mae,
        discrete_internal_node_accuracy=discrete_accuracy,
        discrete_mean_true_state_probability=discrete_mean_true_probability,
        host_internal_node_accuracy=host_node_accuracy,
        host_event_accuracy=host_event_accuracy,
        geographic_internal_node_accuracy=geographic_node_accuracy,
        geographic_event_accuracy=geographic_event_accuracy,
        parameter_row_count=len(report.parameter_recovery_rows),
        threshold_pass_count=sum(
            1 for row in report.threshold_evaluation_rows if row.passed
        ),
        threshold_row_count=len(report.threshold_evaluation_rows),
        workflow_summary_path=workflow_summary_path,
        distance_tree_path=distance_tree_path,
        tree_recovery_path=tree_recovery_path,
        parameter_recovery_path=parameter_recovery_path,
        brownian_fit_summary_path=brownian_fit_summary_path,
        ou_fit_summary_path=ou_fit_summary_path,
        continuous_ancestral_summary_path=continuous_ancestral_summary_path,
        continuous_ancestral_uncertainty_path=continuous_ancestral_uncertainty_path,
        continuous_node_recovery_path=continuous_node_recovery_path,
        discrete_ancestral_summary_path=discrete_ancestral_summary_path,
        discrete_ancestral_probability_path=discrete_ancestral_probability_path,
        discrete_node_recovery_path=discrete_node_recovery_path,
        host_switch_summary_path=host_switch_summary_path,
        host_state_nodes_path=host_state_nodes_path,
        host_switch_branches_path=host_switch_branches_path,
        host_node_recovery_path=host_node_recovery_path,
        host_event_recovery_path=host_event_recovery_path,
        geographic_ancestral_summary_path=geographic_ancestral_summary_path,
        geographic_state_probability_path=geographic_state_probability_path,
        geographic_transition_summary_path=geographic_transition_summary_path,
        geographic_node_recovery_path=geographic_node_recovery_path,
        geographic_event_recovery_path=geographic_event_recovery_path,
        threshold_evaluation_path=threshold_evaluation_path,
    )


def _write_workflow_summary_table(
    path: Path,
    *,
    report: KnownAnswerReferenceWorkflowReport,
    continuous_mae: float,
    discrete_accuracy: float,
    discrete_mean_true_probability: float,
    host_node_accuracy: float,
    host_event_accuracy: float,
    geographic_node_accuracy: float,
    geographic_event_accuracy: float,
) -> Path:
    rows = [
        "\t".join(
            [
                "dataset_id",
                "taxon_count",
                "sequence_length",
                "distance_method",
                "distance_model",
                "rooted_topology_equal",
                "same_unrooted_topology",
                "same_taxa_different_rooting",
                "robinson_foulds_distance",
                "continuous_root_absolute_error",
                "continuous_sigma_squared_absolute_error",
                "ou_alpha_absolute_error",
                "ou_theta_absolute_error",
                "ou_sigma_squared_absolute_error",
                "continuous_internal_node_mean_absolute_error",
                "discrete_internal_node_accuracy",
                "discrete_mean_true_state_probability",
                "host_internal_node_accuracy",
                "host_event_accuracy",
                "geographic_internal_node_accuracy",
                "geographic_event_accuracy",
                "threshold_pass_count",
                "threshold_row_count",
            ]
        ),
        "\t".join(
            [
                report.dataset.dataset_id,
                str(report.dataset.taxon_count),
                str(report.dataset.sequence_length),
                report.dataset.distance_method,
                report.dataset.distance_model,
                str(report.tree_recovery.topology_equal).lower(),
                str(report.tree_recovery.same_unrooted_topology).lower(),
                str(report.tree_recovery.same_taxa_different_rooting).lower(),
                str(report.tree_recovery.robinson_foulds_distance),
                format_number(report.parameter_recovery_rows[0].absolute_error),
                format_number(report.parameter_recovery_rows[1].absolute_error),
                format_number(report.parameter_recovery_rows[2].absolute_error),
                format_number(report.parameter_recovery_rows[3].absolute_error),
                format_number(report.parameter_recovery_rows[4].absolute_error),
                format_number(continuous_mae),
                format_number(discrete_accuracy),
                format_number(discrete_mean_true_probability),
                format_number(host_node_accuracy),
                format_number(host_event_accuracy),
                format_number(geographic_node_accuracy),
                format_number(geographic_event_accuracy),
                str(sum(1 for row in report.threshold_evaluation_rows if row.passed)),
                str(len(report.threshold_evaluation_rows)),
            ]
        ),
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _write_distance_tree(path: Path, newick: str) -> Path:
    path.write_text(f"{newick}\n", encoding="utf-8")
    return path


def _write_tree_recovery_table(
    path: Path,
    report: KnownAnswerReferenceWorkflowReport,
) -> Path:
    topology = report.tree_recovery
    build = report.distance_tree_build
    return write_taxon_rows(
        path,
        columns=[
            "method",
            "model",
            "taxon_count",
            "pair_count",
            "rooted_topology_equal",
            "same_unrooted_topology",
            "same_taxa_different_rooting",
            "robinson_foulds_distance",
            "rooted_robinson_foulds_distance",
            "unrooted_robinson_foulds_distance",
            "normalized_robinson_foulds",
            "rooted_normalized_robinson_foulds",
            "unrooted_normalized_robinson_foulds",
        ],
        rows=[
            {
                "method": build.method,
                "model": build.model,
                "taxon_count": str(build.taxon_count),
                "pair_count": str(build.pair_count),
                "rooted_topology_equal": str(topology.topology_equal).lower(),
                "same_unrooted_topology": str(topology.same_unrooted_topology).lower(),
                "same_taxa_different_rooting": str(
                    topology.same_taxa_different_rooting
                ).lower(),
                "robinson_foulds_distance": str(topology.robinson_foulds_distance),
                "rooted_robinson_foulds_distance": str(
                    topology.rooted_robinson_foulds_distance
                ),
                "unrooted_robinson_foulds_distance": str(
                    topology.unrooted_robinson_foulds_distance
                ),
                "normalized_robinson_foulds": format_number(
                    topology.normalized_robinson_foulds
                ),
                "rooted_normalized_robinson_foulds": format_number(
                    topology.rooted_normalized_robinson_foulds
                ),
                "unrooted_normalized_robinson_foulds": format_number(
                    topology.unrooted_normalized_robinson_foulds
                ),
            }
        ],
    )


def _write_parameter_recovery_table(
    path: Path,
    rows: list[KnownAnswerParameterRecoveryRow],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "parameter",
            "true_value",
            "estimated_value",
            "absolute_error",
            "relative_error",
            "interpretation",
        ],
        rows=[
            {
                "parameter": row.parameter,
                "true_value": format_number(row.true_value),
                "estimated_value": format_number(row.estimated_value),
                "absolute_error": format_number(row.absolute_error),
                "relative_error": format_number(row.relative_error),
                "interpretation": row.interpretation,
            }
            for row in rows
        ],
    )


def _write_continuous_node_recovery_table(
    path: Path,
    rows: list[KnownAnswerContinuousNodeRecoveryRow],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "descendant_taxa",
            "true_value",
            "estimated_value",
            "absolute_error",
            "standard_error",
            "lower_95_interval",
            "upper_95_interval",
            "confidence",
        ],
        rows=[
            {
                "node": row.node,
                "descendant_taxa": ",".join(row.descendant_taxa),
                "true_value": format_number(row.true_value),
                "estimated_value": format_number(row.estimated_value),
                "absolute_error": format_number(row.absolute_error),
                "standard_error": format_number(row.standard_error),
                "lower_95_interval": format_number(row.lower_95_interval),
                "upper_95_interval": format_number(row.upper_95_interval),
                "confidence": format_number(row.confidence),
            }
            for row in rows
        ],
    )


def _write_discrete_node_recovery_table(
    path: Path,
    rows: list[KnownAnswerDiscreteNodeRecoveryRow],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "descendant_taxa",
            "true_state",
            "estimated_state",
            "true_state_probability",
            "confidence",
            "correct",
            "ambiguous",
        ],
        rows=[
            {
                "node": row.node,
                "descendant_taxa": ",".join(row.descendant_taxa),
                "true_state": row.true_state,
                "estimated_state": row.estimated_state,
                "true_state_probability": format_number(row.true_state_probability),
                "confidence": format_number(row.confidence),
                "correct": str(row.correct).lower(),
                "ambiguous": str(row.ambiguous).lower(),
            }
            for row in rows
        ],
    )


def _write_transition_recovery_table(
    path: Path,
    rows: list[KnownAnswerTransitionRecoveryRow],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "parent_node",
            "child_node",
            "true_transition",
            "estimated_transition",
            "true_changed",
            "estimated_changed",
            "true_event_count",
            "estimated_event_count",
            "correct",
        ],
        rows=[
            {
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "true_transition": row.true_transition,
                "estimated_transition": row.estimated_transition,
                "true_changed": str(row.true_changed).lower(),
                "estimated_changed": str(row.estimated_changed).lower(),
                "true_event_count": str(row.true_event_count),
                "estimated_event_count": str(row.estimated_event_count),
                "correct": str(row.correct).lower(),
            }
            for row in rows
        ],
    )


def _write_geographic_summary_table(
    path: Path,
    report: DiscreteStateEvolutionReport,
) -> Path:
    root_estimate = report.estimates[0]
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "model",
            "taxon_count",
            "state_count",
            "transition_count",
            "strongly_supported_transition_count",
            "root_state",
            "root_confidence",
            "warning_count",
        ],
        rows=[
            {
                "trait": report.trait,
                "model": report.model,
                "taxon_count": str(report.taxon_count),
                "state_count": str(len(report.observed_states)),
                "transition_count": str(report.transition_summary.transition_count),
                "strongly_supported_transition_count": str(
                    report.transition_summary.strongly_supported_transition_count
                ),
                "root_state": root_estimate.most_likely_state,
                "root_confidence": format_number(
                    max(root_estimate.state_probabilities.values())
                ),
                "warning_count": str(len(report.warnings)),
            }
        ],
    )


def _write_threshold_evaluation_table(
    path: Path,
    rows: list[KnownAnswerThresholdEvaluationRow],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "metric",
            "comparator",
            "threshold",
            "observed_value",
            "passed",
            "rationale",
        ],
        rows=[
            {
                "metric": row.metric,
                "comparator": row.comparator,
                "threshold": row.threshold,
                "observed_value": row.observed_value,
                "passed": str(row.passed).lower(),
                "rationale": row.rationale,
            }
            for row in rows
        ],
    )
