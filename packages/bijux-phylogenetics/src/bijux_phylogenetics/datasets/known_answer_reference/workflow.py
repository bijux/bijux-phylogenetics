from __future__ import annotations

import csv
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

from bijux_phylogenetics.ancestral.continuous import (
    ContinuousAncestralReport,
    reconstruct_continuous_ancestral_states,
    write_continuous_ancestral_summary_table,
    write_continuous_ancestral_uncertainty_table,
)
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
    write_discrete_ancestral_probability_table,
    write_discrete_ancestral_summary_table,
)
from bijux_phylogenetics.comparative.brownian_trait_evolution import (
    BrownianTraitEvolutionSummaryReport,
    summarize_brownian_trait_evolution,
    write_brownian_trait_evolution_summary_table,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    DiscreteStateEvolutionReport,
    estimate_ancestral_geographic_states,
    write_node_state_probability_table,
    write_transition_summary_table,
)
from bijux_phylogenetics.comparative.ou_trait_evolution import (
    OUTraitEvolutionSummaryReport,
    summarize_ou_trait_evolution,
    write_ou_trait_evolution_summary_table,
)
from bijux_phylogenetics.compare.topology import (
    TreeComparisonReport,
    compare_tree_paths,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.distance import build_distance_tree
from bijux_phylogenetics.host_association import (
    HostSwitchingReport,
    summarize_host_switching,
    write_host_state_node_table,
    write_host_switch_branch_table,
    write_host_switch_summary_table,
)
from bijux_phylogenetics.io.newick import dumps_newick, write_newick

from .export import export_known_answer_reference_dataset
from .models import (
    CONTINUOUS_TRAIT,
    DISCRETE_TRAIT,
    KnownAnswerContinuousNodeRecoveryRow,
    KnownAnswerContinuousNodeTruth,
    KnownAnswerDiscreteNodeRecoveryRow,
    KnownAnswerDiscreteNodeTruth,
    KnownAnswerParameterRecoveryRow,
    KnownAnswerRecoveryThreshold,
    KnownAnswerReferenceDemoResult,
    KnownAnswerReferenceWorkflowBundle,
    KnownAnswerReferenceWorkflowReport,
    KnownAnswerThresholdEvaluationRow,
    KnownAnswerTransitionRecoveryRow,
    KnownAnswerTransitionTruth,
)
from .panel import load_known_answer_reference_dataset


def run_known_answer_reference_workflow() -> KnownAnswerReferenceWorkflowReport:
    """Run the governed recovery workflow over the packaged known-answer panel."""
    dataset = load_known_answer_reference_dataset()
    true_parameters = _load_true_parameter_map(dataset.true_parameters_path)
    continuous_truth = _load_true_continuous_nodes(dataset.true_continuous_nodes_path)
    discrete_truth = _load_true_discrete_nodes(dataset.true_discrete_nodes_path)
    host_truth = _load_true_discrete_nodes(dataset.true_host_nodes_path)
    geographic_truth = _load_true_discrete_nodes(dataset.true_geographic_nodes_path)
    host_event_truth = _load_true_transition_rows(dataset.true_host_switch_events_path)
    geographic_event_truth = _load_true_transition_rows(
        dataset.true_geographic_transition_events_path
    )
    thresholds = _load_recovery_thresholds(dataset.recovery_thresholds_path)

    distance_tree, distance_tree_build = build_distance_tree(
        dataset.alignment_path,
        method=dataset.distance_method,
        model=dataset.distance_model,
    )
    distance_tree_newick = dumps_newick(distance_tree)
    with TemporaryDirectory(prefix="known-answer-reference-") as temporary_root:
        built_tree_path = write_newick(
            Path(temporary_root) / "recovered-distance-tree.nwk",
            distance_tree,
        )
        tree_recovery = compare_tree_paths(built_tree_path, dataset.true_tree_path)

    brownian_fit = summarize_brownian_trait_evolution(
        dataset.true_tree_path,
        dataset.continuous_traits_path,
        trait=CONTINUOUS_TRAIT,
    )
    ou_fit = summarize_ou_trait_evolution(
        dataset.true_tree_path,
        dataset.ou_traits_path,
        trait=CONTINUOUS_TRAIT,
    )
    continuous_ancestral = reconstruct_continuous_ancestral_states(
        dataset.true_tree_path,
        dataset.continuous_traits_path,
        trait=CONTINUOUS_TRAIT,
        model="brownian",
    )
    discrete_ancestral = reconstruct_discrete_ancestral_states(
        dataset.true_tree_path,
        dataset.discrete_traits_path,
        trait=DISCRETE_TRAIT,
        model="equal-rates",
    )
    host_switching = summarize_host_switching(
        dataset.true_tree_path,
        dataset.host_traits_path,
        trait="host_group",
        model="ard",
    )
    geographic_states = estimate_ancestral_geographic_states(
        dataset.true_tree_path,
        dataset.geographic_traits_path,
        trait="region_group",
        model="equal-rates",
    )
    parameter_recovery_rows = _build_parameter_recovery_rows(
        true_parameters=true_parameters,
        brownian_fit=brownian_fit,
        ou_fit=ou_fit,
    )
    continuous_node_recovery_rows = _build_continuous_node_recovery_rows(
        true_nodes=continuous_truth,
        report=continuous_ancestral,
    )
    discrete_node_recovery_rows = _build_discrete_node_recovery_rows(
        true_nodes=discrete_truth,
        report=discrete_ancestral,
    )
    host_node_recovery_rows = _build_host_node_recovery_rows(
        true_nodes=host_truth,
        report=host_switching,
    )
    host_event_recovery_rows = _build_host_event_recovery_rows(
        true_rows=host_event_truth,
        report=host_switching,
    )
    geographic_node_recovery_rows = _build_geographic_node_recovery_rows(
        true_nodes=geographic_truth,
        report=geographic_states,
    )
    geographic_event_recovery_rows = _build_transition_recovery_rows(
        true_rows=geographic_event_truth,
        report=geographic_states,
    )
    threshold_evaluation_rows = _build_threshold_evaluation_rows(
        thresholds=thresholds,
        tree_recovery=tree_recovery,
        parameter_recovery_rows=parameter_recovery_rows,
        discrete_node_recovery_rows=discrete_node_recovery_rows,
        host_node_recovery_rows=host_node_recovery_rows,
        host_event_recovery_rows=host_event_recovery_rows,
        geographic_node_recovery_rows=geographic_node_recovery_rows,
        geographic_event_recovery_rows=geographic_event_recovery_rows,
    )
    return KnownAnswerReferenceWorkflowReport(
        dataset=dataset,
        distance_tree_build=distance_tree_build,
        distance_tree_newick=distance_tree_newick,
        tree_recovery=tree_recovery,
        brownian_fit=brownian_fit,
        ou_fit=ou_fit,
        continuous_ancestral=continuous_ancestral,
        discrete_ancestral=discrete_ancestral,
        host_switching=host_switching,
        geographic_states=geographic_states,
        parameter_recovery_rows=parameter_recovery_rows,
        continuous_node_recovery_rows=continuous_node_recovery_rows,
        discrete_node_recovery_rows=discrete_node_recovery_rows,
        host_node_recovery_rows=host_node_recovery_rows,
        host_event_recovery_rows=host_event_recovery_rows,
        geographic_node_recovery_rows=geographic_node_recovery_rows,
        geographic_event_recovery_rows=geographic_event_recovery_rows,
        threshold_evaluation_rows=threshold_evaluation_rows,
    )


def write_known_answer_reference_workflow_bundle(
    output_root: Path,
    report: KnownAnswerReferenceWorkflowReport,
) -> KnownAnswerReferenceWorkflowBundle:
    """Write the governed recovery outputs for the packaged simulation panel."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    continuous_mae = _mean(
        row.absolute_error for row in report.continuous_node_recovery_rows
    )
    discrete_accuracy = _mean(
        1.0 if row.correct else 0.0 for row in report.discrete_node_recovery_rows
    )
    discrete_mean_true_probability = _mean(
        row.true_state_probability for row in report.discrete_node_recovery_rows
    )
    host_node_accuracy = _mean(
        1.0 if row.correct else 0.0 for row in report.host_node_recovery_rows
    )
    host_event_accuracy = _mean(
        1.0 if row.correct else 0.0 for row in report.host_event_recovery_rows
    )
    geographic_node_accuracy = _mean(
        1.0 if row.correct else 0.0 for row in report.geographic_node_recovery_rows
    )
    geographic_event_accuracy = _mean(
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


def run_known_answer_reference_demo(
    output_root: Path,
) -> KnownAnswerReferenceDemoResult:
    """Materialize the packaged simulation dataset and rerun the recovery outputs."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    report = run_known_answer_reference_workflow()
    dataset_export = export_known_answer_reference_dataset(output_root / "dataset")
    workflow_bundle = write_known_answer_reference_workflow_bundle(
        output_root / "workflow",
        report,
    )
    overview_path = _write_overview(
        output_root / "overview.md", report, workflow_bundle
    )
    return KnownAnswerReferenceDemoResult(
        output_root=output_root,
        dataset=report.dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _build_parameter_recovery_rows(
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


def _build_continuous_node_recovery_rows(
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


def _build_discrete_node_recovery_rows(
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


def _build_geographic_node_recovery_rows(
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


def _build_host_node_recovery_rows(
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


def _build_host_event_recovery_rows(
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
            true_transition=_format_transition(row.source_state, row.target_state),
            estimated_transition=_format_transition(
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
            correct=_transition_recovery_match(
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


def _build_transition_recovery_rows(
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
            true_transition=_format_transition(row.source_state, row.target_state),
            estimated_transition=_format_transition(
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
            correct=_transition_recovery_match(
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


def _build_threshold_evaluation_rows(
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
        "discrete_internal_node_accuracy": _mean(
            1.0 if row.correct else 0.0 for row in discrete_node_recovery_rows
        ),
        "host_internal_node_accuracy": _mean(
            1.0 if row.correct else 0.0 for row in host_node_recovery_rows
        ),
        "host_event_accuracy": _mean(
            1.0 if row.correct else 0.0 for row in host_event_recovery_rows
        ),
        "geographic_internal_node_accuracy": _mean(
            1.0 if row.correct else 0.0 for row in geographic_node_recovery_rows
        ),
        "geographic_event_accuracy": _mean(
            1.0 if row.correct else 0.0 for row in geographic_event_recovery_rows
        ),
    }
    rows: list[KnownAnswerThresholdEvaluationRow] = []
    for threshold in thresholds:
        observed_value = metric_values[threshold.metric]
        passed = _evaluate_threshold(
            comparator=threshold.comparator,
            threshold=threshold.threshold,
            observed_value=observed_value,
        )
        rows.append(
            KnownAnswerThresholdEvaluationRow(
                metric=threshold.metric,
                comparator=threshold.comparator,
                threshold=threshold.threshold,
                observed_value=_format_observed_value(observed_value),
                passed=passed,
                rationale=threshold.rationale,
            )
        )
    return rows


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
                _format_number(report.parameter_recovery_rows[0].absolute_error),
                _format_number(report.parameter_recovery_rows[1].absolute_error),
                _format_number(report.parameter_recovery_rows[2].absolute_error),
                _format_number(report.parameter_recovery_rows[3].absolute_error),
                _format_number(report.parameter_recovery_rows[4].absolute_error),
                _format_number(continuous_mae),
                _format_number(discrete_accuracy),
                _format_number(discrete_mean_true_probability),
                _format_number(host_node_accuracy),
                _format_number(host_event_accuracy),
                _format_number(geographic_node_accuracy),
                _format_number(geographic_event_accuracy),
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
                "normalized_robinson_foulds": _format_number(
                    topology.normalized_robinson_foulds
                ),
                "rooted_normalized_robinson_foulds": _format_number(
                    topology.rooted_normalized_robinson_foulds
                ),
                "unrooted_normalized_robinson_foulds": _format_number(
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
                "true_value": _format_number(row.true_value),
                "estimated_value": _format_number(row.estimated_value),
                "absolute_error": _format_number(row.absolute_error),
                "relative_error": _format_number(row.relative_error),
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
                "true_value": _format_number(row.true_value),
                "estimated_value": _format_number(row.estimated_value),
                "absolute_error": _format_number(row.absolute_error),
                "standard_error": _format_number(row.standard_error),
                "lower_95_interval": _format_number(row.lower_95_interval),
                "upper_95_interval": _format_number(row.upper_95_interval),
                "confidence": _format_number(row.confidence),
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
                "true_state_probability": _format_number(row.true_state_probability),
                "confidence": _format_number(row.confidence),
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
                "root_confidence": _format_number(
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


def _write_overview(
    path: Path,
    report: KnownAnswerReferenceWorkflowReport,
    bundle: KnownAnswerReferenceWorkflowBundle,
) -> Path:
    lines = [
        "# Known-Answer Simulation Demo",
        "",
        f"- dataset id: `{report.dataset.dataset_id}`",
        f"- taxon count: `{report.dataset.taxon_count}`",
        f"- alignment length: `{report.dataset.sequence_length}`",
        f"- distance recovery preserves rooted topology: `{str(bundle.rooted_topology_equal).lower()}`",
        f"- distance recovery preserves unrooted topology: `{str(bundle.same_unrooted_topology).lower()}`",
        f"- continuous internal-node mean absolute error: `{_format_number(bundle.continuous_internal_node_mean_absolute_error)}`",
        f"- discrete internal-node accuracy: `{_format_number(bundle.discrete_internal_node_accuracy)}`",
        f"- host internal-node accuracy: `{_format_number(bundle.host_internal_node_accuracy)}`",
        f"- host event accuracy: `{_format_number(bundle.host_event_accuracy)}`",
        f"- geographic internal-node accuracy: `{_format_number(bundle.geographic_internal_node_accuracy)}`",
        f"- geographic event accuracy: `{_format_number(bundle.geographic_event_accuracy)}`",
        f"- threshold passes: `{bundle.threshold_pass_count}/{bundle.threshold_row_count}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{bundle.workflow_summary_path.name}`",
        f"- recovered distance tree: `{bundle.distance_tree_path.name}`",
        f"- tree recovery ledger: `{bundle.tree_recovery_path.name}`",
        f"- parameter recovery ledger: `{bundle.parameter_recovery_path.name}`",
        f"- OU fit summary: `{bundle.ou_fit_summary_path.name}`",
        f"- continuous node recovery ledger: `{bundle.continuous_node_recovery_path.name}`",
        f"- discrete node recovery ledger: `{bundle.discrete_node_recovery_path.name}`",
        f"- host event recovery ledger: `{bundle.host_event_recovery_path.name}`",
        f"- geographic event recovery ledger: `{bundle.geographic_event_recovery_path.name}`",
        f"- threshold evaluation ledger: `{bundle.threshold_evaluation_path.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _load_true_parameter_map(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["parameter"]: row["value"] for row in reader}


def _load_true_continuous_nodes(path: Path) -> list[KnownAnswerContinuousNodeTruth]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            KnownAnswerContinuousNodeTruth(
                node=row["node"],
                node_name=row["node_name"] or None,
                is_tip=row["is_tip"].strip().lower() == "true",
                descendant_taxa=_split_descendant_taxa(row["descendant_taxa"]),
                true_value=float(row["true_value"]),
            )
            for row in reader
        ]


def _load_true_discrete_nodes(path: Path) -> list[KnownAnswerDiscreteNodeTruth]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            KnownAnswerDiscreteNodeTruth(
                node=row["node"],
                node_name=row["node_name"] or None,
                is_tip=row["is_tip"].strip().lower() == "true",
                descendant_taxa=_split_descendant_taxa(row["descendant_taxa"]),
                true_state=row["true_state"],
            )
            for row in reader
        ]


def _load_true_transition_rows(path: Path) -> list[KnownAnswerTransitionTruth]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            KnownAnswerTransitionTruth(
                parent_node=row["parent_node"],
                child_node=row["child_node"],
                branch_length=float(row["branch_length"]),
                source_state=row["source_state"],
                target_state=row["target_state"],
                changed=row["changed"].strip().lower() == "true",
                event_count=int(row["event_count"]),
            )
            for row in reader
        ]


def _load_recovery_thresholds(path: Path) -> list[KnownAnswerRecoveryThreshold]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            KnownAnswerRecoveryThreshold(
                metric=row["metric"],
                comparator=row["comparator"],
                threshold=row["threshold"],
                rationale=row["rationale"],
            )
            for row in reader
        ]


def _split_descendant_taxa(value: str) -> list[str]:
    if not value:
        return []
    return value.split(",")


def _transition_recovery_match(
    *,
    true_changed: bool,
    estimated_changed: bool,
    true_source: str,
    true_target: str,
    estimated_source: str,
    estimated_target: str,
    true_event_count: int,
    estimated_event_count: int,
) -> bool:
    return (
        true_changed == estimated_changed
        and true_source == estimated_source
        and true_target == estimated_target
        and true_event_count == estimated_event_count
    )


def _format_transition(source_state: str, target_state: str) -> str:
    return f"{source_state}->{target_state}"


def _evaluate_threshold(
    *,
    comparator: str,
    threshold: str,
    observed_value: bool | float,
) -> bool:
    if comparator == "==":
        if isinstance(observed_value, bool):
            return observed_value is (threshold.strip().lower() == "true")
        return float(observed_value) == float(threshold)
    threshold_value = float(threshold)
    numeric_observed = float(observed_value)
    if comparator == "<=":
        return numeric_observed <= threshold_value
    if comparator == ">=":
        return numeric_observed >= threshold_value
    raise ValueError(f"unsupported threshold comparator '{comparator}'")


def _mean(values) -> float:
    materialized = list(values)
    if not materialized:
        return 0.0
    return sum(materialized) / len(materialized)


def _format_number(value: float) -> str:
    return format(value, ".15g")


def _format_observed_value(value: bool | float) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    return _format_number(value)
