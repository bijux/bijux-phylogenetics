from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.comparative.continuous import (
    summarize_brownian_trait_evolution,
    summarize_ou_trait_evolution,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    estimate_ancestral_geographic_states,
)
from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.distance import build_distance_tree
from bijux_phylogenetics.ecology import summarize_host_switching
from bijux_phylogenetics.io.newick import dumps_newick, write_newick

from .models import (
    CONTINUOUS_TRAIT,
    DISCRETE_TRAIT,
    KnownAnswerReferenceWorkflowReport,
)
from .panel import load_known_answer_reference_dataset
from .recovery import (
    build_continuous_node_recovery_rows,
    build_discrete_node_recovery_rows,
    build_geographic_node_recovery_rows,
    build_host_event_recovery_rows,
    build_host_node_recovery_rows,
    build_parameter_recovery_rows,
    build_threshold_evaluation_rows,
    build_transition_recovery_rows,
)
from .truth import (
    load_recovery_thresholds,
    load_true_continuous_nodes,
    load_true_discrete_nodes,
    load_true_parameter_map,
    load_true_transition_rows,
)


def run_known_answer_reference_workflow() -> KnownAnswerReferenceWorkflowReport:
    """Run the governed recovery workflow over the packaged known-answer panel."""
    dataset = load_known_answer_reference_dataset()
    true_parameters = load_true_parameter_map(dataset.true_parameters_path)
    continuous_truth = load_true_continuous_nodes(dataset.true_continuous_nodes_path)
    discrete_truth = load_true_discrete_nodes(dataset.true_discrete_nodes_path)
    host_truth = load_true_discrete_nodes(dataset.true_host_nodes_path)
    geographic_truth = load_true_discrete_nodes(dataset.true_geographic_nodes_path)
    host_event_truth = load_true_transition_rows(dataset.true_host_switch_events_path)
    geographic_event_truth = load_true_transition_rows(
        dataset.true_geographic_transition_events_path
    )
    thresholds = load_recovery_thresholds(dataset.recovery_thresholds_path)

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
    parameter_recovery_rows = build_parameter_recovery_rows(
        true_parameters=true_parameters,
        brownian_fit=brownian_fit,
        ou_fit=ou_fit,
    )
    continuous_node_recovery_rows = build_continuous_node_recovery_rows(
        true_nodes=continuous_truth,
        report=continuous_ancestral,
    )
    discrete_node_recovery_rows = build_discrete_node_recovery_rows(
        true_nodes=discrete_truth,
        report=discrete_ancestral,
    )
    host_node_recovery_rows = build_host_node_recovery_rows(
        true_nodes=host_truth,
        report=host_switching,
    )
    host_event_recovery_rows = build_host_event_recovery_rows(
        true_rows=host_event_truth,
        report=host_switching,
    )
    geographic_node_recovery_rows = build_geographic_node_recovery_rows(
        true_nodes=geographic_truth,
        report=geographic_states,
    )
    geographic_event_recovery_rows = build_transition_recovery_rows(
        true_rows=geographic_event_truth,
        report=geographic_states,
    )
    threshold_evaluation_rows = build_threshold_evaluation_rows(
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
