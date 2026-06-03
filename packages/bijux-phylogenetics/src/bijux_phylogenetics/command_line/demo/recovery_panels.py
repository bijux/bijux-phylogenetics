from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.demo.shared import (
    count_expected_output_entries,
    count_expected_output_files,
    emit_demo_result,
)
from bijux_phylogenetics.datasets.continuous_mode_recovery import (
    run_continuous_mode_recovery_panel_demo,
)
from bijux_phylogenetics.datasets.discrete_mode_recovery import (
    run_discrete_mode_recovery_panel_demo,
)
from bijux_phylogenetics.datasets.known_answer_reference import (
    run_known_answer_reference_demo,
)
from bijux_phylogenetics.datasets.macroevolution_recovery_suite import (
    run_macroevolution_recovery_suite_demo,
)


def add_recovery_demo_commands(demo_subparsers: Any) -> None:
    demo_continuous_mode_recovery = demo_subparsers.add_parser(
        "continuous-mode-recovery-panel",
        help="Materialize the packaged continuous-trait recovery dataset and rerun the governed simulation-recovery outputs.",
    )
    demo_continuous_mode_recovery.add_argument("--out", required=True, type=Path)
    demo_continuous_mode_recovery.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_continuous_mode_recovery)

    demo_discrete_mode_recovery = demo_subparsers.add_parser(
        "discrete-mode-recovery-panel",
        help="Materialize the packaged discrete-trait recovery dataset and rerun the governed simulation-recovery outputs.",
    )
    demo_discrete_mode_recovery.add_argument("--out", required=True, type=Path)
    demo_discrete_mode_recovery.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_discrete_mode_recovery)

    demo_known_answer = demo_subparsers.add_parser(
        "known-answer-reference-panel",
        help="Materialize the packaged known-answer simulation dataset and rerun the governed recovery outputs.",
    )
    demo_known_answer.add_argument("--out", required=True, type=Path)
    demo_known_answer.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_known_answer)

    demo_macroevolution_suite = demo_subparsers.add_parser(
        "macroevolution-recovery-suite",
        help="Materialize the governed macroevolution recovery suite and rerun the bundled component recovery outputs.",
    )
    demo_macroevolution_suite.add_argument("--out", required=True, type=Path)
    demo_macroevolution_suite.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_macroevolution_suite)


def run_recovery_demo_command(args: Any) -> int | None:
    if args.demo_command == "continuous-mode-recovery-panel":
        result = run_continuous_mode_recovery_panel_demo(args.out)
        outputs = [
            result.dataset_export.readme_path,
            result.dataset_export.default_tree_path,
            result.dataset_export.simulation_cases_path,
            result.workflow_bundle.workflow_summary_path,
            result.workflow_bundle.recovery_summary_path,
            result.workflow_bundle.parameter_recovery_path,
            result.workflow_bundle.parameter_comparison_path,
            result.workflow_bundle.model_choice_path,
            result.workflow_bundle.execution_review_path,
            result.workflow_bundle.warning_review_path,
            result.workflow_bundle.geiger_reference_path,
            result.overview_path,
        ]
        return emit_demo_result(
            args,
            outputs=outputs,
            metrics={
                "artifact_count": len(outputs),
                "taxon_count": result.dataset.taxon_count,
                "tree_count": result.dataset.tree_count,
                "case_count": result.dataset.case_count,
                "selection_review_case_count": (
                    result.workflow_bundle.selection_review_case_count
                ),
                "selection_match_count": (result.workflow_bundle.selection_match_count),
                "geiger_selection_match_count": (
                    result.workflow_bundle.geiger_selection_match_count
                ),
                "parameter_pass_count": result.workflow_bundle.parameter_pass_count,
                "parameter_row_count": result.workflow_bundle.parameter_row_count,
                "parameter_comparison_row_count": (
                    result.workflow_bundle.parameter_comparison_row_count
                ),
                "parameter_closer_to_truth_count_bijux": (
                    result.workflow_bundle.parameter_closer_to_truth_count_bijux
                ),
                "parameter_closer_to_truth_count_geiger": (
                    result.workflow_bundle.parameter_closer_to_truth_count_geiger
                ),
                "expected_warning_case_count": (
                    result.workflow_bundle.expected_warning_case_count
                ),
                "expected_warning_present_count": (
                    result.workflow_bundle.expected_warning_present_count
                ),
                "reference_output_count": count_expected_output_files(
                    result.dataset_export.expected_output_root
                ),
            },
            data=result,
            output_root=result.output_root,
        )

    if args.demo_command == "discrete-mode-recovery-panel":
        result = run_discrete_mode_recovery_panel_demo(args.out)
        outputs = [
            result.dataset_export.readme_path,
            result.dataset_export.default_tree_path,
            result.dataset_export.simulation_cases_path,
            result.workflow_bundle.workflow_summary_path,
            result.workflow_bundle.recovery_summary_path,
            result.workflow_bundle.parameter_recovery_path,
            result.workflow_bundle.parameter_comparison_path,
            result.workflow_bundle.rate_recovery_path,
            result.workflow_bundle.rate_comparison_path,
            result.workflow_bundle.model_choice_path,
            result.workflow_bundle.execution_review_path,
            result.workflow_bundle.warning_review_path,
            result.workflow_bundle.geiger_reference_path,
            result.overview_path,
        ]
        return emit_demo_result(
            args,
            outputs=outputs,
            metrics={
                "artifact_count": len(outputs),
                "taxon_count": result.dataset.taxon_count,
                "tree_count": result.dataset.tree_count,
                "case_count": result.dataset.case_count,
                "selection_review_case_count": (
                    result.workflow_bundle.selection_review_case_count
                ),
                "selection_match_count": (result.workflow_bundle.selection_match_count),
                "geiger_selection_match_count": (
                    result.workflow_bundle.geiger_selection_match_count
                ),
                "parameter_pass_count": result.workflow_bundle.parameter_pass_count,
                "governed_parameter_row_count": (
                    result.workflow_bundle.governed_parameter_row_count
                ),
                "parameter_row_count": result.workflow_bundle.parameter_row_count,
                "parameter_comparison_row_count": (
                    result.workflow_bundle.parameter_comparison_row_count
                ),
                "parameter_closer_to_truth_count_bijux": (
                    result.workflow_bundle.parameter_closer_to_truth_count_bijux
                ),
                "parameter_closer_to_truth_count_geiger": (
                    result.workflow_bundle.parameter_closer_to_truth_count_geiger
                ),
                "rate_pass_count": result.workflow_bundle.rate_pass_count,
                "governed_rate_row_count": (
                    result.workflow_bundle.governed_rate_row_count
                ),
                "rate_row_count": result.workflow_bundle.rate_row_count,
                "governed_rate_comparison_row_count": (
                    result.workflow_bundle.governed_rate_comparison_row_count
                ),
                "rate_comparison_row_count": (
                    result.workflow_bundle.rate_comparison_row_count
                ),
                "rate_closer_to_truth_count_bijux": (
                    result.workflow_bundle.rate_closer_to_truth_count_bijux
                ),
                "rate_closer_to_truth_count_geiger": (
                    result.workflow_bundle.rate_closer_to_truth_count_geiger
                ),
                "expected_warning_case_count": (
                    result.workflow_bundle.expected_warning_case_count
                ),
                "expected_warning_present_count": (
                    result.workflow_bundle.expected_warning_present_count
                ),
                "reference_output_count": count_expected_output_files(
                    result.dataset_export.expected_output_root
                ),
            },
            data=result,
            output_root=result.output_root,
        )

    if args.demo_command == "macroevolution-recovery-suite":
        result = run_macroevolution_recovery_suite_demo(args.out)
        outputs = [
            result.dataset_export.readme_path,
            result.workflow_bundle.workflow_summary_path,
            result.workflow_bundle.component_summary_path,
            result.workflow_bundle.requirement_summary_path,
            result.workflow_bundle.sim_char_summary_path,
            result.overview_path,
        ]
        return emit_demo_result(
            args,
            outputs=outputs,
            metrics={
                "artifact_count": len(outputs),
                "component_count": result.workflow_bundle.component_count,
                "geiger_component_count": (
                    result.workflow_bundle.geiger_component_count
                ),
                "case_count": result.workflow_bundle.total_recovery_case_count,
                "geiger_case_count": (
                    result.workflow_bundle.geiger_recovery_case_count
                ),
                "max_taxon_count": result.workflow_bundle.max_taxon_count,
                "selection_review_case_count": (
                    result.workflow_bundle.selection_review_case_count
                ),
                "selection_match_count": (result.workflow_bundle.selection_match_count),
                "geiger_selection_match_count": (
                    result.workflow_bundle.geiger_selection_match_count
                ),
                "governed_value_pass_count": (
                    result.workflow_bundle.governed_value_pass_count
                ),
                "governed_value_row_count": (
                    result.workflow_bundle.governed_value_row_count
                ),
                "governed_comparison_row_count": (
                    result.workflow_bundle.governed_comparison_row_count
                ),
                "truth_threshold_pass_count": (
                    result.workflow_bundle.truth_threshold_pass_count
                ),
                "truth_threshold_row_count": (
                    result.workflow_bundle.truth_threshold_row_count
                ),
                "sim_char_case_count": result.workflow_bundle.sim_char_case_count,
                "sim_char_all_passed": result.workflow_bundle.sim_char_all_passed,
                "requirement_pass_count": (
                    result.workflow_bundle.requirement_pass_count
                ),
                "requirement_row_count": (result.workflow_bundle.requirement_row_count),
                "reference_output_count": count_expected_output_files(
                    result.dataset_export.expected_output_root
                ),
            },
            data=result,
            output_root=result.output_root,
        )

    if args.demo_command != "known-answer-reference-panel":
        return None

    result = run_known_answer_reference_demo(args.out)
    outputs = [
        result.dataset_export.readme_path,
        result.dataset_export.true_tree_path,
        result.dataset_export.alignment_path,
        result.dataset_export.continuous_traits_path,
        result.dataset_export.ou_traits_path,
        result.dataset_export.discrete_traits_path,
        result.dataset_export.host_traits_path,
        result.dataset_export.geographic_traits_path,
        result.dataset_export.true_parameters_path,
        result.dataset_export.true_continuous_nodes_path,
        result.dataset_export.true_ou_nodes_path,
        result.dataset_export.true_discrete_nodes_path,
        result.dataset_export.true_host_nodes_path,
        result.dataset_export.true_geographic_nodes_path,
        result.dataset_export.true_host_switch_events_path,
        result.dataset_export.true_geographic_transition_events_path,
        result.dataset_export.recovery_thresholds_path,
        result.workflow_bundle.workflow_summary_path,
        result.workflow_bundle.distance_tree_path,
        result.workflow_bundle.tree_recovery_path,
        result.workflow_bundle.parameter_recovery_path,
        result.workflow_bundle.brownian_fit_summary_path,
        result.workflow_bundle.ou_fit_summary_path,
        result.workflow_bundle.continuous_ancestral_summary_path,
        result.workflow_bundle.continuous_ancestral_uncertainty_path,
        result.workflow_bundle.continuous_node_recovery_path,
        result.workflow_bundle.discrete_ancestral_summary_path,
        result.workflow_bundle.discrete_ancestral_probability_path,
        result.workflow_bundle.discrete_node_recovery_path,
        result.workflow_bundle.host_switch_summary_path,
        result.workflow_bundle.host_state_nodes_path,
        result.workflow_bundle.host_switch_branches_path,
        result.workflow_bundle.host_node_recovery_path,
        result.workflow_bundle.host_event_recovery_path,
        result.workflow_bundle.geographic_ancestral_summary_path,
        result.workflow_bundle.geographic_state_probability_path,
        result.workflow_bundle.geographic_transition_summary_path,
        result.workflow_bundle.geographic_node_recovery_path,
        result.workflow_bundle.geographic_event_recovery_path,
        result.workflow_bundle.threshold_evaluation_path,
        result.overview_path,
    ]
    return emit_demo_result(
        args,
        outputs=outputs,
        metrics={
            "artifact_count": len(outputs),
            "taxon_count": result.dataset.taxon_count,
            "sequence_length": result.dataset.sequence_length,
            "distance_method": result.dataset.distance_method,
            "distance_model": result.dataset.distance_model,
            "rooted_topology_equal": result.workflow_bundle.rooted_topology_equal,
            "same_unrooted_topology": (result.workflow_bundle.same_unrooted_topology),
            "same_taxa_different_rooting": (
                result.workflow_bundle.same_taxa_different_rooting
            ),
            "robinson_foulds_distance": (
                result.workflow_bundle.robinson_foulds_distance
            ),
            "parameter_row_count": result.workflow_bundle.parameter_row_count,
            "threshold_pass_count": result.workflow_bundle.threshold_pass_count,
            "threshold_row_count": result.workflow_bundle.threshold_row_count,
            "continuous_internal_node_mean_absolute_error": (
                result.workflow_bundle.continuous_internal_node_mean_absolute_error
            ),
            "discrete_internal_node_accuracy": (
                result.workflow_bundle.discrete_internal_node_accuracy
            ),
            "host_internal_node_accuracy": (
                result.workflow_bundle.host_internal_node_accuracy
            ),
            "host_event_accuracy": result.workflow_bundle.host_event_accuracy,
            "geographic_internal_node_accuracy": (
                result.workflow_bundle.geographic_internal_node_accuracy
            ),
            "geographic_event_accuracy": (
                result.workflow_bundle.geographic_event_accuracy
            ),
            "reference_output_count": count_expected_output_entries(
                result.dataset_export.expected_output_root
            ),
        },
        data=result,
        output_root=result.output_root,
    )
