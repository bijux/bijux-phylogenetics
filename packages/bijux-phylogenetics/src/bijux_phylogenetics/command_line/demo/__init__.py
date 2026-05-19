from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.demo.benchmark_panels import (
    add_benchmark_demo_commands,
    run_benchmark_demo_command,
)
from bijux_phylogenetics.command_line.demo.introductory_panels import (
    add_introductory_demo_commands,
    run_introductory_demo_command,
)
from bijux_phylogenetics.command_line.demo.rabies_panels import (
    add_rabies_demo_commands,
    run_rabies_demo_command,
)
from bijux_phylogenetics.command_line.demo.sequence_panels import (
    add_sequence_demo_commands,
    run_sequence_demo_command,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.benchmark import (
    run_real_dataset_macroevolution_benchmark_demo as run_real_dataset_macroevolution_benchmark_demo,
)
from bijux_phylogenetics.datasets import (
    run_rabies_cross_host_geography_panel_demo as run_rabies_cross_host_geography_panel_demo,
)
from bijux_phylogenetics.datasets.continuous_mode_recovery import (
    run_continuous_mode_recovery_panel_demo,
)
from bijux_phylogenetics.datasets.discrete_mode_recovery import (
    run_discrete_mode_recovery_panel_demo,
)
from bijux_phylogenetics.datasets.data_quality_stress import (
    run_catarrhine_data_quality_stress_panel_demo,
)
from bijux_phylogenetics.datasets.known_answer_reference import (
    run_known_answer_reference_demo,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_demo_command(subparsers: Any) -> None:
    demo = subparsers.add_parser(
        get_command_spec("demo").name, help=get_command_spec("demo").summary
    )
    demo_subparsers = demo.add_subparsers(dest="demo_command", required=True)
    add_introductory_demo_commands(demo_subparsers)
    add_benchmark_demo_commands(demo_subparsers)
    add_sequence_demo_commands(demo_subparsers)
    add_rabies_demo_commands(demo_subparsers)
    demo_catarrhine_stress = demo_subparsers.add_parser(
        "catarrhine-data-quality-stress-panel",
        help="Materialize the packaged catarrhine dirty-data stress dataset and rerun the governed audit and cleanup outputs.",
    )
    demo_catarrhine_stress.add_argument("--out", required=True, type=Path)
    demo_catarrhine_stress.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_catarrhine_stress)
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


def run_demo_command(args: Any) -> int:
    introductory_exit_code = run_introductory_demo_command(args)
    if introductory_exit_code is not None:
        return introductory_exit_code
    benchmark_exit_code = run_benchmark_demo_command(args)
    if benchmark_exit_code is not None:
        return benchmark_exit_code
    sequence_exit_code = run_sequence_demo_command(args)
    if sequence_exit_code is not None:
        return sequence_exit_code
    rabies_exit_code = run_rabies_demo_command(args)
    if rabies_exit_code is not None:
        return rabies_exit_code

    if args.demo_command == "catarrhine-data-quality-stress-panel":
        result = run_catarrhine_data_quality_stress_panel_demo(args.out)
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.raw_alignment_path,
                result.dataset_export.raw_sequence_input_path,
                result.dataset_export.raw_coding_sequences_path,
                result.dataset_export.raw_tree_path,
                result.dataset_export.raw_traits_path,
                result.dataset_export.raw_trait_mismatch_path,
                result.workflow_bundle.workflow_summary_path,
                result.workflow_bundle.raw_sequence_findings_path,
                result.workflow_bundle.raw_sequence_repair_path,
                result.workflow_bundle.repaired_sequence_input_path,
                result.workflow_bundle.repaired_sequence_validation_path,
                result.workflow_bundle.coding_sequence_exclusions_path,
                result.workflow_bundle.prepared_coding_sequences_path,
                result.workflow_bundle.raw_trait_linkage_path,
                result.workflow_bundle.trait_duplicates_path,
                result.workflow_bundle.trait_missing_values_path,
                result.workflow_bundle.sequence_outliers_path,
                result.workflow_bundle.tree_issues_path,
                result.workflow_bundle.repair_actions_path,
                result.workflow_bundle.cleaned_traits_path,
                result.workflow_bundle.cleaned_alignment_path,
                result.workflow_bundle.cleaned_tree_path,
                result.workflow_bundle.cleaned_linkage_path,
                result.workflow_bundle.cleaned_validation_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                list(result.dataset_export.expected_output_root.glob("*"))
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "raw_taxon_count": result.workflow_bundle.raw_taxon_count,
                        "cleaned_taxon_count": (
                            result.workflow_bundle.cleaned_taxon_count
                        ),
                        "duplicate_sequence_identifier_count": (
                            result.workflow_bundle.duplicate_sequence_identifier_count
                        ),
                        "illegal_character_count": (
                            result.workflow_bundle.illegal_character_count
                        ),
                        "empty_sequence_count": (
                            result.workflow_bundle.empty_sequence_count
                        ),
                        "raw_sequence_length_outlier_count": (
                            result.workflow_bundle.raw_sequence_length_outlier_count
                        ),
                        "duplicate_trait_taxon_count": (
                            result.workflow_bundle.duplicate_trait_taxon_count
                        ),
                        "missing_trait_value_count": (
                            result.workflow_bundle.missing_trait_value_count
                        ),
                        "sequence_outlier_count": (
                            result.workflow_bundle.sequence_outlier_count
                        ),
                        "tree_zero_length_branch_count": (
                            result.workflow_bundle.tree_zero_length_branch_count
                        ),
                        "tree_negative_branch_count": (
                            result.workflow_bundle.tree_negative_branch_count
                        ),
                        "tree_long_branch_outlier_count": (
                            result.workflow_bundle.tree_long_branch_outlier_count
                        ),
                        "coding_frame_error_count": (
                            result.workflow_bundle.coding_frame_error_count
                        ),
                        "coding_internal_stop_count": (
                            result.workflow_bundle.coding_internal_stop_count
                        ),
                        "raw_trait_missing_from_traits_count": (
                            result.workflow_bundle.raw_trait_missing_from_traits_count
                        ),
                        "raw_trait_extra_taxon_count": (
                            result.workflow_bundle.raw_trait_extra_taxon_count
                        ),
                        "dropped_taxon_count": (
                            result.workflow_bundle.dropped_taxon_count
                        ),
                        "repaired_branch_count": (
                            result.workflow_bundle.repaired_branch_count
                        ),
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "continuous-mode-recovery-panel":
        result = run_continuous_mode_recovery_panel_demo(args.out)
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
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
            ],
        )
        if args.json:
            expected_output_count = len(
                [
                    path
                    for path in result.dataset_export.expected_output_root.rglob("*")
                    if path.is_file()
                ]
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "taxon_count": result.dataset.taxon_count,
                        "tree_count": result.dataset.tree_count,
                        "case_count": result.dataset.case_count,
                        "selection_review_case_count": (
                            result.workflow_bundle.selection_review_case_count
                        ),
                        "selection_match_count": (
                            result.workflow_bundle.selection_match_count
                        ),
                        "geiger_selection_match_count": (
                            result.workflow_bundle.geiger_selection_match_count
                        ),
                        "parameter_pass_count": (
                            result.workflow_bundle.parameter_pass_count
                        ),
                        "parameter_row_count": (
                            result.workflow_bundle.parameter_row_count
                        ),
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
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "discrete-mode-recovery-panel":
        result = run_discrete_mode_recovery_panel_demo(args.out)
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.default_tree_path,
                result.dataset_export.simulation_cases_path,
                result.workflow_bundle.workflow_summary_path,
                result.workflow_bundle.recovery_summary_path,
                result.workflow_bundle.rate_recovery_path,
                result.workflow_bundle.rate_comparison_path,
                result.workflow_bundle.model_choice_path,
                result.workflow_bundle.execution_review_path,
                result.workflow_bundle.warning_review_path,
                result.workflow_bundle.geiger_reference_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                [
                    path
                    for path in result.dataset_export.expected_output_root.rglob("*")
                    if path.is_file()
                ]
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "taxon_count": result.dataset.taxon_count,
                        "tree_count": result.dataset.tree_count,
                        "case_count": result.dataset.case_count,
                        "selection_review_case_count": (
                            result.workflow_bundle.selection_review_case_count
                        ),
                        "selection_match_count": (
                            result.workflow_bundle.selection_match_count
                        ),
                        "geiger_selection_match_count": (
                            result.workflow_bundle.geiger_selection_match_count
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
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "known-answer-reference-panel":
        result = run_known_answer_reference_demo(args.out)
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
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
            ],
        )
        if args.json:
            expected_output_count = len(
                list(result.dataset_export.expected_output_root.glob("*"))
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "taxon_count": result.dataset.taxon_count,
                        "sequence_length": result.dataset.sequence_length,
                        "distance_method": result.dataset.distance_method,
                        "distance_model": result.dataset.distance_model,
                        "rooted_topology_equal": (
                            result.workflow_bundle.rooted_topology_equal
                        ),
                        "same_unrooted_topology": (
                            result.workflow_bundle.same_unrooted_topology
                        ),
                        "same_taxa_different_rooting": (
                            result.workflow_bundle.same_taxa_different_rooting
                        ),
                        "robinson_foulds_distance": (
                            result.workflow_bundle.robinson_foulds_distance
                        ),
                        "parameter_row_count": (
                            result.workflow_bundle.parameter_row_count
                        ),
                        "threshold_pass_count": (
                            result.workflow_bundle.threshold_pass_count
                        ),
                        "threshold_row_count": (
                            result.workflow_bundle.threshold_row_count
                        ),
                        "continuous_internal_node_mean_absolute_error": (
                            result.workflow_bundle.continuous_internal_node_mean_absolute_error
                        ),
                        "discrete_internal_node_accuracy": (
                            result.workflow_bundle.discrete_internal_node_accuracy
                        ),
                        "host_internal_node_accuracy": (
                            result.workflow_bundle.host_internal_node_accuracy
                        ),
                        "host_event_accuracy": (
                            result.workflow_bundle.host_event_accuracy
                        ),
                        "geographic_internal_node_accuracy": (
                            result.workflow_bundle.geographic_internal_node_accuracy
                        ),
                        "geographic_event_accuracy": (
                            result.workflow_bundle.geographic_event_accuracy
                        ),
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    raise NotImplementedError(f"unsupported demo command: {args.demo_command}")
