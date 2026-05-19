from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.report.figure_packages import (
    add_figure_package_report_commands,
    run_figure_package_report_command,
)
from bijux_phylogenetics.command_line.report.input_reports import (
    add_input_report_commands,
    run_input_report_command,
)
from bijux_phylogenetics.command_line.report.methods import (
    add_methods_report_commands,
    run_methods_report_command,
)
from bijux_phylogenetics.command_line.report.publication import (
    add_publication_report_commands,
    run_publication_report_command,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports.service import (
    render_level_one_release_gate_report,
    render_production_scale_readiness_report,
    render_release_truth_report,
    render_workflow_validation_report,
)
from bijux_phylogenetics.reports import (
    write_supplementary_alignment_diagnostics_table,
    write_supplementary_ancestral_state_table,
    write_supplementary_batch_summary_table,
    write_supplementary_clade_support_table,
    write_supplementary_comparative_model_table,
    write_supplementary_diversification_table,
    write_supplementary_model_selection_table,
    write_supplementary_taxon_table,
    write_supplementary_tree_diagnostics_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def _build_production_scale_alignment_classes(
    args: Any,
) -> list[tuple[str, int, int]] | None:
    if args.sequence_counts is None and args.alignment_lengths is None:
        return None
    sequence_counts = args.sequence_counts or []
    alignment_lengths = args.alignment_lengths or []
    if len(sequence_counts) != len(alignment_lengths):
        raise ValueError(
            "report production-scale-readiness requires the same number of --sequence-count and --alignment-length values"
        )
    return [
        (
            f"sequences-{sequence_count}-sites-{alignment_length}",
            sequence_count,
            alignment_length,
        )
        for sequence_count, alignment_length in zip(
            sequence_counts,
            alignment_lengths,
            strict=True,
        )
    ]


def _build_production_scale_tree_set_classes(
    args: Any,
) -> list[tuple[str, int, int]] | None:
    if args.posterior_tree_counts is None and args.tree_set_tip_counts is None:
        return None
    posterior_tree_counts = args.posterior_tree_counts or []
    tree_set_tip_counts = args.tree_set_tip_counts or []
    if len(posterior_tree_counts) != len(tree_set_tip_counts):
        raise ValueError(
            "report production-scale-readiness requires the same number of --posterior-tree-count and --tree-set-tip-count values"
        )
    return [
        (f"trees-{tree_count}-taxa-{tip_count}", tree_count, tip_count)
        for tree_count, tip_count in zip(
            posterior_tree_counts,
            tree_set_tip_counts,
            strict=True,
        )
    ]


def add_report_command(subparsers: Any) -> None:
    report = subparsers.add_parser(
        get_command_spec("report").name, help=get_command_spec("report").summary
    )
    report_subparsers = report.add_subparsers(dest="report_command", required=True)
    add_publication_report_commands(report_subparsers)
    add_methods_report_commands(report_subparsers)
    add_figure_package_report_commands(report_subparsers)
    add_input_report_commands(report_subparsers)

    report_supplementary_taxon_table = report_subparsers.add_parser(
        "supplementary-taxon-table",
        help="Write a supplementary taxon table with IDs, metadata, traits, and exclusion evidence.",
    )
    report_supplementary_taxon_table.add_argument("--tree", required=True, type=Path)
    report_supplementary_taxon_table.add_argument(
        "--metadata", required=True, type=Path
    )
    report_supplementary_taxon_table.add_argument("--traits", required=True, type=Path)
    report_supplementary_taxon_table.add_argument("--alignment", type=Path)
    report_supplementary_taxon_table.add_argument("--filtered-alignment", type=Path)
    report_supplementary_taxon_table.add_argument("--inference-tree", type=Path)
    report_supplementary_taxon_table.add_argument("--reported-taxa", type=Path)
    report_supplementary_taxon_table.add_argument("--tip-dates", type=Path)
    report_supplementary_taxon_table.add_argument("--calibrations", type=Path)
    report_supplementary_taxon_table.add_argument("--out", required=True, type=Path)
    report_supplementary_taxon_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_taxon_table)

    report_supplementary_alignment_table = report_subparsers.add_parser(
        "supplementary-alignment-table",
        help="Write a supplementary alignment diagnostics table with optional filtering status.",
    )
    report_supplementary_alignment_table.add_argument(
        "--alignment", required=True, type=Path
    )
    report_supplementary_alignment_table.add_argument("--filtered-alignment", type=Path)
    report_supplementary_alignment_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_alignment_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_alignment_table)

    report_supplementary_tree_table = report_subparsers.add_parser(
        "supplementary-tree-table",
        help="Write a supplementary tree diagnostics table with topology, support, and warning summaries.",
    )
    report_supplementary_tree_table.add_argument("--tree", required=True, type=Path)
    report_supplementary_tree_table.add_argument("--out", required=True, type=Path)
    report_supplementary_tree_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_tree_table)

    report_supplementary_clade_support_table = report_subparsers.add_parser(
        "supplementary-clade-support-table",
        help="Write a supplementary clade-support table from one reference tree and optional tree-set frequencies.",
    )
    report_supplementary_clade_support_table.add_argument(
        "--tree", required=True, type=Path
    )
    report_supplementary_clade_support_table.add_argument(
        "--comparison-tree-set", type=Path
    )
    report_supplementary_clade_support_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_clade_support_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_clade_support_table)

    report_supplementary_model_selection_table = report_subparsers.add_parser(
        "supplementary-model-selection-table",
        help="Write a supplementary model-selection table from parsed IQ-TREE report artifacts.",
    )
    report_supplementary_model_selection_table.add_argument(
        "--iqtree-report", required=True, type=Path
    )
    report_supplementary_model_selection_table.add_argument("--model-sidecar", type=Path)
    report_supplementary_model_selection_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_model_selection_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_model_selection_table)

    report_supplementary_comparative_model_table = report_subparsers.add_parser(
        "supplementary-comparative-model-table",
        help="Write a supplementary comparative-model table with coefficients, uncertainty, diagnostics, and exclusions.",
    )
    report_supplementary_comparative_model_table.add_argument(
        "--tree", required=True, type=Path
    )
    report_supplementary_comparative_model_table.add_argument(
        "--traits", required=True, type=Path
    )
    report_supplementary_comparative_model_table.add_argument(
        "--formula",
        dest="formulas",
        action="append",
        required=True,
        help="Add one comparative candidate formula. Repeat for each candidate model.",
    )
    report_supplementary_comparative_model_table.add_argument("--taxon-column")
    report_supplementary_comparative_model_table.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    report_supplementary_comparative_model_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_comparative_model_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_comparative_model_table)

    report_supplementary_ancestral_state_table = report_subparsers.add_parser(
        "supplementary-ancestral-state-table",
        help="Write a supplementary ancestral-state table with node estimates, uncertainty, and model settings.",
    )
    report_supplementary_ancestral_state_table.add_argument(
        "--tree", required=True, type=Path
    )
    report_supplementary_ancestral_state_table.add_argument(
        "--traits", required=True, type=Path
    )
    report_supplementary_ancestral_state_table.add_argument("--trait", required=True)
    report_supplementary_ancestral_state_table.add_argument(
        "--reconstruction-kind",
        choices=("continuous", "discrete"),
        required=True,
    )
    report_supplementary_ancestral_state_table.add_argument("--taxon-column")
    report_supplementary_ancestral_state_table.add_argument("--model")
    report_supplementary_ancestral_state_table.add_argument("--estimator")
    report_supplementary_ancestral_state_table.add_argument(
        "--alpha", type=float, default=1.0
    )
    report_supplementary_ancestral_state_table.add_argument(
        "--state-ordering",
        choices=("unordered", "ordered"),
        default="unordered",
    )
    report_supplementary_ancestral_state_table.add_argument(
        "--ordered-state",
        dest="ordered_states",
        action="append",
        help="Add one ordered discrete state. Repeat to define the full state order.",
    )
    report_supplementary_ancestral_state_table.add_argument(
        "--root-prior-mode",
        choices=("equal", "observed-frequency", "fixed"),
        default="equal",
    )
    report_supplementary_ancestral_state_table.add_argument("--fixed-root-state")
    report_supplementary_ancestral_state_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_ancestral_state_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_ancestral_state_table)

    report_supplementary_diversification_table = report_subparsers.add_parser(
        "supplementary-diversification-table",
        help="Write a supplementary diversification table with clade rates, model estimates, sampling correction, and warnings.",
    )
    report_supplementary_diversification_table.add_argument(
        "--tree", required=True, type=Path
    )
    report_supplementary_diversification_table.add_argument("--metadata", type=Path)
    report_supplementary_diversification_table.add_argument("--taxon-column")
    report_supplementary_diversification_table.add_argument("--sampling-column")
    report_supplementary_diversification_table.add_argument(
        "--clade-model",
        choices=("yule", "birth-death"),
        default="birth-death",
    )
    report_supplementary_diversification_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_diversification_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_diversification_table)

    report_supplementary_batch_summary_table = report_subparsers.add_parser(
        "supplementary-batch-summary-table",
        help="Write a supplementary batch summary table from a written workflow bundle.",
    )
    report_supplementary_batch_summary_table.add_argument(
        "--workflow-bundle-root", required=True, type=Path
    )
    report_supplementary_batch_summary_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_batch_summary_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_batch_summary_table)

    report_workflow_validation = report_subparsers.add_parser(
        "workflow-validation",
        help="Render the Level 1 workflow validation fixture report.",
    )
    report_workflow_validation.add_argument("--fixtures-root", type=Path)
    report_workflow_validation.add_argument("--out", required=True, type=Path)
    report_workflow_validation.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_workflow_validation)

    report_release_gate = report_subparsers.add_parser(
        "release-gate",
        help="Render the Level 1 release gate for the checked-in workflow fixtures.",
    )
    report_release_gate.add_argument("--fixtures-root", type=Path)
    report_release_gate.add_argument("--out", required=True, type=Path)
    report_release_gate.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_release_gate)

    report_release_truth = report_subparsers.add_parser(
        "release-truth",
        help="Render one machine-produced release truth report from pytest and workflow evidence.",
    )
    report_release_truth.add_argument(
        "--test-report",
        type=Path,
        action="append",
        required=True,
        help="Path to one pytest JUnit XML report for the full test surface. Repeat to aggregate multiple sessions.",
    )
    report_release_truth.add_argument(
        "--real-engine-test-report",
        type=Path,
        action="append",
        required=True,
        help="Path to one pytest JUnit XML report for real-engine tests. Repeat to aggregate multiple sessions.",
    )
    report_release_truth.add_argument("--fixtures-root", type=Path)
    report_release_truth.add_argument(
        "--stress-tier",
        choices=("small", "heavy"),
        default="small",
        help="Governed stress tier to benchmark during release truth generation.",
    )
    report_release_truth.add_argument(
        "--parity-extended",
        action="store_true",
        help="Include the governed extended reference-parity suite in the release truth report.",
    )
    report_release_truth.add_argument("--out", required=True, type=Path)
    report_release_truth.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_release_truth)

    report_production_scale_readiness = report_subparsers.add_parser(
        "production-scale-readiness",
        help="Render one reviewer-facing production-scale readiness report from governed benchmark evidence.",
    )
    report_production_scale_readiness.add_argument("--replicates", type=int, default=1)
    report_production_scale_readiness.add_argument(
        "--tree-tip-count",
        action="append",
        dest="tree_tip_counts",
        type=int,
        help="Add one large-tree taxon count. Repeat to override the governed tree-size classes.",
    )
    report_production_scale_readiness.add_argument(
        "--sequence-count",
        action="append",
        dest="sequence_counts",
        type=int,
        help="Add one sequence count for the large-alignment classes. Repeat alongside --alignment-length.",
    )
    report_production_scale_readiness.add_argument(
        "--alignment-length",
        action="append",
        dest="alignment_lengths",
        type=int,
        help="Add one aligned-site count for the large-alignment classes. Repeat alongside --sequence-count.",
    )
    report_production_scale_readiness.add_argument(
        "--posterior-tree-count",
        action="append",
        dest="posterior_tree_counts",
        type=int,
        help="Add one posterior tree count for the tree-set classes. Repeat alongside --tree-set-tip-count.",
    )
    report_production_scale_readiness.add_argument(
        "--tree-set-tip-count",
        action="append",
        dest="tree_set_tip_counts",
        type=int,
        help="Add one taxon count for the tree-set classes. Repeat alongside --posterior-tree-count.",
    )
    report_production_scale_readiness.add_argument(
        "--stress-tier",
        action="append",
        dest="stress_tiers",
        choices=("small", "heavy"),
        help="Include one governed stress tier. Repeat to aggregate multiple tiers.",
    )
    report_production_scale_readiness.add_argument("--out", required=True, type=Path)
    report_production_scale_readiness.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_production_scale_readiness)


def run_report_command(args: Any) -> int:
    publication_exit_code = run_publication_report_command(args)
    if publication_exit_code is not None:
        return publication_exit_code
    methods_exit_code = run_methods_report_command(args)
    if methods_exit_code is not None:
        return methods_exit_code
    figure_package_exit_code = run_figure_package_report_command(args)
    if figure_package_exit_code is not None:
        return figure_package_exit_code
    input_report_exit_code = run_input_report_command(args)
    if input_report_exit_code is not None:
        return input_report_exit_code

    parsed_lambda_value: float | str
    if getattr(args, "lambda_value", None) == "estimate":
        parsed_lambda_value = "estimate"
    elif getattr(args, "lambda_value", None) is None:
        parsed_lambda_value = "estimate"
    else:
        parsed_lambda_value = float(args.lambda_value)

    if args.report_command == "supplementary-taxon-table":
        result = write_supplementary_taxon_table(
            args.out,
            tree_path=args.tree,
            metadata_path=args.metadata,
            traits_path=args.traits,
            alignment_path=args.alignment,
            filtered_alignment_path=args.filtered_alignment,
            inference_tree_path=args.inference_tree,
            reported_taxa_path=args.reported_taxa,
            tip_dates_path=args.tip_dates,
            calibration_path=args.calibrations,
        )
        inputs = [args.tree, args.metadata, args.traits]
        if args.alignment is not None:
            inputs.append(args.alignment)
        if args.filtered_alignment is not None:
            inputs.append(args.filtered_alignment)
        if args.inference_tree is not None:
            inputs.append(args.inference_tree)
        if args.reported_taxa is not None:
            inputs.append(args.reported_taxa)
        if args.tip_dates is not None:
            inputs.append(args.tip_dates)
        if args.calibrations is not None:
            inputs.append(args.calibrations)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "analysis_included_count": result.analysis_included_count,
                        "analysis_excluded_count": result.analysis_excluded_count,
                        "reporting_retained_count": result.reporting_retained_count,
                        "reporting_dropped_count": result.reporting_dropped_count,
                        "metadata_column_count": result.metadata_column_count,
                        "trait_column_count": result.trait_column_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-alignment-table":
        result = write_supplementary_alignment_diagnostics_table(
            args.out,
            alignment_path=args.alignment,
            filtered_alignment_path=args.filtered_alignment,
        )
        inputs = [args.alignment]
        if args.filtered_alignment is not None:
            inputs.append(args.filtered_alignment)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "retained_sequence_count": result.retained_sequence_count,
                        "removed_sequence_count": result.removed_sequence_count,
                        "filtered_only_sequence_count": result.filtered_only_sequence_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-tree-table":
        result = write_supplementary_tree_diagnostics_table(
            args.out,
            tree_path=args.tree,
        )
        inputs = [args.tree]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            row = result.rows[0]
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=row.warnings,
                    metrics={
                        "row_count": result.row_count,
                        "tip_count": row.tip_count,
                        "supported_branch_count": row.supported_branch_count,
                        "polytomy_count": row.polytomy_count,
                        "warning_count": row.warning_count,
                        "ultrametric": row.ultrametric,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-clade-support-table":
        result = write_supplementary_clade_support_table(
            args.out,
            tree_path=args.tree,
            comparison_tree_set_path=args.comparison_tree_set,
        )
        inputs = [args.tree]
        if args.comparison_tree_set is not None:
            inputs.append(args.comparison_tree_set)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "supported_clade_count": result.supported_clade_count,
                        "frequency_scored_clade_count": result.frequency_scored_clade_count,
                        "frequency_partial_support_count": result.frequency_partial_support_count,
                        "frequency_absent_clade_count": result.frequency_absent_clade_count,
                        "frequency_unscored_clade_count": result.frequency_unscored_clade_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-model-selection-table":
        result = write_supplementary_model_selection_table(
            args.out,
            iqtree_report_path=args.iqtree_report,
            model_sidecar_path=args.model_sidecar,
        )
        inputs = [args.iqtree_report]
        if args.model_sidecar is not None:
            inputs.append(args.model_sidecar)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "candidate_count": result.candidate_count,
                        "selected_model": result.selected_model,
                        "selected_criterion": result.selected_criterion,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-comparative-model-table":
        result = write_supplementary_comparative_model_table(
            args.out,
            tree_path=args.tree,
            traits_path=args.traits,
            formulas=list(args.formulas),
            taxon_column=args.taxon_column,
            lambda_value=parsed_lambda_value,
        )
        inputs = [args.tree, args.traits]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "model_count": result.model_count,
                        "selected_formula": result.selected_formula,
                        "selected_criterion": result.selected_criterion,
                        "excluded_taxon_count": result.excluded_taxon_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-ancestral-state-table":
        result = write_supplementary_ancestral_state_table(
            args.out,
            tree_path=args.tree,
            traits_path=args.traits,
            trait=args.trait,
            reconstruction_kind=args.reconstruction_kind,
            taxon_column=args.taxon_column,
            model=args.model,
            estimator=args.estimator,
            alpha=args.alpha,
            state_ordering=args.state_ordering,
            ordered_states=args.ordered_states,
            root_prior_mode=args.root_prior_mode,
            fixed_root_state=args.fixed_root_state,
        )
        inputs = [args.tree, args.traits]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "reconstruction_kind": result.reconstruction_kind,
                        "model": result.model,
                        "analysis_taxon_count": result.analysis_taxon_count,
                        "excluded_taxon_count": result.excluded_taxon_count,
                        "unstable_node_count": result.unstable_node_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-diversification-table":
        inputs = [args.tree]
        if args.metadata is not None:
            inputs.append(args.metadata)
        result = write_supplementary_diversification_table(
            args.out,
            tree_path=args.tree,
            metadata_path=args.metadata,
            taxon_column=args.taxon_column,
            sampling_column=args.sampling_column,
            clade_model=args.clade_model,
        )
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "better_model": result.better_model,
                        "clade_model": result.clade_model,
                        "high_clade_count": result.high_clade_count,
                        "low_clade_count": result.low_clade_count,
                        "warning_count": result.warning_count,
                        "sampling_metadata_complete": (
                            result.sampling_metadata_complete
                        ),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-batch-summary-table":
        result = write_supplementary_batch_summary_table(
            args.out,
            workflow_bundle_root=args.workflow_bundle_root,
        )
        inputs = [args.workflow_bundle_root]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "dataset_row_count": result.dataset_row_count,
                        "variant_row_count": result.variant_row_count,
                        "workflow_status": result.workflow_status,
                        "warning_count": result.warning_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "workflow-validation":
        result = render_workflow_validation_report(
            out_path=args.out,
            fixtures_root=args.fixtures_root,
        )
        inputs = [] if args.fixtures_root is None else [args.fixtures_root]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    metrics={
                        "total_fixture_count": result.validation.total_fixture_count,
                        "passed_fixture_count": result.validation.passed_fixture_count,
                        "workflow_count": len(result.validation.workflows),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "release-gate":
        result = render_level_one_release_gate_report(
            out_path=args.out,
            fixtures_root=args.fixtures_root,
        )
        inputs = [] if args.fixtures_root is None else [args.fixtures_root]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=result.release_gate.dataset_warnings,
                    metrics={
                        "decision": result.release_gate.gate.decision,
                        "retained_taxa": len(result.release_gate.gate.retained_taxa),
                        "excluded_taxa": len(result.release_gate.gate.excluded_taxa),
                        "blocked_analysis_count": len(
                            result.release_gate.gate.blocked_analyses
                        ),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "release-truth":
        result = render_release_truth_report(
            out_path=args.out,
            test_report_paths=args.test_report,
            real_engine_test_report_paths=args.real_engine_test_report,
            fixtures_root=args.fixtures_root,
            include_extended_parity=args.parity_extended,
            stress_tier=args.stress_tier,
        )
        inputs = [
            *args.test_report,
            *args.real_engine_test_report,
            *([args.fixtures_root] if args.fixtures_root is not None else []),
        ]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=result.release_truth.known_limitations,
                    metrics={
                        "total_tests": result.release_truth.total_tests.total_tests,
                        "total_tests_passed": result.release_truth.total_tests.passed_tests,
                        "total_tests_failed": result.release_truth.total_tests.failed_tests,
                        "total_tests_skipped": result.release_truth.total_tests.skipped_tests,
                        "real_engine_tests": result.release_truth.real_engine_tests.total_tests,
                        "real_engine_tests_passed": result.release_truth.real_engine_tests.passed_tests,
                        "real_engine_tests_failed": result.release_truth.real_engine_tests.failed_tests,
                        "real_engine_tests_skipped": result.release_truth.real_engine_tests.skipped_tests,
                        "supported_workflow_count": len(
                            result.release_truth.supported_workflows
                        ),
                        "experimental_workflow_count": len(
                            result.release_truth.experimental_workflows
                        ),
                        "flagship_dataset_count": len(
                            result.release_truth.flagship_datasets
                        ),
                        "reference_parity_case_count": result.release_truth.reference_parity.case_count,
                        "stress_workload_count": len(
                            result.release_truth.stress_suite.observations
                        ),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "production-scale-readiness":
        result = render_production_scale_readiness_report(
            out_path=args.out,
            replicates=args.replicates,
            tree_tip_counts=args.tree_tip_counts,
            alignment_size_classes=_build_production_scale_alignment_classes(args),
            tree_set_size_classes=_build_production_scale_tree_set_classes(args),
            stress_tiers=args.stress_tiers,
        )
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[],
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[],
                    outputs=outputs,
                    warnings=result.production_scale_readiness.limitations,
                    metrics=result.machine_manifest["metrics"],
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    raise NotImplementedError(f"unsupported report command: {args.report_command}")
