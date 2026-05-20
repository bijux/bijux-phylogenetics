from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports import (
    write_supplementary_ancestral_state_table,
    write_supplementary_batch_summary_table,
    write_supplementary_comparative_model_table,
    write_supplementary_diversification_table,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .model_selection import (
    add_model_selection_supplementary_table_commands,
    run_model_selection_supplementary_table_command,
)
from .study_inputs import (
    add_study_input_supplementary_table_commands,
    run_study_input_supplementary_table_command,
)
from .tree_review import (
    add_tree_review_supplementary_table_commands,
    run_tree_review_supplementary_table_command,
)


def _parse_lambda_value(value: str | None) -> float | str:
    if value == "estimate" or value is None:
        return "estimate"
    return float(value)


def add_supplementary_table_report_commands(report_subparsers: Any) -> None:
    add_study_input_supplementary_table_commands(report_subparsers)
    add_tree_review_supplementary_table_commands(report_subparsers)
    add_model_selection_supplementary_table_commands(report_subparsers)

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


def run_supplementary_table_report_command(args: Any) -> int | None:
    study_input_result = run_study_input_supplementary_table_command(args)
    if study_input_result is not None:
        return study_input_result

    tree_review_result = run_tree_review_supplementary_table_command(args)
    if tree_review_result is not None:
        return tree_review_result

    model_selection_result = run_model_selection_supplementary_table_command(args)
    if model_selection_result is not None:
        return model_selection_result

    if args.report_command == "supplementary-comparative-model-table":
        result = write_supplementary_comparative_model_table(
            args.out,
            tree_path=args.tree,
            traits_path=args.traits,
            formulas=list(args.formulas),
            taxon_column=args.taxon_column,
            lambda_value=_parse_lambda_value(getattr(args, "lambda_value", None)),
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

    return None
