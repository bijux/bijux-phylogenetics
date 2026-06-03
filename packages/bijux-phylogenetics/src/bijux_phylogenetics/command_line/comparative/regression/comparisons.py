from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.reporting import (
    compare_comparative_results_across_pruning,
    compare_comparative_results_across_trees,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_regression_comparison_commands(comparative_subparsers: Any) -> None:
    comparative_compare_trees = comparative_subparsers.add_parser(
        "compare-trees",
        help="Compare comparative results across two alternative trees.",
    )
    comparative_compare_trees.add_argument("left_tree", type=Path)
    comparative_compare_trees.add_argument("right_tree", type=Path)
    comparative_compare_trees.add_argument("table", type=Path)
    comparative_compare_trees.add_argument("--response")
    comparative_compare_trees.add_argument("--predictors", nargs="+")
    comparative_compare_trees.add_argument("--formula")
    comparative_compare_trees.add_argument("--taxon-column")
    comparative_compare_trees.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_compare_trees.add_argument(
        "--json", action="store_true", help="Emit the comparison as JSON."
    )
    _add_manifest_argument(comparative_compare_trees)

    comparative_compare_pruning = comparative_subparsers.add_parser(
        "compare-pruning",
        help="Compare comparative results before and after explicit pruning.",
    )
    comparative_compare_pruning.add_argument("tree", type=Path)
    comparative_compare_pruning.add_argument("table", type=Path)
    comparative_compare_pruning.add_argument("--response")
    comparative_compare_pruning.add_argument("--predictors", nargs="+")
    comparative_compare_pruning.add_argument("--formula")
    comparative_compare_pruning.add_argument("--drop-taxa", nargs="+")
    comparative_compare_pruning.add_argument("--keep-taxa", nargs="+")
    comparative_compare_pruning.add_argument("--taxon-column")
    comparative_compare_pruning.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_compare_pruning.add_argument(
        "--json", action="store_true", help="Emit the pruning comparison as JSON."
    )
    _add_manifest_argument(comparative_compare_pruning)


def run_regression_comparison_command(
    args: Any,
    *,
    lambda_value: float | str,
) -> int | None:
    if args.comparative_command == "compare-trees":
        report = compare_comparative_results_across_trees(
            args.left_tree,
            args.right_tree,
            args.table,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            taxon_column=args.taxon_column,
            lambda_value=lambda_value,
        )
        outputs = _finalize_outputs(
            args,
            command="comparative",
            inputs=[args.left_tree, args.right_tree, args.table],
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.left_tree, args.right_tree, args.table],
                outputs=outputs,
                metrics={
                    "coefficient_delta_count": len(report.coefficient_deltas),
                    "sign_changed_terms": len(report.sign_changed_terms),
                    "conclusion_changed": report.conclusion_changed,
                    "left_selected_model": report.left_selected_model,
                    "right_selected_model": report.right_selected_model,
                },
                warnings=report.warnings,
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.comparative_command != "compare-pruning":
        return None

    report = compare_comparative_results_across_pruning(
        args.tree,
        args.table,
        response=args.response,
        predictors=list(args.predictors or []),
        formula=args.formula,
        drop_taxa=list(args.drop_taxa or []),
        keep_taxa=list(args.keep_taxa or []),
        taxon_column=args.taxon_column,
        lambda_value=lambda_value,
    )
    outputs = _finalize_outputs(
        args,
        command="comparative",
        inputs=[args.tree, args.table],
    )
    _print_result(
        build_command_result(
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=outputs,
            metrics={
                "baseline_taxa": len(report.baseline_taxa),
                "pruned_taxa": len(report.pruned_taxa),
                "dropped_taxa": len(report.dropped_taxa),
                "sign_changed_terms": len(report.sign_changed_terms),
                "conclusion_changed": report.conclusion_changed,
            },
            warnings=report.warnings,
            data=report,
        ),
        json_output=args.json,
    )
    return 0
