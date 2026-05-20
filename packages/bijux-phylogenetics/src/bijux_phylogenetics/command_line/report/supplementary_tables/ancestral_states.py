from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports import write_supplementary_ancestral_state_table
from bijux_phylogenetics.runtime.results import build_command_result


def add_ancestral_state_supplementary_table_commands(
    report_subparsers: Any,
) -> None:
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


def run_ancestral_state_supplementary_table_command(args: Any) -> int | None:
    if args.report_command != "supplementary-ancestral-state-table":
        return None

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
