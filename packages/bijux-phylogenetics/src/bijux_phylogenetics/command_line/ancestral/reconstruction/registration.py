from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.ancestral.discrete.review import (
    validate_discrete_ancestral_reference_examples,
)
from bijux_phylogenetics.ancestral.comparison import (
    compare_continuous_ancestral_models,
)
from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result

from .continuous import (
    add_continuous_reconstruction_commands,
    run_continuous_reconstruction_command,
)
from .discrete import (
    add_discrete_reconstruction_commands,
    run_discrete_reconstruction_command,
)


def add_reconstruction_ancestral_commands(ancestral_subparsers: Any) -> None:
    add_continuous_reconstruction_commands(ancestral_subparsers)
    add_discrete_reconstruction_commands(ancestral_subparsers)

    ancestral_discrete_reference = ancestral_subparsers.add_parser(
        "discrete-reference",
        help="Validate built-in discrete ancestral reference examples.",
    )
    ancestral_discrete_reference.add_argument(
        "--json",
        action="store_true",
        help="Emit the reference validation report as JSON.",
    )

    ancestral_compare = ancestral_subparsers.add_parser(
        "compare",
        help="Compare two continuous ancestral-state models node by node.",
    )
    ancestral_compare.add_argument("tree", type=Path)
    ancestral_compare.add_argument("table", type=Path)
    ancestral_compare.add_argument("--trait", required=True)
    ancestral_compare.add_argument("--taxon-column")
    ancestral_compare.add_argument(
        "--left-model", choices=("brownian", "ou"), default="brownian"
    )
    ancestral_compare.add_argument(
        "--right-model", choices=("brownian", "ou"), default="ou"
    )
    ancestral_compare.add_argument("--left-alpha", type=float, default=1.0)
    ancestral_compare.add_argument("--right-alpha", type=float, default=1.0)
    ancestral_compare.add_argument(
        "--json", action="store_true", help="Emit the comparison as JSON."
    )
    _add_manifest_argument(ancestral_compare)


def run_reconstruction_ancestral_command(args: Any, *, parser: Any) -> int | None:
    continuous_result = run_continuous_reconstruction_command(args, parser=parser)
    if continuous_result is not None:
        return continuous_result
    discrete_result = run_discrete_reconstruction_command(args, parser=parser)
    if discrete_result is not None:
        return discrete_result

    if args.ancestral_command == "discrete-reference":
        report = validate_discrete_ancestral_reference_examples()
        outputs = _finalize_outputs(args, command="ancestral", inputs=[])
        _print_result(
            build_command_result(
                command="ancestral",
                inputs=[],
                outputs=outputs,
                metrics={
                    "case_count": report.case_count,
                    "external_case_count": report.external_case_count,
                    "all_passed": report.all_passed,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.ancestral_command == "compare":
        report = compare_continuous_ancestral_models(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
            left_model=args.left_model,
            right_model=args.right_model,
            left_alpha=args.left_alpha,
            right_alpha=args.right_alpha,
        )
        outputs = _finalize_outputs(
            args, command="ancestral", inputs=[args.tree, args.table]
        )
        _print_result(
            build_command_result(
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
                metrics={
                    "taxon_count": report.taxon_count,
                    "compared_node_count": len(report.rows),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None
