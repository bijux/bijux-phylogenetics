from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.ancestral.sensitivity import (
    build_ancestral_sensitivity_report,
)
from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _parse_assignment_map,
    _split_csv_values,
    _validate_ancestral_discrete_model_arguments,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_sensitivity_ancestral_commands(ancestral_subparsers: Any) -> None:
    ancestral_sensitivity = ancestral_subparsers.add_parser(
        "sensitivity",
        help="Summarize how ancestral results change across model, tree, pruning, or coding choices.",
    )
    ancestral_sensitivity.add_argument("tree", type=Path)
    ancestral_sensitivity.add_argument("table", type=Path)
    ancestral_sensitivity.add_argument("--trait", required=True)
    ancestral_sensitivity.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_sensitivity.add_argument("--taxon-column")
    ancestral_sensitivity.add_argument("--model")
    ancestral_sensitivity.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_sensitivity.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_sensitivity.add_argument("--alpha", type=float, default=1.0)
    ancestral_sensitivity.add_argument("--compare-model")
    ancestral_sensitivity.add_argument("--compare-tree", type=Path)
    ancestral_sensitivity.add_argument("--drop-taxa", nargs="+")
    ancestral_sensitivity.add_argument(
        "--coding-map",
        help="Comma-delimited KEY=VALUE recoding map for discrete traits.",
    )
    ancestral_sensitivity.add_argument(
        "--json", action="store_true", help="Emit the sensitivity report as JSON."
    )
    _add_manifest_argument(ancestral_sensitivity)


def run_sensitivity_ancestral_command(args: Any, *, parser: Any) -> int | None:
    if args.ancestral_command != "sensitivity":
        return None

    _validate_ancestral_discrete_model_arguments(args, parser)
    resolved_model = args.model or (
        "brownian" if args.kind == "continuous" else "fitch"
    )
    report = build_ancestral_sensitivity_report(
        tree_path=args.tree,
        traits_path=args.table,
        trait=args.trait,
        reconstruction_kind=args.kind,
        model=resolved_model,
        taxon_column=args.taxon_column,
        alpha=args.alpha,
        state_ordering=args.state_ordering,
        ordered_states=_split_csv_values(args.ordered_states) or None,
        compare_tree_path=args.compare_tree,
        compare_model=args.compare_model,
        drop_taxa=args.drop_taxa,
        coding_map=_parse_assignment_map(args.coding_map) or None,
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
                "baseline_node_count": report.baseline_node_count,
                "has_model_sensitivity": report.model_sensitivity is not None,
                "has_tree_sensitivity": report.tree_sensitivity is not None,
                "has_pruning_sensitivity": report.pruning_sensitivity is not None,
                "has_trait_coding_sensitivity": (
                    report.trait_coding_sensitivity is not None
                ),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
