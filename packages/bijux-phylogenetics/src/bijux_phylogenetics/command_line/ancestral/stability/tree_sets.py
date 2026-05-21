from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.ancestral.tree_set import (
    summarize_continuous_ancestral_tree_set,
    summarize_continuous_ancestral_tree_set_report,
    summarize_discrete_ancestral_tree_set,
    summarize_discrete_ancestral_tree_set_report,
    write_ancestral_tree_set_exclusion_table,
    write_ancestral_tree_set_tree_table,
    write_continuous_ancestral_tree_set_clade_table,
    write_continuous_ancestral_tree_set_node_table,
    write_continuous_ancestral_tree_set_summary_table,
    write_discrete_ancestral_tree_set_clade_table,
    write_discrete_ancestral_tree_set_node_table,
    write_discrete_ancestral_tree_set_summary_table,
)
from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _split_csv_values,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_tree_set_stability_commands(ancestral_subparsers: Any) -> None:
    ancestral_tree_set = ancestral_subparsers.add_parser(
        "tree-set",
        help="Summarize ancestral reconstruction stability across a tree set.",
    )
    ancestral_tree_set.add_argument("tree_set", type=Path)
    ancestral_tree_set.add_argument("table", type=Path)
    ancestral_tree_set.add_argument("--trait", required=True)
    ancestral_tree_set.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_tree_set.add_argument("--taxon-column")
    ancestral_tree_set.add_argument(
        "--model",
        choices=(
            "brownian",
            "ou",
            "fitch",
            "equal-rates",
            "symmetric",
            "all-rates-different",
        ),
    )
    ancestral_tree_set.add_argument("--alpha", type=float, default=1.0)
    ancestral_tree_set.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_tree_set.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_tree_set.add_argument("--burnin-fraction", type=float, default=0.0)
    ancestral_tree_set.add_argument("--summary-out", type=Path)
    ancestral_tree_set.add_argument("--trees-out", type=Path)
    ancestral_tree_set.add_argument("--nodes-out", type=Path)
    ancestral_tree_set.add_argument("--clades-out", type=Path)
    ancestral_tree_set.add_argument("--exclusions-out", type=Path)
    ancestral_tree_set.add_argument(
        "--json", action="store_true", help="Emit the tree-set summary as JSON."
    )
    _add_manifest_argument(ancestral_tree_set)


def run_tree_set_stability_command(args: Any, *, parser: Any) -> int | None:
    if args.ancestral_command != "tree-set":
        return None

    if args.kind == "continuous":
        resolved_model = args.model or "brownian"
        if resolved_model not in {"brownian", "ou"}:
            parser.error(
                "continuous ancestral tree-set reconstruction requires model brownian or ou"
            )
        report = summarize_continuous_ancestral_tree_set(
            args.tree_set,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
            model=resolved_model,
            alpha=args.alpha,
            burnin_fraction=args.burnin_fraction,
        )
        summary = summarize_continuous_ancestral_tree_set_report(report)
        outputs: list[Path | str] = []
        if args.summary_out is not None:
            outputs.append(
                write_continuous_ancestral_tree_set_summary_table(
                    args.summary_out,
                    report,
                )
            )
        if args.trees_out is not None:
            outputs.append(write_ancestral_tree_set_tree_table(args.trees_out, report))
        if args.nodes_out is not None:
            outputs.append(
                write_continuous_ancestral_tree_set_node_table(
                    args.nodes_out,
                    report,
                )
            )
        if args.clades_out is not None:
            outputs.append(
                write_continuous_ancestral_tree_set_clade_table(
                    args.clades_out,
                    report,
                )
            )
        if args.exclusions_out is not None:
            outputs.append(
                write_ancestral_tree_set_exclusion_table(
                    args.exclusions_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="ancestral",
            inputs=[args.tree_set, args.table],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="ancestral",
                inputs=[args.tree_set, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "kind": "continuous",
                    "model": report.model,
                    "total_tree_count": report.total_tree_count,
                    "kept_tree_count": report.kept_tree_count,
                    "rooted_topology_count": report.rooted_topology_count,
                    "unrooted_topology_count": report.unrooted_topology_count,
                    "node_row_count": len(report.node_rows),
                    "clade_summary_count": len(report.clade_summaries),
                    "excluded_taxon_count": len(report.exclusions),
                    "unstable_clade_count": summary.unstable_clade_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    resolved_model = args.model or "fitch"
    if resolved_model not in {
        "fitch",
        "equal-rates",
        "symmetric",
        "all-rates-different",
    }:
        parser.error(
            "discrete ancestral tree-set reconstruction requires a discrete model"
        )
    if args.state_ordering == "ordered" and resolved_model == "fitch":
        parser.error(
            "ordered ancestral tree-set discrete reconstruction requires a likelihood model"
        )
    report = summarize_discrete_ancestral_tree_set(
        args.tree_set,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        model=resolved_model,
        state_ordering=args.state_ordering,
        ordered_states=_split_csv_values(args.ordered_states) or None,
        burnin_fraction=args.burnin_fraction,
    )
    summary = summarize_discrete_ancestral_tree_set_report(report)
    outputs: list[Path | str] = []
    if args.summary_out is not None:
        outputs.append(
            write_discrete_ancestral_tree_set_summary_table(
                args.summary_out,
                report,
            )
        )
    if args.trees_out is not None:
        outputs.append(write_ancestral_tree_set_tree_table(args.trees_out, report))
    if args.nodes_out is not None:
        outputs.append(
            write_discrete_ancestral_tree_set_node_table(
                args.nodes_out,
                report,
            )
        )
    if args.clades_out is not None:
        outputs.append(
            write_discrete_ancestral_tree_set_clade_table(
                args.clades_out,
                report,
            )
        )
    if args.exclusions_out is not None:
        outputs.append(
            write_ancestral_tree_set_exclusion_table(
                args.exclusions_out,
                report,
            )
        )
    outputs = _finalize_outputs(
        args,
        command="ancestral",
        inputs=[args.tree_set, args.table],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="ancestral",
            inputs=[args.tree_set, args.table],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "kind": "discrete",
                "model": report.model,
                "total_tree_count": report.total_tree_count,
                "kept_tree_count": report.kept_tree_count,
                "rooted_topology_count": report.rooted_topology_count,
                "unrooted_topology_count": report.unrooted_topology_count,
                "node_row_count": len(report.node_rows),
                "clade_summary_count": len(report.clade_summaries),
                "excluded_taxon_count": len(report.exclusions),
                "unstable_clade_count": summary.unstable_clade_count,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
