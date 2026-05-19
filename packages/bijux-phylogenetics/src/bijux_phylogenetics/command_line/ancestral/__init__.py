from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.ancestral.confidence import (
    build_continuous_ancestral_confidence_rows,
    build_continuous_ancestral_tree_set_confidence_rows,
    build_discrete_ancestral_confidence_rows,
    build_discrete_ancestral_tree_set_confidence_rows,
    summarize_continuous_ancestral_confidence,
    summarize_continuous_ancestral_tree_set_confidence,
    summarize_discrete_ancestral_confidence,
    summarize_discrete_ancestral_tree_set_confidence,
    write_ancestral_confidence_summary_table,
    write_continuous_ancestral_confidence_table,
    write_continuous_ancestral_tree_set_confidence_table,
    write_discrete_ancestral_confidence_table,
    write_discrete_ancestral_tree_set_confidence_table,
)
from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete import (
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.command_line.ancestral.reconstruction import (
    add_reconstruction_ancestral_commands,
    run_reconstruction_ancestral_command,
)
from bijux_phylogenetics.ancestral.irreversible_discrete import (
    summarize_irreversible_discrete_reconstruction,
    summarize_irreversible_discrete_report,
    write_irreversible_discrete_fit_table,
    write_irreversible_discrete_node_table,
    write_irreversible_discrete_summary_table,
    write_irreversible_discrete_transition_table,
)
from bijux_phylogenetics.ancestral.ordered_discrete import (
    summarize_ordered_discrete_reconstruction,
    summarize_ordered_discrete_report,
    write_ordered_discrete_fit_table,
    write_ordered_discrete_node_table,
    write_ordered_discrete_summary_table,
    write_ordered_discrete_transition_table,
)
from bijux_phylogenetics.ancestral.package import build_ancestral_figure_package
from bijux_phylogenetics.ancestral.report_package import (
    build_ancestral_report_package,
)
from bijux_phylogenetics.ancestral.root_sensitivity import (
    summarize_ancestral_root_sensitivity,
    summarize_ancestral_root_sensitivity_report,
    write_ancestral_root_assumption_table,
    write_ancestral_root_sensitivity_node_table,
    write_ancestral_root_sensitivity_summary_table,
)
from bijux_phylogenetics.ancestral.sensitivity import build_ancestral_sensitivity_report
from bijux_phylogenetics.ancestral.service import (
    render_ancestral_state_report,
)
from bijux_phylogenetics.ancestral.transitions import (
    summarize_ancestral_transition_report,
    summarize_ancestral_transition_tree_set,
    summarize_ancestral_transition_tree_set_report,
    summarize_ancestral_transitions,
    write_ancestral_transition_branch_table,
    write_ancestral_transition_count_table,
    write_ancestral_transition_exclusion_table,
    write_ancestral_transition_summary_table,
    write_ancestral_transition_tree_set_branch_table,
    write_ancestral_transition_tree_set_count_table,
    write_ancestral_transition_tree_set_summary_table,
    write_ancestral_transition_tree_set_tree_table,
)
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
from bijux_phylogenetics.ancestral.visualization import (
    render_ancestral_state_visualization,
)
from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _parse_assignment_map,
    _parse_transition_pairs,
    _split_csv_values,
    _validate_ancestral_discrete_model_arguments,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_ancestral_commands(subparsers: Any) -> None:
    ancestral = subparsers.add_parser(
        get_command_spec("ancestral").name,
        help=get_command_spec("ancestral").summary,
    )
    ancestral_subparsers = ancestral.add_subparsers(
        dest="ancestral_command", required=True
    )
    add_reconstruction_ancestral_commands(ancestral_subparsers)
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
    ancestral_confidence = ancestral_subparsers.add_parser(
        "confidence",
        help="Summarize ancestral state confidence on one tree or tree set.",
    )
    ancestral_confidence.add_argument("tree", type=Path)
    ancestral_confidence.add_argument("table", type=Path)
    ancestral_confidence.add_argument("--trait", required=True)
    ancestral_confidence.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_confidence.add_argument("--taxon-column")
    ancestral_confidence.add_argument(
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
    ancestral_confidence.add_argument("--tree-set", action="store_true")
    ancestral_confidence.add_argument("--alpha", type=float, default=1.0)
    ancestral_confidence.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_confidence.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_confidence.add_argument("--burnin-fraction", type=float, default=0.0)
    ancestral_confidence.add_argument("--summary-out", type=Path)
    ancestral_confidence.add_argument("--confidence-out", type=Path)
    ancestral_confidence.add_argument(
        "--json", action="store_true", help="Emit the confidence review as JSON."
    )
    _add_manifest_argument(ancestral_confidence)
    ancestral_root_sensitivity = ancestral_subparsers.add_parser(
        "root-sensitivity",
        help="Summarize how discrete likelihood ancestral reconstructions change under explicit root assumptions.",
    )
    ancestral_root_sensitivity.add_argument("tree", type=Path)
    ancestral_root_sensitivity.add_argument("table", type=Path)
    ancestral_root_sensitivity.add_argument("--trait", required=True)
    ancestral_root_sensitivity.add_argument("--taxon-column")
    ancestral_root_sensitivity.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    ancestral_root_sensitivity.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_root_sensitivity.add_argument(
        "--ordered-states",
        help="Comma-delimited explicit ordered state vocabulary.",
    )
    ancestral_root_sensitivity.add_argument("--fixed-root-state")
    ancestral_root_sensitivity.add_argument("--summary-out", type=Path)
    ancestral_root_sensitivity.add_argument("--assumptions-out", type=Path)
    ancestral_root_sensitivity.add_argument("--nodes-out", type=Path)
    ancestral_root_sensitivity.add_argument(
        "--json", action="store_true", help="Emit the root-sensitivity review as JSON."
    )
    _add_manifest_argument(ancestral_root_sensitivity)
    ancestral_ordered_discrete = ancestral_subparsers.add_parser(
        "ordered-discrete",
        help="Compare ordered and unordered discrete likelihood ancestral reconstructions.",
    )
    ancestral_ordered_discrete.add_argument("tree", type=Path)
    ancestral_ordered_discrete.add_argument("table", type=Path)
    ancestral_ordered_discrete.add_argument("--trait", required=True)
    ancestral_ordered_discrete.add_argument("--taxon-column")
    ancestral_ordered_discrete.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    ancestral_ordered_discrete.add_argument(
        "--ordered-states",
        required=True,
        help="Comma-delimited explicit ordered state vocabulary.",
    )
    ancestral_ordered_discrete.add_argument("--summary-out", type=Path)
    ancestral_ordered_discrete.add_argument("--fits-out", type=Path)
    ancestral_ordered_discrete.add_argument("--nodes-out", type=Path)
    ancestral_ordered_discrete.add_argument("--transitions-out", type=Path)
    ancestral_ordered_discrete.add_argument(
        "--json", action="store_true", help="Emit the ordered discrete review as JSON."
    )
    _add_manifest_argument(ancestral_ordered_discrete)
    ancestral_irreversible_discrete = ancestral_subparsers.add_parser(
        "irreversible-discrete",
        help="Compare constrained and unconstrained discrete ancestral likelihood reconstructions under an allowed transition graph.",
    )
    ancestral_irreversible_discrete.add_argument("tree", type=Path)
    ancestral_irreversible_discrete.add_argument("table", type=Path)
    ancestral_irreversible_discrete.add_argument("--trait", required=True)
    ancestral_irreversible_discrete.add_argument("--taxon-column")
    ancestral_irreversible_discrete.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="all-rates-different",
    )
    ancestral_irreversible_discrete.add_argument(
        "--allowed-transitions",
        required=True,
        help="Comma-delimited directed transition graph in SOURCE->TARGET form.",
    )
    ancestral_irreversible_discrete.add_argument("--summary-out", type=Path)
    ancestral_irreversible_discrete.add_argument("--fits-out", type=Path)
    ancestral_irreversible_discrete.add_argument("--nodes-out", type=Path)
    ancestral_irreversible_discrete.add_argument("--transitions-out", type=Path)
    ancestral_irreversible_discrete.add_argument(
        "--json",
        action="store_true",
        help="Emit the irreversible discrete review as JSON.",
    )
    _add_manifest_argument(ancestral_irreversible_discrete)
    ancestral_transitions = ancestral_subparsers.add_parser(
        "transitions",
        help="Count inferred ancestral transitions on one tree or tree set.",
    )
    ancestral_transitions.add_argument("tree", type=Path)
    ancestral_transitions.add_argument("table", type=Path)
    ancestral_transitions.add_argument("--trait", required=True)
    ancestral_transitions.add_argument("--taxon-column")
    ancestral_transitions.add_argument(
        "--model",
        choices=("fitch", "equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="fitch",
    )
    ancestral_transitions.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_transitions.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_transitions.add_argument("--tree-set", action="store_true")
    ancestral_transitions.add_argument("--burnin-fraction", type=float, default=0.0)
    ancestral_transitions.add_argument("--summary-out", type=Path)
    ancestral_transitions.add_argument("--trees-out", type=Path)
    ancestral_transitions.add_argument("--branches-out", type=Path)
    ancestral_transitions.add_argument("--counts-out", type=Path)
    ancestral_transitions.add_argument("--exclusions-out", type=Path)
    ancestral_transitions.add_argument(
        "--json", action="store_true", help="Emit the transition report as JSON."
    )
    _add_manifest_argument(ancestral_transitions)
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
    ancestral_render = ancestral_subparsers.add_parser(
        "render",
        help="Render a tree annotated with reconstructed ancestral states.",
    )
    ancestral_render.add_argument("tree", type=Path)
    ancestral_render.add_argument("table", type=Path)
    ancestral_render.add_argument("--trait", required=True)
    ancestral_render.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_render.add_argument("--taxon-column")
    ancestral_render.add_argument("--model")
    ancestral_render.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_render.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_render.add_argument("--alpha", type=float, default=1.0)
    ancestral_render.add_argument(
        "--layout", choices=("cladogram", "phylogram", "circular"), default="phylogram"
    )
    ancestral_render.add_argument(
        "--discrete-node-style", choices=("labels", "pies"), default="labels"
    )
    ancestral_render.add_argument(
        "--branch-coloring", choices=("none", "state", "regime"), default="none"
    )
    ancestral_render.add_argument("--out", required=True, type=Path)
    ancestral_render.add_argument(
        "--json", action="store_true", help="Emit the render result as JSON."
    )
    _add_manifest_argument(ancestral_render)
    ancestral_report = ancestral_subparsers.add_parser(
        "report",
        help="Render an HTML report for ancestral-state reconstruction.",
    )
    ancestral_report.add_argument("tree", type=Path)
    ancestral_report.add_argument("table", type=Path)
    ancestral_report.add_argument("--trait", required=True)
    ancestral_report.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_report.add_argument("--taxon-column")
    ancestral_report.add_argument("--model")
    ancestral_report.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_report.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_report.add_argument("--alpha", type=float, default=1.0)
    ancestral_report.add_argument("--compare-model")
    ancestral_report.add_argument("--compare-tree", type=Path)
    ancestral_report.add_argument("--drop-taxa", nargs="+")
    ancestral_report.add_argument(
        "--coding-map",
        help="Comma-delimited KEY=VALUE recoding map for discrete traits.",
    )
    ancestral_report.add_argument("--out", type=Path)
    ancestral_report.add_argument(
        "--out-dir",
        type=Path,
        help="Write a full ancestral reconstruction report package directory.",
    )
    ancestral_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(ancestral_report)
    ancestral_package = ancestral_subparsers.add_parser(
        "package",
        help="Write a publication-ready ancestral-state figure package.",
    )
    ancestral_package.add_argument("tree", type=Path)
    ancestral_package.add_argument("table", type=Path)
    ancestral_package.add_argument("--trait", required=True)
    ancestral_package.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_package.add_argument("--taxon-column")
    ancestral_package.add_argument("--model")
    ancestral_package.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_package.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_package.add_argument("--alpha", type=float, default=1.0)
    ancestral_package.add_argument(
        "--layout", choices=("cladogram", "phylogram", "circular"), default="phylogram"
    )
    ancestral_package.add_argument("--out-dir", required=True, type=Path)
    ancestral_package.add_argument(
        "--json", action="store_true", help="Emit the package build result as JSON."
    )
    _add_manifest_argument(ancestral_package)


def run_ancestral_command(args: Any, *, parser: Any) -> int:
    reconstruction_exit_code = run_reconstruction_ancestral_command(
        args,
        parser=parser,
    )
    if reconstruction_exit_code is not None:
        return reconstruction_exit_code
    if args.ancestral_command == "tree-set":
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
                outputs.append(
                    write_ancestral_tree_set_tree_table(args.trees_out, report)
                )
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
        outputs = []
        if args.summary_out is not None:
            outputs.append(
                write_discrete_ancestral_tree_set_summary_table(
                    args.summary_out,
                    report,
                )
            )
        if args.trees_out is not None:
            outputs.append(
                write_ancestral_tree_set_tree_table(args.trees_out, report)
            )
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
    if args.ancestral_command == "confidence":
        if not args.tree_set and args.burnin_fraction != 0.0:
            parser.error("--burnin-fraction requires --tree-set")
        if args.kind == "continuous":
            resolved_model = args.model or "brownian"
            if resolved_model not in {"brownian", "ou"}:
                parser.error(
                    "continuous ancestral confidence review requires model brownian or ou"
                )
            if args.tree_set:
                report = summarize_continuous_ancestral_tree_set(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=resolved_model,
                    alpha=args.alpha,
                    burnin_fraction=args.burnin_fraction,
                )
                confidence_rows = build_continuous_ancestral_tree_set_confidence_rows(
                    report
                )
                confidence_summary = summarize_continuous_ancestral_tree_set_confidence(
                    report
                )
                outputs: list[Path | str] = []
                if args.summary_out is not None:
                    outputs.append(
                        write_ancestral_confidence_summary_table(
                            args.summary_out,
                            confidence_summary,
                        )
                    )
                if args.confidence_out is not None:
                    outputs.append(
                        write_continuous_ancestral_tree_set_confidence_table(
                            args.confidence_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "kind": "continuous",
                            "source_kind": "tree_set",
                            "model": report.model,
                            "kept_tree_count": report.kept_tree_count,
                            "confidence_row_count": len(confidence_rows),
                            "low_confidence_count": (
                                confidence_summary.low_confidence_count
                            ),
                            "unstable_count": confidence_summary.unstable_count,
                            "high_entropy_count": confidence_summary.high_entropy_count,
                            "top_uncertain_id": confidence_summary.top_uncertain_id,
                        },
                        data={
                            "report": report,
                            "confidence_summary": confidence_summary,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            report = reconstruct_continuous_ancestral_states(
                args.tree,
                args.table,
                trait=args.trait,
                taxon_column=args.taxon_column,
                model=resolved_model,
                alpha=args.alpha,
            )
            confidence_rows = build_continuous_ancestral_confidence_rows(report)
            confidence_summary = summarize_continuous_ancestral_confidence(report)
            outputs = []
            if args.summary_out is not None:
                outputs.append(
                    write_ancestral_confidence_summary_table(
                        args.summary_out,
                        confidence_summary,
                    )
                )
            if args.confidence_out is not None:
                outputs.append(
                    write_continuous_ancestral_confidence_table(
                        args.confidence_out,
                        report,
                    )
                )
            outputs = _finalize_outputs(
                args,
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
            )
            _print_result(
                build_command_result(
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                    warnings=report.warnings,
                    metrics={
                        "kind": "continuous",
                        "source_kind": "tree",
                        "model": report.model,
                        "confidence_row_count": len(confidence_rows),
                        "low_confidence_count": confidence_summary.low_confidence_count,
                        "unstable_count": confidence_summary.unstable_count,
                        "high_entropy_count": confidence_summary.high_entropy_count,
                        "top_uncertain_id": confidence_summary.top_uncertain_id,
                    },
                    data={
                        "report": report,
                        "confidence_summary": confidence_summary,
                    },
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
                "discrete ancestral confidence review requires a discrete model"
            )
        if args.state_ordering == "ordered" and resolved_model == "fitch":
            parser.error(
                "ordered ancestral confidence review requires a likelihood model"
            )
        if args.tree_set:
            report = summarize_discrete_ancestral_tree_set(
                args.tree,
                args.table,
                trait=args.trait,
                taxon_column=args.taxon_column,
                model=resolved_model,
                state_ordering=args.state_ordering,
                ordered_states=_split_csv_values(args.ordered_states) or None,
                burnin_fraction=args.burnin_fraction,
            )
            confidence_rows = build_discrete_ancestral_tree_set_confidence_rows(report)
            confidence_summary = summarize_discrete_ancestral_tree_set_confidence(
                report
            )
            outputs = []
            if args.summary_out is not None:
                outputs.append(
                    write_ancestral_confidence_summary_table(
                        args.summary_out,
                        confidence_summary,
                    )
                )
            if args.confidence_out is not None:
                outputs.append(
                    write_discrete_ancestral_tree_set_confidence_table(
                        args.confidence_out,
                        report,
                    )
                )
            outputs = _finalize_outputs(
                args,
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
            )
            _print_result(
                build_command_result(
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                    warnings=report.warnings,
                    metrics={
                        "kind": "discrete",
                        "source_kind": "tree_set",
                        "model": report.model,
                        "kept_tree_count": report.kept_tree_count,
                        "confidence_row_count": len(confidence_rows),
                        "low_confidence_count": confidence_summary.low_confidence_count,
                        "unstable_count": confidence_summary.unstable_count,
                        "high_entropy_count": confidence_summary.high_entropy_count,
                        "top_uncertain_id": confidence_summary.top_uncertain_id,
                    },
                    data={
                        "report": report,
                        "confidence_summary": confidence_summary,
                    },
                ),
                json_output=args.json,
            )
            return 0
        report = reconstruct_discrete_ancestral_states(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
            model=resolved_model,
            state_ordering=args.state_ordering,
            ordered_states=_split_csv_values(args.ordered_states) or None,
        )
        confidence_rows = build_discrete_ancestral_confidence_rows(report)
        confidence_summary = summarize_discrete_ancestral_confidence(report)
        outputs = []
        if args.summary_out is not None:
            outputs.append(
                write_ancestral_confidence_summary_table(
                    args.summary_out,
                    confidence_summary,
                )
            )
        if args.confidence_out is not None:
            outputs.append(
                write_discrete_ancestral_confidence_table(
                    args.confidence_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="ancestral",
            inputs=[args.tree, args.table],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "kind": "discrete",
                    "source_kind": "tree",
                    "model": report.model,
                    "confidence_row_count": len(confidence_rows),
                    "low_confidence_count": confidence_summary.low_confidence_count,
                    "unstable_count": confidence_summary.unstable_count,
                    "high_entropy_count": confidence_summary.high_entropy_count,
                    "top_uncertain_id": confidence_summary.top_uncertain_id,
                },
                data={
                    "report": report,
                    "confidence_summary": confidence_summary,
                },
            ),
            json_output=args.json,
        )
        return 0
    if args.ancestral_command == "root-sensitivity":
        report = summarize_ancestral_root_sensitivity(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
            model=args.model,
            state_ordering=args.state_ordering,
            ordered_states=_split_csv_values(args.ordered_states) or None,
            fixed_root_state=args.fixed_root_state,
        )
        summary = summarize_ancestral_root_sensitivity_report(report)
        outputs = []
        if args.summary_out is not None:
            outputs.append(
                write_ancestral_root_sensitivity_summary_table(
                    args.summary_out,
                    report,
                )
            )
        if args.assumptions_out is not None:
            outputs.append(
                write_ancestral_root_assumption_table(
                    args.assumptions_out,
                    report,
                )
            )
        if args.nodes_out is not None:
            outputs.append(
                write_ancestral_root_sensitivity_node_table(
                    args.nodes_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="ancestral",
            inputs=[args.tree, args.table],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "model": report.model,
                    "state_ordering": report.state_ordering,
                    "analyzed_taxon_count": report.analyzed_taxon_count,
                    "assumption_count": summary.assumption_count,
                    "compared_node_count": summary.compared_node_count,
                    "state_changed_node_count": summary.state_changed_node_count,
                    "support_changed_node_count": summary.support_changed_node_count,
                    "top_sensitive_node": summary.top_sensitive_node,
                    "fixed_root_state": report.fixed_root_state,
                },
                data={
                    "report": report,
                    "summary": summary,
                },
            ),
            json_output=args.json,
        )
        return 0
    if args.ancestral_command == "ordered-discrete":
        report = summarize_ordered_discrete_reconstruction(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
            model=args.model,
            ordered_states=_split_csv_values(args.ordered_states) or [],
        )
        summary = summarize_ordered_discrete_report(report)
        outputs = []
        if args.summary_out is not None:
            outputs.append(
                write_ordered_discrete_summary_table(
                    args.summary_out,
                    report,
                )
            )
        if args.fits_out is not None:
            outputs.append(
                write_ordered_discrete_fit_table(
                    args.fits_out,
                    report,
                )
            )
        if args.nodes_out is not None:
            outputs.append(
                write_ordered_discrete_node_table(
                    args.nodes_out,
                    report,
                )
            )
        if args.transitions_out is not None:
            outputs.append(
                write_ordered_discrete_transition_table(
                    args.transitions_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="ancestral",
            inputs=[args.tree, args.table],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "model": report.model,
                    "ordered_state_count": len(report.ordered_states),
                    "fit_count": len(report.fit_rows),
                    "differing_node_count": summary.differing_node_count,
                    "ambiguity_change_count": summary.ambiguity_change_count,
                    "restricted_transition_count": (
                        summary.restricted_transition_count
                    ),
                    "preferred_ordering": summary.preferred_ordering,
                },
                data={
                    "report": report,
                    "summary": summary,
                },
            ),
            json_output=args.json,
        )
        return 0
    if args.ancestral_command == "irreversible-discrete":
        report = summarize_irreversible_discrete_reconstruction(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
            model=args.model,
            allowed_transition_pairs=_parse_transition_pairs(args.allowed_transitions),
        )
        summary = summarize_irreversible_discrete_report(report)
        outputs = []
        if args.summary_out is not None:
            outputs.append(
                write_irreversible_discrete_summary_table(
                    args.summary_out,
                    report,
                )
            )
        if args.fits_out is not None:
            outputs.append(
                write_irreversible_discrete_fit_table(
                    args.fits_out,
                    report,
                )
            )
        if args.nodes_out is not None:
            outputs.append(
                write_irreversible_discrete_node_table(
                    args.nodes_out,
                    report,
                )
            )
        if args.transitions_out is not None:
            outputs.append(
                write_irreversible_discrete_transition_table(
                    args.transitions_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="ancestral",
            inputs=[args.tree, args.table],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "model": report.model,
                    "allowed_transition_count": len(report.allowed_transition_pairs),
                    "fit_count": len(report.fit_rows),
                    "differing_node_count": summary.differing_node_count,
                    "ambiguity_change_count": summary.ambiguity_change_count,
                    "forbidden_transition_count": summary.forbidden_transition_count,
                    "preferred_constraint": summary.preferred_constraint,
                },
                data={
                    "report": report,
                    "summary": summary,
                },
            ),
            json_output=args.json,
        )
        return 0
    if args.ancestral_command == "transitions":
        if args.state_ordering == "ordered" and args.model == "fitch":
            parser.error(
                "ordered ancestral transition counting requires a likelihood model"
            )
        if not args.tree_set and args.trees_out is not None:
            parser.error("--trees-out requires --tree-set")
        if not args.tree_set and args.burnin_fraction != 0.0:
            parser.error("--burnin-fraction requires --tree-set")
        if args.tree_set:
            report = summarize_ancestral_transition_tree_set(
                args.tree,
                args.table,
                trait=args.trait,
                taxon_column=args.taxon_column,
                model=args.model,
                state_ordering=args.state_ordering,
                ordered_states=_split_csv_values(args.ordered_states) or None,
                burnin_fraction=args.burnin_fraction,
            )
            summary = summarize_ancestral_transition_tree_set_report(report)
            outputs = []
            if args.summary_out is not None:
                outputs.append(
                    write_ancestral_transition_tree_set_summary_table(
                        args.summary_out,
                        report,
                    )
                )
            if args.trees_out is not None:
                outputs.append(
                    write_ancestral_transition_tree_set_tree_table(
                        args.trees_out,
                        report,
                    )
                )
            if args.branches_out is not None:
                outputs.append(
                    write_ancestral_transition_tree_set_branch_table(
                        args.branches_out,
                        report,
                    )
                )
            if args.counts_out is not None:
                outputs.append(
                    write_ancestral_transition_tree_set_count_table(
                        args.counts_out,
                        report,
                    )
                )
            if args.exclusions_out is not None:
                outputs.append(
                    write_ancestral_transition_exclusion_table(
                        args.exclusions_out,
                        report,
                    )
                )
            outputs = _finalize_outputs(
                args,
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
            )
            _print_result(
                build_command_result(
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                    warnings=report.warnings,
                    metrics={
                        "tree_set": True,
                        "model": report.model,
                        "total_tree_count": report.total_tree_count,
                        "kept_tree_count": report.kept_tree_count,
                        "rooted_topology_count": report.rooted_topology_count,
                        "unrooted_topology_count": report.unrooted_topology_count,
                        "transition_pair_count": len(report.transition_rows),
                        "topology_sensitive_transition_pair_count": (
                            summary.topology_sensitive_transition_pair_count
                        ),
                        "uncertainty_sensitive_transition_pair_count": (
                            summary.uncertainty_sensitive_transition_pair_count
                        ),
                        "excluded_taxon_count": len(report.exclusions),
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        report = summarize_ancestral_transitions(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
            model=args.model,
            state_ordering=args.state_ordering,
            ordered_states=_split_csv_values(args.ordered_states) or None,
        )
        summary = summarize_ancestral_transition_report(report)
        outputs = []
        if args.summary_out is not None:
            outputs.append(
                write_ancestral_transition_summary_table(
                    args.summary_out,
                    report,
                )
            )
        if args.branches_out is not None:
            outputs.append(
                write_ancestral_transition_branch_table(
                    args.branches_out,
                    report,
                )
            )
        if args.counts_out is not None:
            outputs.append(
                write_ancestral_transition_count_table(
                    args.counts_out,
                    report,
                )
            )
        if args.exclusions_out is not None:
            outputs.append(
                write_ancestral_transition_exclusion_table(
                    args.exclusions_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="ancestral",
            inputs=[args.tree, args.table],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "tree_set": False,
                    "model": report.model,
                    "total_branch_count": summary.total_branch_count,
                    "changed_branch_count": summary.changed_branch_count,
                    "certain_change_count": summary.certain_change_count,
                    "uncertain_change_count": summary.uncertain_change_count,
                    "transition_pair_count": len(report.transition_rows),
                    "excluded_taxon_count": len(report.exclusions),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.ancestral_command == "sensitivity":
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
    if args.ancestral_command == "render":
        _validate_ancestral_discrete_model_arguments(args, parser)
        if args.kind == "continuous" and args.branch_coloring == "state":
            parser.error(
                "continuous ancestral rendering does not support branch coloring 'state'"
            )
        if args.kind == "discrete" and args.branch_coloring == "regime":
            parser.error(
                "discrete ancestral rendering does not support branch coloring 'regime'"
            )
        if args.kind == "continuous":
            resolved_model = args.model or "brownian"
            reconstruction = reconstruct_continuous_ancestral_states(
                args.tree,
                args.table,
                trait=args.trait,
                taxon_column=args.taxon_column,
                model=resolved_model,
                alpha=args.alpha,
            )
        else:
            resolved_model = args.model or "fitch"
            reconstruction = reconstruct_discrete_ancestral_states(
                args.tree,
                args.table,
                trait=args.trait,
                taxon_column=args.taxon_column,
                model=resolved_model,
                state_ordering=args.state_ordering,
                ordered_states=_split_csv_values(args.ordered_states) or None,
            )
        result = render_ancestral_state_visualization(
            args.tree,
            reconstruction,
            out_path=args.out,
            layout=args.layout,
            discrete_node_style=args.discrete_node_style,
            branch_coloring=args.branch_coloring,
        )
        rendered_outputs = (
            [result.output_path]
            if result.format == "svg"
            else [result.output_path, result.svg_path]
        )
        outputs = _finalize_outputs(
            args,
            command="ancestral",
            inputs=[args.tree, args.table],
            outputs=rendered_outputs,
        )
        _print_result(
            build_command_result(
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=getattr(reconstruction, "warnings", []),
                metrics={
                    "tip_count": result.tree_render.tip_count,
                    "format": result.format,
                    "layout": result.layout,
                    "rendered_internal_annotation_count": (
                        result.tree_render.rendered_internal_annotation_count
                    ),
                    "rendered_internal_pie_count": (
                        result.tree_render.rendered_internal_pie_count
                    ),
                    "rendered_branch_color_count": (
                        result.tree_render.rendered_branch_color_count
                    ),
                },
                data={
                    "reconstruction": reconstruction,
                    "visualization": result,
                },
            ),
            json_output=args.json,
        )
        return 0
    resolved_model = args.model or ("brownian" if args.kind == "continuous" else "fitch")
    if args.ancestral_command == "package":
        _validate_ancestral_discrete_model_arguments(args, parser)
        result = build_ancestral_figure_package(
            tree_path=args.tree,
            traits_path=args.table,
            trait=args.trait,
            reconstruction_kind=args.kind,
            out_dir=args.out_dir,
            taxon_column=args.taxon_column,
            model=resolved_model,
            alpha=args.alpha,
            state_ordering=args.state_ordering,
            ordered_states=_split_csv_values(args.ordered_states) or None,
            layout=args.layout,
        )
        outputs = _finalize_outputs(
            args,
            command="ancestral",
            inputs=[args.tree, args.table],
            outputs=[
                result.figure_path,
                result.figure_png_path,
                result.figure_html_path,
                result.review_path,
                result.node_table_path,
                result.uncertainty_table_path,
                result.node_review_path,
                result.legend_path,
                result.model_description_path,
                result.caption_path,
                result.manifest_path,
                result.reproducibility_manifest_path,
            ],
        )
        _print_result(
            build_command_result(
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
                metrics={
                    "output_dir": str(result.output_dir),
                    "artifact_count": 12,
                    "publication_ready": result.audit.publication_ready,
                    "internal_state_visible": result.audit.internal_state_visible,
                    "uncertainty_visible": result.audit.uncertainty_visible,
                    "ambiguous_internal_node_count": (
                        result.audit.ambiguous_internal_node_count
                    ),
                    "unstable_internal_node_count": (
                        result.audit.unstable_internal_node_count
                    ),
                    "rendered_internal_annotation_count": (
                        result.audit.rendered_internal_annotation_count
                    ),
                    "rendered_internal_pie_count": (
                        result.audit.rendered_internal_pie_count
                    ),
                },
                data=result,
            ),
            json_output=args.json,
        )
        return 0
    _validate_ancestral_discrete_model_arguments(args, parser)
    if args.out is None and args.out_dir is None:
        parser.error("ancestral report requires --out or --out-dir")
    if args.out_dir is not None:
        result = build_ancestral_report_package(
            tree_path=args.tree,
            traits_path=args.table,
            trait=args.trait,
            reconstruction_kind=args.kind,
            out_dir=args.out_dir,
            taxon_column=args.taxon_column,
            model=resolved_model,
            alpha=args.alpha,
            state_ordering=args.state_ordering,
            ordered_states=_split_csv_values(args.ordered_states) or None,
            compare_model=args.compare_model,
            compare_tree_path=args.compare_tree,
            drop_taxa=args.drop_taxa,
            coding_map=_parse_assignment_map(args.coding_map) or None,
        )
        output_paths: list[Path | str] = [
            result.report_path,
            result.methods_summary_path,
            result.reviewer_audit_checklist_path,
            result.figure_path,
            result.figure_png_path,
            result.figure_html_path,
            result.summary_table_path,
            result.node_table_path,
            result.uncertainty_table_path,
            result.transition_count_table_path,
            result.transition_branch_table_path,
            result.exclusion_table_path,
            result.manifest_path,
        ]
        if args.out is not None:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(
                result.report_path.read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            svg_out = args.out.with_suffix(".svg")
            svg_out.write_text(
                result.figure_path.read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            output_paths.extend([args.out, svg_out])
        outputs = _finalize_outputs(
            args,
            command="ancestral",
            inputs=[args.tree, args.table],
            outputs=output_paths,
        )
        _print_result(
            build_command_result(
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
                metrics={
                    "report_kind": "ancestral-report-package",
                    "reconstruction_kind": result.reconstruction_kind,
                    "output_dir": str(result.output_dir),
                    "artifact_count": 13,
                    "methods_summary_warning_count": (
                        result.methods_summary.warning_count
                    ),
                    "transition_count_row_count": result.machine_manifest["metrics"][
                        "transition_count_row_count"
                    ],
                },
                data=result,
            ),
            json_output=args.json,
        )
        return 0
    if args.out is None:
        raise ValueError("ancestral report rendering requires an explicit output path")
    result = render_ancestral_state_report(
        tree_path=args.tree,
        traits_path=args.table,
        trait=args.trait,
        reconstruction_kind=args.kind,
        out_path=args.out,
        taxon_column=args.taxon_column,
        model=resolved_model,
        alpha=args.alpha,
        state_ordering=args.state_ordering,
        ordered_states=_split_csv_values(args.ordered_states) or None,
        compare_model=args.compare_model,
        compare_tree_path=args.compare_tree,
        drop_taxa=args.drop_taxa,
        coding_map=_parse_assignment_map(args.coding_map) or None,
    )
    outputs = _finalize_outputs(
        args,
        command="ancestral",
        inputs=[args.tree, args.table],
        outputs=[result.output_path, args.out.with_suffix(".svg")],
    )
    _print_result(
        build_command_result(
            command="ancestral",
            inputs=[args.tree, args.table],
            outputs=outputs,
            metrics={
                "report_kind": result.report_kind,
                "reconstruction_kind": result.reconstruction_kind,
            },
            data=result,
        ),
        json_output=args.json,
    )
    return 0
