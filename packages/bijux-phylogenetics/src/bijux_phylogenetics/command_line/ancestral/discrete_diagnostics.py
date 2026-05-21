from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.ancestral.discrete.review import (
    summarize_ancestral_transition_report,
    summarize_ancestral_transition_tree_set,
    summarize_ancestral_transition_tree_set_report,
    summarize_ancestral_transitions,
    summarize_irreversible_discrete_reconstruction,
    summarize_irreversible_discrete_report,
    summarize_ordered_discrete_reconstruction,
    summarize_ordered_discrete_report,
    write_ancestral_transition_branch_table,
    write_ancestral_transition_count_table,
    write_ancestral_transition_exclusion_table,
    write_ancestral_transition_summary_table,
    write_ancestral_transition_tree_set_branch_table,
    write_ancestral_transition_tree_set_count_table,
    write_ancestral_transition_tree_set_summary_table,
    write_ancestral_transition_tree_set_tree_table,
    write_irreversible_discrete_fit_table,
    write_irreversible_discrete_node_table,
    write_irreversible_discrete_summary_table,
    write_irreversible_discrete_transition_table,
    write_ordered_discrete_fit_table,
    write_ordered_discrete_node_table,
    write_ordered_discrete_summary_table,
    write_ordered_discrete_transition_table,
)
from bijux_phylogenetics.ancestral.sensitivity import (
    summarize_ancestral_root_sensitivity,
    summarize_ancestral_root_sensitivity_report,
    write_ancestral_root_assumption_table,
    write_ancestral_root_sensitivity_node_table,
    write_ancestral_root_sensitivity_summary_table,
)
from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _parse_transition_pairs,
    _split_csv_values,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_discrete_diagnostic_ancestral_commands(ancestral_subparsers: Any) -> None:
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
        choices=(
            "fitch",
            "equal-rates",
            "symmetric",
            "all-rates-different",
            "meristic",
        ),
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


def run_discrete_diagnostic_ancestral_command(args: Any, *, parser: Any) -> int | None:
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

    if args.ancestral_command != "transitions":
        return None

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
