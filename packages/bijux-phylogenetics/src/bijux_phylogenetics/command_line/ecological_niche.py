from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.ecology import (
    summarize_niche_transitions,
    write_niche_state_node_table,
    write_niche_transition_branch_table,
    write_niche_transition_clade_table,
    write_niche_transition_count_table,
    write_niche_transition_exclusion_table,
    write_niche_transition_rate_table,
    write_niche_transition_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_ecological_niche_commands(subparsers: Any) -> None:
    ecological_niche = subparsers.add_parser(
        get_command_spec("ecological-niche").name,
        help=get_command_spec("ecological-niche").summary,
    )
    ecological_niche_subparsers = ecological_niche.add_subparsers(
        dest="ecological_niche_command",
        required=True,
    )
    ecological_niche_transitions = ecological_niche_subparsers.add_parser(
        "transitions",
        help="Fit ecological niche transitions, reconstruct ancestral niches, and rank clade-specific shift burden.",
    )
    ecological_niche_transitions.add_argument("tree", type=Path)
    ecological_niche_transitions.add_argument("table", type=Path)
    ecological_niche_transitions.add_argument("--trait", required=True)
    ecological_niche_transitions.add_argument("--taxon-column")
    ecological_niche_transitions.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    ecological_niche_transitions.add_argument("--summary-out", type=Path)
    ecological_niche_transitions.add_argument("--nodes-out", type=Path)
    ecological_niche_transitions.add_argument("--rates-out", type=Path)
    ecological_niche_transitions.add_argument("--branches-out", type=Path)
    ecological_niche_transitions.add_argument("--counts-out", type=Path)
    ecological_niche_transitions.add_argument("--clades-out", type=Path)
    ecological_niche_transitions.add_argument("--exclusions-out", type=Path)
    ecological_niche_transitions.add_argument(
        "--json",
        action="store_true",
        help="Emit the ecological-niche review as JSON.",
    )
    _add_manifest_argument(ecological_niche_transitions)


def run_ecological_niche_command(args: Any) -> int:
    report = summarize_niche_transitions(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        model=args.model,
    )
    outputs: list[Path | str] = []
    if args.summary_out is not None:
        outputs.append(
            write_niche_transition_summary_table(
                args.summary_out,
                report,
            )
        )
    if args.nodes_out is not None:
        outputs.append(
            write_niche_state_node_table(
                args.nodes_out,
                report,
            )
        )
    if args.rates_out is not None:
        outputs.append(
            write_niche_transition_rate_table(
                args.rates_out,
                report,
            )
        )
    if args.branches_out is not None:
        outputs.append(
            write_niche_transition_branch_table(
                args.branches_out,
                report,
            )
        )
    if args.counts_out is not None:
        outputs.append(
            write_niche_transition_count_table(
                args.counts_out,
                report,
            )
        )
    if args.clades_out is not None:
        outputs.append(
            write_niche_transition_clade_table(
                args.clades_out,
                report,
            )
        )
    if args.exclusions_out is not None:
        outputs.append(
            write_niche_transition_exclusion_table(
                args.exclusions_out,
                report,
            )
        )
    outputs = _finalize_outputs(
        args,
        command="ecological-niche",
        inputs=[args.tree, args.table],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="ecological-niche",
            inputs=[args.tree, args.table],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "model": report.model,
                "observed_niche_count": report.summary.observed_niche_count,
                "transition_rate_row_count": report.summary.transition_rate_row_count,
                "changed_branch_count": report.summary.changed_branch_count,
                "certain_transition_count": report.summary.certain_transition_count,
                "uncertain_transition_count": report.summary.uncertain_transition_count,
                "repeated_shift_clade_count": report.summary.repeated_shift_clade_count,
                "excluded_taxon_count": report.summary.excluded_taxon_count,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
