from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.ecology import (
    summarize_host_switching,
    write_host_state_node_table,
    write_host_switch_branch_table,
    write_host_switch_count_table,
    write_host_switch_exclusion_table,
    write_host_switch_fit_table,
    write_host_switch_summary_table,
    write_unsupported_host_switch_claim_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_host_association_commands(subparsers: Any) -> None:
    host_association = subparsers.add_parser(
        get_command_spec("host-association").name,
        help=get_command_spec("host-association").summary,
    )
    host_association_subparsers = host_association.add_subparsers(
        dest="host_association_command",
        required=True,
    )
    host_association_switches = host_association_subparsers.add_parser(
        "switches",
        help="Reconstruct host states, count host switches, and compare constrained host-transition models.",
    )
    host_association_switches.add_argument("tree", type=Path)
    host_association_switches.add_argument("table", type=Path)
    host_association_switches.add_argument("--trait", required=True)
    host_association_switches.add_argument("--taxon-column")
    host_association_switches.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="ard",
    )
    host_association_switches.add_argument(
        "--constraints",
        type=Path,
        help="Optional host-transition constraint ledger with source_host and target_host columns.",
    )
    host_association_switches.add_argument("--summary-out", type=Path)
    host_association_switches.add_argument("--nodes-out", type=Path)
    host_association_switches.add_argument("--branches-out", type=Path)
    host_association_switches.add_argument("--counts-out", type=Path)
    host_association_switches.add_argument("--fits-out", type=Path)
    host_association_switches.add_argument("--unsupported-out", type=Path)
    host_association_switches.add_argument("--exclusions-out", type=Path)
    host_association_switches.add_argument(
        "--json",
        action="store_true",
        help="Emit the host-association review as JSON.",
    )
    _add_manifest_argument(host_association_switches)


def run_host_association_command(args: Any) -> int:
    report = summarize_host_switching(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        model=args.model,
        constraint_path=args.constraints,
    )
    outputs: list[Path | str] = []
    if args.summary_out is not None:
        outputs.append(
            write_host_switch_summary_table(
                args.summary_out,
                report,
            )
        )
    if args.nodes_out is not None:
        outputs.append(
            write_host_state_node_table(
                args.nodes_out,
                report,
            )
        )
    if args.branches_out is not None:
        outputs.append(
            write_host_switch_branch_table(
                args.branches_out,
                report,
            )
        )
    if args.counts_out is not None:
        outputs.append(
            write_host_switch_count_table(
                args.counts_out,
                report,
            )
        )
    if args.fits_out is not None:
        outputs.append(
            write_host_switch_fit_table(
                args.fits_out,
                report,
            )
        )
    if args.unsupported_out is not None:
        outputs.append(
            write_unsupported_host_switch_claim_table(
                args.unsupported_out,
                report,
            )
        )
    if args.exclusions_out is not None:
        outputs.append(
            write_host_switch_exclusion_table(
                args.exclusions_out,
                report,
            )
        )
    inputs = [args.tree, args.table]
    if args.constraints is not None:
        inputs.append(args.constraints)
    outputs = _finalize_outputs(
        args,
        command="host-association",
        inputs=inputs,
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="host-association",
            inputs=inputs,
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "model": report.model,
                "analysis_constraint_mode": report.summary.analysis_constraint_mode,
                "observed_host_count": report.summary.observed_host_count,
                "host_switch_count": report.summary.host_switch_count,
                "certain_host_switch_count": report.summary.certain_host_switch_count,
                "uncertain_host_switch_count": report.summary.uncertain_host_switch_count,
                "preferred_constraint": report.summary.preferred_constraint,
                "unsupported_switch_claim_count": (
                    report.summary.unsupported_switch_claim_count
                ),
                "excluded_taxon_count": report.summary.excluded_taxon_count,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
