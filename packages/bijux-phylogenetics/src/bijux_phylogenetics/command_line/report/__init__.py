from __future__ import annotations

from typing import Any

from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.report.figure_packages import (
    add_figure_package_report_commands,
    run_figure_package_report_command,
)
from bijux_phylogenetics.command_line.report.governance import (
    add_governance_report_commands,
    run_governance_report_command,
)
from bijux_phylogenetics.command_line.report.input_reports import (
    add_input_report_commands,
    run_input_report_command,
)
from bijux_phylogenetics.command_line.report.methods import (
    add_methods_report_commands,
    run_methods_report_command,
)
from bijux_phylogenetics.command_line.report.publication import (
    add_publication_report_commands,
    run_publication_report_command,
)
from bijux_phylogenetics.command_line.report.supplementary_tables import (
    add_supplementary_table_report_commands,
    run_supplementary_table_report_command,
)


def add_report_command(subparsers: Any) -> None:
    report = subparsers.add_parser(
        get_command_spec("report").name, help=get_command_spec("report").summary
    )
    report_subparsers = report.add_subparsers(dest="report_command", required=True)
    add_publication_report_commands(report_subparsers)
    add_methods_report_commands(report_subparsers)
    add_figure_package_report_commands(report_subparsers)
    add_input_report_commands(report_subparsers)
    add_supplementary_table_report_commands(report_subparsers)
    add_governance_report_commands(report_subparsers)


def run_report_command(args: Any) -> int:
    publication_exit_code = run_publication_report_command(args)
    if publication_exit_code is not None:
        return publication_exit_code
    methods_exit_code = run_methods_report_command(args)
    if methods_exit_code is not None:
        return methods_exit_code
    figure_package_exit_code = run_figure_package_report_command(args)
    if figure_package_exit_code is not None:
        return figure_package_exit_code
    input_report_exit_code = run_input_report_command(args)
    if input_report_exit_code is not None:
        return input_report_exit_code
    supplementary_exit_code = run_supplementary_table_report_command(args)
    if supplementary_exit_code is not None:
        return supplementary_exit_code
    governance_exit_code = run_governance_report_command(args)
    if governance_exit_code is not None:
        return governance_exit_code

    raise NotImplementedError(f"unsupported report command: {args.report_command}")
