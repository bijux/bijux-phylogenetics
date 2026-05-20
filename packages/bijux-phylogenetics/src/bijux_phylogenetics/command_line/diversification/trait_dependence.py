from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative import (
    run_trait_dependent_diversification_analysis,
    write_trait_dependent_diversification_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_diversification_trait_dependence_command(
    diversification_subparsers: Any,
) -> None:
    diversification_trait = diversification_subparsers.add_parser(
        "trait-dependent",
        help="Summarize simple trait-linked diversification rates when states form interpretable clades.",
    )
    diversification_trait.add_argument("tree", type=Path)
    diversification_trait.add_argument("table", type=Path)
    diversification_trait.add_argument("--trait", required=True)
    diversification_trait.add_argument("--taxon-column")
    diversification_trait.add_argument(
        "--out",
        type=Path,
        help="Write the trait-dependent diversification table as TSV.",
    )
    diversification_trait.add_argument(
        "--json", action="store_true", help="Emit the trait-dependent report as JSON."
    )
    _add_manifest_argument(diversification_trait)


def run_diversification_trait_dependence_command(args: Any) -> int | None:
    if args.diversification_command != "trait-dependent":
        return None

    report = run_trait_dependent_diversification_analysis(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
    )
    inputs = [args.tree, args.table]
    outputs: list[Path | str] = []
    if args.out is not None:
        outputs.append(write_trait_dependent_diversification_table(args.out, report))
    outputs = _finalize_outputs(
        args,
        command="diversification",
        inputs=inputs,
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="diversification",
            inputs=inputs,
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "state_count": len(report.states),
                "monophyletic_state_count": sum(
                    1 for row in report.states if row.monophyletic
                ),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
