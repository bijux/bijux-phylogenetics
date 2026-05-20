from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.clades.stability import (
    analyze_comparative_clade_stability,
    write_comparative_clade_coefficient_change_table,
    write_comparative_clade_stability_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_clade_stability_commands(comparative_subparsers: Any) -> None:
    comparative_clade_stability = comparative_subparsers.add_parser(
        "clade-stability",
        help="Refit one comparative model after removing each major internal clade.",
    )
    comparative_clade_stability.add_argument("tree", type=Path)
    comparative_clade_stability.add_argument("table", type=Path)
    comparative_clade_stability.add_argument("--response")
    comparative_clade_stability.add_argument("--predictors", nargs="+")
    comparative_clade_stability.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass + habitat'.",
    )
    comparative_clade_stability.add_argument("--taxon-column")
    comparative_clade_stability.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1. Binary-response clade stability requires a numeric value.",
    )
    comparative_clade_stability.add_argument(
        "--clades-out",
        type=Path,
        help="Write the leave-one-clade-out stability summary ledger as TSV or CSV.",
    )
    comparative_clade_stability.add_argument(
        "--terms-out",
        type=Path,
        help="Write the coefficient-delta ledger across clade removals as TSV or CSV.",
    )
    comparative_clade_stability.add_argument(
        "--json", action="store_true", help="Emit the clade-stability report as JSON."
    )
    _add_manifest_argument(comparative_clade_stability)


def run_clade_stability_command(
    args: Any,
    *,
    lambda_value: float | str,
) -> int | None:
    if args.comparative_command != "clade-stability":
        return None

    report = analyze_comparative_clade_stability(
        args.tree,
        args.table,
        response=args.response,
        predictors=list(args.predictors or []),
        formula=args.formula,
        taxon_column=args.taxon_column,
        lambda_value=lambda_value,
    )
    outputs: list[Path | str] = []
    if args.clades_out is not None:
        outputs.append(
            write_comparative_clade_stability_table(
                args.clades_out,
                report,
            )
        )
    if args.terms_out is not None:
        outputs.append(
            write_comparative_clade_coefficient_change_table(
                args.terms_out,
                report,
            )
        )
    outputs = _finalize_outputs(
        args,
        command="comparative",
        inputs=[args.tree, args.table],
        outputs=outputs,
    )
    top_clade = (
        min(
            (
                row
                for row in report.clade_rows
                if row.fit_status == "fit" and row.rank > 0
            ),
            key=lambda row: row.rank,
        ).clade_id
        if any(row.fit_status == "fit" and row.rank > 0 for row in report.clade_rows)
        else None
    )
    _print_result(
        build_command_result(
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "model_family": report.model_family,
                "baseline_taxon_count": len(report.baseline_taxa),
                "baseline_term_count": report.baseline_term_count,
                "candidate_clade_count": report.candidate_clade_count,
                "blocked_clade_count": report.blocked_clade_count,
                "coefficient_change_row_count": len(report.coefficient_rows),
                "top_influential_clade": top_clade,
                "major_clade_fraction": report.major_clade_fraction,
                "minimum_major_clade_size": report.minimum_major_clade_size,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
