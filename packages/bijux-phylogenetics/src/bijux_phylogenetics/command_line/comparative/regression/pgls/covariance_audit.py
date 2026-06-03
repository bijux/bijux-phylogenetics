from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.covariance import (
    summarize_comparative_covariance_audit,
    write_comparative_covariance_audit_candidate_table,
    write_comparative_covariance_audit_excluded_taxa_table,
    write_comparative_covariance_audit_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_covariance_audit_pgls_commands(comparative_subparsers: Any) -> None:
    comparative_covariance_audit = comparative_subparsers.add_parser(
        "covariance-audit",
        help="Inspect the covariance matrix shape and stability for supported comparative analyses.",
    )
    comparative_covariance_audit.add_argument("tree", type=Path)
    comparative_covariance_audit.add_argument("table", type=Path)
    comparative_covariance_audit.add_argument(
        "--analysis",
        choices=("pgls", "brownian-trait", "ou-trait"),
        required=True,
    )
    comparative_covariance_audit.add_argument("--trait")
    comparative_covariance_audit.add_argument("--response")
    comparative_covariance_audit.add_argument("--predictors", nargs="+")
    comparative_covariance_audit.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass * habitat'.",
    )
    comparative_covariance_audit.add_argument("--taxon-column")
    comparative_covariance_audit.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_covariance_audit.add_argument(
        "--alpha",
        default="estimate",
        help="Use 'estimate' or a positive OU alpha value.",
    )
    comparative_covariance_audit.add_argument(
        "--summary-out",
        type=Path,
        help="Write one covariance-audit summary row as TSV or CSV.",
    )
    comparative_covariance_audit.add_argument(
        "--candidates-out",
        type=Path,
        help="Write candidate covariance audit rows as TSV or CSV.",
    )
    comparative_covariance_audit.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write excluded covariance-audit taxa as TSV or CSV.",
    )
    comparative_covariance_audit.add_argument(
        "--json", action="store_true", help="Emit the covariance audit as JSON."
    )
    _add_manifest_argument(comparative_covariance_audit)


def run_covariance_audit_pgls_command(
    args: Any,
    *,
    lambda_value: float | str,
) -> int | None:
    if args.comparative_command != "covariance-audit":
        return None

    resolved_alpha: float | str = (
        "estimate" if args.alpha == "estimate" else float(args.alpha)
    )
    report = summarize_comparative_covariance_audit(
        args.tree,
        args.table,
        analysis=args.analysis,
        trait=args.trait,
        response=args.response,
        predictors=list(args.predictors or []),
        formula=args.formula,
        taxon_column=args.taxon_column,
        lambda_value=lambda_value,
        alpha=resolved_alpha,
    )
    outputs: list[Path | str] = []
    if args.summary_out is not None:
        outputs.append(
            write_comparative_covariance_audit_summary_table(
                args.summary_out,
                report,
            )
        )
    if args.candidates_out is not None:
        outputs.append(
            write_comparative_covariance_audit_candidate_table(
                args.candidates_out,
                report,
            )
        )
    if args.excluded_taxa_out is not None:
        outputs.append(
            write_comparative_covariance_audit_excluded_taxa_table(
                args.excluded_taxa_out,
                report,
            )
        )
    outputs = _finalize_outputs(
        args,
        command="comparative",
        inputs=[args.tree, args.table],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "analysis": report.analysis,
                "covariance_model": report.covariance_model,
                "matrix_dimension": report.matrix_dimension,
                "matrix_rank": report.matrix_rank,
                "condition_number": report.condition_number,
                "fit_strategy": report.fit_strategy,
                "singular": report.singular,
                "near_singular": report.near_singular,
                "matched_taxon_count": len(report.matched_taxa),
                "missing_from_traits_count": len(report.missing_from_traits),
                "extra_trait_taxon_count": len(report.extra_trait_taxa),
                "analysis_taxon_count": len(report.analysis_taxa),
                "duplicate_tree_taxon_count": len(report.duplicate_tree_taxa),
                "duplicate_trait_taxon_count": len(report.duplicate_trait_taxa),
                "zero_length_branch_count": report.zero_length_branch_count,
                "negative_branch_length_count": report.negative_branch_length_count,
                "candidate_row_count": len(report.candidate_rows),
                "blocker_count": len(report.blockers),
                "warning_count": len(report.warnings),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
