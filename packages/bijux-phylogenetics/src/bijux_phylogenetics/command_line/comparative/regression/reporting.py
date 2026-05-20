from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.reporting import (
    build_comparative_method_report,
    build_trait_influence_report,
    write_comparative_method_report,
    write_comparative_methods_summary_text,
)
from bijux_phylogenetics.comparative.reporting.analysis_package import (
    build_comparative_report_package,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_regression_reporting_commands(comparative_subparsers: Any) -> None:
    comparative_report = comparative_subparsers.add_parser(
        "report",
        help="Build an integrated comparative-method report.",
    )
    comparative_report.add_argument("tree", type=Path)
    comparative_report.add_argument("table", type=Path)
    comparative_report.add_argument("--response")
    comparative_report.add_argument("--predictors", nargs="+")
    comparative_report.add_argument("--formula")
    comparative_report.add_argument("--taxon-column")
    comparative_report.add_argument("--out", type=Path)
    comparative_report.add_argument(
        "--methods-summary-out",
        type=Path,
        help="Write reviewer-facing Markdown methods text for the comparative analysis.",
    )
    comparative_report.add_argument(
        "--out-dir",
        type=Path,
        help="Write a full comparative analysis package directory with HTML and reviewer TSV ledgers.",
    )
    comparative_report.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_report.add_argument(
        "--json", action="store_true", help="Emit the comparative report as JSON."
    )
    _add_manifest_argument(comparative_report)

    comparative_influence = comparative_subparsers.add_parser(
        "influence",
        help="Identify predictor terms and taxa driving one comparative result.",
    )
    comparative_influence.add_argument("tree", type=Path)
    comparative_influence.add_argument("table", type=Path)
    comparative_influence.add_argument("--response")
    comparative_influence.add_argument("--predictors", nargs="+")
    comparative_influence.add_argument("--formula")
    comparative_influence.add_argument("--taxon-column")
    comparative_influence.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_influence.add_argument(
        "--json", action="store_true", help="Emit the influence report as JSON."
    )
    _add_manifest_argument(comparative_influence)


def run_regression_reporting_command(
    args: Any,
    *,
    lambda_value: float | str,
) -> int | None:
    if args.comparative_command == "report":
        package_result = None
        if args.out_dir is not None:
            package_result = build_comparative_report_package(
                args.tree,
                args.table,
                out_dir=args.out_dir,
                response=args.response,
                predictors=list(args.predictors or []),
                formula=args.formula,
                taxon_column=args.taxon_column,
                lambda_value=lambda_value,
            )
            report = package_result.report
        else:
            report = build_comparative_method_report(
                args.tree,
                args.table,
                response=args.response,
                predictors=list(args.predictors or []),
                formula=args.formula,
                taxon_column=args.taxon_column,
                lambda_value=lambda_value,
            )
        output_paths: list[Path | str] = []
        if args.out is not None:
            write_comparative_method_report(args.out, report)
            output_paths.append(args.out)
        if args.methods_summary_out is not None:
            write_comparative_methods_summary_text(args.methods_summary_out, report)
            output_paths.append(args.methods_summary_out)
        if package_result is not None:
            output_paths.extend(
                [
                    package_result.report_path,
                    package_result.methods_summary_path,
                    package_result.reviewer_audit_checklist_path,
                    package_result.summary_table_path,
                    package_result.coefficient_table_path,
                    package_result.residual_table_path,
                    package_result.signal_table_path,
                    package_result.model_comparison_table_path,
                    package_result.interpretation_table_path,
                    package_result.audit_table_path,
                    package_result.contrast_table_path,
                    package_result.manifest_path,
                ]
            )
        outputs = _finalize_outputs(
            args,
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=output_paths,
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.tree, args.table],
                outputs=outputs,
                metrics={
                    "taxon_count": report.snapshot.pgls_model.taxon_count,
                    "selected_model": report.snapshot.model_comparison.better_model,
                    "audit_row_count": len(report.snapshot.audit_rows),
                    "excluded_taxa": len(
                        report.snapshot.pgls_inputs.formula_audit.excluded_taxa
                    ),
                    "limitation_count": len(report.snapshot.limitations),
                    "coefficient_count": len(report.snapshot.pgls_model.coefficients),
                    "package_output_count": len(outputs),
                },
                data=report if package_result is None else package_result,
            ),
            json_output=args.json,
        )
        return 0

    if args.comparative_command != "influence":
        return None

    report = build_trait_influence_report(
        args.tree,
        args.table,
        response=args.response,
        predictors=list(args.predictors or []),
        formula=args.formula,
        taxon_column=args.taxon_column,
        lambda_value=lambda_value,
    )
    outputs = _finalize_outputs(
        args,
        command="comparative",
        inputs=[args.tree, args.table],
    )
    _print_result(
        build_command_result(
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=outputs,
            metrics={
                "predictor_count": len(report.predictor_rows),
                "taxon_count": len(report.taxon_rows),
                "top_predictor_terms": len(report.top_predictor_terms),
                "top_taxa": len(report.top_taxa),
                "selected_model": report.selected_model,
            },
            warnings=report.warnings,
            data=report,
        ),
        json_output=args.json,
    )
    return 0
