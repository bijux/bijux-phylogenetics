from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.assessment import (
    run_comparative_sensitivity_analysis,
)
from bijux_phylogenetics.comparative.continuous import compare_brownian_and_ou_models
from bijux_phylogenetics.comparative.reporting.model_comparison_figures import (
    build_comparative_model_figure_package,
)
from bijux_phylogenetics.comparative.validation import (
    audit_comparative_parameter_uncertainty,
    audit_ou_identifiability_reference_examples,
    validate_comparative_reference_examples,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_model_assessment_comparative_evolution_commands(
    comparative_subparsers: Any,
) -> None:
    comparative_compare_models = comparative_subparsers.add_parser(
        "compare-models",
        help="Compare standalone Brownian-motion and OU models for one continuous trait.",
    )
    comparative_compare_models.add_argument("tree", type=Path)
    comparative_compare_models.add_argument("table", type=Path)
    comparative_compare_models.add_argument("--trait", required=True)
    comparative_compare_models.add_argument("--taxon-column")
    comparative_compare_models.add_argument(
        "--json", action="store_true", help="Emit the model comparison as JSON."
    )
    _add_manifest_argument(comparative_compare_models)

    comparative_model_package = comparative_subparsers.add_parser(
        "model-comparison-package",
        help="Build a publication-oriented Brownian versus OU model-comparison figure package.",
    )
    comparative_model_package.add_argument("tree", type=Path)
    comparative_model_package.add_argument("table", type=Path)
    comparative_model_package.add_argument("--trait", required=True)
    comparative_model_package.add_argument("--taxon-column")
    comparative_model_package.add_argument("--out-dir", required=True, type=Path)
    comparative_model_package.add_argument(
        "--json", action="store_true", help="Emit the package build result as JSON."
    )
    _add_manifest_argument(comparative_model_package)

    comparative_validate_reference = comparative_subparsers.add_parser(
        "validate-reference",
        help="Validate built-in Brownian-motion and OU reference examples.",
    )
    comparative_validate_reference.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(comparative_validate_reference)

    comparative_sensitivity = comparative_subparsers.add_parser(
        "sensitivity",
        help="Run leave-one-taxon-out sensitivity for a standalone BM or OU model.",
    )
    comparative_sensitivity.add_argument("tree", type=Path)
    comparative_sensitivity.add_argument("table", type=Path)
    comparative_sensitivity.add_argument("--trait", required=True)
    comparative_sensitivity.add_argument(
        "--model", choices=("brownian", "ou"), required=True
    )
    comparative_sensitivity.add_argument("--taxon-column")
    comparative_sensitivity.add_argument(
        "--json", action="store_true", help="Emit the sensitivity report as JSON."
    )
    _add_manifest_argument(comparative_sensitivity)


def run_model_assessment_comparative_evolution_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    del parser
    if args.comparative_command == "compare-models":
        report = compare_brownian_and_ou_models(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
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
                    "taxon_count": report.taxon_count,
                    "better_model": report.better_model,
                    "model_count": len(report.rows),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.comparative_command == "model-comparison-package":
        result = build_comparative_model_figure_package(
            args.tree,
            args.table,
            trait=args.trait,
            out_dir=args.out_dir,
            taxon_column=args.taxon_column,
        )
        outputs = _finalize_outputs(
            args,
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=[
                result.criteria_figure_path,
                result.likelihood_figure_path,
                result.parameter_figure_path,
                result.fit_figure_path,
                result.criteria_table_path,
                result.likelihood_table_path,
                result.parameter_table_path,
                result.fit_table_path,
                result.legend_path,
                result.caption_path,
                result.review_path,
                result.manifest_path,
                result.reproducibility_manifest_path,
            ],
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.tree, args.table],
                outputs=outputs,
                metrics={
                    "output_dir": str(result.output_dir),
                    "artifact_count": 13,
                    "publication_ready": result.audit.publication_ready,
                    "selected_model": result.audit.selected_model,
                    "support_distinct": result.audit.support_distinct,
                    "aicc_delta": result.audit.aicc_delta,
                    "plotted_model_count": result.audit.plotted_model_count,
                    "rendered_parameter_count": result.audit.rendered_parameter_count,
                    "rendered_fit_row_count": result.audit.rendered_fit_row_count,
                },
                data=result,
            ),
            json_output=args.json,
        )
        return 0

    if args.comparative_command == "validate-reference":
        report = validate_comparative_reference_examples()
        uncertainty_audit = audit_comparative_parameter_uncertainty()
        identifiability_audit = audit_ou_identifiability_reference_examples()
        outputs = _finalize_outputs(args, command="comparative", inputs=[])
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[],
                outputs=outputs,
                metrics={
                    "case_count": len(report.observations),
                    "all_passed": report.all_passed,
                    "interval_audit_passed": (
                        uncertainty_audit.all_reference_estimates_covered
                    ),
                    "identifiability_audit_passed": (
                        identifiability_audit.all_expected_warning_kinds_detected
                    ),
                },
                warnings=[
                    *uncertainty_audit.warnings,
                    *(
                        []
                        if identifiability_audit.all_expected_warning_kinds_detected
                        else [
                            "one or more expected OU warning modes were not detected on the reference fixtures"
                        ]
                    ),
                ],
                data={
                    "reference_validation": report,
                    "parameter_uncertainty_audit": uncertainty_audit,
                    "ou_identifiability_audit": identifiability_audit,
                },
            ),
            json_output=args.json,
        )
        return 0

    if args.comparative_command != "sensitivity":
        return None
    report = run_comparative_sensitivity_analysis(
        args.tree,
        args.table,
        trait=args.trait,
        model=args.model,
        taxon_column=args.taxon_column,
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
                "taxon_count": len(report.rows),
                "model": report.model,
                "influential_taxa": len(report.most_influential_taxa),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
