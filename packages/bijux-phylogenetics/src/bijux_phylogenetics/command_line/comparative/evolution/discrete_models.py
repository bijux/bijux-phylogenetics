from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.discrete_mk import (
    fit_discrete_mk_model,
    write_discrete_mk_rate_table,
    write_discrete_mk_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_discrete_model_comparative_evolution_commands(
    comparative_subparsers: Any,
) -> None:
    comparative_discrete_mk = comparative_subparsers.add_parser(
        "discrete-mk",
        help="Fit one discrete Mk likelihood model for a tip-state trait.",
    )
    comparative_discrete_mk.add_argument("tree", type=Path)
    comparative_discrete_mk.add_argument("table", type=Path)
    comparative_discrete_mk.add_argument("--trait", required=True)
    comparative_discrete_mk.add_argument("--taxon-column")
    comparative_discrete_mk.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
        help="Choose the discrete Mk rate-constraint surface to fit.",
    )
    comparative_discrete_mk.add_argument(
        "--ascertainment",
        choices=("none", "lewis-variable-only"),
        default="none",
        help="Condition the Mk likelihood on the retained character-observation policy.",
    )
    comparative_discrete_mk.add_argument(
        "--transform",
        choices=("lambda", "kappa", "delta", "early-burst"),
        help="Apply one governed branch-length transform before fitting the discrete Mk surface.",
    )
    comparative_discrete_mk.add_argument(
        "--summary-out",
        type=Path,
        help="Write one discrete Mk fit summary ledger as TSV or CSV.",
    )
    comparative_discrete_mk.add_argument(
        "--rates-out",
        type=Path,
        help="Write one fitted directed rate-matrix ledger as TSV or CSV.",
    )
    comparative_discrete_mk.add_argument(
        "--json", action="store_true", help="Emit the discrete Mk report as JSON."
    )
    _add_manifest_argument(comparative_discrete_mk)


def run_discrete_model_comparative_evolution_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    del parser
    if args.comparative_command != "discrete-mk":
        return None
    report = fit_discrete_mk_model(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        model=args.model,
        ascertainment_policy=args.ascertainment,
        transform=args.transform,
    )
    outputs: list[Path | str] = []
    if args.summary_out:
        outputs.append(write_discrete_mk_summary_table(args.summary_out, report))
    if args.rates_out:
        outputs.append(write_discrete_mk_rate_table(args.rates_out, report))
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
            warnings=report.input_audit.warnings,
            metrics={
                "taxon_count": report.taxon_count,
                "model": report.model,
                "ascertainment_policy": report.ascertainment_policy,
                "transform": (
                    None
                    if report.transform_fit is None
                    else report.transform_fit.transform_name
                ),
                "observed_state_count": len(report.input_audit.observed_states),
                "sparse_state_count": len(report.input_audit.sparse_states),
                "pruned_missing_value_taxon_count": len(
                    report.input_audit.pruned_missing_value_taxa
                ),
                "log_likelihood": report.log_likelihood,
                "ascertainment_conditioning_log_probability": (
                    report.ascertainment_conditioning_log_probability
                ),
                "invariant_pattern_log_probability": (
                    report.invariant_pattern_log_probability
                ),
                "parameter_count": report.parameter_count,
                "aic": report.aic,
                "aicc": report.aicc,
                "optimizer_name": report.optimizer_diagnostics.optimizer_name,
                "optimizer_converged": report.optimizer_diagnostics.converged,
                "optimizer_iteration_count": (
                    report.optimizer_diagnostics.iteration_count
                ),
                "optimizer_function_evaluation_count": (
                    report.optimizer_diagnostics.function_evaluation_count
                ),
                "optimizer_hit_lower_parameter_bound": (
                    report.optimizer_diagnostics.hit_lower_parameter_bound
                ),
                "optimizer_hit_upper_parameter_bound": (
                    report.optimizer_diagnostics.hit_upper_parameter_bound
                ),
                "transform_parameter_name": (
                    None
                    if report.transform_fit is None
                    else report.transform_fit.parameter_name
                ),
                "transform_parameter_value": (
                    None
                    if report.transform_fit is None
                    else report.transform_fit.parameter_value
                ),
                "transform_function_evaluation_count": (
                    None
                    if report.transform_fit is None
                    else report.transform_fit.function_evaluation_count
                ),
                "transform_warning_count": (
                    0
                    if report.transform_fit is None
                    else len(report.transform_fit.warnings)
                ),
                "overparameterized": report.overparameterized,
                "transition_rate_count": len(report.transition_rate_rows),
                "baseline_model": (
                    None
                    if report.baseline_comparison is None
                    else report.baseline_comparison.baseline_model
                ),
                "baseline_aic": (
                    None
                    if report.baseline_comparison is None
                    else report.baseline_comparison.baseline_aic
                ),
                "delta_aic": (
                    None
                    if report.baseline_comparison is None
                    else report.baseline_comparison.delta_aic
                ),
                "preferred_model_by_aic": (
                    None
                    if report.baseline_comparison is None
                    else report.baseline_comparison.preferred_model_by_aic
                ),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
