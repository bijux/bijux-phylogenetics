from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.ancestral.comparison import (
    compare_discrete_ancestral_reconstructions,
    write_discrete_ancestral_comparison_table,
)
from bijux_phylogenetics.ancestral.discrete import (
    discrete_ancestral_exclusions,
    reconstruct_discrete_ancestral_states,
    summarize_discrete_ancestral_report,
    write_discrete_ancestral_exclusion_table,
    write_discrete_ancestral_fit_table,
    write_discrete_ancestral_probability_table,
    write_discrete_ancestral_summary_table,
    write_discrete_ancestral_transition_table,
)
from bijux_phylogenetics.ancestral.presentation.report_rendering import (
    write_ancestral_state_table,
)
from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _split_csv_values,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_discrete_reconstruction_commands(ancestral_subparsers: Any) -> None:
    ancestral_discrete = ancestral_subparsers.add_parser(
        "discrete",
        help="Reconstruct ancestral states for a discrete trait.",
    )
    ancestral_discrete.add_argument("tree", type=Path)
    ancestral_discrete.add_argument("table", type=Path)
    ancestral_discrete.add_argument("--trait", required=True)
    ancestral_discrete.add_argument("--taxon-column")
    ancestral_discrete.add_argument(
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
    ancestral_discrete.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_discrete.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_discrete.add_argument(
        "--compare-model",
        choices=(
            "fitch",
            "equal-rates",
            "symmetric",
            "all-rates-different",
            "meristic",
        ),
        help="Optionally compare this reconstruction directly against another discrete model.",
    )
    ancestral_discrete.add_argument(
        "--root-prior-mode",
        choices=("equal", "empirical", "fixed"),
        default="equal",
        help="Likelihood-only root prior policy for discrete ancestral reconstruction.",
    )
    ancestral_discrete.add_argument(
        "--fixed-root-state",
        help="Required when --root-prior-mode fixed; names the forced root state.",
    )
    ancestral_discrete.add_argument("--table-out", type=Path)
    ancestral_discrete.add_argument("--summary-out", type=Path)
    ancestral_discrete.add_argument("--probabilities-out", type=Path)
    ancestral_discrete.add_argument("--transitions-out", type=Path)
    ancestral_discrete.add_argument("--fit-out", type=Path)
    ancestral_discrete.add_argument("--comparison-out", type=Path)
    ancestral_discrete.add_argument("--exclusions-out", type=Path)
    ancestral_discrete.add_argument(
        "--json", action="store_true", help="Emit the reconstruction as JSON."
    )
    _add_manifest_argument(ancestral_discrete)


def run_discrete_reconstruction_command(args: Any, *, parser: Any) -> int | None:
    if args.ancestral_command != "discrete":
        return None

    if args.state_ordering == "ordered" and args.model == "fitch":
        parser.error(
            "ordered ancestral discrete reconstruction requires a likelihood model"
        )
    if args.compare_model == args.model:
        parser.error(
            "discrete ancestral compare-model must differ from the primary model"
        )
    report = reconstruct_discrete_ancestral_states(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        model=args.model,
        state_ordering=args.state_ordering,
        ordered_states=_split_csv_values(args.ordered_states) or None,
        root_prior_mode=args.root_prior_mode,
        fixed_root_state=args.fixed_root_state,
    )
    summary = summarize_discrete_ancestral_report(report)
    exclusions = discrete_ancestral_exclusions(report)
    comparison = (
        compare_discrete_ancestral_reconstructions(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
            left_model=args.model,
            right_model=args.compare_model,
            state_ordering=args.state_ordering,
            ordered_states=_split_csv_values(args.ordered_states) or None,
            root_prior_mode=args.root_prior_mode,
            fixed_root_state=args.fixed_root_state,
        )
        if args.compare_model is not None
        else None
    )
    outputs: list[Path | str] = []
    if args.table_out is not None:
        outputs.append(write_ancestral_state_table(args.table_out, report))
    if args.summary_out is not None:
        outputs.append(write_discrete_ancestral_summary_table(args.summary_out, report))
    if args.probabilities_out is not None:
        outputs.append(
            write_discrete_ancestral_probability_table(args.probabilities_out, report)
        )
    if args.transitions_out is not None:
        outputs.append(
            write_discrete_ancestral_transition_table(args.transitions_out, report)
        )
    if args.fit_out is not None:
        outputs.append(write_discrete_ancestral_fit_table(args.fit_out, report))
    if args.comparison_out is not None and comparison is not None:
        outputs.append(
            write_discrete_ancestral_comparison_table(args.comparison_out, comparison)
        )
    if args.exclusions_out is not None:
        outputs.append(
            write_discrete_ancestral_exclusion_table(args.exclusions_out, report)
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
                "taxon_count": report.taxon_count,
                "estimate_count": len(report.estimates),
                "internal_node_count": summary.internal_node_count,
                "ambiguous_internal_node_count": summary.ambiguous_internal_node_count,
                "excluded_taxon_count": len(exclusions),
                "state_count": len(report.observed_states),
                "minimal_change_count": report.minimal_change_count,
                "parsimonious_root_state_count": report.parsimonious_root_state_count,
                "unstable_node_count": summary.unstable_node_count,
                "log_likelihood": report.log_likelihood,
                "parameter_count": report.parameter_count,
                "aic": report.aic,
                "root_prior_mode": report.root_prior_mode,
                "fixed_root_state": report.fixed_root_state,
                "phytools_rerooting_method_comparable": (
                    report.rerooting_method_compatibility.comparable
                ),
                "optimizer_converged": (
                    None
                    if report.optimizer_diagnostics is None
                    else report.optimizer_diagnostics.converged
                ),
                "optimizer_iteration_count": (
                    None
                    if report.optimizer_diagnostics is None
                    else report.optimizer_diagnostics.iteration_count
                ),
                "optimizer_function_evaluation_count": (
                    None
                    if report.optimizer_diagnostics is None
                    else report.optimizer_diagnostics.function_evaluation_count
                ),
                "overparameterized": report.overparameterized,
                "transition_rate_count": len(report.transition_rate_rows),
                "baseline_model": (
                    None
                    if report.baseline_comparison is None
                    else report.baseline_comparison.baseline_model
                ),
                "baseline_delta_aic": (
                    None
                    if report.baseline_comparison is None
                    else report.baseline_comparison.delta_aic
                ),
                "preferred_model_by_aic": (
                    None
                    if report.baseline_comparison is None
                    else report.baseline_comparison.preferred_model_by_aic
                ),
                "comparison_node_count": len(comparison.rows)
                if comparison is not None
                else 0,
                "comparison_differing_node_count": (
                    comparison.differing_node_count if comparison is not None else 0
                ),
                "model": report.model,
            },
            data={
                "report": report,
                "comparison": comparison,
            },
        ),
        json_output=args.json,
    )
    return 0
