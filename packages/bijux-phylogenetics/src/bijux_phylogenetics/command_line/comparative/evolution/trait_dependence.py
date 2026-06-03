from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.trait_dependence import (
    summarize_correlated_trait_evolution,
    write_correlated_trait_comparison_table,
    write_correlated_trait_exclusion_table,
    write_correlated_trait_observation_table,
    write_correlated_trait_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_trait_dependence_comparative_evolution_commands(
    comparative_subparsers: Any,
) -> None:
    comparative_correlated_traits = comparative_subparsers.add_parser(
        "correlated-traits",
        help="Review correlated evolution between two traits on one tree.",
    )
    comparative_correlated_traits.add_argument("tree", type=Path)
    comparative_correlated_traits.add_argument("table", type=Path)
    comparative_correlated_traits.add_argument("--left-trait", required=True)
    comparative_correlated_traits.add_argument("--right-trait", required=True)
    comparative_correlated_traits.add_argument("--taxon-column")
    comparative_correlated_traits.add_argument(
        "--analysis-kind",
        choices=("auto", "continuous", "binary"),
        default="auto",
        help="Choose auto trait-kind detection or force continuous or binary coupling analysis.",
    )
    comparative_correlated_traits.add_argument(
        "--binary-model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="all-rates-different",
        help="Discrete transition surface for binary-binary correlated-trait review.",
    )
    comparative_correlated_traits.add_argument(
        "--summary-out",
        type=Path,
        help="Write the correlated-trait summary ledger as TSV or CSV.",
    )
    comparative_correlated_traits.add_argument(
        "--comparison-out",
        type=Path,
        help="Write the independent-versus-correlated model comparison ledger as TSV or CSV.",
    )
    comparative_correlated_traits.add_argument(
        "--observations-out",
        type=Path,
        help="Write the contrast or tip-state observation ledger as TSV or CSV.",
    )
    comparative_correlated_traits.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write the explicit excluded-taxa ledger as TSV or CSV.",
    )
    comparative_correlated_traits.add_argument(
        "--json", action="store_true", help="Emit the correlated-trait report as JSON."
    )
    _add_manifest_argument(comparative_correlated_traits)


def run_trait_dependence_comparative_evolution_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    del parser
    if args.comparative_command != "correlated-traits":
        return None
    report = summarize_correlated_trait_evolution(
        args.tree,
        args.table,
        left_trait=args.left_trait,
        right_trait=args.right_trait,
        taxon_column=args.taxon_column,
        analysis_kind=args.analysis_kind,
        binary_model=args.binary_model,
    )
    outputs: list[Path | str] = []
    if args.summary_out is not None:
        outputs.append(write_correlated_trait_summary_table(args.summary_out, report))
    if args.comparison_out is not None:
        outputs.append(
            write_correlated_trait_comparison_table(args.comparison_out, report)
        )
    if args.observations_out is not None:
        outputs.append(
            write_correlated_trait_observation_table(args.observations_out, report)
        )
    if args.excluded_taxa_out is not None:
        outputs.append(
            write_correlated_trait_exclusion_table(args.excluded_taxa_out, report)
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
                "analysis_kind": report.analysis_kind,
                "tree_taxon_count": report.tree_taxon_count,
                "analyzed_taxon_count": len(report.analyzed_taxa),
                "excluded_taxon_count": len(report.excluded_taxa),
                "observation_row_count": len(report.observation_rows),
                "comparison_row_count": len(report.comparison_rows),
                "association_measure_name": report.association_measure_name,
                "association_measure_value": report.association_measure_value,
                "evolutionary_covariance": report.evolutionary_covariance,
                "evolutionary_correlation": report.evolutionary_correlation,
                "better_model": report.better_model,
                "likelihood_ratio_p_value": report.likelihood_ratio_p_value,
                "joint_state_count": len(report.joint_state_counts),
                "warning_count": len(report.warnings),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
