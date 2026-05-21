from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _split_csv_values,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.continuous import (
    summarize_brownian_regime_rates,
    write_brownian_regime_branch_table,
    write_brownian_regime_comparison_table,
    write_brownian_regime_exclusion_table,
    write_brownian_regime_profile_table,
    write_brownian_regime_rate_table,
    write_brownian_regime_summary_table,
)
from bijux_phylogenetics.comparative.traits.regime_mapping import (
    render_trait_regime_map,
    summarize_trait_regime_mapping,
    write_trait_regime_branch_table,
    write_trait_regime_exclusion_table,
    write_trait_regime_node_table,
    write_trait_regime_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_regime_workflow_comparative_evolution_commands(
    comparative_subparsers: Any,
) -> None:
    comparative_brownian_regimes = comparative_subparsers.add_parser(
        "brownian-regimes",
        help="Fit a multi-rate Brownian trait-evolution model from a branch regime map.",
    )
    comparative_brownian_regimes.add_argument("tree", type=Path)
    comparative_brownian_regimes.add_argument("table", type=Path)
    comparative_brownian_regimes.add_argument("regime_map", type=Path)
    comparative_brownian_regimes.add_argument("--trait", required=True)
    comparative_brownian_regimes.add_argument("--taxon-column")
    comparative_brownian_regimes.add_argument("--branch-id-column")
    comparative_brownian_regimes.add_argument("--regime-column", default="regime")
    comparative_brownian_regimes.add_argument(
        "--summary-out",
        type=Path,
        help="Write one overall multi-rate Brownian summary ledger as TSV or CSV.",
    )
    comparative_brownian_regimes.add_argument(
        "--rates-out",
        type=Path,
        help="Write one regime-rate ledger as TSV or CSV.",
    )
    comparative_brownian_regimes.add_argument(
        "--comparison-out",
        type=Path,
        help="Write one single-rate versus multi-rate comparison ledger as TSV or CSV.",
    )
    comparative_brownian_regimes.add_argument(
        "--profile-out",
        type=Path,
        help="Write one conditional regime-rate profile ledger as TSV or CSV.",
    )
    comparative_brownian_regimes.add_argument(
        "--branches-out",
        type=Path,
        help="Write one normalized branch-regime assignment ledger as TSV or CSV.",
    )
    comparative_brownian_regimes.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for the multi-rate Brownian fit as TSV or CSV.",
    )
    comparative_brownian_regimes.add_argument(
        "--json",
        action="store_true",
        help="Emit the multi-rate Brownian model fit as JSON.",
    )
    _add_manifest_argument(comparative_brownian_regimes)

    comparative_regime_map = comparative_subparsers.add_parser(
        "regime-map",
        help="Reconstruct or normalize branch regime assignments for comparative workflows.",
    )
    comparative_regime_map.add_argument("tree", type=Path)
    regime_map_source = comparative_regime_map.add_mutually_exclusive_group(
        required=True
    )
    regime_map_source.add_argument(
        "--table",
        type=Path,
        help="Tip-state table used to reconstruct branch regimes.",
    )
    regime_map_source.add_argument(
        "--regime-map",
        type=Path,
        help="User-provided branch regime table to validate and normalize.",
    )
    comparative_regime_map.add_argument(
        "--trait",
        help="Discrete trait column used when reconstructing regimes from tip states.",
    )
    comparative_regime_map.add_argument("--taxon-column")
    comparative_regime_map.add_argument(
        "--reconstruction-model",
        default="fitch",
        choices=(
            "fitch",
            "equal-rates",
            "symmetric",
            "all-rates-different",
            "meristic",
        ),
    )
    comparative_regime_map.add_argument(
        "--state-ordering",
        default="unordered",
        choices=("unordered", "ordered"),
    )
    comparative_regime_map.add_argument(
        "--ordered-states",
        help="Comma-delimited explicit ordered state vocabulary.",
    )
    comparative_regime_map.add_argument("--branch-id-column")
    comparative_regime_map.add_argument("--regime-column", default="regime")
    comparative_regime_map.add_argument(
        "--summary-out",
        type=Path,
        help="Write one regime-map summary ledger as TSV or CSV.",
    )
    comparative_regime_map.add_argument(
        "--branches-out",
        type=Path,
        help="Write one normalized branch-regime ledger as TSV or CSV.",
    )
    comparative_regime_map.add_argument(
        "--nodes-out",
        type=Path,
        help="Write one node-reconstruction ledger as TSV or CSV.",
    )
    comparative_regime_map.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for tip-state regime mapping as TSV or CSV.",
    )
    comparative_regime_map.add_argument(
        "--svg-out",
        type=Path,
        help="Render one SVG regime-map figure.",
    )
    comparative_regime_map.add_argument(
        "--layout",
        default="cladogram",
        choices=("cladogram", "phylogram", "circular"),
    )
    comparative_regime_map.add_argument(
        "--json",
        action="store_true",
        help="Emit the regime-map report as JSON.",
    )
    _add_manifest_argument(comparative_regime_map)


def run_regime_workflow_comparative_evolution_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    if args.comparative_command == "brownian-regimes":
        report = summarize_brownian_regime_rates(
            args.tree,
            args.table,
            args.regime_map,
            trait=args.trait,
            taxon_column=args.taxon_column,
            branch_id_column=args.branch_id_column,
            regime_column=args.regime_column,
        )
        if args.summary_out:
            write_brownian_regime_summary_table(args.summary_out, report)
        if args.rates_out:
            write_brownian_regime_rate_table(args.rates_out, report)
        if args.comparison_out:
            write_brownian_regime_comparison_table(args.comparison_out, report)
        if args.profile_out:
            write_brownian_regime_profile_table(args.profile_out, report)
        if args.branches_out:
            write_brownian_regime_branch_table(args.branches_out, report)
        if args.excluded_taxa_out:
            write_brownian_regime_exclusion_table(args.excluded_taxa_out, report)
        outputs = _finalize_outputs(
            args,
            command="comparative",
            inputs=[args.tree, args.table, args.regime_map],
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.tree, args.table, args.regime_map],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "tree_taxon_count": report.tree_taxon_count,
                    "analyzed_taxon_count": report.analyzed_taxon_count,
                    "excluded_taxon_count": len(report.excluded_taxa),
                    "regime_count": len(report.regime_rows),
                    "root_state": report.root_state,
                    "log_likelihood": report.log_likelihood,
                    "aic": report.aic,
                    "aicc": report.aicc,
                    "better_model": report.better_model,
                    "likelihood_ratio_statistic": report.likelihood_ratio_statistic,
                    "likelihood_ratio_p_value": report.likelihood_ratio_p_value,
                    "identifiability_warning_count": len(
                        report.identifiability_warnings
                    ),
                    "profile_row_count": len(report.profile_rows),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.comparative_command != "regime-map":
        return None
    if args.table and not args.trait:
        parser.error("--trait is required when reconstructing regimes from --table")
    report = summarize_trait_regime_mapping(
        args.tree,
        tip_states_path=args.table,
        regime_map_path=args.regime_map,
        trait=args.trait,
        taxon_column=args.taxon_column,
        reconstruction_model=args.reconstruction_model,
        state_ordering=args.state_ordering,
        ordered_states=_split_csv_values(args.ordered_states) or None,
        branch_id_column=args.branch_id_column,
        regime_column=args.regime_column,
    )
    render = None
    if args.summary_out:
        write_trait_regime_summary_table(args.summary_out, report)
    if args.branches_out:
        write_trait_regime_branch_table(args.branches_out, report)
    if args.nodes_out:
        write_trait_regime_node_table(args.nodes_out, report)
    if args.excluded_taxa_out:
        write_trait_regime_exclusion_table(args.excluded_taxa_out, report)
    if args.svg_out:
        render = render_trait_regime_map(
            report,
            out_path=args.svg_out,
            layout=args.layout,
        )
    inputs = [args.tree, args.table or args.regime_map]
    outputs = _finalize_outputs(
        args,
        command="comparative",
        inputs=inputs,
    )
    _print_result(
        build_command_result(
            command="comparative",
            inputs=inputs,
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "source_kind": report.source_kind,
                "tree_taxon_count": report.tree_taxon_count,
                "analyzed_taxon_count": report.analyzed_taxon_count,
                "excluded_taxon_count": len(report.excluded_taxa),
                "regime_count": len(report.observed_regimes),
                "branch_count": len(report.branch_rows),
                "node_count": len(report.node_rows),
                "ambiguous_branch_count": report.ambiguous_branch_count,
                "rendered_internal_annotation_count": (
                    0 if render is None else render.rendered_internal_annotation_count
                ),
                "rendered_categorical_trait_count": (
                    0 if render is None else render.rendered_categorical_trait_count
                ),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
