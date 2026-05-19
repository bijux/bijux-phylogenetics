from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _split_csv_values,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.command_line.comparative_evolution.discrete_models import (
    add_discrete_model_comparative_evolution_commands,
    run_discrete_model_comparative_evolution_command,
)
from bijux_phylogenetics.command_line.comparative_evolution.continuous_models import (
    add_continuous_model_comparative_evolution_commands,
    run_continuous_model_comparative_evolution_command,
)
from bijux_phylogenetics.command_line.comparative_evolution.trait_dependence import (
    add_trait_dependence_comparative_evolution_commands,
    run_trait_dependence_comparative_evolution_command,
)
from bijux_phylogenetics.comparative.brownian_regime_rates import (
    summarize_brownian_regime_rates,
    write_brownian_regime_branch_table,
    write_brownian_regime_comparison_table,
    write_brownian_regime_exclusion_table,
    write_brownian_regime_profile_table,
    write_brownian_regime_rate_table,
    write_brownian_regime_summary_table,
)
from bijux_phylogenetics.comparative.disparity_through_time import (
    render_disparity_through_time_svg,
    summarize_continuous_clade_disparity,
    summarize_disparity_through_time,
    write_continuous_clade_disparity_table,
    write_continuous_clade_disparity_summary_table,
    write_disparity_through_time_bin_table,
    write_disparity_through_time_curve_table,
    write_disparity_through_time_exclusion_table,
    write_disparity_through_time_summary_table,
)
from bijux_phylogenetics.comparative.models import (
    audit_comparative_parameter_uncertainty,
    audit_ou_identifiability_reference_examples,
    compare_brownian_and_ou_models,
    run_comparative_sensitivity_analysis,
    validate_comparative_reference_examples,
)
from bijux_phylogenetics.comparative.model_comparison_package import (
    build_comparative_model_figure_package,
)
from bijux_phylogenetics.comparative.trait_regime_mapping import (
    render_trait_regime_map,
    summarize_trait_regime_mapping,
    write_trait_regime_branch_table,
    write_trait_regime_exclusion_table,
    write_trait_regime_node_table,
    write_trait_regime_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_comparative_evolution_commands(comparative_subparsers: Any) -> None:
    add_discrete_model_comparative_evolution_commands(comparative_subparsers)
    add_trait_dependence_comparative_evolution_commands(comparative_subparsers)
    add_continuous_model_comparative_evolution_commands(comparative_subparsers)

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
        choices=("fitch", "equal-rates", "symmetric", "all-rates-different", "meristic"),
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

    comparative_dtt = comparative_subparsers.add_parser(
        "dtt",
        help="Summarize geiger-style disparity through time for one continuous trait matrix.",
    )
    comparative_dtt.add_argument("tree", type=Path)
    comparative_dtt.add_argument("table", type=Path)
    comparative_dtt.add_argument(
        "--traits",
        required=True,
        help="Comma-delimited continuous trait columns used as the disparity matrix.",
    )
    comparative_dtt.add_argument("--taxon-column")
    comparative_dtt.add_argument(
        "--time-bin-count",
        type=int,
        help="Optional equal-width time-bin count used to summarize the raw DTT curve.",
    )
    comparative_dtt.add_argument(
        "--summary-out",
        type=Path,
        help="Write one disparity-through-time summary ledger as TSV or CSV.",
    )
    comparative_dtt.add_argument(
        "--curve-out",
        type=Path,
        help="Write one raw disparity-through-time curve ledger as TSV or CSV.",
    )
    comparative_dtt.add_argument(
        "--clades-out",
        type=Path,
        help="Write one internal-clade disparity ledger as TSV or CSV.",
    )
    comparative_dtt.add_argument(
        "--bins-out",
        type=Path,
        help="Write one equal-width time-bin disparity ledger as TSV or CSV.",
    )
    comparative_dtt.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for the DTT fit as TSV or CSV.",
    )
    comparative_dtt.add_argument(
        "--svg-out",
        type=Path,
        help="Render one SVG disparity-through-time figure.",
    )
    comparative_dtt.add_argument(
        "--json",
        action="store_true",
        help="Emit the disparity-through-time report as JSON.",
    )
    _add_manifest_argument(comparative_dtt)

    comparative_disparity = comparative_subparsers.add_parser(
        "disparity",
        help="Summarize geiger-style continuous clade disparity for one rooted tree.",
    )
    comparative_disparity.add_argument("tree", type=Path)
    comparative_disparity.add_argument("table", type=Path)
    comparative_disparity.add_argument(
        "--traits",
        required=True,
        help="Comma-delimited continuous trait columns used as the disparity matrix.",
    )
    comparative_disparity.add_argument("--taxon-column")
    comparative_disparity.add_argument(
        "--summary-out",
        type=Path,
        help="Write one clade disparity summary ledger as TSV or CSV.",
    )
    comparative_disparity.add_argument(
        "--clades-out",
        type=Path,
        help="Write one internal-clade disparity ledger as TSV or CSV.",
    )
    comparative_disparity.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for clade disparity as TSV or CSV.",
    )
    comparative_disparity.add_argument(
        "--json",
        action="store_true",
        help="Emit the clade disparity report as JSON.",
    )
    _add_manifest_argument(comparative_disparity)

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


def run_comparative_evolution_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    discrete_model_exit_code = run_discrete_model_comparative_evolution_command(
        args,
        parser=parser,
    )
    if discrete_model_exit_code is not None:
        return discrete_model_exit_code

    trait_dependence_exit_code = run_trait_dependence_comparative_evolution_command(
        args,
        parser=parser,
    )
    if trait_dependence_exit_code is not None:
        return trait_dependence_exit_code

    continuous_model_exit_code = run_continuous_model_comparative_evolution_command(
        args,
        parser=parser,
    )
    if continuous_model_exit_code is not None:
        return continuous_model_exit_code

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

    if args.comparative_command == "regime-map":
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
                        0
                        if render is None
                        else render.rendered_internal_annotation_count
                    ),
                    "rendered_categorical_trait_count": (
                        0
                        if render is None
                        else render.rendered_categorical_trait_count
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.comparative_command == "dtt":
        report = summarize_disparity_through_time(
            args.tree,
            args.table,
            trait_columns=_split_csv_values(args.traits),
            taxon_column=args.taxon_column,
            time_bin_count=args.time_bin_count,
        )
        rendered_point_count = 0
        if args.summary_out:
            write_disparity_through_time_summary_table(args.summary_out, report)
        if args.curve_out:
            write_disparity_through_time_curve_table(args.curve_out, report)
        if args.clades_out:
            write_continuous_clade_disparity_table(args.clades_out, report)
        if args.bins_out:
            write_disparity_through_time_bin_table(args.bins_out, report)
        if args.excluded_taxa_out:
            write_disparity_through_time_exclusion_table(
                args.excluded_taxa_out,
                report,
            )
        if args.svg_out:
            rendered_point_count = render_disparity_through_time_svg(
                args.svg_out,
                report,
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
                warnings=report.warnings,
                metrics={
                    "tree_taxon_count": report.tree_taxon_count,
                    "analyzed_taxon_count": report.analyzed_taxon_count,
                    "excluded_taxon_count": len(report.excluded_taxa),
                    "trait_column_count": len(report.trait_columns),
                    "curve_point_count": len(report.curve_rows),
                    "time_bin_count": len(report.time_bin_rows),
                    "root_disparity": report.root_disparity,
                    "relative_scaling_applied": report.relative_scaling_applied,
                    "rendered_point_count": rendered_point_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.comparative_command == "disparity":
        report = summarize_continuous_clade_disparity(
            args.tree,
            args.table,
            trait_columns=_split_csv_values(args.traits),
            taxon_column=args.taxon_column,
        )
        if args.summary_out:
            write_continuous_clade_disparity_summary_table(args.summary_out, report)
        if args.clades_out:
            write_continuous_clade_disparity_table(args.clades_out, report)
        if args.excluded_taxa_out:
            write_disparity_through_time_exclusion_table(
                args.excluded_taxa_out,
                report,
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
                warnings=report.warnings,
                metrics={
                    "tree_taxon_count": report.tree_taxon_count,
                    "analyzed_taxon_count": report.analyzed_taxon_count,
                    "excluded_taxon_count": len(report.excluded_taxa),
                    "trait_column_count": len(report.trait_columns),
                    "clade_count": len(report.clade_rows),
                    "root_disparity": report.root_disparity,
                    "minimum_clade_disparity": report.minimum_clade_disparity,
                    "maximum_clade_disparity": report.maximum_clade_disparity,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

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

    if args.comparative_command == "sensitivity":
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

    return None
