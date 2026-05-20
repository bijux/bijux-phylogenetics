from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative import (
    build_diversification_figure_package,
    build_diversification_method_report,
    compare_diversification_models,
    detect_diversification_outlier_clades,
    estimate_diversification_rate,
    render_diversification_report,
    run_trait_dependent_diversification_analysis,
    summarize_geiger_birth_death_exclusion,
    summarize_medusa_exclusion,
    write_clade_diversification_table,
    write_diversification_methods_summary_text,
    write_trait_dependent_diversification_table,
)
from bijux_phylogenetics.runtime.errors import DiversificationAnalysisError
from bijux_phylogenetics.runtime.results import build_command_result

from .inputs import tree_and_metadata_inputs, tree_metadata_traits_inputs
from .inspection import (
    add_diversification_inspection_commands,
    run_diversification_inspection_command,
)


def add_diversification_commands(subparsers: Any) -> None:
    diversification = subparsers.add_parser(
        get_command_spec("diversification").name,
        help=get_command_spec("diversification").summary,
    )
    diversification_subparsers = diversification.add_subparsers(
        dest="diversification_command",
        required=True,
    )
    add_diversification_inspection_commands(diversification_subparsers)

    diversification_estimate = diversification_subparsers.add_parser(
        "estimate",
        help="Estimate a simple Yule or birth-death diversification model.",
    )
    diversification_estimate.add_argument("tree", type=Path)
    diversification_estimate.add_argument("--metadata", type=Path)
    diversification_estimate.add_argument("--taxon-column")
    diversification_estimate.add_argument("--sampling-column")
    diversification_estimate.add_argument(
        "--model", choices=("yule", "birth-death"), default="birth-death"
    )
    diversification_estimate.add_argument(
        "--json", action="store_true", help="Emit the diversification estimate as JSON."
    )
    _add_manifest_argument(diversification_estimate)

    diversification_compare = diversification_subparsers.add_parser(
        "compare-models",
        help="Compare Yule and birth-death diversification fits.",
    )
    diversification_compare.add_argument("tree", type=Path)
    diversification_compare.add_argument("--metadata", type=Path)
    diversification_compare.add_argument("--taxon-column")
    diversification_compare.add_argument("--sampling-column")
    diversification_compare.add_argument(
        "--json", action="store_true", help="Emit the model comparison as JSON."
    )
    _add_manifest_argument(diversification_compare)

    diversification_clades = diversification_subparsers.add_parser(
        "clades",
        help="Detect clades with unusually high or low diversification.",
    )
    diversification_clades.add_argument("tree", type=Path)
    diversification_clades.add_argument(
        "--model", choices=("yule", "birth-death"), default="birth-death"
    )
    diversification_clades.add_argument("--min-tip-count", type=int, default=2)
    diversification_clades.add_argument(
        "--out", type=Path, help="Write the clade diversification table as TSV."
    )
    diversification_clades.add_argument(
        "--json", action="store_true", help="Emit the clade scan report as JSON."
    )
    _add_manifest_argument(diversification_clades)

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

    diversification_package = diversification_subparsers.add_parser(
        "package",
        help="Build a publication-oriented diversification figure package.",
    )
    diversification_package.add_argument("tree", type=Path)
    diversification_package.add_argument("--metadata", type=Path)
    diversification_package.add_argument("--taxon-column")
    diversification_package.add_argument("--sampling-column")
    diversification_package.add_argument(
        "--model", choices=("yule", "birth-death"), default="birth-death"
    )
    diversification_package.add_argument("--min-tip-count", type=int, default=2)
    diversification_package.add_argument("--out-dir", required=True, type=Path)
    diversification_package.add_argument(
        "--json", action="store_true", help="Emit the package build result as JSON."
    )
    _add_manifest_argument(diversification_package)

    diversification_report = diversification_subparsers.add_parser(
        "report",
        help="Render an HTML diversification and macroevolution report.",
    )
    diversification_report.add_argument("tree", type=Path)
    diversification_report.add_argument("--metadata", type=Path)
    diversification_report.add_argument("--taxon-column")
    diversification_report.add_argument("--sampling-column")
    diversification_report.add_argument("--traits", type=Path)
    diversification_report.add_argument("--trait")
    diversification_report.add_argument("--out", required=True, type=Path)
    diversification_report.add_argument(
        "--methods-summary-out",
        type=Path,
        help="Write reviewer-facing Markdown methods text for the diversification analysis.",
    )
    diversification_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(diversification_report)

    diversification_methods_summary = diversification_subparsers.add_parser(
        "methods-summary",
        help="Write reviewer-facing Markdown methods text for one diversification analysis.",
    )
    diversification_methods_summary.add_argument("tree", type=Path)
    diversification_methods_summary.add_argument("--metadata", type=Path)
    diversification_methods_summary.add_argument("--taxon-column")
    diversification_methods_summary.add_argument("--sampling-column")
    diversification_methods_summary.add_argument("--traits", type=Path)
    diversification_methods_summary.add_argument("--trait")
    diversification_methods_summary.add_argument(
        "--estimate-model",
        choices=("yule", "birth-death"),
        default="birth-death",
    )
    diversification_methods_summary.add_argument(
        "--clade-model",
        choices=("yule", "birth-death"),
        default="birth-death",
    )
    diversification_methods_summary.add_argument(
        "--min-tip-count",
        type=int,
        default=2,
        help="Minimum clade size included in the clade outlier review.",
    )
    diversification_methods_summary.add_argument("--out", required=True, type=Path)
    diversification_methods_summary.add_argument(
        "--json", action="store_true", help="Emit the methods summary metrics as JSON."
    )
    _add_manifest_argument(diversification_methods_summary)

    diversification_medusa = diversification_subparsers.add_parser(
        "medusa",
        help="Explain the explicit exclusion boundary for geiger::medusa parity.",
    )
    diversification_medusa.add_argument("tree", type=Path)
    diversification_medusa.add_argument("--metadata", type=Path)
    diversification_medusa.add_argument("--taxon-column")
    diversification_medusa.add_argument("--sampling-column")
    diversification_medusa.add_argument(
        "--json", action="store_true", help="Emit the MEDUSA exclusion as JSON."
    )
    _add_manifest_argument(diversification_medusa)

    diversification_bd_ms = diversification_subparsers.add_parser(
        "bd-ms",
        help="Explain the explicit exclusion boundary for geiger::bd.ms birth-death parity.",
    )
    diversification_bd_ms.add_argument("tree", type=Path)
    diversification_bd_ms.add_argument("--metadata", type=Path)
    diversification_bd_ms.add_argument("--taxon-column")
    diversification_bd_ms.add_argument("--sampling-column")
    diversification_bd_ms.add_argument(
        "--json",
        action="store_true",
        help="Emit the birth-death exclusion as JSON.",
    )
    _add_manifest_argument(diversification_bd_ms)


def _run_estimate(args: Any) -> int:
    inputs = tree_and_metadata_inputs(args)
    report = estimate_diversification_rate(
        args.tree,
        metadata_path=args.metadata,
        taxon_column=args.taxon_column,
        sampling_column=args.sampling_column,
        model=args.model,
    )
    outputs = _finalize_outputs(args, command="diversification", inputs=inputs)
    _print_result(
        build_command_result(
            command="diversification",
            inputs=inputs,
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "model": report.model,
                "sampling_fraction": report.sampling_fraction,
                "net_diversification_rate": report.net_diversification_rate,
                "aic": report.aic,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def _run_compare_models(args: Any) -> int:
    inputs = tree_and_metadata_inputs(args)
    report = compare_diversification_models(
        args.tree,
        metadata_path=args.metadata,
        taxon_column=args.taxon_column,
        sampling_column=args.sampling_column,
    )
    outputs = _finalize_outputs(args, command="diversification", inputs=inputs)
    _print_result(
        build_command_result(
            command="diversification",
            inputs=inputs,
            outputs=outputs,
            metrics={
                "better_model": report.better_model,
                "model_count": len(report.rows),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def _run_clades(args: Any) -> int:
    report = detect_diversification_outlier_clades(
        args.tree,
        min_tip_count=args.min_tip_count,
        model=args.model,
    )
    outputs: list[Path | str] = []
    if args.out is not None:
        outputs.append(write_clade_diversification_table(args.out, report))
    outputs = _finalize_outputs(
        args,
        command="diversification",
        inputs=[args.tree],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="diversification",
            inputs=[args.tree],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "global_rate": report.global_rate,
                "high_clade_count": len(report.high_diversification_clades),
                "low_clade_count": len(report.low_diversification_clades),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def _run_trait_dependent(args: Any) -> int:
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


def _run_package(args: Any) -> int:
    inputs = tree_and_metadata_inputs(args)
    result = build_diversification_figure_package(
        args.tree,
        out_dir=args.out_dir,
        metadata_path=args.metadata,
        taxon_column=args.taxon_column,
        sampling_column=args.sampling_column,
        min_tip_count=args.min_tip_count,
        model=args.model,
    )
    outputs = _finalize_outputs(
        args,
        command="diversification",
        inputs=inputs,
        outputs=[
            result.lineage_figure_path,
            result.clade_figure_path,
            result.model_figure_path,
            result.lineage_table_path,
            result.clade_table_path,
            result.model_table_path,
            result.legend_path,
            result.caption_path,
            result.methods_summary_path,
            result.review_path,
            result.manifest_path,
            result.reproducibility_manifest_path,
        ],
    )
    warnings = [] if result.sampling_report is None else list(result.sampling_report.warnings)
    _print_result(
        build_command_result(
            command="diversification",
            inputs=inputs,
            outputs=outputs,
            warnings=warnings,
            metrics={
                "publication_ready": result.audit.publication_ready,
                "sampling_metadata_complete": result.audit.sampling_metadata_complete,
                "plotted_ltt_point_count": result.audit.plotted_ltt_point_count,
                "plotted_clade_count": result.audit.plotted_clade_count,
                "highlighted_outlier_count": result.audit.highlighted_outlier_count,
                "plotted_model_count": result.audit.plotted_model_count,
                "better_model": result.audit.better_model,
                "methods_summary_warning_count": result.methods_summary.warning_count,
            },
            data=result.machine_manifest,
        ),
        json_output=args.json,
    )
    return 0


def _run_methods_summary(args: Any) -> int:
    inputs = tree_metadata_traits_inputs(args)
    report = build_diversification_method_report(
        args.tree,
        metadata_path=args.metadata,
        taxon_column=args.taxon_column,
        sampling_column=args.sampling_column,
        traits_path=args.traits,
        trait=args.trait,
        estimate_model=args.estimate_model,
        clade_model=args.clade_model,
        clade_min_tip_count=args.min_tip_count,
    )
    result = write_diversification_methods_summary_text(args.out, report)
    outputs = _finalize_outputs(
        args,
        command="diversification",
        inputs=inputs,
        outputs=[result.output_path],
    )
    warnings = [] if result.report.sampling_report is None else list(result.report.sampling_report.warnings)
    _print_result(
        build_command_result(
            command="diversification",
            inputs=inputs,
            outputs=outputs,
            warnings=warnings,
            metrics={
                "warning_count": result.warning_count,
                "better_model": result.better_model,
                "sampling_metadata_complete": result.sampling_metadata_complete,
                "clade_observation_count": result.clade_observation_count,
                "trait_state_count": (
                    0
                    if result.report.trait_report is None
                    else len(result.report.trait_report.states)
                ),
            },
            data=result,
        ),
        json_output=args.json,
    )
    return 0


def _raise_medusa_exclusion(args: Any) -> None:
    report = summarize_medusa_exclusion(
        args.tree,
        metadata_path=args.metadata,
        taxon_column=args.taxon_column,
        sampling_column=args.sampling_column,
    )
    raise DiversificationAnalysisError(
        report.exclusion_reason,
        code="diversification_medusa_explicitly_excluded",
        details={
            "failure_reason": report.exclusion_code,
            "supported_surfaces": report.supported_surfaces,
            "missing_surfaces": report.missing_surfaces,
            "tip_count": report.validation.tip_count,
            "rooted": report.validation.rooted,
            "ultrametric": report.validation.ultrametric,
            "sampling_metadata_complete": (
                None if report.sampling_report is None else report.sampling_report.complete
            ),
        },
    )


def _raise_bd_ms_exclusion(args: Any) -> None:
    report = summarize_geiger_birth_death_exclusion(
        args.tree,
        metadata_path=args.metadata,
        taxon_column=args.taxon_column,
        sampling_column=args.sampling_column,
    )
    raise DiversificationAnalysisError(
        report.exclusion_reason,
        code="diversification_birth_death_explicitly_excluded",
        details={
            "failure_reason": report.exclusion_code,
            "geiger_reference_surface": report.geiger_reference_surface,
            "geiger_reference_arguments": report.geiger_reference_arguments,
            "owned_surface": report.owned_surface,
            "tip_count": report.validation.tip_count,
            "rooted": report.validation.rooted,
            "ultrametric": report.validation.ultrametric,
            "sampling_metadata_complete": (
                None if report.sampling_report is None else report.sampling_report.complete
            ),
        },
    )


def _run_report(args: Any) -> int:
    inputs = tree_metadata_traits_inputs(args)
    result = render_diversification_report(
        tree_path=args.tree,
        out_path=args.out,
        metadata_path=args.metadata,
        taxon_column=args.taxon_column,
        sampling_column=args.sampling_column,
        traits_path=args.traits,
        trait=args.trait,
        methods_summary_path=args.methods_summary_out,
    )
    output_paths: list[Path | str] = [result.output_path]
    if result.methods_summary_path is not None:
        output_paths.append(result.methods_summary_path)
    outputs = _finalize_outputs(
        args,
        command="diversification",
        inputs=inputs,
        outputs=output_paths,
    )
    _print_result(
        build_command_result(
            command="diversification",
            inputs=inputs,
            outputs=outputs,
            metrics={
                "report_kind": result.report_kind,
                "methods_summary_warning_count": result.methods_summary_warning_count,
                "better_model": result.report.model_comparison.better_model,
            },
            data=result,
        ),
        json_output=args.json,
    )
    return 0


def run_diversification_command(args: Any) -> int:
    inspection_exit_code = run_diversification_inspection_command(args)
    if inspection_exit_code is not None:
        return inspection_exit_code

    if args.diversification_command == "estimate":
        return _run_estimate(args)
    if args.diversification_command == "compare-models":
        return _run_compare_models(args)
    if args.diversification_command == "clades":
        return _run_clades(args)
    if args.diversification_command == "trait-dependent":
        return _run_trait_dependent(args)
    if args.diversification_command == "package":
        return _run_package(args)
    if args.diversification_command == "methods-summary":
        return _run_methods_summary(args)
    if args.diversification_command == "medusa":
        _raise_medusa_exclusion(args)
    if args.diversification_command == "bd-ms":
        _raise_bd_ms_exclusion(args)
    return _run_report(args)
