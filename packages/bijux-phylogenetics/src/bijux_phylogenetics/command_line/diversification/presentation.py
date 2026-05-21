from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative import (
    build_diversification_figure_package,
    build_diversification_method_report,
    render_diversification_report,
    write_diversification_methods_summary_text,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .inputs import tree_and_metadata_inputs, tree_metadata_traits_inputs


def add_diversification_presentation_commands(
    diversification_subparsers: Any,
) -> None:
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


def run_diversification_presentation_command(args: Any) -> int | None:
    if args.diversification_command == "package":
        return _run_package(args)
    if args.diversification_command == "methods-summary":
        return _run_methods_summary(args)
    if args.diversification_command == "report":
        return _run_report(args)
    return None


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
    warnings = (
        [] if result.sampling_report is None else list(result.sampling_report.warnings)
    )
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
    warnings = (
        []
        if result.report.sampling_report is None
        else list(result.report.sampling_report.warnings)
    )
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
