from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import build_time_tree_figure_package
from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _split_csv_values,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.render.annotated_trait_tree_package import (
    build_annotated_trait_tree_package,
)
from bijux_phylogenetics.reports import build_alignment_figure_package
from bijux_phylogenetics.runtime.results import build_command_result


def add_figure_package_report_commands(report_subparsers: Any) -> None:
    report_trait_tree_package = report_subparsers.add_parser(
        "trait-tree-package",
        help="Build a publication-oriented annotated trait tree package with coverage and reviewer audits.",
    )
    report_trait_tree_package.add_argument("tree", type=Path)
    report_trait_tree_package.add_argument("--metadata", type=Path)
    report_trait_tree_package.add_argument("--traits", type=Path)
    report_trait_tree_package.add_argument("--taxon-column")
    report_trait_tree_package.add_argument("--label-column")
    report_trait_tree_package.add_argument("--categorical-column")
    report_trait_tree_package.add_argument("--continuous-column")
    report_trait_tree_package.add_argument("--metadata-strip-columns")
    report_trait_tree_package.add_argument("--heatmap-columns")
    report_trait_tree_package.add_argument(
        "--layout", choices=("cladogram", "phylogram", "circular"), default="phylogram"
    )
    report_trait_tree_package.add_argument("--support-labels", action="store_true")
    report_trait_tree_package.add_argument("--out-dir", required=True, type=Path)
    report_trait_tree_package.add_argument(
        "--json", action="store_true", help="Emit the package build result as JSON."
    )
    _add_manifest_argument(report_trait_tree_package)

    report_time_tree_package = report_subparsers.add_parser(
        "time-tree-package",
        help="Build a publication-oriented time-tree package with node-age labels and HPD intervals.",
    )
    report_time_tree_package.add_argument("posterior_trees", type=Path)
    report_time_tree_package.add_argument(
        "--source-format",
        choices=("generic", "beast", "mrbayes"),
        default="generic",
    )
    report_time_tree_package.add_argument("--burnin-fraction", type=float, default=0.25)
    report_time_tree_package.add_argument("--metadata", type=Path)
    report_time_tree_package.add_argument("--label-column")
    report_time_tree_package.add_argument("--taxon-column")
    report_time_tree_package.add_argument("--tip-dates", type=Path)
    report_time_tree_package.add_argument("--calibrations", type=Path)
    report_time_tree_package.add_argument("--alignment", type=Path)
    report_time_tree_package.add_argument(
        "--title",
        default="Bijux Time Tree Figure",
        help="Reviewer-facing title for the time-tree figure package.",
    )
    report_time_tree_package.add_argument("--out-dir", required=True, type=Path)
    report_time_tree_package.add_argument(
        "--json", action="store_true", help="Emit the package build result as JSON."
    )
    _add_manifest_argument(report_time_tree_package)

    report_alignment_package = report_subparsers.add_parser(
        "alignment-package",
        help="Build a publication-oriented alignment quality figure package.",
    )
    report_alignment_package.add_argument("alignment", type=Path)
    report_alignment_package.add_argument("--out-dir", required=True, type=Path)
    report_alignment_package.add_argument("--maximum-site-bins", type=int, default=120)
    report_alignment_package.add_argument("--window-size", type=int, default=30)
    report_alignment_package.add_argument("--step-size", type=int, default=10)
    report_alignment_package.add_argument("--panel-row-limit", type=int, default=12)
    report_alignment_package.add_argument(
        "--json", action="store_true", help="Emit the package build result as JSON."
    )
    _add_manifest_argument(report_alignment_package)


def run_figure_package_report_command(args: Any) -> int | None:
    if args.report_command == "trait-tree-package":
        result = build_annotated_trait_tree_package(
            args.tree,
            out_dir=args.out_dir,
            metadata_path=args.metadata,
            traits_path=args.traits,
            taxon_column=args.taxon_column,
            label_column=args.label_column,
            categorical_column=args.categorical_column,
            continuous_column=args.continuous_column,
            metadata_strip_columns=_split_csv_values(args.metadata_strip_columns),
            heatmap_columns=_split_csv_values(args.heatmap_columns),
            layout=args.layout,
            show_support_values=args.support_labels,
        )
        inputs = [args.tree]
        if args.metadata is not None:
            inputs.append(args.metadata)
        if args.traits is not None:
            inputs.append(args.traits)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[
                result.review_path,
                result.figure_package.figure_path,
                result.figure_package.caption_path,
                result.figure_package.legend_path,
                result.coverage_path,
                result.summary_path,
                result.manifest_path,
                result.reproducibility_manifest_path,
            ],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=result.audit.limitations,
                    metrics={
                        "publication_ready": result.audit.publication_ready,
                        "required_surface_count": result.audit.required_surface_count,
                        "complete_surface_count": result.audit.complete_surface_count,
                        "missing_surface_count": result.audit.missing_surface_count,
                        "visible_tip_count": result.figure_package.render.visible_tip_count,
                        "legend_entry_count": result.audit.legend_entry_count,
                        "caption_ready": result.audit.caption_ready,
                        "legible": result.audit.legible,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_dir)
        return 0

    if args.report_command == "time-tree-package":
        result = build_time_tree_figure_package(
            args.posterior_trees,
            out_dir=args.out_dir,
            source_format=args.source_format,
            burnin_fraction=args.burnin_fraction,
            metadata_path=args.metadata,
            label_column=args.label_column,
            taxon_column=args.taxon_column,
            tip_dates_path=args.tip_dates,
            calibration_path=args.calibrations,
            alignment_path=args.alignment,
            title=args.title,
        )
        inputs = [args.posterior_trees]
        if args.metadata is not None:
            inputs.append(args.metadata)
        if args.tip_dates is not None:
            inputs.append(args.tip_dates)
        if args.calibrations is not None:
            inputs.append(args.calibrations)
        if args.alignment is not None:
            inputs.append(args.alignment)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[
                result.figure_path,
                result.retained_tree_set_path,
                result.mcc_tree_path,
                result.interval_table_path,
                result.legend_path,
                result.caption_path,
                result.review_path,
                result.manifest_path,
                result.reproducibility_manifest_path,
            ],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=result.audit.limitations,
                    metrics={
                        "publication_ready": result.audit.publication_ready,
                        "retained_tree_count": result.retained_tree_count,
                        "root_age": result.render.root_age,
                        "rendered_interval_count": result.render.rendered_interval_count,
                        "expected_interval_count": result.audit.expected_interval_count,
                        "ultrametric": result.audit.ultrametric,
                        "readiness_decision": result.audit.readiness_decision,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_dir)
        return 0

    if args.report_command == "alignment-package":
        result = build_alignment_figure_package(
            args.alignment,
            out_dir=args.out_dir,
            maximum_site_bins=args.maximum_site_bins,
            window_size=args.window_size,
            step_size=args.step_size,
            panel_row_limit=args.panel_row_limit,
        )
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[args.alignment],
            outputs=[
                result.heatmap_figure_path,
                result.site_summary_figure_path,
                result.sequence_panel_figure_path,
                result.heatmap_table_path,
                result.window_table_path,
                result.ranking_table_path,
                result.legend_path,
                result.caption_path,
                result.review_path,
                result.manifest_path,
                result.reproducibility_manifest_path,
                result.reviewer_audit_checklist_path,
            ],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[args.alignment],
                    outputs=outputs,
                    warnings=result.audit.limitations,
                    metrics={
                        "publication_ready": result.audit.publication_ready,
                        "quality_score": result.audit.quality_score,
                        "suspicious_alignment": result.audit.suspicious_alignment,
                        "heatmap_row_count": result.audit.heatmap_row_count,
                        "heatmap_bin_count": result.audit.heatmap_bin_count,
                        "plotted_window_count": result.audit.plotted_window_count,
                        "plotted_sequence_count": result.audit.plotted_sequence_count,
                        "reviewer_audit_item_count": len(
                            result.reviewer_audit_checklist.items
                        ),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_dir)
        return 0

    return None
