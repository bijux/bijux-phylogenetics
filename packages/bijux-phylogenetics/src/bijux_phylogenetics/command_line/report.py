from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _split_csv_values,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.provenance.method_tiers import method_tier_metrics
from bijux_phylogenetics.bayesian import build_time_tree_figure_package
from bijux_phylogenetics.render.trait_tree_package import (
    build_annotated_trait_tree_package,
)
from bijux_phylogenetics.reports.service import (
    render_alignment_report,
    render_dataset_report,
    render_level_one_release_gate_report,
    render_phylo_inputs_report,
    render_production_scale_readiness_report,
    render_release_truth_report,
    render_taxon_report,
    render_tree_report,
    render_workflow_validation_report,
)
from bijux_phylogenetics.reports import (
    build_alignment_figure_package,
    write_supplementary_alignment_diagnostics_table,
    write_supplementary_ancestral_state_table,
    write_supplementary_clade_support_table,
    write_supplementary_comparative_model_table,
    write_supplementary_model_selection_table,
    write_supplementary_taxon_table,
    write_supplementary_tree_diagnostics_table,
)
from bijux_phylogenetics.reports.tree_package import build_tree_report_package
from bijux_phylogenetics.runtime.results import build_command_result


def _build_production_scale_alignment_classes(
    args: Any,
) -> list[tuple[str, int, int]] | None:
    if args.sequence_counts is None and args.alignment_lengths is None:
        return None
    sequence_counts = args.sequence_counts or []
    alignment_lengths = args.alignment_lengths or []
    if len(sequence_counts) != len(alignment_lengths):
        raise ValueError(
            "report production-scale-readiness requires the same number of --sequence-count and --alignment-length values"
        )
    return [
        (
            f"sequences-{sequence_count}-sites-{alignment_length}",
            sequence_count,
            alignment_length,
        )
        for sequence_count, alignment_length in zip(
            sequence_counts,
            alignment_lengths,
            strict=True,
        )
    ]


def _build_production_scale_tree_set_classes(
    args: Any,
) -> list[tuple[str, int, int]] | None:
    if args.posterior_tree_counts is None and args.tree_set_tip_counts is None:
        return None
    posterior_tree_counts = args.posterior_tree_counts or []
    tree_set_tip_counts = args.tree_set_tip_counts or []
    if len(posterior_tree_counts) != len(tree_set_tip_counts):
        raise ValueError(
            "report production-scale-readiness requires the same number of --posterior-tree-count and --tree-set-tip-count values"
        )
    return [
        (f"trees-{tree_count}-taxa-{tip_count}", tree_count, tip_count)
        for tree_count, tip_count in zip(
            posterior_tree_counts,
            tree_set_tip_counts,
            strict=True,
        )
    ]


def add_report_command(subparsers: Any) -> None:
    report = subparsers.add_parser(
        get_command_spec("report").name, help=get_command_spec("report").summary
    )
    report_subparsers = report.add_subparsers(dest="report_command", required=True)

    report_tree = report_subparsers.add_parser(
        "tree", help="Render a deterministic single-tree HTML report."
    )
    report_tree.add_argument("tree", type=Path)
    report_tree.add_argument("--out", required=True, type=Path)
    report_tree.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_tree)

    report_tree_package = report_subparsers.add_parser(
        "tree-package",
        help="Build a full tree report package with figure and TSV ledgers.",
    )
    report_tree_package.add_argument("tree", type=Path)
    report_tree_package.add_argument("--out-dir", required=True, type=Path)
    report_tree_package.add_argument(
        "--json", action="store_true", help="Emit the package build result as JSON."
    )
    _add_manifest_argument(report_tree_package)

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

    report_alignment = report_subparsers.add_parser(
        "alignment", help="Render an alignment-only HTML diagnostic report."
    )
    report_alignment.add_argument("--alignment", required=True, type=Path)
    report_alignment.add_argument("--out", required=True, type=Path)
    report_alignment.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_alignment)

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

    report_dataset = report_subparsers.add_parser(
        "dataset", help="Render a tree plus table dataset HTML report."
    )
    report_dataset.add_argument("--tree", required=True, type=Path)
    report_dataset.add_argument("--metadata", required=True, type=Path)
    report_dataset.add_argument("--traits", type=Path)
    report_dataset.add_argument("--alignment", type=Path)
    report_dataset.add_argument("--tip-dates", type=Path)
    report_dataset.add_argument("--calibrations", type=Path)
    report_dataset.add_argument("--out", required=True, type=Path)
    report_dataset.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_dataset)

    report_phylo_inputs = report_subparsers.add_parser(
        "phylo-inputs",
        help="Render a tree plus alignment HTML input report.",
    )
    report_phylo_inputs.add_argument("--tree", required=True, type=Path)
    report_phylo_inputs.add_argument("--alignment", required=True, type=Path)
    report_phylo_inputs.add_argument("--out", required=True, type=Path)
    report_phylo_inputs.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_phylo_inputs)

    report_taxonomy = report_subparsers.add_parser(
        "taxonomy", help="Render a reviewer-facing taxon audit HTML report."
    )
    report_taxonomy.add_argument("--tree", required=True, type=Path)
    report_taxonomy.add_argument("--synonym-table", type=Path)
    report_taxonomy.add_argument("--metadata", type=Path)
    report_taxonomy.add_argument("--traits", type=Path)
    report_taxonomy.add_argument("--alignment", type=Path)
    report_taxonomy.add_argument("--filtered-alignment", type=Path)
    report_taxonomy.add_argument("--inference-tree", type=Path)
    report_taxonomy.add_argument("--reported-taxa", type=Path)
    report_taxonomy.add_argument("--out", required=True, type=Path)
    report_taxonomy.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_taxonomy)

    report_supplementary_taxon_table = report_subparsers.add_parser(
        "supplementary-taxon-table",
        help="Write a supplementary taxon table with IDs, metadata, traits, and exclusion evidence.",
    )
    report_supplementary_taxon_table.add_argument("--tree", required=True, type=Path)
    report_supplementary_taxon_table.add_argument(
        "--metadata", required=True, type=Path
    )
    report_supplementary_taxon_table.add_argument("--traits", required=True, type=Path)
    report_supplementary_taxon_table.add_argument("--alignment", type=Path)
    report_supplementary_taxon_table.add_argument("--filtered-alignment", type=Path)
    report_supplementary_taxon_table.add_argument("--inference-tree", type=Path)
    report_supplementary_taxon_table.add_argument("--reported-taxa", type=Path)
    report_supplementary_taxon_table.add_argument("--tip-dates", type=Path)
    report_supplementary_taxon_table.add_argument("--calibrations", type=Path)
    report_supplementary_taxon_table.add_argument("--out", required=True, type=Path)
    report_supplementary_taxon_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_taxon_table)

    report_supplementary_alignment_table = report_subparsers.add_parser(
        "supplementary-alignment-table",
        help="Write a supplementary alignment diagnostics table with optional filtering status.",
    )
    report_supplementary_alignment_table.add_argument(
        "--alignment", required=True, type=Path
    )
    report_supplementary_alignment_table.add_argument("--filtered-alignment", type=Path)
    report_supplementary_alignment_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_alignment_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_alignment_table)

    report_supplementary_tree_table = report_subparsers.add_parser(
        "supplementary-tree-table",
        help="Write a supplementary tree diagnostics table with topology, support, and warning summaries.",
    )
    report_supplementary_tree_table.add_argument("--tree", required=True, type=Path)
    report_supplementary_tree_table.add_argument("--out", required=True, type=Path)
    report_supplementary_tree_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_tree_table)

    report_supplementary_clade_support_table = report_subparsers.add_parser(
        "supplementary-clade-support-table",
        help="Write a supplementary clade-support table from one reference tree and optional tree-set frequencies.",
    )
    report_supplementary_clade_support_table.add_argument(
        "--tree", required=True, type=Path
    )
    report_supplementary_clade_support_table.add_argument(
        "--comparison-tree-set", type=Path
    )
    report_supplementary_clade_support_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_clade_support_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_clade_support_table)

    report_supplementary_model_selection_table = report_subparsers.add_parser(
        "supplementary-model-selection-table",
        help="Write a supplementary model-selection table from parsed IQ-TREE report artifacts.",
    )
    report_supplementary_model_selection_table.add_argument(
        "--iqtree-report", required=True, type=Path
    )
    report_supplementary_model_selection_table.add_argument("--model-sidecar", type=Path)
    report_supplementary_model_selection_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_model_selection_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_model_selection_table)

    report_supplementary_comparative_model_table = report_subparsers.add_parser(
        "supplementary-comparative-model-table",
        help="Write a supplementary comparative-model table with coefficients, uncertainty, diagnostics, and exclusions.",
    )
    report_supplementary_comparative_model_table.add_argument(
        "--tree", required=True, type=Path
    )
    report_supplementary_comparative_model_table.add_argument(
        "--traits", required=True, type=Path
    )
    report_supplementary_comparative_model_table.add_argument(
        "--formula",
        dest="formulas",
        action="append",
        required=True,
        help="Add one comparative candidate formula. Repeat for each candidate model.",
    )
    report_supplementary_comparative_model_table.add_argument("--taxon-column")
    report_supplementary_comparative_model_table.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    report_supplementary_comparative_model_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_comparative_model_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_comparative_model_table)

    report_supplementary_ancestral_state_table = report_subparsers.add_parser(
        "supplementary-ancestral-state-table",
        help="Write a supplementary ancestral-state table with node estimates, uncertainty, and model settings.",
    )
    report_supplementary_ancestral_state_table.add_argument(
        "--tree", required=True, type=Path
    )
    report_supplementary_ancestral_state_table.add_argument(
        "--traits", required=True, type=Path
    )
    report_supplementary_ancestral_state_table.add_argument("--trait", required=True)
    report_supplementary_ancestral_state_table.add_argument(
        "--reconstruction-kind",
        choices=("continuous", "discrete"),
        required=True,
    )
    report_supplementary_ancestral_state_table.add_argument("--taxon-column")
    report_supplementary_ancestral_state_table.add_argument("--model")
    report_supplementary_ancestral_state_table.add_argument("--estimator")
    report_supplementary_ancestral_state_table.add_argument(
        "--alpha", type=float, default=1.0
    )
    report_supplementary_ancestral_state_table.add_argument(
        "--state-ordering",
        choices=("unordered", "ordered"),
        default="unordered",
    )
    report_supplementary_ancestral_state_table.add_argument(
        "--ordered-state",
        dest="ordered_states",
        action="append",
        help="Add one ordered discrete state. Repeat to define the full state order.",
    )
    report_supplementary_ancestral_state_table.add_argument(
        "--root-prior-mode",
        choices=("equal", "observed-frequency", "fixed"),
        default="equal",
    )
    report_supplementary_ancestral_state_table.add_argument("--fixed-root-state")
    report_supplementary_ancestral_state_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_ancestral_state_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_ancestral_state_table)

    report_workflow_validation = report_subparsers.add_parser(
        "workflow-validation",
        help="Render the Level 1 workflow validation fixture report.",
    )
    report_workflow_validation.add_argument("--fixtures-root", type=Path)
    report_workflow_validation.add_argument("--out", required=True, type=Path)
    report_workflow_validation.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_workflow_validation)

    report_release_gate = report_subparsers.add_parser(
        "release-gate",
        help="Render the Level 1 release gate for the checked-in workflow fixtures.",
    )
    report_release_gate.add_argument("--fixtures-root", type=Path)
    report_release_gate.add_argument("--out", required=True, type=Path)
    report_release_gate.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_release_gate)

    report_release_truth = report_subparsers.add_parser(
        "release-truth",
        help="Render one machine-produced release truth report from pytest and workflow evidence.",
    )
    report_release_truth.add_argument(
        "--test-report",
        type=Path,
        action="append",
        required=True,
        help="Path to one pytest JUnit XML report for the full test surface. Repeat to aggregate multiple sessions.",
    )
    report_release_truth.add_argument(
        "--real-engine-test-report",
        type=Path,
        action="append",
        required=True,
        help="Path to one pytest JUnit XML report for real-engine tests. Repeat to aggregate multiple sessions.",
    )
    report_release_truth.add_argument("--fixtures-root", type=Path)
    report_release_truth.add_argument(
        "--stress-tier",
        choices=("small", "heavy"),
        default="small",
        help="Governed stress tier to benchmark during release truth generation.",
    )
    report_release_truth.add_argument(
        "--parity-extended",
        action="store_true",
        help="Include the governed extended reference-parity suite in the release truth report.",
    )
    report_release_truth.add_argument("--out", required=True, type=Path)
    report_release_truth.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_release_truth)

    report_production_scale_readiness = report_subparsers.add_parser(
        "production-scale-readiness",
        help="Render one reviewer-facing production-scale readiness report from governed benchmark evidence.",
    )
    report_production_scale_readiness.add_argument("--replicates", type=int, default=1)
    report_production_scale_readiness.add_argument(
        "--tree-tip-count",
        action="append",
        dest="tree_tip_counts",
        type=int,
        help="Add one large-tree taxon count. Repeat to override the governed tree-size classes.",
    )
    report_production_scale_readiness.add_argument(
        "--sequence-count",
        action="append",
        dest="sequence_counts",
        type=int,
        help="Add one sequence count for the large-alignment classes. Repeat alongside --alignment-length.",
    )
    report_production_scale_readiness.add_argument(
        "--alignment-length",
        action="append",
        dest="alignment_lengths",
        type=int,
        help="Add one aligned-site count for the large-alignment classes. Repeat alongside --sequence-count.",
    )
    report_production_scale_readiness.add_argument(
        "--posterior-tree-count",
        action="append",
        dest="posterior_tree_counts",
        type=int,
        help="Add one posterior tree count for the tree-set classes. Repeat alongside --tree-set-tip-count.",
    )
    report_production_scale_readiness.add_argument(
        "--tree-set-tip-count",
        action="append",
        dest="tree_set_tip_counts",
        type=int,
        help="Add one taxon count for the tree-set classes. Repeat alongside --posterior-tree-count.",
    )
    report_production_scale_readiness.add_argument(
        "--stress-tier",
        action="append",
        dest="stress_tiers",
        choices=("small", "heavy"),
        help="Include one governed stress tier. Repeat to aggregate multiple tiers.",
    )
    report_production_scale_readiness.add_argument("--out", required=True, type=Path)
    report_production_scale_readiness.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_production_scale_readiness)


def run_report_command(args: Any) -> int:
    parsed_lambda_value: float | str
    if getattr(args, "lambda_value", None) == "estimate":
        parsed_lambda_value = "estimate"
    elif getattr(args, "lambda_value", None) is None:
        parsed_lambda_value = "estimate"
    else:
        parsed_lambda_value = float(args.lambda_value)

    if args.report_command == "tree":
        result = render_tree_report(tree_path=args.tree, out_path=args.out)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[args.tree],
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[args.tree],
                    outputs=outputs,
                    warnings=result.validation.warnings + result.inspection.warnings,
                    metrics={"tip_count": result.inspection.tip_count},
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "tree-package":
        result = build_tree_report_package(args.tree, out_dir=args.out_dir)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[args.tree],
            outputs=[
                result.report_path,
                result.figure_path,
                result.support_table_path,
                result.clade_table_path,
                result.branch_stats_path,
                result.manifest_path,
            ],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[args.tree],
                    outputs=outputs,
                    warnings=result.validation.warnings + result.inspection.warnings,
                    metrics={
                        "tip_count": result.inspection.tip_count,
                        "supported_branch_count": sum(
                            1 for row in result.support_rows if row.support is not None
                        ),
                        "rendered_support_count": result.figure.rendered_support_count,
                        "long_outlier_count": result.branch_stats.long_outlier_count,
                        **method_tier_metrics(result.method_tier),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_dir)
        return 0

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

    if args.report_command == "alignment":
        result = render_alignment_report(
            alignment_path=args.alignment, out_path=args.out
        )
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[args.alignment],
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[args.alignment],
                    outputs=outputs,
                    warnings=result.alignment_forensic.warnings,
                    metrics={
                        "sequence_count": result.alignment.sequence_count,
                        "alignment_length": result.alignment.alignment_length,
                        "quality_score": result.alignment_quality.quality_score,
                        "warning_count": len(result.alignment_forensic.warnings),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
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
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_dir)
        return 0

    if args.report_command == "dataset":
        result = render_dataset_report(
            tree_path=args.tree,
            metadata_path=args.metadata,
            traits_path=args.traits,
            alignment_path=args.alignment,
            tip_dates_path=args.tip_dates,
            calibration_path=args.calibrations,
            out_path=args.out,
        )
        inputs = [args.tree, args.metadata]
        if args.traits is not None:
            inputs.append(args.traits)
        if args.alignment is not None:
            inputs.append(args.alignment)
        if args.tip_dates is not None:
            inputs.append(args.tip_dates)
        if args.calibrations is not None:
            inputs.append(args.calibrations)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=result.validation.warnings + result.inspection.warnings,
                    metrics={
                        "tip_count": result.inspection.tip_count,
                        "linked_taxa": result.metadata_linkage.linked_taxa,
                        "readiness_decision": None
                        if result.dataset_audit is None
                        else result.dataset_audit.readiness_decision,
                        "excluded_taxa": 0
                        if result.dataset_audit is None
                        else len(result.dataset_audit.exclusion_table.rows),
                        "blocked_analysis_count": 0
                        if result.dataset_audit is None
                        else len(result.dataset_audit.blocked_analyses),
                        "risky_analysis_count": 0
                        if result.dataset_audit is None
                        else sum(
                            1
                            for row in result.dataset_audit.analysis_decisions
                            if row.decision == "risky"
                        ),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "phylo-inputs":
        result = render_phylo_inputs_report(
            tree_path=args.tree,
            alignment_path=args.alignment,
            out_path=args.out,
        )
        inputs = [args.tree, args.alignment]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=result.validation.warnings + result.inspection.warnings,
                    metrics={
                        "tip_count": result.inspection.tip_count,
                        "alignment_length": result.alignment.alignment_length,
                        "linked_taxa": result.alignment_linkage.linked_taxa,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "taxonomy":
        result = render_taxon_report(
            tree_path=args.tree,
            synonym_table_path=args.synonym_table,
            metadata_path=args.metadata,
            traits_path=args.traits,
            alignment_path=args.alignment,
            filtered_alignment_path=args.filtered_alignment,
            inference_tree_path=args.inference_tree,
            reported_taxa_path=args.reported_taxa,
            out_path=args.out,
        )
        inputs = [args.tree, *([args.synonym_table] if args.synonym_table is not None else [])]
        if args.metadata is not None:
            inputs.append(args.metadata)
        if args.traits is not None:
            inputs.append(args.traits)
        if args.alignment is not None:
            inputs.append(args.alignment)
        if args.filtered_alignment is not None:
            inputs.append(args.filtered_alignment)
        if args.inference_tree is not None:
            inputs.append(args.inference_tree)
        if args.reported_taxa is not None:
            inputs.append(args.reported_taxa)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=result.taxon_audit.warnings,
                    metrics={
                        "tree_tip_count": result.taxon_audit.tree_tip_count,
                        "status": result.taxon_audit.status,
                        "conflict_count": len(result.taxon_audit.mapping_conflicts.rows),
                        "crosswalk_rows": 0
                        if result.taxon_crosswalk is None
                        else len(result.taxon_crosswalk.rows),
                        "excluded_taxa": 0
                        if result.taxon_exclusions is None
                        else len(result.taxon_exclusions.rows),
                        "loss_stage_count": 0
                        if result.taxon_workflow_loss is None
                        else len(result.taxon_workflow_loss.loss_stage_counts),
                        "unstable_taxa": 0
                        if result.taxon_stability is None
                        else len(result.taxon_stability.unstable_taxa),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-taxon-table":
        result = write_supplementary_taxon_table(
            args.out,
            tree_path=args.tree,
            metadata_path=args.metadata,
            traits_path=args.traits,
            alignment_path=args.alignment,
            filtered_alignment_path=args.filtered_alignment,
            inference_tree_path=args.inference_tree,
            reported_taxa_path=args.reported_taxa,
            tip_dates_path=args.tip_dates,
            calibration_path=args.calibrations,
        )
        inputs = [args.tree, args.metadata, args.traits]
        if args.alignment is not None:
            inputs.append(args.alignment)
        if args.filtered_alignment is not None:
            inputs.append(args.filtered_alignment)
        if args.inference_tree is not None:
            inputs.append(args.inference_tree)
        if args.reported_taxa is not None:
            inputs.append(args.reported_taxa)
        if args.tip_dates is not None:
            inputs.append(args.tip_dates)
        if args.calibrations is not None:
            inputs.append(args.calibrations)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "analysis_included_count": result.analysis_included_count,
                        "analysis_excluded_count": result.analysis_excluded_count,
                        "reporting_retained_count": result.reporting_retained_count,
                        "reporting_dropped_count": result.reporting_dropped_count,
                        "metadata_column_count": result.metadata_column_count,
                        "trait_column_count": result.trait_column_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-alignment-table":
        result = write_supplementary_alignment_diagnostics_table(
            args.out,
            alignment_path=args.alignment,
            filtered_alignment_path=args.filtered_alignment,
        )
        inputs = [args.alignment]
        if args.filtered_alignment is not None:
            inputs.append(args.filtered_alignment)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "retained_sequence_count": result.retained_sequence_count,
                        "removed_sequence_count": result.removed_sequence_count,
                        "filtered_only_sequence_count": result.filtered_only_sequence_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-tree-table":
        result = write_supplementary_tree_diagnostics_table(
            args.out,
            tree_path=args.tree,
        )
        inputs = [args.tree]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            row = result.rows[0]
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=row.warnings,
                    metrics={
                        "row_count": result.row_count,
                        "tip_count": row.tip_count,
                        "supported_branch_count": row.supported_branch_count,
                        "polytomy_count": row.polytomy_count,
                        "warning_count": row.warning_count,
                        "ultrametric": row.ultrametric,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-clade-support-table":
        result = write_supplementary_clade_support_table(
            args.out,
            tree_path=args.tree,
            comparison_tree_set_path=args.comparison_tree_set,
        )
        inputs = [args.tree]
        if args.comparison_tree_set is not None:
            inputs.append(args.comparison_tree_set)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "supported_clade_count": result.supported_clade_count,
                        "frequency_scored_clade_count": result.frequency_scored_clade_count,
                        "frequency_partial_support_count": result.frequency_partial_support_count,
                        "frequency_absent_clade_count": result.frequency_absent_clade_count,
                        "frequency_unscored_clade_count": result.frequency_unscored_clade_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-model-selection-table":
        result = write_supplementary_model_selection_table(
            args.out,
            iqtree_report_path=args.iqtree_report,
            model_sidecar_path=args.model_sidecar,
        )
        inputs = [args.iqtree_report]
        if args.model_sidecar is not None:
            inputs.append(args.model_sidecar)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "candidate_count": result.candidate_count,
                        "selected_model": result.selected_model,
                        "selected_criterion": result.selected_criterion,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-comparative-model-table":
        result = write_supplementary_comparative_model_table(
            args.out,
            tree_path=args.tree,
            traits_path=args.traits,
            formulas=list(args.formulas),
            taxon_column=args.taxon_column,
            lambda_value=parsed_lambda_value,
        )
        inputs = [args.tree, args.traits]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "model_count": result.model_count,
                        "selected_formula": result.selected_formula,
                        "selected_criterion": result.selected_criterion,
                        "excluded_taxon_count": result.excluded_taxon_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-ancestral-state-table":
        result = write_supplementary_ancestral_state_table(
            args.out,
            tree_path=args.tree,
            traits_path=args.traits,
            trait=args.trait,
            reconstruction_kind=args.reconstruction_kind,
            taxon_column=args.taxon_column,
            model=args.model,
            estimator=args.estimator,
            alpha=args.alpha,
            state_ordering=args.state_ordering,
            ordered_states=args.ordered_states,
            root_prior_mode=args.root_prior_mode,
            fixed_root_state=args.fixed_root_state,
        )
        inputs = [args.tree, args.traits]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "reconstruction_kind": result.reconstruction_kind,
                        "model": result.model,
                        "analysis_taxon_count": result.analysis_taxon_count,
                        "excluded_taxon_count": result.excluded_taxon_count,
                        "unstable_node_count": result.unstable_node_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "workflow-validation":
        result = render_workflow_validation_report(
            out_path=args.out,
            fixtures_root=args.fixtures_root,
        )
        inputs = [] if args.fixtures_root is None else [args.fixtures_root]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    metrics={
                        "total_fixture_count": result.validation.total_fixture_count,
                        "passed_fixture_count": result.validation.passed_fixture_count,
                        "workflow_count": len(result.validation.workflows),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "release-gate":
        result = render_level_one_release_gate_report(
            out_path=args.out,
            fixtures_root=args.fixtures_root,
        )
        inputs = [] if args.fixtures_root is None else [args.fixtures_root]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=result.release_gate.dataset_warnings,
                    metrics={
                        "decision": result.release_gate.gate.decision,
                        "retained_taxa": len(result.release_gate.gate.retained_taxa),
                        "excluded_taxa": len(result.release_gate.gate.excluded_taxa),
                        "blocked_analysis_count": len(
                            result.release_gate.gate.blocked_analyses
                        ),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "release-truth":
        result = render_release_truth_report(
            out_path=args.out,
            test_report_paths=args.test_report,
            real_engine_test_report_paths=args.real_engine_test_report,
            fixtures_root=args.fixtures_root,
            include_extended_parity=args.parity_extended,
            stress_tier=args.stress_tier,
        )
        inputs = [
            *args.test_report,
            *args.real_engine_test_report,
            *([args.fixtures_root] if args.fixtures_root is not None else []),
        ]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=result.release_truth.known_limitations,
                    metrics={
                        "total_tests": result.release_truth.total_tests.total_tests,
                        "total_tests_passed": result.release_truth.total_tests.passed_tests,
                        "total_tests_failed": result.release_truth.total_tests.failed_tests,
                        "total_tests_skipped": result.release_truth.total_tests.skipped_tests,
                        "real_engine_tests": result.release_truth.real_engine_tests.total_tests,
                        "real_engine_tests_passed": result.release_truth.real_engine_tests.passed_tests,
                        "real_engine_tests_failed": result.release_truth.real_engine_tests.failed_tests,
                        "real_engine_tests_skipped": result.release_truth.real_engine_tests.skipped_tests,
                        "supported_workflow_count": len(
                            result.release_truth.supported_workflows
                        ),
                        "experimental_workflow_count": len(
                            result.release_truth.experimental_workflows
                        ),
                        "flagship_dataset_count": len(
                            result.release_truth.flagship_datasets
                        ),
                        "reference_parity_case_count": result.release_truth.reference_parity.case_count,
                        "stress_workload_count": len(
                            result.release_truth.stress_suite.observations
                        ),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "production-scale-readiness":
        result = render_production_scale_readiness_report(
            out_path=args.out,
            replicates=args.replicates,
            tree_tip_counts=args.tree_tip_counts,
            alignment_size_classes=_build_production_scale_alignment_classes(args),
            tree_set_size_classes=_build_production_scale_tree_set_classes(args),
            stress_tiers=args.stress_tiers,
        )
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[],
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[],
                    outputs=outputs,
                    warnings=result.production_scale_readiness.limitations,
                    metrics=result.machine_manifest["metrics"],
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    raise NotImplementedError(f"unsupported report command: {args.report_command}")
