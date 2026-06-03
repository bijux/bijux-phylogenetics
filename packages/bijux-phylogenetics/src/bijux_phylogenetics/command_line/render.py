from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _build_annotation_strips,
    _build_numeric_trait_map,
    _build_string_trait_map,
    _split_csv_values,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.render.tree_figure_package import build_tree_figure_package
from bijux_phylogenetics.render.tree_svg import (
    audit_support_label_rendering,
    render_tree_svg,
)
from bijux_phylogenetics.runtime.errors import MetadataJoinError
from bijux_phylogenetics.runtime.results import build_command_result


def add_render_command(subparsers: Any) -> None:
    render = subparsers.add_parser(
        get_command_spec("render").name, help=get_command_spec("render").summary
    )
    render.add_argument("tree", type=Path)
    render.add_argument("--metadata", type=Path)
    render.add_argument("--traits", type=Path)
    render.add_argument("--taxon-column")
    render.add_argument("--label-column")
    render.add_argument(
        "--layout", choices=["cladogram", "phylogram", "circular"], default="cladogram"
    )
    render.add_argument("--support-labels", action="store_true")
    render.add_argument("--categorical-column")
    render.add_argument("--continuous-column")
    render.add_argument("--metadata-strip-columns")
    render.add_argument("--heatmap-columns")
    render.add_argument("--collapse-clades")
    render.add_argument("--package-dir", type=Path)
    render.add_argument("--out", required=True, type=Path)
    render.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(render)


def run_render_command(args: Any) -> int:
    metadata_table = (
        load_taxon_table(args.metadata, taxon_column=args.taxon_column)
        if args.metadata is not None
        else None
    )
    traits_table = (
        load_taxon_table(args.traits, taxon_column=args.taxon_column)
        if args.traits is not None
        else None
    )
    labels: dict[str, str] | None = None
    if metadata_table is not None and args.label_column is not None:
        if args.label_column not in metadata_table.columns:
            raise MetadataJoinError(
                f"metadata table does not contain label column '{args.label_column}'"
            )
        labels = {
            row[metadata_table.taxon_column]: row[args.label_column]
            for row in metadata_table.rows
            if row[args.label_column]
        }
    categorical_traits = (
        _build_string_trait_map(traits_table, args.categorical_column)
        if traits_table is not None and args.categorical_column is not None
        else None
    )
    continuous_traits = (
        _build_numeric_trait_map(traits_table, args.continuous_column)
        if traits_table is not None and args.continuous_column is not None
        else None
    )
    metadata_strips = (
        _build_annotation_strips(
            metadata_table, _split_csv_values(args.metadata_strip_columns)
        )
        if metadata_table is not None
        else []
    )
    heatmap_columns = (
        _build_annotation_strips(traits_table, _split_csv_values(args.heatmap_columns))
        if traits_table is not None
        else []
    )
    collapsed_clades = _split_csv_values(args.collapse_clades)
    support_audit = (
        audit_support_label_rendering(args.tree) if args.support_labels else None
    )
    result = render_tree_svg(
        args.tree,
        out_path=args.out,
        labels=labels,
        layout=args.layout,
        show_support_values=args.support_labels
        and (support_audit.validated if support_audit is not None else False),
        categorical_traits=categorical_traits,
        continuous_traits=continuous_traits,
        metadata_strips=metadata_strips,
        heatmap_columns=heatmap_columns,
        collapsed_clades=collapsed_clades,
        validated_support_labels={}
        if support_audit is None
        else support_audit.labels_by_node,
        support_validation_warnings=[]
        if support_audit is None
        else support_audit.warnings,
    )
    inputs = [args.tree]
    if args.metadata is not None:
        inputs.append(args.metadata)
    if args.traits is not None:
        inputs.append(args.traits)
    outputs = [result.output_path]
    package_result = None
    if args.package_dir is not None:
        package_result = build_tree_figure_package(
            args.tree,
            out_dir=args.package_dir,
            labels=labels,
            layout=args.layout,
            show_support_values=args.support_labels,
            categorical_traits=categorical_traits,
            continuous_traits=continuous_traits,
            metadata_strips=metadata_strips,
            heatmap_columns=heatmap_columns,
            collapsed_clades=collapsed_clades,
        )
        outputs.append(package_result.output_dir)
    outputs = _finalize_outputs(args, command="render", inputs=inputs, outputs=outputs)
    if args.json:
        _print_result(
            build_command_result(
                command="render",
                inputs=inputs,
                outputs=outputs,
                warnings=result.missing_metadata_labels
                + ([] if support_audit is None else support_audit.warnings),
                metrics={
                    "tip_count": result.tip_count,
                    "visible_tip_count": result.visible_tip_count,
                    "label_count": result.label_count,
                    "rendered_support_count": result.rendered_support_count,
                    "rendered_categorical_trait_count": result.rendered_categorical_trait_count,
                    "rendered_continuous_trait_count": result.rendered_continuous_trait_count,
                    "rendered_metadata_strip_count": result.rendered_metadata_strip_count,
                    "rendered_heatmap_column_count": result.rendered_heatmap_column_count,
                    "collapsed_clade_count": result.collapsed_clade_count,
                    "figure_package_legible": None
                    if package_result is None
                    else package_result.legibility_audit.legible,
                    "figure_package_legend_entry_count": 0
                    if package_result is None
                    else len(package_result.legend_entries),
                    "figure_package_caption_ready": None
                    if package_result is None
                    else package_result.caption_draft.caption_ready,
                },
                data={
                    "render": result,
                    "figure_package_dir": package_result.output_dir
                    if package_result is not None
                    else None,
                    "figure_package_audit": None
                    if package_result is None
                    else package_result.audit,
                    "figure_package_legibility_audit": None
                    if package_result is None
                    else package_result.legibility_audit,
                    "figure_package_caption_draft": None
                    if package_result is None
                    else package_result.caption_draft,
                    "figure_package_legend_entries": None
                    if package_result is None
                    else package_result.legend_entries,
                    "support_audit": support_audit,
                },
            ),
            json_output=True,
        )
        return 0
    print(result.output_path)
    return 0
