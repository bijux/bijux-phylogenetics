from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete import (
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.ancestral.presentation import (
    build_ancestral_figure_package,
    build_ancestral_report_package,
    render_ancestral_state_visualization,
)
from bijux_phylogenetics.ancestral.presentation.report_rendering import (
    render_ancestral_state_report,
)
from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _parse_assignment_map,
    _split_csv_values,
    _validate_ancestral_discrete_model_arguments,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_presentation_ancestral_commands(ancestral_subparsers: Any) -> None:
    ancestral_render = ancestral_subparsers.add_parser(
        "render",
        help="Render a tree annotated with reconstructed ancestral states.",
    )
    ancestral_render.add_argument("tree", type=Path)
    ancestral_render.add_argument("table", type=Path)
    ancestral_render.add_argument("--trait", required=True)
    ancestral_render.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_render.add_argument("--taxon-column")
    ancestral_render.add_argument("--model")
    ancestral_render.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_render.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_render.add_argument("--alpha", type=float, default=1.0)
    ancestral_render.add_argument(
        "--layout", choices=("cladogram", "phylogram", "circular"), default="phylogram"
    )
    ancestral_render.add_argument(
        "--discrete-node-style", choices=("labels", "pies"), default="labels"
    )
    ancestral_render.add_argument(
        "--branch-coloring", choices=("none", "state", "regime"), default="none"
    )
    ancestral_render.add_argument("--out", required=True, type=Path)
    ancestral_render.add_argument(
        "--json", action="store_true", help="Emit the render result as JSON."
    )
    _add_manifest_argument(ancestral_render)

    ancestral_report = ancestral_subparsers.add_parser(
        "report",
        help="Render an HTML report for ancestral-state reconstruction.",
    )
    ancestral_report.add_argument("tree", type=Path)
    ancestral_report.add_argument("table", type=Path)
    ancestral_report.add_argument("--trait", required=True)
    ancestral_report.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_report.add_argument("--taxon-column")
    ancestral_report.add_argument("--model")
    ancestral_report.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_report.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_report.add_argument("--alpha", type=float, default=1.0)
    ancestral_report.add_argument("--compare-model")
    ancestral_report.add_argument("--compare-tree", type=Path)
    ancestral_report.add_argument("--drop-taxa", nargs="+")
    ancestral_report.add_argument(
        "--coding-map",
        help="Comma-delimited KEY=VALUE recoding map for discrete traits.",
    )
    ancestral_report.add_argument("--out", type=Path)
    ancestral_report.add_argument(
        "--out-dir",
        type=Path,
        help="Write a full ancestral reconstruction report package directory.",
    )
    ancestral_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(ancestral_report)

    ancestral_package = ancestral_subparsers.add_parser(
        "package",
        help="Write a publication-ready ancestral-state figure package.",
    )
    ancestral_package.add_argument("tree", type=Path)
    ancestral_package.add_argument("table", type=Path)
    ancestral_package.add_argument("--trait", required=True)
    ancestral_package.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_package.add_argument("--taxon-column")
    ancestral_package.add_argument("--model")
    ancestral_package.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_package.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_package.add_argument("--alpha", type=float, default=1.0)
    ancestral_package.add_argument(
        "--layout", choices=("cladogram", "phylogram", "circular"), default="phylogram"
    )
    ancestral_package.add_argument("--out-dir", required=True, type=Path)
    ancestral_package.add_argument(
        "--json", action="store_true", help="Emit the package build result as JSON."
    )
    _add_manifest_argument(ancestral_package)


def run_presentation_ancestral_command(args: Any, *, parser: Any) -> int | None:
    if args.ancestral_command == "render":
        _validate_ancestral_discrete_model_arguments(args, parser)
        if args.kind == "continuous" and args.branch_coloring == "state":
            parser.error(
                "continuous ancestral rendering does not support branch coloring 'state'"
            )
        if args.kind == "discrete" and args.branch_coloring == "regime":
            parser.error(
                "discrete ancestral rendering does not support branch coloring 'regime'"
            )
        if args.kind == "continuous":
            resolved_model = args.model or "brownian"
            reconstruction = reconstruct_continuous_ancestral_states(
                args.tree,
                args.table,
                trait=args.trait,
                taxon_column=args.taxon_column,
                model=resolved_model,
                alpha=args.alpha,
            )
        else:
            resolved_model = args.model or "fitch"
            reconstruction = reconstruct_discrete_ancestral_states(
                args.tree,
                args.table,
                trait=args.trait,
                taxon_column=args.taxon_column,
                model=resolved_model,
                state_ordering=args.state_ordering,
                ordered_states=_split_csv_values(args.ordered_states) or None,
            )
        result = render_ancestral_state_visualization(
            args.tree,
            reconstruction,
            out_path=args.out,
            layout=args.layout,
            discrete_node_style=args.discrete_node_style,
            branch_coloring=args.branch_coloring,
        )
        rendered_outputs = (
            [result.output_path]
            if result.format == "svg"
            else [result.output_path, result.svg_path]
        )
        outputs = _finalize_outputs(
            args,
            command="ancestral",
            inputs=[args.tree, args.table],
            outputs=rendered_outputs,
        )
        _print_result(
            build_command_result(
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=getattr(reconstruction, "warnings", []),
                metrics={
                    "tip_count": result.tree_render.tip_count,
                    "format": result.format,
                    "layout": result.layout,
                    "rendered_internal_annotation_count": (
                        result.tree_render.rendered_internal_annotation_count
                    ),
                    "rendered_internal_pie_count": (
                        result.tree_render.rendered_internal_pie_count
                    ),
                    "rendered_branch_color_count": (
                        result.tree_render.rendered_branch_color_count
                    ),
                },
                data={
                    "reconstruction": reconstruction,
                    "visualization": result,
                },
            ),
            json_output=args.json,
        )
        return 0

    if args.ancestral_command not in {"package", "report"}:
        return None

    resolved_model = args.model or (
        "brownian" if args.kind == "continuous" else "fitch"
    )
    if args.ancestral_command == "package":
        _validate_ancestral_discrete_model_arguments(args, parser)
        result = build_ancestral_figure_package(
            tree_path=args.tree,
            traits_path=args.table,
            trait=args.trait,
            reconstruction_kind=args.kind,
            out_dir=args.out_dir,
            taxon_column=args.taxon_column,
            model=resolved_model,
            alpha=args.alpha,
            state_ordering=args.state_ordering,
            ordered_states=_split_csv_values(args.ordered_states) or None,
            layout=args.layout,
        )
        outputs = _finalize_outputs(
            args,
            command="ancestral",
            inputs=[args.tree, args.table],
            outputs=[
                result.figure_path,
                result.figure_png_path,
                result.figure_html_path,
                result.review_path,
                result.node_table_path,
                result.uncertainty_table_path,
                result.node_review_path,
                result.legend_path,
                result.model_description_path,
                result.caption_path,
                result.manifest_path,
                result.reproducibility_manifest_path,
            ],
        )
        _print_result(
            build_command_result(
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
                metrics={
                    "output_dir": str(result.output_dir),
                    "artifact_count": 12,
                    "publication_ready": result.audit.publication_ready,
                    "internal_state_visible": result.audit.internal_state_visible,
                    "uncertainty_visible": result.audit.uncertainty_visible,
                    "ambiguous_internal_node_count": (
                        result.audit.ambiguous_internal_node_count
                    ),
                    "unstable_internal_node_count": (
                        result.audit.unstable_internal_node_count
                    ),
                    "rendered_internal_annotation_count": (
                        result.audit.rendered_internal_annotation_count
                    ),
                    "rendered_internal_pie_count": (
                        result.audit.rendered_internal_pie_count
                    ),
                },
                data=result,
            ),
            json_output=args.json,
        )
        return 0

    _validate_ancestral_discrete_model_arguments(args, parser)
    if args.out is None and args.out_dir is None:
        parser.error("ancestral report requires --out or --out-dir")
    if args.out_dir is not None:
        result = build_ancestral_report_package(
            tree_path=args.tree,
            traits_path=args.table,
            trait=args.trait,
            reconstruction_kind=args.kind,
            out_dir=args.out_dir,
            taxon_column=args.taxon_column,
            model=resolved_model,
            alpha=args.alpha,
            state_ordering=args.state_ordering,
            ordered_states=_split_csv_values(args.ordered_states) or None,
            compare_model=args.compare_model,
            compare_tree_path=args.compare_tree,
            drop_taxa=args.drop_taxa,
            coding_map=_parse_assignment_map(args.coding_map) or None,
        )
        output_paths: list[Path | str] = [
            result.report_path,
            result.methods_summary_path,
            result.reviewer_audit_checklist_path,
            result.figure_path,
            result.figure_png_path,
            result.figure_html_path,
            result.summary_table_path,
            result.node_table_path,
            result.uncertainty_table_path,
            result.transition_count_table_path,
            result.transition_branch_table_path,
            result.exclusion_table_path,
            result.manifest_path,
        ]
        if args.out is not None:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(
                result.report_path.read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            svg_out = args.out.with_suffix(".svg")
            svg_out.write_text(
                result.figure_path.read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            output_paths.extend([args.out, svg_out])
        outputs = _finalize_outputs(
            args,
            command="ancestral",
            inputs=[args.tree, args.table],
            outputs=output_paths,
        )
        _print_result(
            build_command_result(
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=outputs,
                metrics={
                    "report_kind": "ancestral-report-package",
                    "reconstruction_kind": result.reconstruction_kind,
                    "output_dir": str(result.output_dir),
                    "artifact_count": 13,
                    "methods_summary_warning_count": (
                        result.methods_summary.warning_count
                    ),
                    "transition_count_row_count": result.machine_manifest["metrics"][
                        "transition_count_row_count"
                    ],
                },
                data=result,
            ),
            json_output=args.json,
        )
        return 0

    if args.out is None:
        raise ValueError("ancestral report rendering requires an explicit output path")
    result = render_ancestral_state_report(
        tree_path=args.tree,
        traits_path=args.table,
        trait=args.trait,
        reconstruction_kind=args.kind,
        out_path=args.out,
        taxon_column=args.taxon_column,
        model=resolved_model,
        alpha=args.alpha,
        state_ordering=args.state_ordering,
        ordered_states=_split_csv_values(args.ordered_states) or None,
        compare_model=args.compare_model,
        compare_tree_path=args.compare_tree,
        drop_taxa=args.drop_taxa,
        coding_map=_parse_assignment_map(args.coding_map) or None,
    )
    outputs = _finalize_outputs(
        args,
        command="ancestral",
        inputs=[args.tree, args.table],
        outputs=[result.output_path, args.out.with_suffix(".svg")],
    )
    _print_result(
        build_command_result(
            command="ancestral",
            inputs=[args.tree, args.table],
            outputs=outputs,
            metrics={
                "report_kind": result.report_kind,
                "reconstruction_kind": result.reconstruction_kind,
            },
            data=result,
        ),
        json_output=args.json,
    )
    return 0
