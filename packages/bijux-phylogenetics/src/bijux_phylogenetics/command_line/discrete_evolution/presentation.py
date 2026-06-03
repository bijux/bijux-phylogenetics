from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.discrete_evolution import (
    estimate_ancestral_geographic_states,
    render_discrete_state_evolution_report,
    render_tree_with_geographic_states,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .shared import COMMAND_NAME, allowed_states, model_inputs, ordered_states


def add_presentation_commands(discrete_evolution_subparsers: Any) -> None:
    discrete_render = discrete_evolution_subparsers.add_parser(
        "render",
        help="Render a tree annotated with reconstructed geographic or other discrete states.",
    )
    discrete_render.add_argument("tree", type=Path)
    discrete_render.add_argument("table", type=Path)
    discrete_render.add_argument("--trait", required=True)
    discrete_render.add_argument("--taxon-column")
    discrete_render.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    discrete_render.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    discrete_render.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    discrete_render.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_render.add_argument(
        "--layout", choices=("cladogram", "phylogram", "circular"), default="phylogram"
    )
    discrete_render.add_argument("--out", required=True, type=Path)
    discrete_render.add_argument(
        "--json", action="store_true", help="Emit the render result as JSON."
    )
    _add_manifest_argument(discrete_render)

    discrete_report = discrete_evolution_subparsers.add_parser(
        "report",
        help="Render an HTML report for one discrete-state evolution analysis.",
    )
    discrete_report.add_argument("tree", type=Path)
    discrete_report.add_argument("table", type=Path)
    discrete_report.add_argument("--trait", required=True)
    discrete_report.add_argument("--taxon-column")
    discrete_report.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    discrete_report.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    discrete_report.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    discrete_report.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_report.add_argument(
        "--compare-model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
    )
    discrete_report.add_argument("--out", required=True, type=Path)
    discrete_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(discrete_report)


def run_presentation_command(args: Any) -> int | None:
    if args.discrete_evolution_command == "render":
        return _run_render(args)
    if args.discrete_evolution_command == "report":
        return _run_report(args)
    return None


def _run_render(args: Any) -> int:
    report = estimate_ancestral_geographic_states(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        model=args.model,
        allowed_states=allowed_states(args),
        state_ordering=args.state_ordering,
        ordered_states=ordered_states(args),
    )
    result = render_tree_with_geographic_states(
        args.tree,
        report,
        out_path=args.out,
        layout=args.layout,
    )
    outputs = _finalize_outputs(
        args,
        command=COMMAND_NAME,
        inputs=model_inputs(args),
        outputs=[result.output_path],
    )
    _print_result(
        build_command_result(
            command=COMMAND_NAME,
            inputs=model_inputs(args),
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "tip_count": result.tip_count,
                "rendered_internal_annotation_count": (
                    result.rendered_internal_annotation_count
                ),
                "layout": result.layout,
                "model": report.model,
                "state_ordering": report.state_ordering,
            },
            data={"reconstruction": report, "render": result},
        ),
        json_output=args.json,
    )
    return 0


def _run_report(args: Any) -> int:
    result = render_discrete_state_evolution_report(
        tree_path=args.tree,
        traits_path=args.table,
        trait=args.trait,
        out_path=args.out,
        taxon_column=args.taxon_column,
        model=args.model,
        allowed_states=allowed_states(args),
        state_ordering=args.state_ordering,
        ordered_states=ordered_states(args),
        compare_model=args.compare_model,
    )
    outputs = _finalize_outputs(
        args,
        command=COMMAND_NAME,
        inputs=model_inputs(args),
        outputs=[result.output_path, args.out.with_suffix(".svg")],
    )
    _print_result(
        build_command_result(
            command=COMMAND_NAME,
            inputs=model_inputs(args),
            outputs=outputs,
            metrics={
                "report_kind": result.report_kind,
                "model": result.model,
                "state_ordering": args.state_ordering,
            },
            data=result,
        ),
        json_output=args.json,
    )
    return 0
