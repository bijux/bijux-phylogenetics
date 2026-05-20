from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.discrete_evolution import (
    compare_discrete_state_models,
    count_discrete_stochastic_map_transitions,
    estimate_ancestral_geographic_states,
    load_stochastic_map_collection,
    render_discrete_state_evolution_report,
    render_stochastic_map_density_artifact,
    render_tree_with_geographic_states,
    simulate_discrete_stochastic_maps,
    summarize_discrete_stochastic_map_density,
    summarize_discrete_stochastic_maps,
    write_discrete_model_comparison_table,
    write_node_state_probability_table,
    write_stochastic_map_aggregate_transition_matrix,
    write_stochastic_map_branch_occupancy_table,
    write_stochastic_map_branch_probability_table,
    write_stochastic_map_branch_transition_count_table,
    write_stochastic_map_collection,
    write_stochastic_map_density_branch_table,
    write_stochastic_map_density_slice_table,
    write_stochastic_map_event_table,
    write_stochastic_map_segment_table,
    write_stochastic_map_state_time_table,
    write_stochastic_map_summary_table,
    write_stochastic_map_transition_count_matrix,
    write_transition_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .shared import COMMAND_NAME, allowed_states, model_inputs, ordered_states, render_density_outputs
from .validation import add_validation_commands, run_validation_command


def add_discrete_evolution_commands(subparsers: Any) -> None:
    discrete_evolution = subparsers.add_parser(
        get_command_spec("discrete-evolution").name,
        help=get_command_spec("discrete-evolution").summary,
    )
    discrete_evolution_subparsers = discrete_evolution.add_subparsers(
        dest="discrete_evolution_command",
        required=True,
    )
    add_validation_commands(discrete_evolution_subparsers)

    discrete_model = discrete_evolution_subparsers.add_parser(
        "model",
        help="Run one discrete-state transition model and export node or branch summaries.",
    )
    discrete_model.add_argument("tree", type=Path)
    discrete_model.add_argument("table", type=Path)
    discrete_model.add_argument("--trait", required=True)
    discrete_model.add_argument("--taxon-column")
    discrete_model.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    discrete_model.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    discrete_model.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    discrete_model.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_model.add_argument(
        "--node-table-out", type=Path, help="Write node-state probabilities as TSV."
    )
    discrete_model.add_argument(
        "--transitions-out", type=Path, help="Write branch transition summaries as TSV."
    )
    discrete_model.add_argument(
        "--json", action="store_true", help="Emit the model report as JSON."
    )
    _add_manifest_argument(discrete_model)

    discrete_compare = discrete_evolution_subparsers.add_parser(
        "compare-models",
        help="Compare two supported discrete-state evolution models node by node.",
    )
    discrete_compare.add_argument("tree", type=Path)
    discrete_compare.add_argument("table", type=Path)
    discrete_compare.add_argument("--trait", required=True)
    discrete_compare.add_argument("--taxon-column")
    discrete_compare.add_argument(
        "--left-model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    discrete_compare.add_argument(
        "--right-model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="all-rates-different",
    )
    discrete_compare.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    discrete_compare.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    discrete_compare.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_compare.add_argument(
        "--table-out", type=Path, help="Write node-wise model differences as TSV."
    )
    discrete_compare.add_argument(
        "--json", action="store_true", help="Emit the comparison report as JSON."
    )
    _add_manifest_argument(discrete_compare)

    discrete_stochastic = discrete_evolution_subparsers.add_parser(
        "stochastic-map",
        help="Generate seeded stochastic character maps from a fitted discrete-state CTMC.",
    )
    discrete_stochastic.add_argument("tree", type=Path)
    discrete_stochastic.add_argument("table", type=Path)
    discrete_stochastic.add_argument("--trait", required=True)
    discrete_stochastic.add_argument("--taxon-column")
    discrete_stochastic.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    discrete_stochastic.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    discrete_stochastic.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    discrete_stochastic.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_stochastic.add_argument("--replicates", type=int, default=100)
    discrete_stochastic.add_argument("--seed", type=int, default=0)
    discrete_stochastic.add_argument(
        "--collection-out", type=Path, help="Write stochastic maps as JSON."
    )
    discrete_stochastic.add_argument(
        "--summary-out", type=Path, help="Write stochastic-map summary as TSV."
    )
    discrete_stochastic.add_argument(
        "--state-times-out",
        type=Path,
        help="Write per-state time-in-state summaries as TSV.",
    )
    discrete_stochastic.add_argument(
        "--branch-occupancy-out",
        type=Path,
        help="Write per-branch state-occupancy summaries as TSV.",
    )
    discrete_stochastic.add_argument(
        "--count-matrix-out",
        type=Path,
        help="Write one per-replicate transition count matrix as TSV.",
    )
    discrete_stochastic.add_argument(
        "--aggregate-matrix-out",
        type=Path,
        help="Write one aggregate mean transition matrix as TSV.",
    )
    discrete_stochastic.add_argument(
        "--branch-transition-out",
        type=Path,
        help="Write per-branch transition-count summaries as TSV.",
    )
    discrete_stochastic.add_argument(
        "--segments-out",
        type=Path,
        help="Write flat branch-state segment rows as TSV.",
    )
    discrete_stochastic.add_argument(
        "--events-out",
        type=Path,
        help="Write flat stochastic transition-event rows as TSV.",
    )
    discrete_stochastic.add_argument(
        "--focal-state",
        help="Resolve one focal state for density slices and density rendering. When omitted, binary collections default to the second state in the fitted state order.",
    )
    discrete_stochastic.add_argument(
        "--density-resolution",
        type=int,
        default=100,
        help="Set the branch-slice resolution used for density summaries.",
    )
    discrete_stochastic.add_argument(
        "--branch-probabilities-out",
        type=Path,
        help="Write per-branch state-probability summaries as TSV.",
    )
    discrete_stochastic.add_argument(
        "--density-branches-out",
        type=Path,
        help="Write per-branch focal-state density summaries as TSV.",
    )
    discrete_stochastic.add_argument(
        "--density-slices-out",
        type=Path,
        help="Write flat branch-slice density rows as TSV.",
    )
    discrete_stochastic.add_argument(
        "--density-figure-out",
        type=Path,
        help="Write one branch-colored density artifact as .svg or .html.",
    )
    discrete_stochastic.add_argument(
        "--layout",
        choices=("phylogram", "cladogram", "circular"),
        default="phylogram",
        help="Choose the layout for any density artifact written from this command.",
    )
    discrete_stochastic.add_argument(
        "--json",
        action="store_true",
        help="Emit the stochastic-map collection as JSON.",
    )
    _add_manifest_argument(discrete_stochastic)

    discrete_summarize_maps = discrete_evolution_subparsers.add_parser(
        "summarize-maps",
        help="Summarize a previously written stochastic-map collection.",
    )
    discrete_summarize_maps.add_argument("input_path", type=Path)
    discrete_summarize_maps.add_argument(
        "--summary-out", type=Path, help="Write stochastic-map summary as TSV."
    )
    discrete_summarize_maps.add_argument(
        "--state-times-out",
        type=Path,
        help="Write per-state time-in-state summaries as TSV.",
    )
    discrete_summarize_maps.add_argument(
        "--branch-occupancy-out",
        type=Path,
        help="Write per-branch state-occupancy summaries as TSV.",
    )
    discrete_summarize_maps.add_argument(
        "--json", action="store_true", help="Emit the stochastic-map summary as JSON."
    )
    _add_manifest_argument(discrete_summarize_maps)

    discrete_count_maps = discrete_evolution_subparsers.add_parser(
        "count-maps",
        help="Count directional transitions in a previously written stochastic-map collection.",
    )
    discrete_count_maps.add_argument("input_path", type=Path)
    discrete_count_maps.add_argument(
        "--count-matrix-out",
        type=Path,
        help="Write one per-replicate transition count matrix as TSV.",
    )
    discrete_count_maps.add_argument(
        "--aggregate-matrix-out",
        type=Path,
        help="Write one aggregate mean transition matrix as TSV.",
    )
    discrete_count_maps.add_argument(
        "--branch-transition-out",
        type=Path,
        help="Write per-branch transition-count summaries as TSV.",
    )
    discrete_count_maps.add_argument(
        "--events-out",
        type=Path,
        help="Write flat stochastic transition-event rows as TSV.",
    )
    discrete_count_maps.add_argument(
        "--json",
        action="store_true",
        help="Emit the stochastic-map count report as JSON.",
    )
    _add_manifest_argument(discrete_count_maps)

    discrete_density_maps = discrete_evolution_subparsers.add_parser(
        "density-maps",
        help="Summarize posterior density over a previously written stochastic-map collection.",
    )
    discrete_density_maps.add_argument("input_path", type=Path)
    discrete_density_maps.add_argument(
        "--focal-state",
        help="Resolve one focal state for density slices and density rendering. When omitted, binary collections default to the second state in the fitted state order.",
    )
    discrete_density_maps.add_argument(
        "--resolution",
        type=int,
        default=100,
        help="Set the branch-slice resolution used for density summaries.",
    )
    discrete_density_maps.add_argument(
        "--branch-probabilities-out",
        type=Path,
        help="Write per-branch state-probability summaries as TSV.",
    )
    discrete_density_maps.add_argument(
        "--density-branches-out",
        type=Path,
        help="Write per-branch focal-state density summaries as TSV.",
    )
    discrete_density_maps.add_argument(
        "--density-slices-out",
        type=Path,
        help="Write flat branch-slice density rows as TSV.",
    )
    discrete_density_maps.add_argument(
        "--out",
        type=Path,
        help="Write one branch-colored density artifact as .svg or .html.",
    )
    discrete_density_maps.add_argument(
        "--layout",
        choices=("phylogram", "cladogram", "circular"),
        default="phylogram",
    )
    discrete_density_maps.add_argument(
        "--json",
        action="store_true",
        help="Emit the stochastic-map density report as JSON.",
    )
    _add_manifest_argument(discrete_density_maps)

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


def _run_model(args: Any) -> int:
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
    outputs: list[Path | str] = []
    if args.node_table_out is not None:
        outputs.append(write_node_state_probability_table(args.node_table_out, report))
    if args.transitions_out is not None:
        outputs.append(write_transition_summary_table(args.transitions_out, report))
    outputs = _finalize_outputs(
        args,
        command="discrete-evolution",
        inputs=model_inputs(args),
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="discrete-evolution",
            inputs=model_inputs(args),
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "taxon_count": report.taxon_count,
                "observed_state_count": len(report.observed_states),
                "transition_count": report.transition_summary.transition_count,
                "strongly_supported_transition_count": (
                    report.transition_summary.strongly_supported_transition_count
                ),
                "model": report.model,
                "state_ordering": report.state_ordering,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def _run_stochastic_map(args: Any) -> int:
    report = simulate_discrete_stochastic_maps(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        model=args.model,
        allowed_states=allowed_states(args),
        state_ordering=args.state_ordering,
        ordered_states=ordered_states(args),
        replicates=args.replicates,
        seed=args.seed,
    )
    count_report = count_discrete_stochastic_map_transitions(report)
    density_report = None
    density_result = None
    if any(
        output is not None
        for output in (
            args.branch_probabilities_out,
            args.density_branches_out,
            args.density_slices_out,
            args.density_figure_out,
        )
    ):
        density_report = summarize_discrete_stochastic_map_density(
            report,
            resolution=args.density_resolution,
            focal_state=args.focal_state,
        )
    outputs: list[Path | str] = []
    if args.collection_out is not None:
        outputs.append(write_stochastic_map_collection(args.collection_out, report))
    if args.summary_out is not None:
        outputs.append(write_stochastic_map_summary_table(args.summary_out, report.summary))
    if args.state_times_out is not None:
        outputs.append(
            write_stochastic_map_state_time_table(args.state_times_out, report.summary)
        )
    if args.branch_occupancy_out is not None:
        outputs.append(
            write_stochastic_map_branch_occupancy_table(
                args.branch_occupancy_out,
                report.summary,
            )
        )
    if args.count_matrix_out is not None:
        outputs.append(
            write_stochastic_map_transition_count_matrix(
                args.count_matrix_out,
                count_report,
            )
        )
    if args.aggregate_matrix_out is not None:
        outputs.append(
            write_stochastic_map_aggregate_transition_matrix(
                args.aggregate_matrix_out,
                count_report,
            )
        )
    if args.branch_transition_out is not None:
        outputs.append(
            write_stochastic_map_branch_transition_count_table(
                args.branch_transition_out,
                count_report,
            )
        )
    if args.segments_out is not None:
        outputs.append(write_stochastic_map_segment_table(args.segments_out, report))
    if args.events_out is not None:
        outputs.append(write_stochastic_map_event_table(args.events_out, report))
    if density_report is not None and args.branch_probabilities_out is not None:
        outputs.append(
            write_stochastic_map_branch_probability_table(
                args.branch_probabilities_out,
                density_report,
            )
        )
    if density_report is not None and args.density_branches_out is not None:
        outputs.append(
            write_stochastic_map_density_branch_table(
                args.density_branches_out,
                density_report,
            )
        )
    if density_report is not None and args.density_slices_out is not None:
        outputs.append(
            write_stochastic_map_density_slice_table(
                args.density_slices_out,
                density_report,
            )
        )
    if density_report is not None and args.density_figure_out is not None:
        density_result = render_stochastic_map_density_artifact(
            density_report,
            tree_path=report.tree_path,
            out_path=args.density_figure_out,
            layout=args.layout,
        )
        render_density_outputs(outputs, density_result)
    outputs = _finalize_outputs(
        args,
        command="discrete-evolution",
        inputs=model_inputs(args),
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="discrete-evolution",
            inputs=model_inputs(args),
            outputs=outputs,
            warnings=list(
                dict.fromkeys(
                    [
                        *report.warnings,
                        *count_report.warnings,
                        *([] if density_report is None else density_report.warnings),
                    ]
                )
            ),
            metrics={
                "requested_replicate_count": report.replicates,
                "successful_replicate_count": report.summary.replicate_count,
                "simulation_failure_count": report.summary.simulation_failure_count,
                "mean_total_transition_count": count_report.mean_total_transition_count,
                "branch_state_row_count": len(report.summary.branch_occupancy_rows),
                "count_matrix_row_count": len(count_report.matrix_rows),
                "branch_transition_row_count": len(count_report.branch_rows),
                "conditioned_on_node_estimates": report.conditioned_on_node_estimates,
                "parameter_count": report.fit_audit.parameter_count,
                "optimizer_converged": report.fit_audit.optimizer_converged,
                "overparameterized": report.fit_audit.overparameterized,
                "preferred_model_by_aic": report.fit_audit.preferred_model_by_aic,
                "fit_warning_count": len(report.fit_audit.warnings),
                "density_branch_row_count": (
                    0 if density_report is None else len(density_report.branch_rows)
                ),
                "density_slice_row_count": (
                    0 if density_report is None else len(density_report.density_rows)
                ),
                "rendered_branch_color_count": (
                    0
                    if density_result is None
                    else density_result.rendered_branch_color_count
                ),
            },
            data={
                "collection": report,
                "counts": count_report,
                "density": density_report,
                "density_render": density_result,
            },
        ),
        json_output=args.json,
    )
    return 0


def _run_summarize_maps(args: Any) -> int:
    collection = load_stochastic_map_collection(args.input_path)
    report = summarize_discrete_stochastic_maps(collection)
    outputs: list[Path | str] = []
    if args.summary_out is not None:
        outputs.append(write_stochastic_map_summary_table(args.summary_out, report))
    if args.state_times_out is not None:
        outputs.append(write_stochastic_map_state_time_table(args.state_times_out, report))
    if args.branch_occupancy_out is not None:
        outputs.append(
            write_stochastic_map_branch_occupancy_table(args.branch_occupancy_out, report)
        )
    outputs = _finalize_outputs(
        args,
        command="discrete-evolution",
        inputs=[args.input_path],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="discrete-evolution",
            inputs=[args.input_path],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "replicate_count": report.replicate_count,
                "mean_total_transition_count": report.mean_total_transition_count,
                "simulation_failure_count": report.simulation_failure_count,
                "branch_state_row_count": len(report.branch_occupancy_rows),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def _run_count_maps(args: Any) -> int:
    collection = load_stochastic_map_collection(args.input_path)
    report = count_discrete_stochastic_map_transitions(collection)
    outputs: list[Path | str] = []
    if args.count_matrix_out is not None:
        outputs.append(
            write_stochastic_map_transition_count_matrix(
                args.count_matrix_out,
                report,
            )
        )
    if args.aggregate_matrix_out is not None:
        outputs.append(
            write_stochastic_map_aggregate_transition_matrix(
                args.aggregate_matrix_out,
                report,
            )
        )
    if args.branch_transition_out is not None:
        outputs.append(
            write_stochastic_map_branch_transition_count_table(
                args.branch_transition_out,
                report,
            )
        )
    if args.events_out is not None:
        outputs.append(write_stochastic_map_event_table(args.events_out, collection))
    outputs = _finalize_outputs(
        args,
        command="discrete-evolution",
        inputs=[args.input_path],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="discrete-evolution",
            inputs=[args.input_path],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "replicate_count": report.replicate_count,
                "mean_total_transition_count": report.mean_total_transition_count,
                "count_matrix_row_count": len(report.matrix_rows),
                "branch_transition_row_count": len(report.branch_rows),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def _run_density_maps(args: Any) -> int:
    collection = load_stochastic_map_collection(args.input_path)
    report = summarize_discrete_stochastic_map_density(
        collection,
        resolution=args.resolution,
        focal_state=args.focal_state,
    )
    outputs: list[Path | str] = []
    density_result = None
    if args.branch_probabilities_out is not None:
        outputs.append(
            write_stochastic_map_branch_probability_table(
                args.branch_probabilities_out,
                report,
            )
        )
    if args.density_branches_out is not None:
        outputs.append(
            write_stochastic_map_density_branch_table(
                args.density_branches_out,
                report,
            )
        )
    if args.density_slices_out is not None:
        outputs.append(
            write_stochastic_map_density_slice_table(
                args.density_slices_out,
                report,
            )
        )
    if args.out is not None:
        density_result = render_stochastic_map_density_artifact(
            report,
            tree_path=collection.tree_path,
            out_path=args.out,
            layout=args.layout,
        )
        render_density_outputs(outputs, density_result)
    outputs = _finalize_outputs(
        args,
        command="discrete-evolution",
        inputs=[args.input_path],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="discrete-evolution",
            inputs=[args.input_path],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "replicate_count": report.replicate_count,
                "resolution": report.resolution,
                "branch_probability_row_count": len(report.branch_state_rows),
                "density_branch_row_count": len(report.branch_rows),
                "density_slice_row_count": len(report.density_rows),
                "focal_state": report.focal_state,
                "baseline_state": report.baseline_state,
                "rendered_branch_color_count": (
                    0
                    if density_result is None
                    else density_result.rendered_branch_color_count
                ),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


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
        command="discrete-evolution",
        inputs=model_inputs(args),
        outputs=[result.output_path],
    )
    _print_result(
        build_command_result(
            command="discrete-evolution",
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
        command="discrete-evolution",
        inputs=model_inputs(args),
        outputs=[result.output_path, args.out.with_suffix(".svg")],
    )
    _print_result(
        build_command_result(
            command="discrete-evolution",
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


def _run_compare_models(args: Any) -> int:
    comparison = compare_discrete_state_models(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        left_model=args.left_model,
        right_model=args.right_model,
        allowed_states=allowed_states(args),
        state_ordering=args.state_ordering,
        ordered_states=ordered_states(args),
    )
    outputs: list[Path | str] = []
    if args.table_out is not None:
        outputs.append(write_discrete_model_comparison_table(args.table_out, comparison))
    outputs = _finalize_outputs(
        args,
        command="discrete-evolution",
        inputs=model_inputs(args),
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="discrete-evolution",
            inputs=model_inputs(args),
            outputs=outputs,
            metrics={
                "better_model": comparison.better_model,
                "model_count": len(comparison.rows),
                "differing_node_count": sum(
                    1 for row in comparison.node_differences if row.differs
                ),
                "state_ordering": args.state_ordering,
            },
            data=comparison,
        ),
        json_output=args.json,
    )
    return 0


def run_discrete_evolution_command(args: Any) -> int:
    validation_result = run_validation_command(args)
    if validation_result is not None:
        return validation_result
    if args.discrete_evolution_command == "model":
        return _run_model(args)
    if args.discrete_evolution_command == "stochastic-map":
        return _run_stochastic_map(args)
    if args.discrete_evolution_command == "summarize-maps":
        return _run_summarize_maps(args)
    if args.discrete_evolution_command == "count-maps":
        return _run_count_maps(args)
    if args.discrete_evolution_command == "density-maps":
        return _run_density_maps(args)
    if args.discrete_evolution_command == "render":
        return _run_render(args)
    if args.discrete_evolution_command == "report":
        return _run_report(args)
    return _run_compare_models(args)
