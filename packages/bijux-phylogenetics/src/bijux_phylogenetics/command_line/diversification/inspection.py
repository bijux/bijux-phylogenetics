from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative import (
    compute_diversification_gamma_statistic,
    compute_lineage_through_time_curve,
    detect_incomplete_taxon_sampling_metadata,
    write_diversification_gamma_statistic_table,
    write_lineage_through_time_table,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .inputs import tree_and_metadata_inputs


def add_diversification_inspection_commands(diversification_subparsers: Any) -> None:
    diversification_ltt = diversification_subparsers.add_parser(
        "ltt",
        help="Compute a lineage-through-time curve for one rooted ultrametric tree.",
    )
    diversification_ltt.add_argument("tree", type=Path)
    diversification_ltt.add_argument(
        "--out", type=Path, help="Write the lineage-through-time table as TSV."
    )
    diversification_ltt.add_argument(
        "--json", action="store_true", help="Emit the LTT report as JSON."
    )
    _add_manifest_argument(diversification_ltt)

    diversification_sampling = diversification_subparsers.add_parser(
        "sampling",
        help="Inspect taxon sampling-fraction metadata against the tree tips.",
    )
    diversification_sampling.add_argument("tree", type=Path)
    diversification_sampling.add_argument("table", type=Path)
    diversification_sampling.add_argument("--taxon-column")
    diversification_sampling.add_argument("--sampling-column")
    diversification_sampling.add_argument(
        "--json", action="store_true", help="Emit the sampling report as JSON."
    )
    _add_manifest_argument(diversification_sampling)

    diversification_gamma = diversification_subparsers.add_parser(
        "gamma-stat",
        help="Compute the Pybus-Harvey diversification gamma statistic.",
    )
    diversification_gamma.add_argument("tree", type=Path)
    diversification_gamma.add_argument("--metadata", type=Path)
    diversification_gamma.add_argument("--taxon-column")
    diversification_gamma.add_argument("--sampling-column")
    diversification_gamma.add_argument(
        "--out",
        type=Path,
        help="Write the diversification gamma-statistic table as TSV.",
    )
    diversification_gamma.add_argument(
        "--json",
        action="store_true",
        help="Emit the diversification gamma-statistic report as JSON.",
    )
    _add_manifest_argument(diversification_gamma)


def run_diversification_inspection_command(args: Any) -> int | None:
    if args.diversification_command == "ltt":
        return _run_ltt(args)
    if args.diversification_command == "sampling":
        return _run_sampling(args)
    if args.diversification_command == "gamma-stat":
        return _run_gamma_stat(args)
    return None


def _run_ltt(args: Any) -> int:
    report = compute_lineage_through_time_curve(args.tree)
    outputs: list[Path | str] = []
    if args.out is not None:
        outputs.append(write_lineage_through_time_table(args.out, report))
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
            metrics={
                "tip_count": report.tip_count,
                "root_age": report.root_age,
                "point_count": len(report.points),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def _run_sampling(args: Any) -> int:
    report = detect_incomplete_taxon_sampling_metadata(
        args.tree,
        args.table,
        taxon_column=args.taxon_column,
        sampling_column=args.sampling_column,
    )
    inputs = [args.tree, args.table]
    outputs = _finalize_outputs(args, command="diversification", inputs=inputs)
    _print_result(
        build_command_result(
            command="diversification",
            inputs=inputs,
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "complete": report.complete,
                "matched_taxon_count": len(report.matched_taxa),
                "missing_taxon_count": len(report.missing_taxa),
                "invalid_row_count": len(report.invalid_rows),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def _run_gamma_stat(args: Any) -> int:
    inputs = tree_and_metadata_inputs(args)
    report = compute_diversification_gamma_statistic(
        args.tree,
        metadata_path=args.metadata,
        taxon_column=args.taxon_column,
        sampling_column=args.sampling_column,
    )
    outputs: list[Path | str] = []
    if args.out is not None:
        outputs.append(write_diversification_gamma_statistic_table(args.out, report))
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
                "tip_count": report.tip_count,
                "branching_time_count": report.branching_time_count,
                "gamma_statistic": report.gamma_statistic,
                "sampling_fraction": report.sampling_fraction,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
