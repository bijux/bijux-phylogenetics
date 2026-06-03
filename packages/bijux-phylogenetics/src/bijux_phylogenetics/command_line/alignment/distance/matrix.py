from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.distance import (
    compute_pairwise_genetic_distance_matrix,
    write_genetic_distance_component_table,
    write_genetic_distance_matrix,
    write_genetic_distance_parameter_table,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .shared import (
    add_ambiguity_policy_option,
    add_distance_model_option,
    add_gap_handling_option,
)


def add_distance_matrix_command(alignment_subparsers: Any) -> None:
    alignment_distance = alignment_subparsers.add_parser(
        "distance-matrix",
        help="Compute a pairwise DNA genetic distance matrix.",
    )
    alignment_distance.add_argument("alignment", type=Path)
    add_distance_model_option(alignment_distance)
    add_gap_handling_option(alignment_distance)
    add_ambiguity_policy_option(alignment_distance)
    alignment_distance.add_argument("--out", type=Path, help="Write the matrix as TSV.")
    alignment_distance.add_argument(
        "--components-out",
        type=Path,
        help="Write pairwise mismatch, transition, and transversion components as TSV.",
    )
    alignment_distance.add_argument(
        "--parameters-out",
        type=Path,
        help="Write alignment-wide distance-model parameters as TSV.",
    )
    alignment_distance.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_distance)


def run_distance_matrix_command(args: Any) -> int | None:
    if args.alignment_command != "distance-matrix":
        return None

    report = compute_pairwise_genetic_distance_matrix(
        args.alignment,
        model=args.model,
        gap_handling=args.gap_handling,
        ambiguity_policy=args.ambiguity_policy,
    )
    outputs: list[Path | str] = []
    if args.out is not None:
        outputs.append(write_genetic_distance_matrix(args.out, report))
    if args.components_out is not None:
        outputs.append(
            write_genetic_distance_component_table(
                args.components_out,
                report,
            )
        )
    if args.parameters_out is not None:
        outputs.append(
            write_genetic_distance_parameter_table(
                args.parameters_out,
                report,
            )
        )
    warnings = list(report.warnings)
    if any(
        pair.saturated
        for pair in report.pairs
        if pair.left_identifier != pair.right_identifier
    ):
        warnings.append(
            "one or more pairwise distances are saturated or non-finite under the selected model"
        )
    outputs = _finalize_outputs(
        args,
        command="alignment",
        inputs=[args.alignment],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="alignment",
            inputs=[args.alignment],
            outputs=outputs,
            warnings=warnings,
            metrics={
                "sequence_count": len(report.identifiers),
                "pair_count": len(report.pairs),
                "model": report.model,
                "gap_handling": report.gap_handling,
                "ambiguity_policy": report.ambiguity_policy,
                "alphabet": report.inferred_alphabet,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
