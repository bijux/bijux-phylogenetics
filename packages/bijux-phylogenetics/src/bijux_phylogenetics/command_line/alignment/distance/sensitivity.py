from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.distance import (
    compare_distance_gap_policies,
    compare_distance_models,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .shared import (
    add_ambiguity_policy_option,
    add_distance_model_option,
    add_gap_handling_option,
)


def add_distance_sensitivity_commands(alignment_subparsers: Any) -> None:
    alignment_distance_models = alignment_subparsers.add_parser(
        "distance-models",
        help="Compare all supported distance models on the same alignment.",
    )
    alignment_distance_models.add_argument("alignment", type=Path)
    add_gap_handling_option(alignment_distance_models)
    add_ambiguity_policy_option(alignment_distance_models)
    alignment_distance_models.add_argument(
        "--json", action="store_true", help="Emit the model comparison as JSON."
    )
    _add_manifest_argument(alignment_distance_models)

    alignment_distance_gap = alignment_subparsers.add_parser(
        "distance-gap-sensitivity",
        help="Compare pairwise versus complete deletion for the same distance workflow.",
    )
    alignment_distance_gap.add_argument("alignment", type=Path)
    add_distance_model_option(alignment_distance_gap)
    add_ambiguity_policy_option(alignment_distance_gap)
    alignment_distance_gap.add_argument(
        "--json", action="store_true", help="Emit the gap-policy sensitivity as JSON."
    )
    _add_manifest_argument(alignment_distance_gap)


def run_distance_sensitivity_command(args: Any) -> int | None:
    if args.alignment_command == "distance-models":
        report = compare_distance_models(
            args.alignment,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "model_count": len(report.rows),
                    "alphabet": report.inferred_alphabet,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.alignment_command == "distance-gap-sensitivity":
        report = compare_distance_gap_policies(
            args.alignment,
            model=args.model,
            ambiguity_policy=args.ambiguity_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "changed_pair_count": report.changed_pair_count,
                    "pair_count": report.pair_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None
