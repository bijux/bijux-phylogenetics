from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.distance import (
    assess_distance_method_assumptions,
    inspect_distance_matrix_quality,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .shared import (
    add_ambiguity_policy_option,
    add_distance_model_option,
    add_gap_handling_option,
)


def add_distance_diagnostic_commands(alignment_subparsers: Any) -> None:
    alignment_distance_quality = alignment_subparsers.add_parser(
        "distance-quality",
        help="Inspect saturation, divergence, and low-information risks in a computed distance matrix.",
    )
    alignment_distance_quality.add_argument("alignment", type=Path)
    add_distance_model_option(alignment_distance_quality)
    add_gap_handling_option(alignment_distance_quality)
    add_ambiguity_policy_option(alignment_distance_quality)
    alignment_distance_quality.add_argument(
        "--json", action="store_true", help="Emit the diagnostics as JSON."
    )
    _add_manifest_argument(alignment_distance_quality)

    alignment_distance_suitability = alignment_subparsers.add_parser(
        "distance-suitability",
        help="Emit the explicit suitability decision for distance-method use on one alignment.",
    )
    alignment_distance_suitability.add_argument("alignment", type=Path)
    add_distance_model_option(alignment_distance_suitability)
    add_gap_handling_option(alignment_distance_suitability)
    add_ambiguity_policy_option(alignment_distance_suitability)
    alignment_distance_suitability.add_argument(
        "--json", action="store_true", help="Emit the suitability decision as JSON."
    )
    _add_manifest_argument(alignment_distance_suitability)

    alignment_distance_assumptions = alignment_subparsers.add_parser(
        "distance-assumptions",
        help="Audit NJ and UPGMA assumptions, including UPGMA ultrametric compatibility.",
    )
    alignment_distance_assumptions.add_argument("alignment", type=Path)
    add_distance_model_option(alignment_distance_assumptions)
    add_gap_handling_option(alignment_distance_assumptions)
    add_ambiguity_policy_option(alignment_distance_assumptions)
    alignment_distance_assumptions.add_argument(
        "--json", action="store_true", help="Emit the assumption audit as JSON."
    )
    _add_manifest_argument(alignment_distance_assumptions)


def run_distance_diagnostic_command(args: Any) -> int | None:
    if args.alignment_command == "distance-quality":
        report = inspect_distance_matrix_quality(
            args.alignment,
            model=args.model,
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
                    "taxon_count": report.taxon_count,
                    "saturated_pair_count": len(report.saturated_pairs),
                    "low_information_pair_count": len(report.low_information_pairs),
                    "decision": report.method_assessment.decision,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.alignment_command == "distance-suitability":
        report = inspect_distance_matrix_quality(
            args.alignment,
            model=args.model,
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
                warnings=report.method_assessment.reasons,
                metrics={
                    "decision": report.method_assessment.decision,
                    "reason_count": len(report.method_assessment.reasons),
                },
                data=report.method_assessment,
            ),
            json_output=args.json,
        )
        return 0

    if args.alignment_command == "distance-assumptions":
        report = assess_distance_method_assumptions(
            args.alignment,
            model=args.model,
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
                    "taxon_count": report.taxon_count,
                    "ultrametric_compatible": report.ultrametric_compatible,
                    "upgma_violation_count": len(report.upgma_ultrametric_violations),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None
