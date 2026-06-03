from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative import (
    compare_diversification_models,
    estimate_diversification_rate,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .inputs import tree_and_metadata_inputs


def add_diversification_modeling_commands(diversification_subparsers: Any) -> None:
    diversification_estimate = diversification_subparsers.add_parser(
        "estimate",
        help="Estimate a simple Yule or birth-death diversification model.",
    )
    diversification_estimate.add_argument("tree", type=Path)
    diversification_estimate.add_argument("--metadata", type=Path)
    diversification_estimate.add_argument("--taxon-column")
    diversification_estimate.add_argument("--sampling-column")
    diversification_estimate.add_argument(
        "--model", choices=("yule", "birth-death"), default="birth-death"
    )
    diversification_estimate.add_argument(
        "--json", action="store_true", help="Emit the diversification estimate as JSON."
    )
    _add_manifest_argument(diversification_estimate)

    diversification_compare = diversification_subparsers.add_parser(
        "compare-models",
        help="Compare Yule and birth-death diversification fits.",
    )
    diversification_compare.add_argument("tree", type=Path)
    diversification_compare.add_argument("--metadata", type=Path)
    diversification_compare.add_argument("--taxon-column")
    diversification_compare.add_argument("--sampling-column")
    diversification_compare.add_argument(
        "--json", action="store_true", help="Emit the model comparison as JSON."
    )
    _add_manifest_argument(diversification_compare)


def run_diversification_modeling_command(args: Any) -> int | None:
    if args.diversification_command == "estimate":
        return _run_estimate(args)
    if args.diversification_command == "compare-models":
        return _run_compare_models(args)
    return None


def _run_estimate(args: Any) -> int:
    inputs = tree_and_metadata_inputs(args)
    report = estimate_diversification_rate(
        args.tree,
        metadata_path=args.metadata,
        taxon_column=args.taxon_column,
        sampling_column=args.sampling_column,
        model=args.model,
    )
    outputs = _finalize_outputs(args, command="diversification", inputs=inputs)
    _print_result(
        build_command_result(
            command="diversification",
            inputs=inputs,
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "model": report.model,
                "sampling_fraction": report.sampling_fraction,
                "net_diversification_rate": report.net_diversification_rate,
                "aic": report.aic,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def _run_compare_models(args: Any) -> int:
    inputs = tree_and_metadata_inputs(args)
    report = compare_diversification_models(
        args.tree,
        metadata_path=args.metadata,
        taxon_column=args.taxon_column,
        sampling_column=args.sampling_column,
    )
    outputs = _finalize_outputs(args, command="diversification", inputs=inputs)
    _print_result(
        build_command_result(
            command="diversification",
            inputs=inputs,
            outputs=outputs,
            metrics={
                "better_model": report.better_model,
                "model_count": len(report.rows),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
