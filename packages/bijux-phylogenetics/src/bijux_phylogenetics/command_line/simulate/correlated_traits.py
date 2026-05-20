from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _parse_float_csv_row,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.simulation import (
    simulate_correlated_brownian_trait_collection,
    write_correlated_continuous_trait_collection_summary_table,
    write_correlated_continuous_trait_collection_table,
)


def add_simulate_correlated_trait_commands(simulate_subparsers: Any) -> None:
    simulate_brownian_correlated = simulate_subparsers.add_parser(
        "traits-brownian-correlated",
        help="Simulate correlated continuous tip traits under multivariate Brownian motion.",
    )
    simulate_brownian_correlated.add_argument("tree", type=Path)
    simulate_brownian_correlated.add_argument(
        "--trait",
        action="append",
        required=True,
        help="One trait name. Repeat to define the multivariate trait order.",
    )
    simulate_brownian_correlated.add_argument(
        "--root-state",
        action="append",
        type=float,
        default=[],
        help="One root state value aligned to the declared trait order.",
    )
    covariance_group = simulate_brownian_correlated.add_mutually_exclusive_group(
        required=True
    )
    covariance_group.add_argument(
        "--covariance-row",
        action="append",
        dest="covariance_rows",
        help="One comma-separated covariance-matrix row. Repeat to build the full matrix.",
    )
    covariance_group.add_argument(
        "--correlation-row",
        action="append",
        dest="correlation_rows",
        help="One comma-separated correlation-matrix row. Repeat to build the full matrix.",
    )
    simulate_brownian_correlated.add_argument(
        "--trait-standard-deviation",
        action="append",
        type=float,
        default=[],
        help="One trait standard deviation aligned to the declared trait order when using --correlation-row.",
    )
    simulate_brownian_correlated.add_argument("--replicates", type=int, default=128)
    simulate_brownian_correlated.add_argument("--seed", type=int, default=1)
    simulate_brownian_correlated.add_argument("--out", required=True, type=Path)
    simulate_brownian_correlated.add_argument("--summary-out", type=Path)
    simulate_brownian_correlated.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_brownian_correlated)


def run_simulate_correlated_trait_command(args: Any, *, parser: Any) -> int | None:
    if args.simulate_command != "traits-brownian-correlated":
        return None
    if args.covariance_rows is None and args.correlation_rows is None:
        parser.error(
            "correlated Brownian simulation requires either --covariance-row or --correlation-row"
        )
    if args.covariance_rows is not None and args.trait_standard_deviation:
        parser.error(
            "trait standard deviations can only be used with --correlation-row"
        )
    if args.correlation_rows is not None and not args.trait_standard_deviation:
        parser.error(
            "correlated Brownian simulation requires --trait-standard-deviation with --correlation-row"
        )

    report = simulate_correlated_brownian_trait_collection(
        args.tree,
        trait_names=args.trait,
        evolutionary_covariance_matrix=(
            None
            if args.covariance_rows is None
            else [_parse_float_csv_row(raw_row) for raw_row in args.covariance_rows]
        ),
        evolutionary_correlation_matrix=(
            None
            if args.correlation_rows is None
            else [_parse_float_csv_row(raw_row) for raw_row in args.correlation_rows]
        ),
        trait_standard_deviations=(
            None if not args.trait_standard_deviation else args.trait_standard_deviation
        ),
        root_states=None if not args.root_state else args.root_state,
        replicates=args.replicates,
        seed=args.seed,
    )
    outputs: list[Path | str] = [
        write_correlated_continuous_trait_collection_table(args.out, report)
    ]
    if args.summary_out is not None:
        outputs.append(
            write_correlated_continuous_trait_collection_summary_table(
                args.summary_out,
                report,
            )
        )
    outputs = _finalize_outputs(
        args,
        command="simulate",
        inputs=[args.tree],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="simulate",
            inputs=[args.tree],
            outputs=outputs,
            metrics={
                "tip_count": report.tip_count,
                "trait_count": len(report.trait_names),
                "replicate_count": report.replicate_count,
                "summary_row_count": len(report.rows),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
