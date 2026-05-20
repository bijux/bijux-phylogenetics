from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.simulation import (
    simulate_brownian_traits,
    simulate_early_burst_traits,
    simulate_ou_traits,
    simulate_speciational_traits,
    write_continuous_trait_table,
)


def add_simulate_continuous_trait_commands(simulate_subparsers: Any) -> None:
    simulate_brownian = simulate_subparsers.add_parser(
        "traits-brownian",
        help="Simulate a continuous tip trait under Brownian motion.",
    )
    simulate_brownian.add_argument("tree", type=Path)
    simulate_brownian.add_argument("--root-state", type=float, default=0.0)
    simulate_brownian.add_argument("--sigma", type=float)
    simulate_brownian.add_argument("--sigma-squared", type=float)
    simulate_brownian.add_argument("--seed", type=int, default=1)
    simulate_brownian.add_argument("--out", required=True, type=Path)
    simulate_brownian.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_brownian)

    simulate_speciational = simulate_subparsers.add_parser(
        "traits-speciational",
        help="Simulate a continuous tip trait under geiger-style speciational Brownian motion.",
    )
    simulate_speciational.add_argument("tree", type=Path)
    simulate_speciational.add_argument("--root-state", type=float, default=0.0)
    simulate_speciational.add_argument("--sigma", type=float)
    simulate_speciational.add_argument("--sigma-squared", type=float)
    simulate_speciational.add_argument("--seed", type=int, default=1)
    simulate_speciational.add_argument("--out", required=True, type=Path)
    simulate_speciational.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_speciational)

    simulate_ou = simulate_subparsers.add_parser(
        "traits-ou",
        help="Simulate a continuous tip trait under an OU process.",
    )
    simulate_ou.add_argument("tree", type=Path)
    simulate_ou.add_argument("--root-state", type=float, default=0.0)
    simulate_ou.add_argument("--sigma", type=float, default=1.0)
    simulate_ou.add_argument("--alpha", type=float, default=1.0)
    simulate_ou.add_argument("--theta", type=float, default=0.0)
    simulate_ou.add_argument("--seed", type=int, default=1)
    simulate_ou.add_argument("--out", required=True, type=Path)
    simulate_ou.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_ou)

    simulate_early_burst = simulate_subparsers.add_parser(
        "traits-early-burst",
        help="Simulate a continuous tip trait under an early-burst branch-rate process.",
    )
    simulate_early_burst.add_argument("tree", type=Path)
    simulate_early_burst.add_argument("--root-state", type=float, default=0.0)
    simulate_early_burst.add_argument("--sigma", type=float, default=1.0)
    simulate_early_burst.add_argument("--rate-change", type=float, default=1.0)
    simulate_early_burst.add_argument("--seed", type=int, default=1)
    simulate_early_burst.add_argument("--out", required=True, type=Path)
    simulate_early_burst.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_early_burst)


def run_simulate_continuous_trait_command(args: Any) -> int | None:
    if args.simulate_command == "traits-brownian":
        report = simulate_brownian_traits(
            args.tree,
            root_state=args.root_state,
            sigma=args.sigma,
            sigma_squared=args.sigma_squared,
            seed=args.seed,
        )
        _print_continuous_trait_result(
            args,
            report=report,
            metrics={
                "tip_count": report.tip_count,
                "trait_count": len(report.traits),
                "sigma_squared": report.sigma_squared,
            },
        )
        return 0
    if args.simulate_command == "traits-speciational":
        report = simulate_speciational_traits(
            args.tree,
            root_state=args.root_state,
            sigma=args.sigma,
            sigma_squared=args.sigma_squared,
            seed=args.seed,
        )
        _print_continuous_trait_result(
            args,
            report=report,
            metrics={
                "tip_count": report.tip_count,
                "trait_count": len(report.traits),
                "sigma_squared": report.sigma_squared,
            },
        )
        return 0
    if args.simulate_command == "traits-ou":
        report = simulate_ou_traits(
            args.tree,
            root_state=args.root_state,
            sigma=args.sigma,
            alpha=args.alpha,
            theta=args.theta,
            seed=args.seed,
        )
        _print_continuous_trait_result(
            args,
            report=report,
            metrics={
                "tip_count": report.tip_count,
                "trait_count": len(report.traits),
            },
        )
        return 0
    if args.simulate_command != "traits-early-burst":
        return None

    report = simulate_early_burst_traits(
        args.tree,
        root_state=args.root_state,
        sigma=args.sigma,
        rate_change=args.rate_change,
        seed=args.seed,
    )
    _print_continuous_trait_result(
        args,
        report=report,
        metrics={
            "tip_count": report.tip_count,
            "trait_count": len(report.traits),
            "rate_change": report.rate_change,
        },
    )
    return 0


def _print_continuous_trait_result(
    args: Any,
    *,
    report: Any,
    metrics: dict[str, Any],
) -> None:
    output_path = write_continuous_trait_table(args.out, report)
    outputs = _finalize_outputs(
        args, command="simulate", inputs=[args.tree], outputs=[output_path]
    )
    _print_result(
        build_command_result(
            command="simulate",
            inputs=[args.tree],
            outputs=outputs,
            metrics=metrics,
            data=report,
        ),
        json_output=args.json,
    )
