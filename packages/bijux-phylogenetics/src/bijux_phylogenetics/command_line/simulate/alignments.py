from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.simulation import (
    simulate_dna_alignment,
    simulate_protein_alignment,
    write_simulated_alignment,
)


def add_simulate_alignment_commands(simulate_subparsers: Any) -> None:
    simulate_dna = simulate_subparsers.add_parser(
        "alignment-dna",
        help="Simulate a DNA alignment along a rooted tree.",
    )
    simulate_dna.add_argument("tree", type=Path)
    simulate_dna.add_argument("--sequence-length", type=int, required=True)
    simulate_dna.add_argument("--substitution-rate", type=float, default=1.0)
    simulate_dna.add_argument("--seed", type=int, default=1)
    simulate_dna.add_argument("--out", required=True, type=Path)
    simulate_dna.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_dna)

    simulate_protein = simulate_subparsers.add_parser(
        "alignment-protein",
        help="Simulate a protein alignment along a rooted tree.",
    )
    simulate_protein.add_argument("tree", type=Path)
    simulate_protein.add_argument("--sequence-length", type=int, required=True)
    simulate_protein.add_argument("--substitution-rate", type=float, default=1.0)
    simulate_protein.add_argument("--seed", type=int, default=1)
    simulate_protein.add_argument("--out", required=True, type=Path)
    simulate_protein.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_protein)


def run_simulate_alignment_command(args: Any) -> int | None:
    if args.simulate_command == "alignment-dna":
        report = simulate_dna_alignment(
            args.tree,
            sequence_length=args.sequence_length,
            substitution_rate=args.substitution_rate,
            seed=args.seed,
        )
        _print_alignment_result(args, report=report)
        return 0
    if args.simulate_command != "alignment-protein":
        return None

    report = simulate_protein_alignment(
        args.tree,
        sequence_length=args.sequence_length,
        substitution_rate=args.substitution_rate,
        seed=args.seed,
    )
    _print_alignment_result(args, report=report)
    return 0


def _print_alignment_result(args: Any, *, report: Any) -> None:
    output_path = write_simulated_alignment(args.out, report)
    outputs = _finalize_outputs(
        args, command="simulate", inputs=[args.tree], outputs=[output_path]
    )
    _print_result(
        build_command_result(
            command="simulate",
            inputs=[args.tree],
            outputs=outputs,
            metrics={
                "tip_count": report.tip_count,
                "sequence_length": report.sequence_length,
            },
            data=report,
        ),
        json_output=args.json,
    )
