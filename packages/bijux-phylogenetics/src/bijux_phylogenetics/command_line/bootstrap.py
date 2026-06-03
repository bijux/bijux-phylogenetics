from __future__ import annotations

import argparse
import sys
from typing import Any

__all__ = ["build_parser", "main", "run_command"]

from bijux_phylogenetics import __version__
from bijux_phylogenetics.command_line.adapters import (
    add_adapter_commands,
    run_adapter_command,
)
from bijux_phylogenetics.command_line.alignment import (
    add_alignment_commands,
    run_alignment_command,
)
from bijux_phylogenetics.command_line.ancestral import (
    add_ancestral_commands,
    run_ancestral_command,
)
from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _json_requested,
)
from bijux_phylogenetics.command_line.benchmark import (
    add_benchmark_commands,
    run_benchmark_command,
)
from bijux_phylogenetics.command_line.biogeography import (
    add_biogeography_commands,
    run_biogeography_command,
)
from bijux_phylogenetics.command_line.comparative import (
    add_comparative_commands,
    run_comparative_command,
)
from bijux_phylogenetics.command_line.compare import (
    add_compare_command,
    run_compare_command,
)
from bijux_phylogenetics.command_line.demo import (
    add_demo_command,
    run_demo_command,
)
from bijux_phylogenetics.command_line.diagnose import (
    add_diagnose_command,
    run_diagnose_command,
)
from bijux_phylogenetics.command_line.discrete_evolution import (
    add_discrete_evolution_commands,
    run_discrete_evolution_command,
)
from bijux_phylogenetics.command_line.distance import (
    add_distance_commands,
    run_distance_command,
)
from bijux_phylogenetics.command_line.diversification import (
    add_diversification_commands,
    run_diversification_command,
)
from bijux_phylogenetics.command_line.ecological_niche import (
    add_ecological_niche_commands,
    run_ecological_niche_command,
)
from bijux_phylogenetics.command_line.engines import (
    add_phylo_commands,
    run_phylo_command,
)
from bijux_phylogenetics.command_line.evidence import (
    add_evidence_command,
    run_evidence_command,
)
from bijux_phylogenetics.command_line.host_association import (
    add_host_association_commands,
    run_host_association_command,
)
from bijux_phylogenetics.command_line.output import (
    _print_commands,
    _print_result,
)
from bijux_phylogenetics.command_line.parity import (
    add_parity_command,
    run_parity_command,
)
from bijux_phylogenetics.command_line.phylogeography import (
    add_phylogeography_commands,
    run_phylogeography_command,
)
from bijux_phylogenetics.command_line.prune import (
    add_prune_command,
    run_prune_command,
)
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.render import (
    add_render_command,
    run_render_command,
)
from bijux_phylogenetics.command_line.report import (
    add_report_command,
    run_report_command,
)
from bijux_phylogenetics.command_line.routing import (
    _command_inputs,
    _finalize_outputs,
)
from bijux_phylogenetics.command_line.simulate import (
    add_simulate_command,
    run_simulate_command,
)
from bijux_phylogenetics.command_line.study_inputs.annotation import (
    add_annotate_command,
    run_annotate_command,
)
from bijux_phylogenetics.command_line.study_inputs.metadata import (
    add_metadata_commands,
    run_metadata_command,
)
from bijux_phylogenetics.command_line.study_inputs.traits import (
    add_traits_commands,
    run_traits_command,
)
from bijux_phylogenetics.command_line.taxonomy import (
    add_taxonomy_commands,
    run_taxonomy_command,
)
from bijux_phylogenetics.command_line.tree import (
    add_topology_commands,
    add_tree_inspection_commands,
    add_tree_normalization_commands,
    add_tree_set_commands,
    run_inspect_command,
    run_normalize_command,
    run_normalize_taxa_command,
    run_topology_command,
    run_tree_set_command,
    run_validate_command,
)
from bijux_phylogenetics.core.environment import inspect_environment
from bijux_phylogenetics.runtime.errors import PhylogeneticsError
from bijux_phylogenetics.runtime.results import build_command_result, build_error_result


def build_parser() -> argparse.ArgumentParser:
    """Build the repository CLI parser."""
    parser = argparse.ArgumentParser(prog="bijux-phylogenetics")
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    commands = subparsers.add_parser(
        "commands", help="List the registered command taxonomy."
    )
    commands.add_argument("--format", choices=("text", "json"), default="text")

    env = subparsers.add_parser(
        get_command_spec("env").name, help=get_command_spec("env").summary
    )
    env_subparsers = env.add_subparsers(dest="env_command", required=True)
    env_inspect = env_subparsers.add_parser(
        "inspect", help="Inspect runtime dependency availability."
    )
    env_inspect.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(env_inspect)

    add_phylo_commands(subparsers)

    add_metadata_commands(subparsers)

    add_traits_commands(subparsers)

    add_prune_command(subparsers)

    add_alignment_commands(subparsers)

    add_comparative_commands(subparsers)

    add_ancestral_commands(subparsers)

    add_biogeography_commands(subparsers)
    add_host_association_commands(subparsers)
    add_ecological_niche_commands(subparsers)
    add_phylogeography_commands(subparsers)
    add_discrete_evolution_commands(subparsers)
    add_diversification_commands(subparsers)

    add_distance_commands(subparsers)
    add_tree_set_commands(subparsers)

    add_simulate_command(subparsers)

    add_benchmark_commands(subparsers)

    add_parity_command(subparsers)

    add_tree_inspection_commands(subparsers)

    add_tree_normalization_commands(subparsers)

    add_taxonomy_commands(subparsers)

    add_topology_commands(subparsers)

    add_compare_command(subparsers)

    add_annotate_command(subparsers)

    add_diagnose_command(subparsers)

    add_render_command(subparsers)

    add_evidence_command(subparsers)

    add_report_command(subparsers)

    add_demo_command(subparsers)

    add_adapter_commands(subparsers)

    return parser


def run_command(args: Any, *, parser: argparse.ArgumentParser) -> int:
    """Run the selected command."""
    try:
        if args.command == "commands":
            _print_commands(output_format=args.format)
            return 0
        if args.command == "env":
            report = inspect_environment()
            outputs = _finalize_outputs(args, command="env", inputs=[])
            _print_result(
                build_command_result(
                    command="env",
                    inputs=[],
                    outputs=outputs,
                    metrics={"dependency_count": len(report.dependencies)},
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "phylo":
            return run_phylo_command(args)
        if args.command == "metadata":
            return run_metadata_command(args)
        if args.command == "traits":
            return run_traits_command(args)
        if args.command == "validate":
            return run_validate_command(args)
        if args.command == "prune":
            return run_prune_command(args)
        if args.command == "alignment":
            alignment_exit_code = run_alignment_command(args)
            if alignment_exit_code is not None:
                return alignment_exit_code
        if args.command == "comparative":
            return run_comparative_command(args, parser=parser)
        if args.command == "ancestral":
            return run_ancestral_command(args, parser=parser)
        if args.command == "biogeography":
            return run_biogeography_command(args)
        if args.command == "host-association":
            return run_host_association_command(args)
        if args.command == "ecological-niche":
            return run_ecological_niche_command(args)
        if args.command == "phylogeography":
            return run_phylogeography_command(args)
        if args.command == "discrete-evolution":
            return run_discrete_evolution_command(args)
        if args.command == "diversification":
            return run_diversification_command(args)
        if args.command == "distance":
            return run_distance_command(args)
        if args.command == "tree-set":
            return run_tree_set_command(args)
        if args.command == "simulate":
            return run_simulate_command(args, parser=parser)
        if args.command == "benchmark":
            return run_benchmark_command(args)
        if args.command == "parity":
            return run_parity_command(args)
        if args.command == "inspect":
            return run_inspect_command(args)
        if args.command == "normalize":
            return run_normalize_command(args)
        if args.command == "normalize-taxa":
            return run_normalize_taxa_command(args)
        if args.command == "taxonomy":
            return run_taxonomy_command(args)
        if args.command == "topology":
            return run_topology_command(args)
        if args.command == "diagnose":
            return run_diagnose_command(args, parser=parser)
        if args.command == "compare":
            return run_compare_command(args, parser=parser)
        if args.command == "annotate":
            return run_annotate_command(args)
        if args.command == "render":
            return run_render_command(args)
        if args.command == "evidence":
            return run_evidence_command(args)
        if args.command == "demo":
            return run_demo_command(args)
        if args.command == "report":
            return run_report_command(args)
        adapter_exit_code = run_adapter_command(args)
        if adapter_exit_code is not None:
            return adapter_exit_code
    except PhylogeneticsError as error:
        if _json_requested(args):
            _print_result(
                build_error_result(
                    command=args.command, inputs=_command_inputs(args), error=error
                ),
                json_output=True,
            )
            return 2
        parser.exit(status=2, message=f"{error.code}: {error.message}\n")
    except FileNotFoundError as error:
        parser.exit(status=2, message=f"{error}\n")
    except ValueError as error:
        parser.exit(status=2, message=f"{error}\n")
    except NotImplementedError as error:
        parser.exit(status=2, message=f"{error}\n")
    except Exception as error:  # pragma: no cover - defensive CLI guard
        parser.exit(status=1, message=f"unexpected error: {error}\n")

    parser.print_help(sys.stderr)
    return 2


def main(argv: list[str] | None = None) -> int:
    """Run the canonical phylogenetics command-line interface."""
    parser = build_parser()
    args = parser.parse_args(argv)
    args._argv = list(argv) if argv is not None else list(sys.argv[1:])
    return run_command(args, parser=parser)
