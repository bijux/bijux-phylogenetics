from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.comparative import (
    summarize_geiger_birth_death_exclusion,
    summarize_medusa_exclusion,
)
from bijux_phylogenetics.runtime.errors import DiversificationAnalysisError

from .clades import (
    add_diversification_clade_command,
    run_diversification_clade_command,
)
from .inspection import (
    add_diversification_inspection_commands,
    run_diversification_inspection_command,
)
from .modeling import (
    add_diversification_modeling_commands,
    run_diversification_modeling_command,
)
from .presentation import (
    add_diversification_presentation_commands,
    run_diversification_presentation_command,
)
from .trait_dependence import (
    add_diversification_trait_dependence_command,
    run_diversification_trait_dependence_command,
)


def add_diversification_commands(subparsers: Any) -> None:
    diversification = subparsers.add_parser(
        get_command_spec("diversification").name,
        help=get_command_spec("diversification").summary,
    )
    diversification_subparsers = diversification.add_subparsers(
        dest="diversification_command",
        required=True,
    )
    add_diversification_inspection_commands(diversification_subparsers)
    add_diversification_modeling_commands(diversification_subparsers)
    add_diversification_clade_command(diversification_subparsers)
    add_diversification_trait_dependence_command(diversification_subparsers)
    add_diversification_presentation_commands(diversification_subparsers)
    _add_medusa_exclusion_command(diversification_subparsers)
    _add_birth_death_exclusion_command(diversification_subparsers)


def run_diversification_command(args: Any) -> int:
    inspection_exit_code = run_diversification_inspection_command(args)
    if inspection_exit_code is not None:
        return inspection_exit_code

    modeling_exit_code = run_diversification_modeling_command(args)
    if modeling_exit_code is not None:
        return modeling_exit_code

    clade_exit_code = run_diversification_clade_command(args)
    if clade_exit_code is not None:
        return clade_exit_code

    trait_dependence_exit_code = run_diversification_trait_dependence_command(args)
    if trait_dependence_exit_code is not None:
        return trait_dependence_exit_code

    presentation_exit_code = run_diversification_presentation_command(args)
    if presentation_exit_code is not None:
        return presentation_exit_code

    if args.diversification_command == "medusa":
        _raise_medusa_exclusion(args)
    if args.diversification_command == "bd-ms":
        _raise_birth_death_exclusion(args)

    raise NotImplementedError(
        f"unsupported diversification command: {args.diversification_command}"
    )


def _add_medusa_exclusion_command(diversification_subparsers: Any) -> None:
    diversification_medusa = diversification_subparsers.add_parser(
        "medusa",
        help="Explain the explicit exclusion boundary for geiger::medusa parity.",
    )
    diversification_medusa.add_argument("tree", type=Path)
    diversification_medusa.add_argument("--metadata", type=Path)
    diversification_medusa.add_argument("--taxon-column")
    diversification_medusa.add_argument("--sampling-column")
    diversification_medusa.add_argument(
        "--json", action="store_true", help="Emit the MEDUSA exclusion as JSON."
    )
    _add_manifest_argument(diversification_medusa)


def _add_birth_death_exclusion_command(diversification_subparsers: Any) -> None:
    diversification_bd_ms = diversification_subparsers.add_parser(
        "bd-ms",
        help="Explain the explicit exclusion boundary for geiger::bd.ms birth-death parity.",
    )
    diversification_bd_ms.add_argument("tree", type=Path)
    diversification_bd_ms.add_argument("--metadata", type=Path)
    diversification_bd_ms.add_argument("--taxon-column")
    diversification_bd_ms.add_argument("--sampling-column")
    diversification_bd_ms.add_argument(
        "--json",
        action="store_true",
        help="Emit the birth-death exclusion as JSON.",
    )
    _add_manifest_argument(diversification_bd_ms)


def _raise_medusa_exclusion(args: Any) -> None:
    report = summarize_medusa_exclusion(
        args.tree,
        metadata_path=args.metadata,
        taxon_column=args.taxon_column,
        sampling_column=args.sampling_column,
    )
    raise DiversificationAnalysisError(
        report.exclusion_reason,
        code="diversification_medusa_explicitly_excluded",
        details={
            "failure_reason": report.exclusion_code,
            "supported_surfaces": report.supported_surfaces,
            "missing_surfaces": report.missing_surfaces,
            "tip_count": report.validation.tip_count,
            "rooted": report.validation.rooted,
            "ultrametric": report.validation.ultrametric,
            "sampling_metadata_complete": (
                None if report.sampling_report is None else report.sampling_report.complete
            ),
        },
    )


def _raise_birth_death_exclusion(args: Any) -> None:
    report = summarize_geiger_birth_death_exclusion(
        args.tree,
        metadata_path=args.metadata,
        taxon_column=args.taxon_column,
        sampling_column=args.sampling_column,
    )
    raise DiversificationAnalysisError(
        report.exclusion_reason,
        code="diversification_birth_death_explicitly_excluded",
        details={
            "failure_reason": report.exclusion_code,
            "geiger_reference_surface": report.geiger_reference_surface,
            "geiger_reference_arguments": report.geiger_reference_arguments,
            "owned_surface": report.owned_surface,
            "tip_count": report.validation.tip_count,
            "rooted": report.validation.rooted,
            "ultrametric": report.validation.ultrametric,
            "sampling_metadata_complete": (
                None if report.sampling_report is None else report.sampling_report.complete
            ),
        },
    )


__all__ = [
    "add_diversification_commands",
    "run_diversification_command",
]
