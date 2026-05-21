from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.comparative import (
    summarize_geiger_birth_death_exclusion,
    summarize_medusa_exclusion,
)
from bijux_phylogenetics.runtime.errors import DiversificationAnalysisError


def add_diversification_exclusion_commands(diversification_subparsers: Any) -> None:
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


def run_diversification_exclusion_command(args: Any) -> int | None:
    if args.diversification_command == "medusa":
        _raise_medusa_exclusion(args)
    if args.diversification_command == "bd-ms":
        _raise_birth_death_exclusion(args)
    return None


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
                None
                if report.sampling_report is None
                else report.sampling_report.complete
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
                None
                if report.sampling_report is None
                else report.sampling_report.complete
            ),
        },
    )
