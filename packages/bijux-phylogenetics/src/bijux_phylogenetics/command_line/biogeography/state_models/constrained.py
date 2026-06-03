from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result

from ..shared import command_line_api


def add_constrained_geography_command(biogeography_subparsers: Any) -> None:
    biogeography_constrained = biogeography_subparsers.add_parser(
        "constrained",
        help="Compare constrained and unconstrained geographic fits under an explicit region adjacency matrix.",
    )
    biogeography_constrained.add_argument("tree", type=Path)
    biogeography_constrained.add_argument("table", type=Path)
    biogeography_constrained.add_argument("adjacency", type=Path)
    biogeography_constrained.add_argument("--trait", required=True)
    biogeography_constrained.add_argument("--taxon-column")
    biogeography_constrained.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="ard",
    )
    biogeography_constrained.add_argument("--summary-out", type=Path)
    biogeography_constrained.add_argument("--fits-out", type=Path)
    biogeography_constrained.add_argument("--transitions-out", type=Path)
    biogeography_constrained.add_argument("--unsupported-out", type=Path)
    biogeography_constrained.add_argument("--exclusions-out", type=Path)
    biogeography_constrained.add_argument(
        "--json", action="store_true", help="Emit the biogeography review as JSON."
    )
    _add_manifest_argument(biogeography_constrained)


def run_constrained_geography_command(args: Any) -> int | None:
    if args.biogeography_command != "constrained":
        return None

    cli_api = command_line_api()
    report = cli_api.summarize_constrained_geographic_model(
        args.tree,
        args.table,
        args.adjacency,
        trait=args.trait,
        taxon_column=args.taxon_column,
        model=args.model,
    )
    summary = cli_api.summarize_constrained_geographic_report(report)
    outputs: list[Path | str] = []
    if args.summary_out is not None:
        outputs.append(
            cli_api.write_constrained_geographic_summary_table(
                args.summary_out,
                report,
            )
        )
    if args.fits_out is not None:
        outputs.append(
            cli_api.write_constrained_geographic_fit_table(
                args.fits_out,
                report,
            )
        )
    if args.transitions_out is not None:
        outputs.append(
            cli_api.write_constrained_geographic_transition_table(
                args.transitions_out,
                report,
            )
        )
    if args.unsupported_out is not None:
        outputs.append(
            cli_api.write_unsupported_geographic_transition_claim_table(
                args.unsupported_out,
                report,
            )
        )
    if args.exclusions_out is not None:
        outputs.append(
            cli_api.write_constrained_geographic_exclusion_table(
                args.exclusions_out,
                report,
            )
        )
    outputs = _finalize_outputs(
        args,
        command="biogeography",
        inputs=[args.tree, args.table, args.adjacency],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="biogeography",
            inputs=[args.tree, args.table, args.adjacency],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "model": report.model,
                "allowed_transition_count": summary.allowed_transition_count,
                "forbidden_transition_count": summary.forbidden_transition_count,
                "unsupported_transition_claim_count": (
                    summary.unsupported_transition_claim_count
                ),
                "preferred_constraint": summary.preferred_constraint,
                "excluded_taxon_count": summary.excluded_taxon_count,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
