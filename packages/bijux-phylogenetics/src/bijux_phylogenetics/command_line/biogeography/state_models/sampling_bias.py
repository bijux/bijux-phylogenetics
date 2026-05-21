from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result

from ..shared import command_line_api


def add_sampling_bias_geography_command(biogeography_subparsers: Any) -> None:
    biogeography_sampling_bias = biogeography_subparsers.add_parser(
        "sampling-bias",
        help="Reweight geographic region sampling and compare weighted versus unweighted ancestral conclusions.",
    )
    biogeography_sampling_bias.add_argument("tree", type=Path)
    biogeography_sampling_bias.add_argument("table", type=Path)
    biogeography_sampling_bias.add_argument("--trait", required=True)
    biogeography_sampling_bias.add_argument("--taxon-column")
    biogeography_sampling_bias.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    biogeography_sampling_bias.add_argument(
        "--allowed-regions",
        help="Comma-delimited explicit region vocabulary.",
    )
    biogeography_sampling_bias.add_argument(
        "--weights",
        type=Path,
        help="Optional region-weight table with explicit region and weight columns.",
    )
    biogeography_sampling_bias.add_argument(
        "--region-column",
        default="region",
        help="Region column in the optional region-weight table.",
    )
    biogeography_sampling_bias.add_argument(
        "--weight-column",
        default="weight",
        help="Numeric weight column in the optional region-weight table.",
    )
    biogeography_sampling_bias.add_argument("--summary-out", type=Path)
    biogeography_sampling_bias.add_argument("--regions-out", type=Path)
    biogeography_sampling_bias.add_argument("--nodes-out", type=Path)
    biogeography_sampling_bias.add_argument("--transitions-out", type=Path)
    biogeography_sampling_bias.add_argument("--exclusions-out", type=Path)
    biogeography_sampling_bias.add_argument(
        "--json", action="store_true", help="Emit the biogeography review as JSON."
    )
    _add_manifest_argument(biogeography_sampling_bias)


def run_sampling_bias_geography_command(args: Any) -> int | None:
    if args.biogeography_command != "sampling-bias":
        return None

    cli_api = command_line_api()
    report = cli_api.summarize_geographic_sampling_bias(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        model=args.model,
        allowed_regions=cli_api._split_csv_values(args.allowed_regions) or None,
        weights_path=args.weights,
        region_column=args.region_column,
        weight_column=args.weight_column,
    )
    outputs: list[Path | str] = []
    if args.summary_out is not None:
        outputs.append(
            cli_api.write_geographic_sampling_bias_summary_table(
                args.summary_out,
                report,
            )
        )
    if args.regions_out is not None:
        outputs.append(
            cli_api.write_geographic_sampling_count_table(
                args.regions_out,
                report,
            )
        )
    if args.nodes_out is not None:
        outputs.append(
            cli_api.write_geographic_sampling_bias_node_table(
                args.nodes_out,
                report,
            )
        )
    if args.transitions_out is not None:
        outputs.append(
            cli_api.write_geographic_sampling_bias_transition_table(
                args.transitions_out,
                report,
            )
        )
    if args.exclusions_out is not None:
        outputs.append(
            cli_api.write_geographic_sampling_bias_exclusion_table(
                args.exclusions_out,
                report,
            )
        )
    inputs = (
        [args.tree, args.table]
        if args.weights is None
        else [args.tree, args.table, args.weights]
    )
    outputs = _finalize_outputs(
        args,
        command="biogeography",
        inputs=inputs,
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="biogeography",
            inputs=inputs,
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "model": report.summary.model,
                "weighting_mode": report.summary.weighting_mode,
                "region_dominated": report.summary.region_dominated,
                "dominant_region": report.summary.dominant_region,
                "dominant_region_fraction": (report.summary.dominant_region_fraction),
                "root_region_changed": report.summary.root_region_changed,
                "changed_internal_node_count": (
                    report.summary.changed_internal_node_count
                ),
                "changed_transition_count": report.summary.changed_transition_count,
                "excluded_taxon_count": report.summary.excluded_taxon_count,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
