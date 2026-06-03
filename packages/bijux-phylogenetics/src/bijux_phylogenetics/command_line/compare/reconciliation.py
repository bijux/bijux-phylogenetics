from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.compare.topology import (
    reconcile_duplication_loss_transfer,
    write_duplication_loss_transfer_event_table,
    write_duplication_loss_transfer_taxon_map_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def run_compare_duplication_loss_transfer_command(
    args: Any, *, parser: argparse.ArgumentParser
) -> int | None:
    if args.left != "duplication-loss-transfer":
        return None
    if args.third is None:
        parser.exit(
            status=2,
            message=(
                "compare duplication-loss-transfer requires a species tree and a gene tree\n"
            ),
        )

    species_tree_path = Path(args.right)
    gene_tree_path = Path(args.third)
    report = reconcile_duplication_loss_transfer(
        species_tree_path,
        gene_tree_path,
        taxon_map_path=args.taxon_map,
        duplication_cost=args.duplication_cost,
        loss_cost=args.loss_cost,
        transfer_cost=args.transfer_cost,
    )
    output_paths: list[Path | str] = []
    if args.out is not None:
        output_paths.append(
            write_duplication_loss_transfer_event_table(
                args.out,
                species_tree_path,
                gene_tree_path,
                taxon_map_path=args.taxon_map,
                duplication_cost=args.duplication_cost,
                loss_cost=args.loss_cost,
                transfer_cost=args.transfer_cost,
            )
        )
    if args.mapping_out is not None:
        output_paths.append(
            write_duplication_loss_transfer_taxon_map_table(
                args.mapping_out,
                species_tree_path,
                gene_tree_path,
                taxon_map_path=args.taxon_map,
                duplication_cost=args.duplication_cost,
                loss_cost=args.loss_cost,
                transfer_cost=args.transfer_cost,
            )
        )
    inputs: list[Path | str] = [species_tree_path, gene_tree_path]
    if args.taxon_map is not None:
        inputs.append(args.taxon_map)
    outputs = _finalize_outputs(
        args,
        command="compare",
        inputs=inputs,
        outputs=output_paths,
    )
    _print_result(
        build_command_result(
            command="compare",
            inputs=inputs,
            outputs=outputs,
            metrics={
                "observed_species_taxa": len(report.observed_species_taxa),
                "species_only_taxa": len(report.species_only_taxa),
                "gene_tip_count": report.gene_tip_count,
                "reconciliation_score": report.reconciliation_score,
                "duplication_event_count": report.duplication_event_count,
                "loss_event_count": report.loss_event_count,
                "transfer_event_count": report.transfer_event_count,
                "speciation_event_count": report.speciation_event_count,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
