from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.compare.topology import (
    count_deep_coalescences,
    write_deep_coalescence_branch_table,
    write_deep_coalescence_taxon_map_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def run_compare_deep_coalescence_command(
    args: Any, *, parser: argparse.ArgumentParser
) -> int | None:
    if args.left != "deep-coalescence":
        return None
    if args.third is None:
        parser.exit(
            status=2,
            message="compare deep-coalescence requires a species tree and a gene tree\n",
        )

    species_tree_path = Path(args.right)
    gene_tree_path = Path(args.third)
    report = count_deep_coalescences(
        species_tree_path,
        gene_tree_path,
        taxon_map_path=args.taxon_map,
    )
    output_paths: list[Path | str] = []
    if args.out is not None:
        output_paths.append(
            write_deep_coalescence_branch_table(
                args.out,
                species_tree_path,
                gene_tree_path,
                taxon_map_path=args.taxon_map,
            )
        )
    if args.mapping_out is not None:
        output_paths.append(
            write_deep_coalescence_taxon_map_table(
                args.mapping_out,
                species_tree_path,
                gene_tree_path,
                taxon_map_path=args.taxon_map,
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
                "deep_coalescence_total": report.deep_coalescence_total,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
