from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import (
    join_table_to_taxa,
    load_taxon_table,
)
from bijux_phylogenetics.diagnostics.validation.structure import _load_tree
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.phylo.alignment import AlignmentSummary

from .models import TableLinkageReport


def summarise_alignment_path(path: Path) -> AlignmentSummary:
    """Expose FASTA alignment summary for external callers."""
    return summarise_fasta(path)


def annotate_tree_against_table(
    tree_path: Path,
    table_path: Path,
    *,
    taxon_column: str | None = None,
) -> TableLinkageReport:
    """Summarise how a TSV table links against tree tips."""
    tree = _load_tree(tree_path)
    table = load_taxon_table(table_path, taxon_column=taxon_column)
    join = join_table_to_taxa(tree.tip_names, table_path, taxon_column=taxon_column)
    annotated_taxa = [row.taxon for row in join.joined_rows if row.matched]
    return TableLinkageReport(
        tree_path=tree_path,
        table_path=table_path,
        tree_taxa=tree.tip_count,
        table_rows=table.row_count,
        linked_taxa=len(annotated_taxa),
        missing_from_table=join.missing_from_metadata,
        extra_table_entries=join.extra_metadata_taxa,
        index_column=table.index_column,
        annotated_taxa=annotated_taxa,
        joined_rows=join.joined_rows,
    )


def write_annotation_report(path: Path, report: TableLinkageReport) -> Path:
    """Write a linkage report to a deterministic JSON artifact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path
