from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from bijux_phylogenetics.core.alignment import AlignmentSummary
from bijux_phylogenetics.core.metadata import load_taxon_table
from bijux_phylogenetics.diagnostics.validation import TreeInspectionReport, TreeValidationReport, inspect_tree_path, validate_tree_path
from bijux_phylogenetics.diagnostics.validation import _load_tree
from bijux_phylogenetics.io.fasta import summarise_fasta
from bijux_phylogenetics.render.html import write_html_report


@dataclass(slots=True)
class TableLinkageReport:
    tree_path: Path
    table_path: Path
    tree_taxa: int
    table_rows: int
    linked_taxa: int
    missing_from_table: list[str]
    extra_table_entries: list[str]
    index_column: str
    annotated_taxa: list[str]


@dataclass(slots=True)
class ReportBuildResult:
    output_path: Path
    validation: TreeValidationReport
    inspection: TreeInspectionReport
    metadata_linkage: TableLinkageReport | None
    traits_linkage: TableLinkageReport | None
    alignment: AlignmentSummary | None


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
    table = load_taxon_table(table_path, taxon_column=taxon_column)
    full_tip_names = set(_load_tree(tree_path).tip_names)
    missing = sorted(full_tip_names - table.indexed_values)
    extras = sorted(table.indexed_values - full_tip_names)
    annotated_taxa = sorted(full_tip_names & table.indexed_values)
    linked = len(annotated_taxa)
    return TableLinkageReport(
        tree_path=tree_path,
        table_path=table_path,
        tree_taxa=len(full_tip_names),
        table_rows=table.row_count,
        linked_taxa=linked,
        missing_from_table=missing,
        extra_table_entries=extras,
        index_column=table.index_column,
        annotated_taxa=annotated_taxa,
    )


def write_annotation_report(path: Path, report: TableLinkageReport) -> Path:
    """Write a linkage report to a deterministic JSON artifact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(report), default=str, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _section(name: str, payload: object) -> tuple[str, str]:
    return name, json.dumps(payload, default=str, indent=2, sort_keys=True)


def render_phylogenetics_report(
    *,
    tree_path: Path,
    out_path: Path,
    alignment_path: Path | None = None,
    traits_path: Path | None = None,
    metadata_path: Path | None = None,
) -> ReportBuildResult:
    """Build an HTML report around a tree and optional evidence tables."""
    validation = validate_tree_path(tree_path)
    inspection = inspect_tree_path(tree_path)
    alignment = summarise_fasta(alignment_path) if alignment_path else None
    traits_linkage = annotate_tree_against_table(tree_path, traits_path) if traits_path else None
    metadata_linkage = annotate_tree_against_table(tree_path, metadata_path) if metadata_path else None

    sections = [
        _section("tree-validation", asdict(validation)),
        _section("tree-inspection", asdict(inspection)),
    ]
    if alignment is not None:
        sections.append(_section("alignment-summary", asdict(alignment)))
    if traits_linkage is not None:
        sections.append(_section("traits-linkage", asdict(traits_linkage)))
    if metadata_linkage is not None:
        sections.append(_section("metadata-linkage", asdict(metadata_linkage)))

    write_html_report(title="Bijux Phylogenetics Report", sections=sections, out_path=out_path)
    return ReportBuildResult(
        output_path=out_path,
        validation=validation,
        inspection=inspection,
        metadata_linkage=metadata_linkage,
        traits_linkage=traits_linkage,
        alignment=alignment,
    )
