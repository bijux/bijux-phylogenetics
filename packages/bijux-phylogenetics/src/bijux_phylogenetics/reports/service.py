from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from bijux_phylogenetics.core.alignment import AlignmentLinkageReport, AlignmentSummary
from bijux_phylogenetics.core.dataset import DatasetReadinessSummary, summarize_dataset_readiness
from bijux_phylogenetics.core.metadata import MetadataJoinRow, join_table_to_taxa, load_taxon_table
from bijux_phylogenetics.core.traits import TraitMissingValueReport, detect_missing_trait_values
from bijux_phylogenetics.diagnostics.validation import TreeInspectionReport, TreeValidationReport, inspect_tree_path, validate_tree_path
from bijux_phylogenetics.diagnostics.validation import _load_tree
from bijux_phylogenetics.io.fasta import link_alignment_to_tree, summarise_fasta
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
    joined_rows: list[MetadataJoinRow]


@dataclass(slots=True)
class ReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    validation: TreeValidationReport
    inspection: TreeInspectionReport
    metadata_linkage: TableLinkageReport | None
    traits_linkage: TableLinkageReport | None
    trait_missing_values: TraitMissingValueReport | None
    alignment: AlignmentSummary | None
    alignment_linkage: AlignmentLinkageReport | None
    dataset_readiness: DatasetReadinessSummary | None
    machine_manifest: dict[str, object]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
    path.write_text(json.dumps(asdict(report), default=str, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _section(name: str, payload: object) -> tuple[str, str]:
    return name, json.dumps(payload, default=str, indent=2, sort_keys=True)


def _build_machine_manifest(
    *,
    report_kind: str,
    title: str,
    input_paths: list[Path],
    sections: list[tuple[str, str]],
    inspection: TreeInspectionReport,
) -> dict[str, object]:
    return {
        "report_kind": report_kind,
        "title": title,
        "input_paths": [str(path) for path in input_paths],
        "input_checksums": {str(path): _sha256(path) for path in input_paths},
        "sections": [name for name, _ in sections],
        "metrics": {
            "tip_count": inspection.tip_count,
            "node_count": inspection.node_count,
            "clade_count": inspection.clade_count,
        },
    }


def render_tree_report(*, tree_path: Path, out_path: Path) -> ReportBuildResult:
    """Build the explicit single-tree report contract."""
    validation = validate_tree_path(tree_path)
    inspection = inspect_tree_path(tree_path)
    title = "Bijux Tree Report"
    sections = [
        _section("tree-validation", asdict(validation)),
        _section("tree-inspection", asdict(inspection)),
    ]
    machine_manifest = _build_machine_manifest(
        report_kind="tree",
        title=title,
        input_paths=[tree_path],
        sections=sections,
        inspection=inspection,
    )
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return ReportBuildResult(
        output_path=out_path,
        report_kind="tree",
        title=title,
        validation=validation,
        inspection=inspection,
        metadata_linkage=None,
        traits_linkage=None,
        trait_missing_values=None,
        alignment=None,
        alignment_linkage=None,
        dataset_readiness=None,
        machine_manifest=machine_manifest,
    )


def render_dataset_report(
    *,
    tree_path: Path,
    metadata_path: Path,
    out_path: Path,
    traits_path: Path | None = None,
) -> ReportBuildResult:
    """Build the explicit tree plus table dataset report contract."""
    validation = validate_tree_path(tree_path)
    inspection = inspect_tree_path(tree_path)
    metadata_linkage = annotate_tree_against_table(tree_path, metadata_path)
    traits_linkage = annotate_tree_against_table(tree_path, traits_path) if traits_path is not None else None
    trait_missing_values = detect_missing_trait_values(traits_path) if traits_path is not None else None
    dataset_readiness = (
        summarize_dataset_readiness(tree_path, metadata_path, traits_path)
        if traits_path is not None
        else None
    )
    title = "Bijux Dataset Report"
    sections = [
        _section("tree-validation", asdict(validation)),
        _section("tree-inspection", asdict(inspection)),
        _section("metadata-linkage", asdict(metadata_linkage)),
    ]
    if traits_linkage is not None:
        sections.append(_section("traits-linkage", asdict(traits_linkage)))
    if trait_missing_values is not None:
        sections.append(_section("trait-missing-values", asdict(trait_missing_values)))
    if dataset_readiness is not None:
        sections.append(_section("dataset-readiness", asdict(dataset_readiness)))
    input_paths = [tree_path, metadata_path]
    if traits_path is not None:
        input_paths.append(traits_path)
    machine_manifest = _build_machine_manifest(
        report_kind="dataset",
        title=title,
        input_paths=input_paths,
        sections=sections,
        inspection=inspection,
    )
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return ReportBuildResult(
        output_path=out_path,
        report_kind="dataset",
        title=title,
        validation=validation,
        inspection=inspection,
        metadata_linkage=metadata_linkage,
        traits_linkage=traits_linkage,
        trait_missing_values=trait_missing_values,
        alignment=None,
        alignment_linkage=None,
        dataset_readiness=dataset_readiness,
        machine_manifest=machine_manifest,
    )


def render_phylo_inputs_report(
    *,
    tree_path: Path,
    alignment_path: Path,
    out_path: Path,
) -> ReportBuildResult:
    """Build the explicit tree plus alignment input report contract."""
    validation = validate_tree_path(tree_path)
    inspection = inspect_tree_path(tree_path)
    alignment = summarise_fasta(alignment_path)
    alignment_linkage = link_alignment_to_tree(tree_path, alignment_path)
    title = "Bijux Phylo Inputs Report"
    sections = [
        _section("tree-validation", asdict(validation)),
        _section("tree-inspection", asdict(inspection)),
        _section("alignment-summary", asdict(alignment)),
        _section("alignment-linkage", asdict(alignment_linkage)),
    ]
    machine_manifest = _build_machine_manifest(
        report_kind="phylo-inputs",
        title=title,
        input_paths=[tree_path, alignment_path],
        sections=sections,
        inspection=inspection,
    )
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return ReportBuildResult(
        output_path=out_path,
        report_kind="phylo-inputs",
        title=title,
        validation=validation,
        inspection=inspection,
        metadata_linkage=None,
        traits_linkage=None,
        trait_missing_values=None,
        alignment=alignment,
        alignment_linkage=alignment_linkage,
        dataset_readiness=None,
        machine_manifest=machine_manifest,
    )


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
    trait_missing_values = detect_missing_trait_values(traits_path) if traits_path else None
    metadata_linkage = annotate_tree_against_table(tree_path, metadata_path) if metadata_path else None
    dataset_readiness = (
        summarize_dataset_readiness(tree_path, metadata_path, traits_path)
        if traits_path and metadata_path
        else None
    )

    sections = [
        _section("tree-validation", asdict(validation)),
        _section("tree-inspection", asdict(inspection)),
    ]
    if alignment is not None:
        sections.append(_section("alignment-summary", asdict(alignment)))
    if traits_linkage is not None:
        sections.append(_section("traits-linkage", asdict(traits_linkage)))
    if trait_missing_values is not None:
        sections.append(_section("trait-missing-values", asdict(trait_missing_values)))
    if metadata_linkage is not None:
        sections.append(_section("metadata-linkage", asdict(metadata_linkage)))
    if dataset_readiness is not None:
        sections.append(_section("dataset-readiness", asdict(dataset_readiness)))

    title = "Bijux Phylogenetics Report"
    input_paths = [tree_path]
    if alignment_path is not None:
        input_paths.append(alignment_path)
    if traits_path is not None:
        input_paths.append(traits_path)
    if metadata_path is not None:
        input_paths.append(metadata_path)
    machine_manifest = _build_machine_manifest(
        report_kind="phylogenetics",
        title=title,
        input_paths=input_paths,
        sections=sections,
        inspection=inspection,
    )
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return ReportBuildResult(
        output_path=out_path,
        report_kind="phylogenetics",
        title=title,
        validation=validation,
        inspection=inspection,
        metadata_linkage=metadata_linkage,
        traits_linkage=traits_linkage,
        trait_missing_values=trait_missing_values,
        alignment=alignment,
        alignment_linkage=None,
        dataset_readiness=dataset_readiness,
        machine_manifest=machine_manifest,
    )
