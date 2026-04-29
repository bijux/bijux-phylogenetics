from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from bijux_phylogenetics.core.alignment import (
    AlignmentForensicReport,
    AlignmentLinkageReport,
    AlignmentQualityReport,
    AlignmentSummary,
    CodingAlignmentDiagnostics,
    SequenceIdentityMatrix,
)
from bijux_phylogenetics.core.dataset import DatasetAuditReport, DatasetReadinessSummary, audit_dataset_inputs, summarize_dataset_readiness
from bijux_phylogenetics.core.metadata import MetadataJoinRow, join_table_to_taxa, load_taxon_table
from bijux_phylogenetics.core.traits import TraitMissingValueReport, detect_missing_trait_values
from bijux_phylogenetics.distance import (
    assess_distance_method_assumptions,
    assess_imported_distance_method_assumptions,
    build_distance_tree,
    build_tree_from_imported_distance_matrix,
    compare_distance_tree_topologies,
    compute_pairwise_genetic_distance_matrix,
    inspect_distance_matrix_quality,
    validate_distance_reference_examples,
    validate_imported_distance_matrix,
)
from bijux_phylogenetics.diagnostics.validation import (
    TreeForensicReport,
    TreeInspectionReport,
    TreeValidationReport,
    forensic_tree_path,
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.diagnostics.validation import _load_tree
from bijux_phylogenetics.io.fasta import (
    build_alignment_forensic_report,
    build_alignment_quality_report,
    compute_pairwise_sequence_identity_matrix,
    inspect_coding_alignment,
    link_alignment_to_tree,
    summarise_fasta,
)
from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.tree_set import (
    cluster_trees_by_topology,
    compare_posterior_tree_sets,
    compute_clade_frequency_table,
    compute_consensus_tree,
    compute_tree_distance_matrix,
    detect_unstable_clades,
    detect_unstable_taxa,
    load_tree_set,
)


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
    forensic: TreeForensicReport
    metadata_linkage: TableLinkageReport | None
    traits_linkage: TableLinkageReport | None
    trait_missing_values: TraitMissingValueReport | None
    alignment: AlignmentSummary | None
    alignment_quality: AlignmentQualityReport | None
    alignment_forensic: AlignmentForensicReport | None
    alignment_coding: CodingAlignmentDiagnostics | None
    alignment_identity_matrix: SequenceIdentityMatrix | None
    alignment_linkage: AlignmentLinkageReport | None
    dataset_readiness: DatasetReadinessSummary | None
    dataset_audit: DatasetAuditReport | None
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class DistanceReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    source_path: Path
    source_kind: str
    method_limitations: list[str]
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class TreeUncertaintyReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    source_path: Path
    tree_count: int
    rooted_topology_count: int
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class TreeSetComparisonReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    left_path: Path
    right_path: Path
    shared_rooted_topology_count: int
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


def distance_method_limitations() -> list[str]:
    """Explain why distance-based tree building is approximate."""
    return [
        "distance methods collapse site-by-site sequence evidence into pairwise summaries before tree building",
        "different evolutionary histories can yield similar pairwise distances, so topology is not uniquely identified by the matrix alone",
        "UPGMA additionally assumes an ultrametric clock-like process and can misplace taxa when rates vary across lineages",
        "Neighbor-Joining is often useful for quick structure, but it is still a summary approximation rather than a full likelihood inference",
    ]


def render_tree_report(*, tree_path: Path, out_path: Path) -> ReportBuildResult:
    """Build the explicit single-tree report contract."""
    validation = validate_tree_path(tree_path)
    inspection = inspect_tree_path(tree_path)
    forensic = forensic_tree_path(tree_path)
    title = "Bijux Tree Report"
    sections = [
        _section("tree-validation", asdict(validation)),
        _section("tree-inspection", asdict(inspection)),
        _section("tree-forensic", asdict(forensic)),
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
        forensic=forensic,
        metadata_linkage=None,
        traits_linkage=None,
        trait_missing_values=None,
        alignment=None,
        alignment_quality=None,
        alignment_forensic=None,
        alignment_coding=None,
        alignment_identity_matrix=None,
        alignment_linkage=None,
        dataset_readiness=None,
        dataset_audit=None,
        machine_manifest=machine_manifest,
    )


def render_dataset_report(
    *,
    tree_path: Path,
    metadata_path: Path,
    out_path: Path,
    traits_path: Path | None = None,
    alignment_path: Path | None = None,
    tip_dates_path: Path | None = None,
    calibration_path: Path | None = None,
) -> ReportBuildResult:
    """Build the explicit tree plus table dataset report contract."""
    validation = validate_tree_path(tree_path)
    inspection = inspect_tree_path(tree_path)
    forensic = forensic_tree_path(tree_path)
    metadata_linkage = annotate_tree_against_table(tree_path, metadata_path)
    traits_linkage = annotate_tree_against_table(tree_path, traits_path) if traits_path is not None else None
    trait_missing_values = detect_missing_trait_values(traits_path) if traits_path is not None else None
    dataset_readiness = (
        summarize_dataset_readiness(tree_path, metadata_path, traits_path)
        if traits_path is not None
        else None
    )
    dataset_audit = (
        audit_dataset_inputs(
            tree_path,
            metadata_path,
            traits_path,
            alignment_path=alignment_path,
            tip_dates_path=tip_dates_path,
            calibration_path=calibration_path,
        )
        if traits_path is not None
        else None
    )
    title = "Bijux Dataset Report"
    sections = [
        _section("tree-validation", asdict(validation)),
        _section("tree-inspection", asdict(inspection)),
        _section("tree-forensic", asdict(forensic)),
        _section("metadata-linkage", asdict(metadata_linkage)),
    ]
    if traits_linkage is not None:
        sections.append(_section("traits-linkage", asdict(traits_linkage)))
    if trait_missing_values is not None:
        sections.append(_section("trait-missing-values", asdict(trait_missing_values)))
    if dataset_readiness is not None:
        sections.append(_section("dataset-readiness", asdict(dataset_readiness)))
    if dataset_audit is not None:
        sections.append(_section("dataset-audit", asdict(dataset_audit)))
        sections.append(_section("dataset-findings", [asdict(row) for row in dataset_audit.findings]))
        sections.append(_section("dataset-analysis-decisions", [asdict(row) for row in dataset_audit.analysis_decisions]))
        sections.append(_section("dataset-readiness-levels", [asdict(row) for row in dataset_audit.readiness_levels]))
        sections.append(_section("dataset-crosswalk", asdict(dataset_audit.crosswalk)))
        sections.append(_section("dataset-completeness", asdict(dataset_audit.completeness_matrix)))
        sections.append(_section("dataset-exclusions", asdict(dataset_audit.exclusion_table)))
        sections.append(_section("dataset-ordering", asdict(dataset_audit.ordering_audit)))
        sections.append(_section("dataset-pruning", [asdict(row) for row in dataset_audit.pruning_steps]))
        sections.append(_section("dataset-group-imbalance", [asdict(row) for row in dataset_audit.group_imbalance_warnings]))
    input_paths = [tree_path, metadata_path]
    if traits_path is not None:
        input_paths.append(traits_path)
    if alignment_path is not None:
        input_paths.append(alignment_path)
    if tip_dates_path is not None:
        input_paths.append(tip_dates_path)
    if calibration_path is not None:
        input_paths.append(calibration_path)
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
        forensic=forensic,
        metadata_linkage=metadata_linkage,
        traits_linkage=traits_linkage,
        trait_missing_values=trait_missing_values,
        alignment=None,
        alignment_quality=None,
        alignment_forensic=None,
        alignment_coding=None,
        alignment_identity_matrix=None,
        alignment_linkage=None,
        dataset_readiness=dataset_readiness,
        dataset_audit=dataset_audit,
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
    forensic = forensic_tree_path(tree_path)
    alignment = summarise_fasta(alignment_path)
    alignment_quality = build_alignment_quality_report(alignment_path)
    alignment_forensic = build_alignment_forensic_report(alignment_path)
    alignment_coding = (
        inspect_coding_alignment(alignment_path)
        if alignment.inferred_alphabet in {"dna", "rna"}
        else None
    )
    alignment_identity_matrix = compute_pairwise_sequence_identity_matrix(alignment_path)
    alignment_linkage = link_alignment_to_tree(tree_path, alignment_path)
    title = "Bijux Phylo Inputs Report"
    sections = [
        _section("tree-validation", asdict(validation)),
        _section("tree-inspection", asdict(inspection)),
        _section("tree-forensic", asdict(forensic)),
        _section("alignment-summary", asdict(alignment)),
        _section("alignment-quality", asdict(alignment_quality)),
        _section("alignment-forensic", asdict(alignment_forensic)),
        *([_section("alignment-coding", asdict(alignment_coding))] if alignment_coding is not None else []),
        _section("alignment-identity-matrix", asdict(alignment_identity_matrix)),
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
        forensic=forensic,
        metadata_linkage=None,
        traits_linkage=None,
        trait_missing_values=None,
        alignment=alignment,
        alignment_quality=alignment_quality,
        alignment_forensic=alignment_forensic,
        alignment_coding=alignment_coding,
        alignment_identity_matrix=alignment_identity_matrix,
        alignment_linkage=alignment_linkage,
        dataset_readiness=None,
        dataset_audit=None,
        machine_manifest=machine_manifest,
    )


def render_distance_report(
    *,
    out_path: Path,
    alignment_path: Path | None = None,
    matrix_path: Path | None = None,
) -> DistanceReportBuildResult:
    """Build a deterministic HTML report for computed or imported distance analysis."""
    if (alignment_path is None) == (matrix_path is None):
        raise ValueError("render_distance_report requires exactly one of alignment_path or matrix_path")

    method_limitations = distance_method_limitations()
    if alignment_path is not None:
        matrix = compute_pairwise_genetic_distance_matrix(alignment_path)
        quality = inspect_distance_matrix_quality(alignment_path)
        assumptions = assess_distance_method_assumptions(alignment_path)
        reference_validation = validate_distance_reference_examples()
        nj_tree, _ = build_distance_tree(alignment_path, method="neighbor-joining")
        upgma_tree, _ = build_distance_tree(alignment_path, method="upgma")
        comparison = compare_distance_tree_topologies(alignment_path)
        title = "Bijux Distance Analysis Report"
        sections = [
            _section("computed-distance-matrix", asdict(matrix)),
            _section("distance-quality", asdict(quality)),
            _section("distance-method-assumptions", asdict(assumptions)),
            _section("distance-reference-validation", asdict(reference_validation)),
            _section("neighbor-joining-tree", {"newick": dumps_newick(nj_tree)}),
            _section("upgma-tree", {"newick": dumps_newick(upgma_tree)}),
            _section("distance-tree-comparison", asdict(comparison)),
            _section("distance-method-limitations", method_limitations),
        ]
        machine_manifest = {
            "report_kind": "distance-analysis",
            "source_kind": "alignment",
            "source_path": str(alignment_path),
            "method_limitations": method_limitations,
            "sections": [name for name, _ in sections],
        }
        write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
        return DistanceReportBuildResult(
            output_path=out_path,
            report_kind="distance-analysis",
            title=title,
            source_path=alignment_path,
            source_kind="alignment",
            method_limitations=method_limitations,
            machine_manifest=machine_manifest,
        )

    validation = validate_imported_distance_matrix(matrix_path)
    assumptions = assess_imported_distance_method_assumptions(matrix_path)
    title = "Bijux Imported Distance Report"
    sections = [
        _section("imported-distance-matrix-validation", asdict(validation)),
        _section("distance-method-assumptions", asdict(assumptions)),
        _section("distance-method-limitations", method_limitations),
    ]
    if validation.complete and validation.symmetric and validation.zero_diagonal and validation.nonnegative:
        nj_tree, _ = build_tree_from_imported_distance_matrix(matrix_path, method="neighbor-joining")
        upgma_tree, _ = build_tree_from_imported_distance_matrix(matrix_path, method="upgma")
        sections.extend(
            [
                _section("neighbor-joining-tree", {"newick": dumps_newick(nj_tree)}),
                _section("upgma-tree", {"newick": dumps_newick(upgma_tree)}),
            ]
        )
    machine_manifest = {
        "report_kind": "distance-analysis",
        "source_kind": "imported-distance-matrix",
        "source_path": str(matrix_path),
        "method_limitations": method_limitations,
        "sections": [name for name, _ in sections],
    }
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return DistanceReportBuildResult(
        output_path=out_path,
        report_kind="distance-analysis",
        title=title,
        source_path=matrix_path,
        source_kind="imported-distance-matrix",
        method_limitations=method_limitations,
        machine_manifest=machine_manifest,
    )


def render_tree_uncertainty_report(*, tree_set_path: Path, out_path: Path) -> TreeUncertaintyReportBuildResult:
    """Build a deterministic HTML report for consensus and uncertainty across a tree set."""
    summary = load_tree_set(tree_set_path)
    consensus_tree, consensus = compute_consensus_tree(tree_set_path)
    clade_frequencies = compute_clade_frequency_table(tree_set_path)
    distances = compute_tree_distance_matrix(tree_set_path)
    clusters = cluster_trees_by_topology(tree_set_path)
    unstable_taxa = detect_unstable_taxa(tree_set_path)
    unstable_clades = detect_unstable_clades(tree_set_path)
    title = "Bijux Tree Uncertainty Report"
    sections = [
        _section("tree-set-summary", asdict(summary)),
        _section("consensus-tree", {"newick": dumps_newick(consensus_tree), "report": asdict(consensus)}),
        _section("clade-frequencies", asdict(clade_frequencies)),
        _section("pairwise-tree-distances", asdict(distances)),
        _section("topology-clusters", asdict(clusters)),
        _section("unstable-taxa", asdict(unstable_taxa)),
        _section("unstable-clades", asdict(unstable_clades)),
    ]
    machine_manifest = {
        "report_kind": "tree-uncertainty",
        "title": title,
        "source_path": str(tree_set_path),
        "input_checksum": _sha256(tree_set_path),
        "tree_count": summary.tree_count,
        "rooted_topology_count": summary.rooted_topology_count,
        "sections": [name for name, _ in sections],
    }
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return TreeUncertaintyReportBuildResult(
        output_path=out_path,
        report_kind="tree-uncertainty",
        title=title,
        source_path=tree_set_path,
        tree_count=summary.tree_count,
        rooted_topology_count=summary.rooted_topology_count,
        machine_manifest=machine_manifest,
    )


def render_tree_set_comparison_report(
    *,
    left_tree_set_path: Path,
    right_tree_set_path: Path,
    out_path: Path,
) -> TreeSetComparisonReportBuildResult:
    """Render an HTML comparison report for two tree sets."""
    comparison = compare_posterior_tree_sets(left_tree_set_path, right_tree_set_path)
    left_summary = load_tree_set(left_tree_set_path)
    right_summary = load_tree_set(right_tree_set_path)
    left_clusters = cluster_trees_by_topology(left_tree_set_path)
    right_clusters = cluster_trees_by_topology(right_tree_set_path)
    left_unstable_taxa = detect_unstable_taxa(left_tree_set_path)
    right_unstable_taxa = detect_unstable_taxa(right_tree_set_path)
    left_unstable_clades = detect_unstable_clades(left_tree_set_path)
    right_unstable_clades = detect_unstable_clades(right_tree_set_path)
    sections = [
        _section("tree-set-comparison", asdict(comparison)),
        _section("left-tree-set-summary", asdict(left_summary)),
        _section("right-tree-set-summary", asdict(right_summary)),
        _section("left-topology-clusters", asdict(left_clusters)),
        _section("right-topology-clusters", asdict(right_clusters)),
        _section("left-unstable-taxa", asdict(left_unstable_taxa)),
        _section("right-unstable-taxa", asdict(right_unstable_taxa)),
        _section("left-unstable-clades", asdict(left_unstable_clades)),
        _section("right-unstable-clades", asdict(right_unstable_clades)),
    ]
    title = "Bijux Tree-Set Comparison Report"
    machine_manifest = {
        "report_kind": "tree-set-comparison",
        "title": title,
        "left_path": str(left_tree_set_path),
        "right_path": str(right_tree_set_path),
        "shared_rooted_topology_count": comparison.shared_rooted_topology_count,
        "sections": [name for name, _ in sections],
    }
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return TreeSetComparisonReportBuildResult(
        output_path=out_path,
        report_kind="tree-set-comparison",
        title=title,
        left_path=left_tree_set_path,
        right_path=right_tree_set_path,
        shared_rooted_topology_count=comparison.shared_rooted_topology_count,
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
    forensic = forensic_tree_path(tree_path)
    alignment = summarise_fasta(alignment_path) if alignment_path else None
    alignment_quality = build_alignment_quality_report(alignment_path) if alignment_path else None
    alignment_forensic = build_alignment_forensic_report(alignment_path) if alignment_path else None
    alignment_coding = (
        inspect_coding_alignment(alignment_path)
        if alignment_path is not None and alignment is not None and alignment.inferred_alphabet in {"dna", "rna"}
        else None
    )
    alignment_identity_matrix = (
        compute_pairwise_sequence_identity_matrix(alignment_path)
        if alignment_path is not None
        else None
    )
    traits_linkage = annotate_tree_against_table(tree_path, traits_path) if traits_path else None
    trait_missing_values = detect_missing_trait_values(traits_path) if traits_path else None
    metadata_linkage = annotate_tree_against_table(tree_path, metadata_path) if metadata_path else None
    dataset_readiness = (
        summarize_dataset_readiness(tree_path, metadata_path, traits_path)
        if traits_path and metadata_path
        else None
    )
    dataset_audit = (
        audit_dataset_inputs(
            tree_path,
            metadata_path,
            traits_path,
            alignment_path=alignment_path,
        )
        if traits_path and metadata_path
        else None
    )

    sections = [
        _section("tree-validation", asdict(validation)),
        _section("tree-inspection", asdict(inspection)),
        _section("tree-forensic", asdict(forensic)),
    ]
    if alignment is not None:
        sections.append(_section("alignment-summary", asdict(alignment)))
    if alignment_quality is not None:
        sections.append(_section("alignment-quality", asdict(alignment_quality)))
    if alignment_forensic is not None:
        sections.append(_section("alignment-forensic", asdict(alignment_forensic)))
    if alignment_coding is not None:
        sections.append(_section("alignment-coding", asdict(alignment_coding)))
    if alignment_identity_matrix is not None:
        sections.append(_section("alignment-identity-matrix", asdict(alignment_identity_matrix)))
    if traits_linkage is not None:
        sections.append(_section("traits-linkage", asdict(traits_linkage)))
    if trait_missing_values is not None:
        sections.append(_section("trait-missing-values", asdict(trait_missing_values)))
    if metadata_linkage is not None:
        sections.append(_section("metadata-linkage", asdict(metadata_linkage)))
    if dataset_readiness is not None:
        sections.append(_section("dataset-readiness", asdict(dataset_readiness)))
    if dataset_audit is not None:
        sections.append(_section("dataset-audit", asdict(dataset_audit)))
        sections.append(_section("dataset-crosswalk", asdict(dataset_audit.crosswalk)))
        sections.append(_section("dataset-completeness", asdict(dataset_audit.completeness_matrix)))

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
        forensic=forensic,
        metadata_linkage=metadata_linkage,
        traits_linkage=traits_linkage,
        trait_missing_values=trait_missing_values,
        alignment=alignment,
        alignment_quality=alignment_quality,
        alignment_forensic=alignment_forensic,
        alignment_coding=alignment_coding,
        alignment_identity_matrix=alignment_identity_matrix,
        alignment_linkage=None,
        dataset_readiness=dataset_readiness,
        dataset_audit=dataset_audit,
        machine_manifest=machine_manifest,
    )
