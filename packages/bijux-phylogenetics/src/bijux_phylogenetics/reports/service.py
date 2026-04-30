from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
import csv

from bijux_phylogenetics.core.alignment import (
    AlignmentForensicReport,
    AlignmentLinkageReport,
    AlignmentLowInformationReport,
    AlignmentQualityReport,
    AlignmentSummary,
    AlignmentAmbiguousColumnReport,
    CodingAlignmentDiagnostics,
    DuplicateSequencePolicyReport,
    SequenceIdentityMatrix,
    SequenceQualityRankingReport,
)
from bijux_phylogenetics.core.dataset import (
    DatasetAuditReport,
    DatasetCrosswalkReport,
    DatasetReadinessSummary,
    audit_dataset_inputs,
    build_dataset_crosswalk,
    summarize_dataset_readiness,
)
from bijux_phylogenetics.core.metadata import MetadataJoinRow, join_table_to_taxa, load_taxon_table
from bijux_phylogenetics.core.taxonomy import TaxonAuditReport, build_taxon_audit_report
from bijux_phylogenetics.core.taxon_workflows import TaxonWorkflowLossReport, build_taxon_workflow_loss_report
from bijux_phylogenetics.core.traits import TraitMissingValueReport, detect_missing_trait_values
from bijux_phylogenetics.reference_validation import (
    CoreWorkflowValidationReport,
    LevelOneReleaseGateReport,
    build_core_workflow_validation_report,
    build_level_one_release_gate_report,
)
from bijux_phylogenetics.distance import (
    assess_distance_method_assumptions,
    assess_distance_method_maturity,
    assess_imported_distance_method_assumptions,
    build_distance_method_report,
    build_distance_tree,
    build_tree_from_imported_distance_matrix,
    compare_distance_tree_topologies,
    compute_pairwise_genetic_distance_matrix,
    inspect_imported_distance_matrix_quality,
    inspect_distance_matrix_quality,
    validate_distance_reference_examples,
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
    assess_alignment_low_information,
    build_ambiguous_alignment_column_report,
    build_alignment_forensic_report,
    build_alignment_quality_report,
    build_duplicate_sequence_policy_report,
    build_sequence_quality_ranking,
    compute_pairwise_sequence_identity_matrix,
    inspect_coding_alignment,
    list_alignment_filter_profiles,
    link_alignment_to_tree,
    summarise_fasta,
)
from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.tree_set import (
    assess_tree_set_maturity,
    assess_tree_set_storage_risk,
    assess_tree_set_thinning_sensitivity,
    benchmark_tree_set_uncertainty,
    cluster_trees_by_topology,
    compare_consensus_thresholds,
    compare_posterior_topological_diversity,
    compare_posterior_tree_sets,
    compute_clade_frequency_table,
    compute_consensus_tree,
    compute_tree_distance_matrix,
    detect_posterior_topology_multimodality,
    detect_unstable_clades,
    detect_unstable_taxa,
    load_tree_set,
    summarize_clade_credibility_conflicts,
    summarize_uncertainty_aware_conclusions,
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
    machine_manifest_path: Path
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
    alignment_low_information: AlignmentLowInformationReport | None
    alignment_duplicate_policy: DuplicateSequencePolicyReport | None
    alignment_ambiguous_columns: AlignmentAmbiguousColumnReport | None
    alignment_sequence_ranking: SequenceQualityRankingReport | None
    alignment_coding: CodingAlignmentDiagnostics | None
    alignment_identity_matrix: SequenceIdentityMatrix | None
    alignment_linkage: AlignmentLinkageReport | None
    dataset_readiness: DatasetReadinessSummary | None
    dataset_audit: DatasetAuditReport | None
    input_ledger: list["ReportInputLedgerEntry"]
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


@dataclass(slots=True)
class AlignmentReportBuildResult:
    output_path: Path
    machine_manifest_path: Path
    report_kind: str
    title: str
    alignment: AlignmentSummary
    alignment_quality: AlignmentQualityReport
    alignment_forensic: AlignmentForensicReport
    alignment_low_information: AlignmentLowInformationReport
    alignment_duplicate_policy: DuplicateSequencePolicyReport
    alignment_ambiguous_columns: AlignmentAmbiguousColumnReport
    alignment_sequence_ranking: SequenceQualityRankingReport
    alignment_coding: CodingAlignmentDiagnostics | None
    alignment_identity_matrix: SequenceIdentityMatrix
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class TaxonReportBuildResult:
    output_path: Path
    machine_manifest_path: Path
    report_kind: str
    title: str
    tree_path: Path
    taxon_audit: TaxonAuditReport
    taxon_crosswalk: DatasetCrosswalkReport | None
    taxon_workflow_loss: TaxonWorkflowLossReport | None
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class WorkflowValidationReportBuildResult:
    output_path: Path
    machine_manifest_path: Path
    report_kind: str
    title: str
    validation: CoreWorkflowValidationReport
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class ReleaseGateReportBuildResult:
    output_path: Path
    machine_manifest_path: Path
    report_kind: str
    title: str
    release_gate: LevelOneReleaseGateReport
    machine_manifest: dict[str, object]


@dataclass(frozen=True, slots=True)
class ReportInputLedgerEntry:
    path: Path
    role: str
    checksum: str
    taxa_count: int
    usage: list[str]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _dataset_surface_taxa_count(path: Path, role: str) -> int:
    if role == "tree":
        return _load_tree(path).tip_count
    if role == "alignment":
        return summarise_fasta(path).sequence_count
    if role in {"metadata", "traits", "tip_dates"}:
        return load_taxon_table(path).row_count
    if role == "calibrations":
        delimiter = "," if path.suffix.lower() == ".csv" else "\t"
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            return sum(1 for _ in reader)
    raise ValueError(f"unsupported dataset ledger role: {role}")


def _build_input_ledger(entries: list[tuple[Path, str, list[str]]]) -> list[ReportInputLedgerEntry]:
    return [
        ReportInputLedgerEntry(
            path=path,
            role=role,
            checksum=_sha256(path),
            taxa_count=_dataset_surface_taxa_count(path, role),
            usage=usage,
        )
        for path, role, usage in entries
    ]


def _serialize_input_ledger(entries: list[ReportInputLedgerEntry]) -> list[dict[str, object]]:
    return [
        {
            "path": str(entry.path),
            "role": entry.role,
            "checksum": entry.checksum,
            "taxa_count": entry.taxa_count,
            "usage": entry.usage,
        }
        for entry in entries
    ]


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


def _write_machine_manifest(path: Path, payload: dict[str, object]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, default=str, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


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


def _report_sidecar_path(out_path: Path) -> Path:
    return out_path.with_suffix(".json")


def _report_summary_and_limitations(
    *,
    report_kind: str,
    validation: TreeValidationReport,
    inspection: TreeInspectionReport,
    forensic: TreeForensicReport,
    dataset_audit: DatasetAuditReport | None = None,
    alignment_forensic: AlignmentForensicReport | None = None,
) -> tuple[list[str], list[str]]:
    summary = [
        f"tree validity decision: {validation.validity_decision}",
        f"tree quality score: {inspection.tree_quality_score}",
        (
            "tree is currently safe for publication-oriented use"
            if forensic.safe_for_publication
            else "tree still carries publication-facing risks that should be reviewed"
        ),
    ]
    limitations = list(forensic.warnings)
    if report_kind == "dataset" and dataset_audit is not None:
        summary.append(f"dataset readiness decision: {dataset_audit.readiness_decision}")
        summary.append(
            f"blocked analyses: {len(dataset_audit.blocked_analyses)}, risky analyses: "
            f"{sum(1 for row in dataset_audit.analysis_decisions if row.decision == 'risky')}"
        )
        limitations.extend(finding.message for finding in dataset_audit.findings if finding.severity in {"warning", "blocker"})
    if report_kind == "phylo-inputs" and alignment_forensic is not None:
        safe_methods = sum(
            1
            for flag in (
                alignment_forensic.safe_for_distance_analysis,
                alignment_forensic.safe_for_maximum_likelihood,
                alignment_forensic.safe_for_bayesian_inference,
                alignment_forensic.safe_for_coding_analysis,
                alignment_forensic.safe_for_publication,
            )
            if flag
        )
        summary.append(f"alignment safe-for flags passed: {safe_methods}/5")
        limitations.extend(alignment_forensic.limitations)
    if inspection.mixed_support_scales:
        limitations.append("support labels originate from mixed scales and should be interpreted only after normalization audit")
    return summary, sorted(set(limitations))


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
    reviewer_summary, limitations = _report_summary_and_limitations(
        report_kind="tree",
        validation=validation,
        inspection=inspection,
        forensic=forensic,
    )
    sections = [
        _section("reviewer-summary", reviewer_summary),
        _section("tree-validation", asdict(validation)),
        _section("tree-inspection", asdict(inspection)),
        _section("tree-forensic", asdict(forensic)),
        _section("limitations", limitations),
    ]
    machine_manifest = _build_machine_manifest(
        report_kind="tree",
        title=title,
        input_paths=[tree_path],
        sections=sections,
        inspection=inspection,
    )
    machine_manifest["reviewer_summary"] = reviewer_summary
    machine_manifest["limitations"] = limitations
    input_ledger = _build_input_ledger([(tree_path, "tree", ["tree_validation", "tree_inspection", "tree_forensic"])])
    machine_manifest["input_ledger"] = _serialize_input_ledger(input_ledger)
    machine_manifest_path = _write_machine_manifest(_report_sidecar_path(out_path), machine_manifest)
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return ReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
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
        alignment_low_information=None,
        alignment_duplicate_policy=None,
        alignment_ambiguous_columns=None,
        alignment_sequence_ranking=None,
        alignment_coding=None,
        alignment_identity_matrix=None,
        alignment_linkage=None,
        dataset_readiness=None,
        dataset_audit=None,
        input_ledger=input_ledger,
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
    reviewer_summary, limitations = _report_summary_and_limitations(
        report_kind="dataset",
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        dataset_audit=dataset_audit,
    )
    sections = [
        _section("reviewer-summary", reviewer_summary),
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
        sections.append(_section("dataset-mismatches", asdict(dataset_audit.mismatch_report)))
        sections.append(_section("dataset-risk-score", asdict(dataset_audit.risk_score)))
        sections.append(_section("dataset-minimal-fix-plan", asdict(dataset_audit.minimal_fix_plan)))
        sections.append(_section("dataset-reviewer-checklist", asdict(dataset_audit.reviewer_checklist)))
        sections.append(_section("dataset-ordering", asdict(dataset_audit.ordering_audit)))
        sections.append(_section("dataset-pruning", [asdict(row) for row in dataset_audit.pruning_steps]))
        sections.append(_section("dataset-group-imbalance", [asdict(row) for row in dataset_audit.group_imbalance_warnings]))
    input_ledger_entries: list[tuple[Path, str, list[str]]] = [
        (tree_path, "tree", ["tree_validation", "tree_inspection", "tree_forensic", "dataset_audit"]),
        (metadata_path, "metadata", ["metadata_linkage", "dataset_audit"]),
    ]
    if traits_path is not None:
        input_ledger_entries.append((traits_path, "traits", ["traits_linkage", "trait_missing_values", "dataset_audit"]))
    if alignment_path is not None:
        input_ledger_entries.append((alignment_path, "alignment", ["alignment_forensic", "dataset_audit"]))
    if tip_dates_path is not None:
        input_ledger_entries.append((tip_dates_path, "tip_dates", ["dataset_audit"]))
    if calibration_path is not None:
        input_ledger_entries.append((calibration_path, "calibrations", ["dataset_audit"]))
    input_ledger = _build_input_ledger(input_ledger_entries)
    sections.append(_section("dataset-input-ledger", _serialize_input_ledger(input_ledger)))
    sections.append(_section("limitations", limitations))
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
    machine_manifest["reviewer_summary"] = reviewer_summary
    machine_manifest["limitations"] = limitations
    machine_manifest["input_ledger"] = _serialize_input_ledger(input_ledger)
    machine_manifest_path = _write_machine_manifest(_report_sidecar_path(out_path), machine_manifest)
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return ReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
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
        alignment_low_information=None,
        alignment_duplicate_policy=None,
        alignment_ambiguous_columns=None,
        alignment_sequence_ranking=None,
        alignment_coding=None,
        alignment_identity_matrix=None,
        alignment_linkage=None,
        dataset_readiness=dataset_readiness,
        dataset_audit=dataset_audit,
        input_ledger=input_ledger,
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
    alignment_low_information = assess_alignment_low_information(alignment_path)
    alignment_duplicate_policy = build_duplicate_sequence_policy_report(alignment_path)
    alignment_ambiguous_columns = build_ambiguous_alignment_column_report(alignment_path)
    alignment_sequence_ranking = build_sequence_quality_ranking(alignment_path)
    alignment_coding = (
        inspect_coding_alignment(alignment_path)
        if alignment.inferred_alphabet in {"dna", "rna"}
        else None
    )
    alignment_identity_matrix = compute_pairwise_sequence_identity_matrix(alignment_path)
    alignment_linkage = link_alignment_to_tree(tree_path, alignment_path)
    title = "Bijux Phylo Inputs Report"
    reviewer_summary, limitations = _report_summary_and_limitations(
        report_kind="phylo-inputs",
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        alignment_forensic=alignment_forensic,
    )
    sections = [
        _section("reviewer-summary", reviewer_summary),
        _section("tree-validation", asdict(validation)),
        _section("tree-inspection", asdict(inspection)),
        _section("tree-forensic", asdict(forensic)),
        _section("alignment-summary", asdict(alignment)),
        _section("alignment-quality", asdict(alignment_quality)),
        _section("alignment-low-information", asdict(alignment_low_information)),
        _section("alignment-duplicate-policy", asdict(alignment_duplicate_policy)),
        _section("alignment-ambiguous-columns", asdict(alignment_ambiguous_columns)),
        _section("alignment-sequence-ranking", asdict(alignment_sequence_ranking)),
        _section("alignment-forensic", asdict(alignment_forensic)),
        *([_section("alignment-coding", asdict(alignment_coding))] if alignment_coding is not None else []),
        _section("alignment-identity-matrix", asdict(alignment_identity_matrix)),
        _section("alignment-linkage", asdict(alignment_linkage)),
        _section("limitations", limitations),
    ]
    input_ledger = _build_input_ledger([
        (tree_path, "tree", ["tree_validation", "tree_inspection", "tree_forensic", "alignment_linkage"]),
        (alignment_path, "alignment", [
            "alignment_summary",
            "alignment_quality",
            "alignment_low_information",
            "alignment_duplicate_policy",
            "alignment_ambiguous_columns",
            "alignment_sequence_ranking",
            "alignment_forensic",
            "alignment_identity_matrix",
            "alignment_linkage",
        ]),
    ])
    machine_manifest = _build_machine_manifest(
        report_kind="phylo-inputs",
        title=title,
        input_paths=[tree_path, alignment_path],
        sections=sections,
        inspection=inspection,
    )
    machine_manifest["reviewer_summary"] = reviewer_summary
    machine_manifest["limitations"] = limitations
    machine_manifest["input_ledger"] = _serialize_input_ledger(input_ledger)
    machine_manifest_path = _write_machine_manifest(_report_sidecar_path(out_path), machine_manifest)
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return ReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
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
        alignment_low_information=alignment_low_information,
        alignment_duplicate_policy=alignment_duplicate_policy,
        alignment_ambiguous_columns=alignment_ambiguous_columns,
        alignment_sequence_ranking=alignment_sequence_ranking,
        alignment_coding=alignment_coding,
        alignment_identity_matrix=alignment_identity_matrix,
        alignment_linkage=alignment_linkage,
        dataset_readiness=None,
        dataset_audit=None,
        input_ledger=input_ledger,
        machine_manifest=machine_manifest,
    )


def render_alignment_report(*, alignment_path: Path, out_path: Path) -> AlignmentReportBuildResult:
    """Build a reviewer-facing alignment-only report."""
    alignment = summarise_fasta(alignment_path)
    alignment_quality = build_alignment_quality_report(alignment_path)
    alignment_low_information = assess_alignment_low_information(alignment_path)
    alignment_duplicate_policy = build_duplicate_sequence_policy_report(alignment_path)
    alignment_ambiguous_columns = build_ambiguous_alignment_column_report(alignment_path)
    alignment_sequence_ranking = build_sequence_quality_ranking(alignment_path)
    alignment_forensic = build_alignment_forensic_report(alignment_path)
    alignment_coding = inspect_coding_alignment(alignment_path) if alignment.inferred_alphabet in {"dna", "rna"} else None
    alignment_identity_matrix = compute_pairwise_sequence_identity_matrix(alignment_path)
    title = "Bijux Alignment Report"
    reviewer_summary = [
        f"alignment quality score: {alignment_quality.quality_score}",
        (
            "alignment remains suitable for at least one inference family"
            if any(
                (
                    alignment_forensic.safe_for_distance_analysis,
                    alignment_forensic.safe_for_maximum_likelihood,
                    alignment_forensic.safe_for_bayesian_inference,
                    alignment_forensic.safe_for_coding_analysis,
                )
            )
            else "alignment is currently blocked for the main inference families reviewed here"
        ),
        f"reviewer-facing warnings: {len(alignment_forensic.warnings)}",
    ]
    limitations = sorted(
        dict.fromkeys(
            [
                *alignment_forensic.limitations,
                *alignment_forensic.warnings,
            ]
        )
    )
    sections = [
        _section("reviewer-summary", reviewer_summary),
        _section("alignment-summary", asdict(alignment)),
        _section("alignment-quality", asdict(alignment_quality)),
        _section("alignment-readiness", asdict(alignment_forensic.readiness)),
        _section("alignment-low-information", asdict(alignment_low_information)),
        _section("alignment-duplicate-policy", asdict(alignment_duplicate_policy)),
        _section("alignment-ambiguous-columns", asdict(alignment_ambiguous_columns)),
        _section("alignment-sequence-ranking", asdict(alignment_sequence_ranking)),
        _section("alignment-filter-profiles", [asdict(profile) for profile in list_alignment_filter_profiles()]),
        _section("alignment-suspicious-windows", {
            "over_aligned_regions": [asdict(row) for row in alignment_forensic.over_aligned_regions],
            "under_aligned_regions": [asdict(row) for row in alignment_forensic.under_aligned_regions],
        }),
        _section("alignment-forensic", asdict(alignment_forensic)),
        *([_section("alignment-coding", asdict(alignment_coding))] if alignment_coding is not None else []),
        _section("alignment-identity-matrix", asdict(alignment_identity_matrix)),
        _section("limitations", limitations),
    ]
    machine_manifest = {
        "report_kind": "alignment",
        "title": title,
        "input_paths": [str(alignment_path)],
        "input_checksums": {str(alignment_path): _sha256(alignment_path)},
        "sections": [name for name, _ in sections],
        "metrics": {
            "sequence_count": alignment.sequence_count,
            "alignment_length": alignment.alignment_length,
            "quality_score": alignment_quality.quality_score,
        },
        "reviewer_summary": reviewer_summary,
        "limitations": limitations,
    }
    machine_manifest_path = _write_machine_manifest(_report_sidecar_path(out_path), machine_manifest)
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return AlignmentReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="alignment",
        title=title,
        alignment=alignment,
        alignment_quality=alignment_quality,
        alignment_forensic=alignment_forensic,
        alignment_low_information=alignment_low_information,
        alignment_duplicate_policy=alignment_duplicate_policy,
        alignment_ambiguous_columns=alignment_ambiguous_columns,
        alignment_sequence_ranking=alignment_sequence_ranking,
        alignment_coding=alignment_coding,
        alignment_identity_matrix=alignment_identity_matrix,
        machine_manifest=machine_manifest,
    )


def render_taxon_report(
    *,
    tree_path: Path,
    out_path: Path,
    synonym_table_path: Path | None = None,
    metadata_path: Path | None = None,
    traits_path: Path | None = None,
    alignment_path: Path | None = None,
    filtered_alignment_path: Path | None = None,
    inference_tree_path: Path | None = None,
    reported_taxa_path: Path | None = None,
) -> TaxonReportBuildResult:
    """Build a reviewer-facing taxon audit report."""
    tree = _load_tree(tree_path)
    audit = build_taxon_audit_report(tree, synonym_table_path=synonym_table_path)
    taxon_crosswalk = (
        None
        if metadata_path is None or traits_path is None
        else build_dataset_crosswalk(
            tree_path,
            metadata_path,
            traits_path,
            alignment_path=alignment_path,
        )
    )
    taxon_workflow_loss = (
        None
        if metadata_path is None or traits_path is None
        else build_taxon_workflow_loss_report(
            tree_path,
            metadata_path,
            traits_path,
            alignment_path=alignment_path,
            filtered_alignment_path=filtered_alignment_path,
            inference_tree_path=inference_tree_path,
            reported_taxa_path=reported_taxa_path,
        )
    )
    title = "Bijux Taxon Audit Report"
    reviewer_summary = [
        f"taxon audit status: {audit.status}",
        f"tree tip count: {audit.tree_tip_count}",
        *audit.summary,
    ]
    limitations = sorted(dict.fromkeys(audit.warnings))
    sections = [
        _section("reviewer-summary", reviewer_summary),
        _section("taxon-audit", asdict(audit)),
        _section("taxon-identity", asdict(audit.identity_audit)),
        _section("taxon-safety", asdict(audit.safety_report)),
        _section("taxon-namespaces", asdict(audit.namespace_report)),
        _section("taxon-rank-consistency", asdict(audit.rank_consistency)),
        *([_section("taxon-synonyms", asdict(audit.synonym_audit))] if audit.synonym_audit is not None else []),
        _section("taxon-duplicate-identities", asdict(audit.duplicate_identities)),
        _section("taxon-mapping-conflicts", asdict(audit.mapping_conflicts)),
        *([_section("taxon-accepted-names", asdict(audit.accepted_name_export))] if audit.accepted_name_export is not None else []),
        *([_section("taxon-crosswalk", asdict(taxon_crosswalk))] if taxon_crosswalk is not None else []),
        *([_section("taxon-loss", asdict(taxon_workflow_loss))] if taxon_workflow_loss is not None else []),
        *(
            [
                _section(
                    "taxon-exclusions",
                    [
                        {
                            "taxon": row.taxon,
                            "first_loss_stage": row.first_loss_stage,
                            "loss_events": [asdict(event) for event in row.loss_events],
                        }
                        for row in taxon_workflow_loss.rows
                        if row.loss_events
                    ],
                )
            ]
            if taxon_workflow_loss is not None
            else []
        ),
        _section("limitations", limitations),
    ]
    input_paths = [
        tree_path,
        *([synonym_table_path] if synonym_table_path is not None else []),
        *([metadata_path] if metadata_path is not None else []),
        *([traits_path] if traits_path is not None else []),
        *([alignment_path] if alignment_path is not None else []),
        *([filtered_alignment_path] if filtered_alignment_path is not None else []),
        *([inference_tree_path] if inference_tree_path is not None else []),
        *([reported_taxa_path] if reported_taxa_path is not None else []),
    ]
    machine_manifest = {
        "report_kind": "taxonomy",
        "title": title,
        "input_paths": [str(path) for path in input_paths],
        "input_checksums": {str(path): _sha256(path) for path in input_paths},
        "sections": [name for name, _ in sections],
        "metrics": {
            "tree_tip_count": audit.tree_tip_count,
            "warning_count": len(audit.warnings),
            "conflict_count": len(audit.mapping_conflicts.rows),
        },
        "reviewer_summary": reviewer_summary,
        "limitations": limitations,
    }
    machine_manifest_path = _write_machine_manifest(_report_sidecar_path(out_path), machine_manifest)
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return TaxonReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="taxonomy",
        title=title,
        tree_path=tree_path,
        taxon_audit=audit,
        taxon_crosswalk=taxon_crosswalk,
        taxon_workflow_loss=taxon_workflow_loss,
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
        report = build_distance_method_report(alignment_path)
        title = "Bijux Distance Analysis Report"
        sections = [
            _section("computed-distance-matrix", asdict(report.matrix)),
            _section("distance-quality", asdict(report.quality)),
            _section("distance-method-assumptions", asdict(report.assumptions)),
            _section("distance-reference-validation", asdict(report.reference_validation)),
            _section("neighbor-joining-tree", {"newick": report.built_tree_newick}),
            _section("upgma-tree", {"newick": report.alternative_tree_newick}),
            _section("distance-tree-comparison", asdict(report.topology_comparison)),
            _section("distance-bootstrap-summary", asdict(report.bootstrap_summary)),
            _section("distance-model-comparison", asdict(report.model_comparison)),
            _section("distance-gap-policy-sensitivity", asdict(report.gap_policy_sensitivity)),
            _section("distance-maturity-gate", asdict(report.maturity_gate)),
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

    quality = inspect_imported_distance_matrix_quality(matrix_path)
    assumptions = assess_imported_distance_method_assumptions(matrix_path)
    title = "Bijux Imported Distance Report"
    sections = [
        _section("imported-distance-matrix-quality", asdict(quality)),
        _section("distance-method-assumptions", asdict(assumptions)),
        _section("distance-method-limitations", method_limitations),
    ]
    validation = quality.validation
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
    diversity = compare_posterior_topological_diversity(tree_set_path, tree_set_path)
    multimodality = detect_posterior_topology_multimodality(tree_set_path)
    unstable_taxa = detect_unstable_taxa(tree_set_path)
    unstable_clades = detect_unstable_clades(tree_set_path)
    clade_conflicts = summarize_clade_credibility_conflicts(tree_set_path)
    conclusion_summary = summarize_uncertainty_aware_conclusions(tree_set_path)
    storage_risk = assess_tree_set_storage_risk(tree_set_path)
    thinning_sensitivity = assess_tree_set_thinning_sensitivity(tree_set_path)
    consensus_sensitivity = compare_consensus_thresholds(tree_set_path)
    maturity = assess_tree_set_maturity(tree_set_path)
    benchmark = benchmark_tree_set_uncertainty(
        tree_counts=[summary.tree_count],
        taxon_counts=[max(len(summary.shared_taxa), 2)],
    )
    title = "Bijux Tree Uncertainty Report"
    sections = [
        _section("tree-set-summary", asdict(summary)),
        _section("consensus-tree", {"newick": dumps_newick(consensus_tree), "report": asdict(consensus)}),
        _section("clade-frequencies", asdict(clade_frequencies)),
        _section("pairwise-tree-distances", asdict(distances)),
        _section("topology-clusters", asdict(clusters)),
        _section("topological-diversity", asdict(diversity)),
        _section("topology-multimodality", asdict(multimodality)),
        _section("unstable-taxa", asdict(unstable_taxa)),
        _section("unstable-clades", asdict(unstable_clades)),
        _section("clade-credibility-conflicts", asdict(clade_conflicts)),
        _section("uncertainty-aware-conclusions", asdict(conclusion_summary)),
        _section("storage-risk", asdict(storage_risk)),
        _section("thinning-sensitivity", asdict(thinning_sensitivity)),
        _section("consensus-threshold-sensitivity", asdict(consensus_sensitivity)),
        _section("tree-set-benchmark", asdict(benchmark)),
        _section("maturity-gate", asdict(maturity)),
    ]
    core_sections = sections[:11]
    supplemental_sections = sections[11:]
    machine_manifest = {
        "report_kind": "tree-uncertainty",
        "title": title,
        "source_path": str(tree_set_path),
        "input_checksum": _sha256(tree_set_path),
        "tree_count": summary.tree_count,
        "rooted_topology_count": summary.rooted_topology_count,
        "sections": [name for name, _ in core_sections],
        "supplemental_sections": [name for name, _ in supplemental_sections],
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
    diversity = compare_posterior_topological_diversity(left_tree_set_path, right_tree_set_path)
    left_multimodality = detect_posterior_topology_multimodality(left_tree_set_path)
    right_multimodality = detect_posterior_topology_multimodality(right_tree_set_path)
    left_unstable_taxa = detect_unstable_taxa(left_tree_set_path)
    right_unstable_taxa = detect_unstable_taxa(right_tree_set_path)
    left_unstable_clades = detect_unstable_clades(left_tree_set_path)
    right_unstable_clades = detect_unstable_clades(right_tree_set_path)
    left_conflicts = summarize_clade_credibility_conflicts(left_tree_set_path)
    right_conflicts = summarize_clade_credibility_conflicts(right_tree_set_path)
    left_conclusions = summarize_uncertainty_aware_conclusions(left_tree_set_path)
    right_conclusions = summarize_uncertainty_aware_conclusions(right_tree_set_path)
    sections = [
        _section("tree-set-comparison", asdict(comparison)),
        _section("topological-diversity-comparison", asdict(diversity)),
        _section("left-tree-set-summary", asdict(left_summary)),
        _section("right-tree-set-summary", asdict(right_summary)),
        _section("left-topology-clusters", asdict(left_clusters)),
        _section("right-topology-clusters", asdict(right_clusters)),
        _section("left-topology-multimodality", asdict(left_multimodality)),
        _section("right-topology-multimodality", asdict(right_multimodality)),
        _section("left-unstable-taxa", asdict(left_unstable_taxa)),
        _section("right-unstable-taxa", asdict(right_unstable_taxa)),
        _section("left-unstable-clades", asdict(left_unstable_clades)),
        _section("right-unstable-clades", asdict(right_unstable_clades)),
        _section("left-clade-credibility-conflicts", asdict(left_conflicts)),
        _section("right-clade-credibility-conflicts", asdict(right_conflicts)),
        _section("left-uncertainty-aware-conclusions", asdict(left_conclusions)),
        _section("right-uncertainty-aware-conclusions", asdict(right_conclusions)),
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


def render_workflow_validation_report(
    *,
    out_path: Path,
    fixtures_root: Path | None = None,
) -> WorkflowValidationReportBuildResult:
    """Render the Level 1 workflow validation fixture report."""
    validation = build_core_workflow_validation_report(fixtures_root=fixtures_root)
    title = "Bijux Core Workflow Validation Report"
    reviewer_summary = [
        f"fixture checks passed: {validation.passed_fixture_count}/{validation.total_fixture_count}",
        f"validated workflow surfaces: {len(validation.workflows)}",
        f"known failure-gallery cases: {len(validation.failure_gallery)}",
    ]
    sections = [
        _section("reviewer-summary", reviewer_summary),
        _section(
            "validation-overview",
            {
                "total_fixture_count": validation.total_fixture_count,
                "passed_fixture_count": validation.passed_fixture_count,
                "failed_fixture_count": validation.failed_fixture_count,
            },
        ),
        _section("validation-suites", [asdict(suite) for suite in validation.suites]),
        _section("workflow-coverage", [asdict(row) for row in validation.workflows]),
        _section("failure-gallery", [asdict(row) for row in validation.failure_gallery]),
        _section("maturity-classification", [asdict(row) for row in validation.maturity_classifications]),
        _section("limitations", validation.limitations),
    ]
    fixture_paths = sorted(
        {
            path
            for suite in validation.suites
            for fixture in suite.fixtures
            for path in fixture.fixture_paths
        }
    )
    machine_manifest = {
        "report_kind": "workflow-validation",
        "title": title,
        "input_paths": [str(path) for path in fixture_paths],
        "input_checksums": {str(path): _sha256(path) for path in fixture_paths if path.exists()},
        "sections": [name for name, _ in sections],
        "metrics": {
            "total_fixture_count": validation.total_fixture_count,
            "passed_fixture_count": validation.passed_fixture_count,
            "workflow_count": len(validation.workflows),
        },
        "reviewer_summary": reviewer_summary,
        "limitations": validation.limitations,
    }
    machine_manifest_path = _write_machine_manifest(_report_sidecar_path(out_path), machine_manifest)
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return WorkflowValidationReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="workflow-validation",
        title=title,
        validation=validation,
        machine_manifest=machine_manifest,
    )


def render_level_one_release_gate_report(
    *,
    out_path: Path,
    fixtures_root: Path | None = None,
) -> ReleaseGateReportBuildResult:
    """Render the Level 1 release gate for the checked-in workflow fixtures."""
    release_gate = build_level_one_release_gate_report(fixtures_root=fixtures_root)
    title = "Bijux Level 1 Release Gate Report"
    reviewer_summary = [
        f"gate decision: {release_gate.gate.decision}",
        f"dataset readiness: {release_gate.dataset_readiness_decision}",
        f"retained taxa: {len(release_gate.gate.retained_taxa)}, excluded taxa: {len(release_gate.gate.excluded_taxa)}",
    ]
    sections = [
        _section("reviewer-summary", reviewer_summary),
        _section("gate-decision", asdict(release_gate.gate)),
        _section(
            "dataset-readiness",
            {
                "decision": release_gate.dataset_readiness_decision,
                "blockers": release_gate.dataset_blockers,
                "warnings": release_gate.dataset_warnings,
            },
        ),
        _section(
            "taxon-loss-traceability",
            {
                "first_loss_stage": release_gate.taxon_first_loss_stage,
                "exclusion_causes": release_gate.exclusion_causes,
            },
        ),
        _section("workflow-validation", asdict(release_gate.validation)),
        _section("limitations", release_gate.validation.limitations),
    ]
    fixture_paths = sorted(
        {
            path
            for suite in release_gate.validation.suites
            for fixture in suite.fixtures
            for path in fixture.fixture_paths
        }
    )
    machine_manifest = {
        "report_kind": "release-gate",
        "title": title,
        "input_paths": [str(path) for path in fixture_paths],
        "input_checksums": {str(path): _sha256(path) for path in fixture_paths if path.exists()},
        "sections": [name for name, _ in sections],
        "metrics": {
            "retained_taxa": len(release_gate.gate.retained_taxa),
            "excluded_taxa": len(release_gate.gate.excluded_taxa),
            "blocked_analysis_count": len(release_gate.gate.blocked_analyses),
        },
        "reviewer_summary": reviewer_summary,
        "limitations": release_gate.validation.limitations,
    }
    machine_manifest_path = _write_machine_manifest(_report_sidecar_path(out_path), machine_manifest)
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return ReleaseGateReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="release-gate",
        title=title,
        release_gate=release_gate,
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
    reviewer_summary, limitations = _report_summary_and_limitations(
        report_kind="dataset" if dataset_audit is not None else "tree",
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        dataset_audit=dataset_audit,
        alignment_forensic=alignment_forensic,
    )

    sections = [
        _section("reviewer-summary", reviewer_summary),
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
    sections.append(_section("limitations", limitations))

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
    machine_manifest["reviewer_summary"] = reviewer_summary
    machine_manifest["limitations"] = limitations
    input_ledger: list[ReportInputLedgerEntry] = []
    machine_manifest["input_ledger"] = _serialize_input_ledger(input_ledger)
    machine_manifest_path = _write_machine_manifest(_report_sidecar_path(out_path), machine_manifest)
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return ReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
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
        alignment_low_information=None,
        alignment_duplicate_policy=None,
        alignment_ambiguous_columns=None,
        alignment_sequence_ranking=None,
        alignment_coding=alignment_coding,
        alignment_identity_matrix=alignment_identity_matrix,
        alignment_linkage=None,
        dataset_readiness=dataset_readiness,
        dataset_audit=dataset_audit,
        input_ledger=input_ledger,
        machine_manifest=machine_manifest,
    )
