from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.alignment import (
    AlignmentLowInformationReport,
    AlignmentQualityReport,
    AlignmentSummary,
    SequenceQualityRankingReport,
    SequenceQualityRankingRow,
    SequenceUncertaintyProfile,
)
from bijux_phylogenetics.diagnostics.validation import (
    TreeForensicReport,
    TreeInspectionReport,
    TreeValidationReport,
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.core.dataset import (
    DatasetAuditReport,
    DatasetCompletenessRow,
    DatasetCrosswalkRow,
    DatasetExclusionRow,
    audit_dataset_inputs,
)
from bijux_phylogenetics.core.metadata import TaxonTable, load_taxon_table, write_taxon_rows
from bijux_phylogenetics.core.taxon_workflows import (
    TaxonWorkflowLossReport,
    TaxonWorkflowLossRow,
    build_taxon_workflow_loss_report,
)
from bijux_phylogenetics.io.fasta import (
    assess_alignment_low_information,
    build_alignment_quality_report,
    build_sequence_quality_ranking,
    compare_alignment_summaries,
    summarise_fasta,
)
from bijux_phylogenetics.reports.tree_package import (
    TreeBranchStatisticsRow,
    TreeSupportRow,
    summarize_tree_branch_statistics,
    summarize_tree_support,
)
from bijux_phylogenetics.trees import (
    analyze_branch_length_distribution,
    extract_tree_clades,
)


@dataclass(frozen=True, slots=True)
class SupplementaryTaxonTableRow:
    """One supplementary reviewer-facing row for a taxon across dataset surfaces."""

    taxon: str
    tree_tip_id: str | None
    alignment_id: str | None
    metadata_id: str | None
    trait_id: str | None
    tip_date_id: str | None
    geography_source: str | None
    calibration_targets: list[str]
    external_taxonomy_ids: dict[str, str]
    analysis_status: str
    analysis_exclusion_reason: str | None
    analysis_exclusion_causes: list[str]
    analysis_first_failed_surface: str | None
    affected_analyses: list[str]
    reporting_status: str
    reporting_loss_reason: str | None
    workflow_first_loss_stage: str | None
    workflow_loss_stages: list[str]
    workflow_loss_reasons: list[str]
    in_tree: bool
    in_alignment: bool
    in_metadata: bool
    in_traits: bool
    in_tip_dates: bool
    in_geography: bool
    in_calibrations: bool
    metadata_values: dict[str, str]
    trait_values: dict[str, str]


@dataclass(slots=True)
class SupplementaryTaxonTableResult:
    output_path: Path
    row_count: int
    analysis_included_count: int
    analysis_excluded_count: int
    reporting_retained_count: int
    reporting_dropped_count: int
    metadata_column_count: int
    trait_column_count: int
    columns: list[str]
    rows: list[SupplementaryTaxonTableRow]


@dataclass(frozen=True, slots=True)
class SupplementaryAlignmentDiagnosticsRow:
    """One reviewer-facing alignment diagnostics row for one sequence identifier."""

    sequence_id: str
    original_sequence_present: bool
    filtered_sequence_present: bool
    filtering_status: str
    filtering_reason: str | None
    original_missing_fraction: float | None
    original_gap_fraction: float | None
    original_ambiguity_fraction: float | None
    original_quality_score: float | None
    duplicate_status: str | None
    composition_outlier: bool | None
    original_alignment_length: int | None
    original_sequence_count: int | None
    original_missing_data_fraction: float | None
    original_gap_fraction_alignment: float | None
    original_variable_site_count: int | None
    original_parsimony_informative_site_count: int | None
    original_suspicious_alignment: bool | None
    original_low_information: bool | None
    original_low_information_reasons: list[str]
    filtered_alignment_length: int | None
    filtered_sequence_count: int | None
    filtered_missing_data_fraction: float | None
    filtered_gap_fraction_alignment: float | None
    filtered_variable_site_count: int | None
    filtered_parsimony_informative_site_count: int | None
    filtered_low_information: bool | None
    filtered_low_information_reasons: list[str]


@dataclass(slots=True)
class SupplementaryAlignmentDiagnosticsTableResult:
    output_path: Path
    row_count: int
    retained_sequence_count: int
    removed_sequence_count: int
    filtered_only_sequence_count: int
    columns: list[str]
    rows: list[SupplementaryAlignmentDiagnosticsRow]


@dataclass(frozen=True, slots=True)
class SupplementaryTreeDiagnosticsRow:
    """One reviewer-facing tree diagnostics row for one tree source."""

    tree_source: str
    source_format: str
    tip_count: int
    internal_node_count: int
    edge_count: int
    clade_count: int
    topology_shape: str
    is_binary: bool
    star_like: bool
    comb_like: bool
    polytomy_count: int
    polytomy_nodes: list[str]
    rooted: bool
    root_state_classification: str
    root_state_suspicious: bool
    branch_length_status: str
    has_complete_branch_lengths: bool
    total_branch_length: float
    minimum_branch_length: float | None
    maximum_branch_length: float | None
    mean_branch_length: float | None
    median_branch_length: float | None
    positive_branch_median: float | None
    missing_branch_count: int
    zero_length_branch_count: int
    negative_branch_count: int
    long_branch_outlier_count: int
    short_branch_outlier_count: int
    supported_branch_count: int
    strong_support_branch_count: int
    moderate_support_branch_count: int
    weak_support_branch_count: int
    missing_support_branch_count: int
    support_value_range_warnings: list[str]
    ultrametric: bool | None
    min_root_to_tip: float | None
    max_root_to_tip: float | None
    tree_diameter: float | None
    tree_quality_score: float
    safe_for_topology_comparison: bool
    safe_for_time_tree_analysis: bool
    safe_for_comparative_methods: bool
    safe_for_visualization: bool
    safe_for_publication: bool
    warning_count: int
    warnings: list[str]


@dataclass(slots=True)
class SupplementaryTreeDiagnosticsTableResult:
    output_path: Path
    row_count: int
    columns: list[str]
    rows: list[SupplementaryTreeDiagnosticsRow]


def _row_lookup(table: TaxonTable) -> dict[str, dict[str, str]]:
    return {row[table.taxon_column]: row for row in table.rows}


def _stringify_list(values: list[str]) -> str:
    return "|".join(values)


def _stringify_mapping(values: dict[str, str]) -> str:
    return "|".join(f"{key}={value}" for key, value in sorted(values.items()))


def _table_delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def _build_dynamic_columns(
    metadata_table: TaxonTable,
    traits_table: TaxonTable,
) -> tuple[list[str], list[str]]:
    metadata_columns = [
        f"metadata_{column}"
        for column in metadata_table.columns
        if column != metadata_table.taxon_column
    ]
    trait_columns = [
        f"trait_{column}"
        for column in traits_table.columns
        if column != traits_table.taxon_column
    ]
    return metadata_columns, trait_columns


def _table_columns(metadata_table: TaxonTable, traits_table: TaxonTable) -> list[str]:
    metadata_columns, trait_columns = _build_dynamic_columns(
        metadata_table, traits_table
    )
    return [
        "taxon",
        "tree_tip_id",
        "alignment_id",
        "metadata_id",
        "trait_id",
        "tip_date_id",
        "geography_source",
        "calibration_targets",
        "external_taxonomy_ids",
        "analysis_status",
        "analysis_exclusion_reason",
        "analysis_exclusion_causes",
        "analysis_first_failed_surface",
        "affected_analyses",
        "reporting_status",
        "reporting_loss_reason",
        "workflow_first_loss_stage",
        "workflow_loss_stages",
        "workflow_loss_reasons",
        "in_tree",
        "in_alignment",
        "in_metadata",
        "in_traits",
        "in_tip_dates",
        "in_geography",
        "in_calibrations",
        *metadata_columns,
        *trait_columns,
    ]


def _alignment_table_columns() -> list[str]:
    return [
        "sequence_id",
        "original_sequence_present",
        "filtered_sequence_present",
        "filtering_status",
        "filtering_reason",
        "original_missing_fraction",
        "original_gap_fraction",
        "original_ambiguity_fraction",
        "original_quality_score",
        "duplicate_status",
        "composition_outlier",
        "original_alignment_length",
        "original_sequence_count",
        "original_missing_data_fraction",
        "original_gap_fraction_alignment",
        "original_variable_site_count",
        "original_parsimony_informative_site_count",
        "original_suspicious_alignment",
        "original_low_information",
        "original_low_information_reasons",
        "filtered_alignment_length",
        "filtered_sequence_count",
        "filtered_missing_data_fraction",
        "filtered_gap_fraction_alignment",
        "filtered_variable_site_count",
        "filtered_parsimony_informative_site_count",
        "filtered_low_information",
        "filtered_low_information_reasons",
    ]


def _tree_table_columns() -> list[str]:
    return [
        "tree_source",
        "source_format",
        "tip_count",
        "internal_node_count",
        "edge_count",
        "clade_count",
        "topology_shape",
        "is_binary",
        "star_like",
        "comb_like",
        "polytomy_count",
        "polytomy_nodes",
        "rooted",
        "root_state_classification",
        "root_state_suspicious",
        "branch_length_status",
        "has_complete_branch_lengths",
        "total_branch_length",
        "minimum_branch_length",
        "maximum_branch_length",
        "mean_branch_length",
        "median_branch_length",
        "positive_branch_median",
        "missing_branch_count",
        "zero_length_branch_count",
        "negative_branch_count",
        "long_branch_outlier_count",
        "short_branch_outlier_count",
        "supported_branch_count",
        "strong_support_branch_count",
        "moderate_support_branch_count",
        "weak_support_branch_count",
        "missing_support_branch_count",
        "support_value_range_warnings",
        "ultrametric",
        "min_root_to_tip",
        "max_root_to_tip",
        "tree_diameter",
        "tree_quality_score",
        "safe_for_topology_comparison",
        "safe_for_time_tree_analysis",
        "safe_for_comparative_methods",
        "safe_for_visualization",
        "safe_for_publication",
        "warning_count",
        "warnings",
    ]


def _serialize_row(
    row: SupplementaryTaxonTableRow,
    *,
    metadata_table: TaxonTable,
    traits_table: TaxonTable,
) -> dict[str, object]:
    serialized: dict[str, object] = {
        "taxon": row.taxon,
        "tree_tip_id": row.tree_tip_id or "",
        "alignment_id": row.alignment_id or "",
        "metadata_id": row.metadata_id or "",
        "trait_id": row.trait_id or "",
        "tip_date_id": row.tip_date_id or "",
        "geography_source": row.geography_source or "",
        "calibration_targets": _stringify_list(row.calibration_targets),
        "external_taxonomy_ids": _stringify_mapping(row.external_taxonomy_ids),
        "analysis_status": row.analysis_status,
        "analysis_exclusion_reason": row.analysis_exclusion_reason or "",
        "analysis_exclusion_causes": _stringify_list(row.analysis_exclusion_causes),
        "analysis_first_failed_surface": row.analysis_first_failed_surface or "",
        "affected_analyses": _stringify_list(row.affected_analyses),
        "reporting_status": row.reporting_status,
        "reporting_loss_reason": row.reporting_loss_reason or "",
        "workflow_first_loss_stage": row.workflow_first_loss_stage or "",
        "workflow_loss_stages": _stringify_list(row.workflow_loss_stages),
        "workflow_loss_reasons": _stringify_list(row.workflow_loss_reasons),
        "in_tree": row.in_tree,
        "in_alignment": row.in_alignment,
        "in_metadata": row.in_metadata,
        "in_traits": row.in_traits,
        "in_tip_dates": row.in_tip_dates,
        "in_geography": row.in_geography,
        "in_calibrations": row.in_calibrations,
    }
    for column in metadata_table.columns:
        if column == metadata_table.taxon_column:
            continue
        serialized[f"metadata_{column}"] = row.metadata_values.get(column, "")
    for column in traits_table.columns:
        if column == traits_table.taxon_column:
            continue
        serialized[f"trait_{column}"] = row.trait_values.get(column, "")
    return serialized


def _serialize_alignment_row(
    row: SupplementaryAlignmentDiagnosticsRow,
) -> dict[str, object]:
    return {
        "sequence_id": row.sequence_id,
        "original_sequence_present": row.original_sequence_present,
        "filtered_sequence_present": row.filtered_sequence_present,
        "filtering_status": row.filtering_status,
        "filtering_reason": row.filtering_reason or "",
        "original_missing_fraction": ""
        if row.original_missing_fraction is None
        else row.original_missing_fraction,
        "original_gap_fraction": ""
        if row.original_gap_fraction is None
        else row.original_gap_fraction,
        "original_ambiguity_fraction": ""
        if row.original_ambiguity_fraction is None
        else row.original_ambiguity_fraction,
        "original_quality_score": ""
        if row.original_quality_score is None
        else row.original_quality_score,
        "duplicate_status": row.duplicate_status or "",
        "composition_outlier": ""
        if row.composition_outlier is None
        else row.composition_outlier,
        "original_alignment_length": ""
        if row.original_alignment_length is None
        else row.original_alignment_length,
        "original_sequence_count": ""
        if row.original_sequence_count is None
        else row.original_sequence_count,
        "original_missing_data_fraction": ""
        if row.original_missing_data_fraction is None
        else row.original_missing_data_fraction,
        "original_gap_fraction_alignment": ""
        if row.original_gap_fraction_alignment is None
        else row.original_gap_fraction_alignment,
        "original_variable_site_count": ""
        if row.original_variable_site_count is None
        else row.original_variable_site_count,
        "original_parsimony_informative_site_count": ""
        if row.original_parsimony_informative_site_count is None
        else row.original_parsimony_informative_site_count,
        "original_suspicious_alignment": ""
        if row.original_suspicious_alignment is None
        else row.original_suspicious_alignment,
        "original_low_information": ""
        if row.original_low_information is None
        else row.original_low_information,
        "original_low_information_reasons": _stringify_list(
            row.original_low_information_reasons
        ),
        "filtered_alignment_length": ""
        if row.filtered_alignment_length is None
        else row.filtered_alignment_length,
        "filtered_sequence_count": ""
        if row.filtered_sequence_count is None
        else row.filtered_sequence_count,
        "filtered_missing_data_fraction": ""
        if row.filtered_missing_data_fraction is None
        else row.filtered_missing_data_fraction,
        "filtered_gap_fraction_alignment": ""
        if row.filtered_gap_fraction_alignment is None
        else row.filtered_gap_fraction_alignment,
        "filtered_variable_site_count": ""
        if row.filtered_variable_site_count is None
        else row.filtered_variable_site_count,
        "filtered_parsimony_informative_site_count": ""
        if row.filtered_parsimony_informative_site_count is None
        else row.filtered_parsimony_informative_site_count,
        "filtered_low_information": ""
        if row.filtered_low_information is None
        else row.filtered_low_information,
        "filtered_low_information_reasons": _stringify_list(
            row.filtered_low_information_reasons
        ),
    }


def _serialize_tree_row(
    row: SupplementaryTreeDiagnosticsRow,
) -> dict[str, object]:
    return {
        "tree_source": row.tree_source,
        "source_format": row.source_format,
        "tip_count": row.tip_count,
        "internal_node_count": row.internal_node_count,
        "edge_count": row.edge_count,
        "clade_count": row.clade_count,
        "topology_shape": row.topology_shape,
        "is_binary": row.is_binary,
        "star_like": row.star_like,
        "comb_like": row.comb_like,
        "polytomy_count": row.polytomy_count,
        "polytomy_nodes": _stringify_list(row.polytomy_nodes),
        "rooted": row.rooted,
        "root_state_classification": row.root_state_classification,
        "root_state_suspicious": row.root_state_suspicious,
        "branch_length_status": row.branch_length_status,
        "has_complete_branch_lengths": row.has_complete_branch_lengths,
        "total_branch_length": row.total_branch_length,
        "minimum_branch_length": ""
        if row.minimum_branch_length is None
        else row.minimum_branch_length,
        "maximum_branch_length": ""
        if row.maximum_branch_length is None
        else row.maximum_branch_length,
        "mean_branch_length": ""
        if row.mean_branch_length is None
        else row.mean_branch_length,
        "median_branch_length": ""
        if row.median_branch_length is None
        else row.median_branch_length,
        "positive_branch_median": ""
        if row.positive_branch_median is None
        else row.positive_branch_median,
        "missing_branch_count": row.missing_branch_count,
        "zero_length_branch_count": row.zero_length_branch_count,
        "negative_branch_count": row.negative_branch_count,
        "long_branch_outlier_count": row.long_branch_outlier_count,
        "short_branch_outlier_count": row.short_branch_outlier_count,
        "supported_branch_count": row.supported_branch_count,
        "strong_support_branch_count": row.strong_support_branch_count,
        "moderate_support_branch_count": row.moderate_support_branch_count,
        "weak_support_branch_count": row.weak_support_branch_count,
        "missing_support_branch_count": row.missing_support_branch_count,
        "support_value_range_warnings": _stringify_list(
            row.support_value_range_warnings
        ),
        "ultrametric": "" if row.ultrametric is None else row.ultrametric,
        "min_root_to_tip": ""
        if row.min_root_to_tip is None
        else row.min_root_to_tip,
        "max_root_to_tip": ""
        if row.max_root_to_tip is None
        else row.max_root_to_tip,
        "tree_diameter": "" if row.tree_diameter is None else row.tree_diameter,
        "tree_quality_score": row.tree_quality_score,
        "safe_for_topology_comparison": row.safe_for_topology_comparison,
        "safe_for_time_tree_analysis": row.safe_for_time_tree_analysis,
        "safe_for_comparative_methods": row.safe_for_comparative_methods,
        "safe_for_visualization": row.safe_for_visualization,
        "safe_for_publication": row.safe_for_publication,
        "warning_count": row.warning_count,
        "warnings": _stringify_list(row.warnings),
    }


def _build_supplementary_row(
    *,
    taxon: str,
    analysis_taxa: set[str],
    crosswalk: DatasetCrosswalkRow | None,
    completeness: DatasetCompletenessRow | None,
    exclusion: DatasetExclusionRow | None,
    workflow_loss: TaxonWorkflowLossRow | None,
    metadata_values: dict[str, str],
    trait_values: dict[str, str],
) -> SupplementaryTaxonTableRow:
    workflow_events = [] if workflow_loss is None else workflow_loss.loss_events
    analysis_exclusion_reason = (
        None if exclusion is None or not exclusion.causes else exclusion.causes[0]
    )
    reporting_loss_reason = (
        None
        if not workflow_events or workflow_loss is None or workflow_loss.retained_for_reporting
        else f"{workflow_events[0].stage}:{workflow_events[0].reason}"
    )
    return SupplementaryTaxonTableRow(
        taxon=taxon,
        tree_tip_id=None if crosswalk is None else crosswalk.tree_tip,
        alignment_id=None if crosswalk is None else crosswalk.alignment_id,
        metadata_id=None if crosswalk is None else crosswalk.metadata_id,
        trait_id=None if crosswalk is None else crosswalk.trait_id,
        tip_date_id=None if crosswalk is None else crosswalk.tip_date_id,
        geography_source=None if crosswalk is None else crosswalk.geography_source,
        calibration_targets=[] if crosswalk is None else crosswalk.calibration_targets,
        external_taxonomy_ids={}
        if crosswalk is None
        else crosswalk.external_taxonomy_ids,
        analysis_status="included" if taxon in analysis_taxa else "excluded",
        analysis_exclusion_reason=analysis_exclusion_reason,
        analysis_exclusion_causes=[] if exclusion is None else exclusion.causes,
        analysis_first_failed_surface=None
        if exclusion is None
        else exclusion.first_failed_surface,
        affected_analyses=[] if exclusion is None else exclusion.affected_analyses,
        reporting_status="retained"
        if workflow_loss is None or workflow_loss.retained_for_reporting
        else "dropped",
        reporting_loss_reason=reporting_loss_reason,
        workflow_first_loss_stage=None
        if workflow_loss is None
        else workflow_loss.first_loss_stage,
        workflow_loss_stages=[event.stage for event in workflow_events],
        workflow_loss_reasons=[event.reason for event in workflow_events],
        in_tree=False if completeness is None else completeness.in_tree,
        in_alignment=False if completeness is None else completeness.in_alignment,
        in_metadata=False if completeness is None else completeness.in_metadata,
        in_traits=False if completeness is None else completeness.in_traits,
        in_tip_dates=False if completeness is None else completeness.in_tip_dates,
        in_geography=False if completeness is None else completeness.in_geography,
        in_calibrations=False
        if completeness is None
        else completeness.in_calibrations,
        metadata_values=metadata_values,
        trait_values=trait_values,
    )


def _sequence_uncertainty_lookup(
    summary: AlignmentSummary,
) -> dict[str, SequenceUncertaintyProfile]:
    return {row.identifier: row for row in summary.per_sequence_uncertainty}


def _sequence_quality_lookup(
    ranking: SequenceQualityRankingReport,
) -> dict[str, SequenceQualityRankingRow]:
    return {row.identifier: row for row in ranking.rows}


def _filtering_status(
    *,
    sequence_id: str,
    original_ids: set[str],
    filtered_ids: set[str] | None,
) -> tuple[str, str | None]:
    if filtered_ids is None:
        return "not_requested", None
    if sequence_id in original_ids and sequence_id in filtered_ids:
        return "retained_after_filtering", None
    if sequence_id in original_ids and sequence_id not in filtered_ids:
        return "removed_during_filtering", "absent_from_filtered_alignment"
    return "only_in_filtered_alignment", "absent_from_original_alignment"


def _alignment_sequence_order(
    original_summary: AlignmentSummary,
    filtered_summary: AlignmentSummary | None,
) -> list[str]:
    ordered = list(original_summary.ids)
    if filtered_summary is None:
        return ordered
    filtered_only = [
        identifier for identifier in filtered_summary.ids if identifier not in set(ordered)
    ]
    return [*ordered, *filtered_only]


def _build_alignment_row(
    *,
    sequence_id: str,
    original_summary: AlignmentSummary,
    original_quality: AlignmentQualityReport,
    original_low_information: AlignmentLowInformationReport,
    original_quality_lookup: dict[str, SequenceQualityRankingRow],
    original_uncertainty_lookup: dict[str, SequenceUncertaintyProfile],
    filtered_summary: AlignmentSummary | None,
    filtered_low_information: AlignmentLowInformationReport | None,
) -> SupplementaryAlignmentDiagnosticsRow:
    original_ids = set(original_summary.ids)
    filtered_ids = None if filtered_summary is None else set(filtered_summary.ids)
    filtering_status, filtering_reason = _filtering_status(
        sequence_id=sequence_id,
        original_ids=original_ids,
        filtered_ids=filtered_ids,
    )
    original_uncertainty = original_uncertainty_lookup.get(sequence_id)
    original_ranking = original_quality_lookup.get(sequence_id)
    return SupplementaryAlignmentDiagnosticsRow(
        sequence_id=sequence_id,
        original_sequence_present=sequence_id in original_ids,
        filtered_sequence_present=False
        if filtered_ids is None
        else sequence_id in filtered_ids,
        filtering_status=filtering_status,
        filtering_reason=filtering_reason,
        original_missing_fraction=None
        if original_uncertainty is None
        else original_uncertainty.missing_fraction,
        original_gap_fraction=None
        if original_uncertainty is None
        else original_uncertainty.gap_fraction,
        original_ambiguity_fraction=None
        if original_uncertainty is None
        else original_uncertainty.ambiguity_fraction,
        original_quality_score=None
        if original_ranking is None
        else original_ranking.score,
        duplicate_status=None if original_ranking is None else original_ranking.duplicate_status,
        composition_outlier=None
        if original_ranking is None
        else original_ranking.composition_outlier,
        original_alignment_length=original_summary.alignment_length,
        original_sequence_count=original_summary.sequence_count,
        original_missing_data_fraction=original_summary.missing_data_fraction,
        original_gap_fraction_alignment=original_summary.gap_fraction,
        original_variable_site_count=original_summary.variable_site_count,
        original_parsimony_informative_site_count=original_summary.parsimony_informative_site_count,
        original_suspicious_alignment=original_quality.suspicious_alignment,
        original_low_information=original_low_information.low_information,
        original_low_information_reasons=original_low_information.reasons,
        filtered_alignment_length=None
        if filtered_summary is None
        else filtered_summary.alignment_length,
        filtered_sequence_count=None
        if filtered_summary is None
        else filtered_summary.sequence_count,
        filtered_missing_data_fraction=None
        if filtered_summary is None
        else filtered_summary.missing_data_fraction,
        filtered_gap_fraction_alignment=None
        if filtered_summary is None
        else filtered_summary.gap_fraction,
        filtered_variable_site_count=None
        if filtered_summary is None
        else filtered_summary.variable_site_count,
        filtered_parsimony_informative_site_count=None
        if filtered_summary is None
        else filtered_summary.parsimony_informative_site_count,
        filtered_low_information=None
        if filtered_low_information is None
        else filtered_low_information.low_information,
        filtered_low_information_reasons=[]
        if filtered_low_information is None
        else filtered_low_information.reasons,
    )


def _write_alignment_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryAlignmentDiagnosticsRow],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter=_table_delimiter(path),
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(_serialize_alignment_row(row))
    return path


def _support_counts(rows: list[TreeSupportRow]) -> dict[str, int]:
    counts = {"strong": 0, "moderate": 0, "weak": 0, "missing": 0}
    for row in rows:
        counts[row.support_class] = counts.get(row.support_class, 0) + 1
    return counts


def _topology_shape(inspection: TreeInspectionReport) -> str:
    if inspection.star_like:
        return "star"
    if inspection.comb_like:
        return "comb"
    if inspection.is_binary:
        return "binary"
    if inspection.polytomy_count:
        return "polytomy"
    return "mixed"


def _tree_warning_ledger(
    *,
    validation: TreeValidationReport,
    inspection: TreeInspectionReport,
    forensic: TreeForensicReport,
) -> list[str]:
    return sorted(
        dict.fromkeys(
            [
                *validation.warnings,
                *inspection.warnings,
                *forensic.warnings,
            ]
        )
    )


def _build_tree_forensic_review(
    *,
    tree_path: Path,
    validation: TreeValidationReport,
    inspection: TreeInspectionReport,
) -> TreeForensicReport:
    context_lookup = {
        context.context: context for context in validation.branch_length_contexts
    }
    safe_for_topology_comparison = (
        validation.syntax_valid
        and not validation.duplicate_taxa
        and validation.missing_taxa == 0
    )
    safe_for_time_tree_analysis = (
        context_lookup["time_tree"].allowed and validation.biologically_safe
    )
    safe_for_comparative_methods = (
        context_lookup["comparative_methods"].allowed and validation.biologically_safe
    )
    safe_for_visualization = validation.syntax_valid
    safe_for_publication = (
        validation.biologically_safe and not inspection.internal_label_conflicts
    )
    warnings = sorted(dict.fromkeys([*validation.warnings, *inspection.warnings]))
    return TreeForensicReport(
        path=tree_path,
        source_format=validation.source_format,
        syntax_valid=validation.syntax_valid,
        biologically_safe=validation.biologically_safe,
        validity_decision=validation.validity_decision,
        integrity_issues=validation.integrity_issues,
        findings=validation.warning_details,
        root_state_confidence=validation.root_state_confidence,
        branch_length_contexts=validation.branch_length_contexts,
        branch_length_repair_suggestions=validation.branch_length_repair_suggestions,
        internal_label_conflicts=validation.internal_label_conflicts,
        stable_node_identities=validation.stable_node_identities,
        unsafe_external_labels=validation.unsafe_external_labels,
        taxon_identity_audit=validation.taxon_identity_audit,
        safe_for_topology_comparison=safe_for_topology_comparison,
        safe_for_time_tree_analysis=safe_for_time_tree_analysis,
        safe_for_comparative_methods=safe_for_comparative_methods,
        safe_for_visualization=safe_for_visualization,
        safe_for_publication=safe_for_publication,
        warnings=warnings,
    )


def _build_tree_row(
    *,
    tree_path: Path,
    validation: TreeValidationReport,
    inspection: TreeInspectionReport,
    forensic: TreeForensicReport,
    support_rows: list[TreeSupportRow],
    branch_stats: TreeBranchStatisticsRow,
) -> SupplementaryTreeDiagnosticsRow:
    warnings = _tree_warning_ledger(
        validation=validation,
        inspection=inspection,
        forensic=forensic,
    )
    support_counts = _support_counts(support_rows)
    supported_branch_count = sum(
        1 for row in support_rows if row.support is not None
    )
    return SupplementaryTreeDiagnosticsRow(
        tree_source=str(tree_path),
        source_format=inspection.source_format,
        tip_count=inspection.tip_count,
        internal_node_count=inspection.internal_node_count,
        edge_count=inspection.edge_count,
        clade_count=inspection.clade_count,
        topology_shape=_topology_shape(inspection),
        is_binary=inspection.is_binary,
        star_like=inspection.star_like,
        comb_like=inspection.comb_like,
        polytomy_count=inspection.polytomy_count,
        polytomy_nodes=inspection.polytomy_nodes,
        rooted=inspection.rooted,
        root_state_classification=inspection.root_state_confidence.classification,
        root_state_suspicious=inspection.root_state_confidence.suspicious_placement,
        branch_length_status=inspection.branch_length_status,
        has_complete_branch_lengths=validation.has_complete_branch_lengths,
        total_branch_length=inspection.total_branch_length,
        minimum_branch_length=branch_stats.minimum_branch_length,
        maximum_branch_length=branch_stats.maximum_branch_length,
        mean_branch_length=branch_stats.mean_branch_length,
        median_branch_length=branch_stats.median_branch_length,
        positive_branch_median=branch_stats.positive_branch_median,
        missing_branch_count=branch_stats.missing_branch_count,
        zero_length_branch_count=branch_stats.zero_length_branch_count,
        negative_branch_count=branch_stats.negative_branch_count,
        long_branch_outlier_count=branch_stats.long_outlier_count,
        short_branch_outlier_count=branch_stats.short_outlier_count,
        supported_branch_count=supported_branch_count,
        strong_support_branch_count=support_counts["strong"],
        moderate_support_branch_count=support_counts["moderate"],
        weak_support_branch_count=support_counts["weak"],
        missing_support_branch_count=support_counts["missing"],
        support_value_range_warnings=inspection.suspicious_support_value_ranges,
        ultrametric=inspection.is_ultrametric,
        min_root_to_tip=inspection.min_root_to_tip,
        max_root_to_tip=inspection.max_root_to_tip,
        tree_diameter=inspection.tree_diameter,
        tree_quality_score=inspection.tree_quality_score,
        safe_for_topology_comparison=forensic.safe_for_topology_comparison,
        safe_for_time_tree_analysis=forensic.safe_for_time_tree_analysis,
        safe_for_comparative_methods=forensic.safe_for_comparative_methods,
        safe_for_visualization=forensic.safe_for_visualization,
        safe_for_publication=forensic.safe_for_publication,
        warning_count=len(warnings),
        warnings=warnings,
    )


def _write_tree_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryTreeDiagnosticsRow],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter=_table_delimiter(path),
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(_serialize_tree_row(row))
    return path


def _build_taxon_rows(
    *,
    dataset_audit: DatasetAuditReport,
    workflow_loss: TaxonWorkflowLossReport,
    metadata_table: TaxonTable,
    traits_table: TaxonTable,
) -> list[SupplementaryTaxonTableRow]:
    analysis_taxa = set(dataset_audit.analysis_taxa)
    crosswalk_by_taxon = {
        row.taxon: row for row in dataset_audit.crosswalk.rows
    }
    completeness_by_taxon = {
        row.taxon: row for row in dataset_audit.completeness_matrix.rows
    }
    exclusion_by_taxon = {
        row.taxon: row for row in dataset_audit.exclusion_table.rows
    }
    workflow_loss_by_taxon = {row.taxon: row for row in workflow_loss.rows}
    metadata_by_taxon = _row_lookup(metadata_table)
    traits_by_taxon = _row_lookup(traits_table)
    taxa = sorted(
        set(crosswalk_by_taxon)
        | set(completeness_by_taxon)
        | set(exclusion_by_taxon)
        | set(workflow_loss_by_taxon)
        | set(metadata_by_taxon)
        | set(traits_by_taxon)
    )
    return [
        _build_supplementary_row(
            taxon=taxon,
            analysis_taxa=analysis_taxa,
            crosswalk=crosswalk_by_taxon.get(taxon),
            completeness=completeness_by_taxon.get(taxon),
            exclusion=exclusion_by_taxon.get(taxon),
            workflow_loss=workflow_loss_by_taxon.get(taxon),
            metadata_values=metadata_by_taxon.get(taxon, {}),
            trait_values=traits_by_taxon.get(taxon, {}),
        )
        for taxon in taxa
    ]


def write_supplementary_taxon_table(
    path: Path,
    *,
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    alignment_path: Path | None = None,
    filtered_alignment_path: Path | None = None,
    inference_tree_path: Path | None = None,
    reported_taxa_path: Path | None = None,
    tip_dates_path: Path | None = None,
    calibration_path: Path | None = None,
) -> SupplementaryTaxonTableResult:
    """Write one supplementary taxon table with IDs, annotations, and exclusion evidence."""
    dataset_audit = audit_dataset_inputs(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
    )
    workflow_loss = build_taxon_workflow_loss_report(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        filtered_alignment_path=filtered_alignment_path,
        inference_tree_path=inference_tree_path,
        reported_taxa_path=reported_taxa_path,
    )
    metadata_table = load_taxon_table(metadata_path)
    traits_table = load_taxon_table(traits_path)
    rows = _build_taxon_rows(
        dataset_audit=dataset_audit,
        workflow_loss=workflow_loss,
        metadata_table=metadata_table,
        traits_table=traits_table,
    )
    columns = _table_columns(metadata_table, traits_table)
    write_taxon_rows(
        path,
        columns=columns,
        rows=[
            _serialize_row(
                row,
                metadata_table=metadata_table,
                traits_table=traits_table,
            )
            for row in rows
        ],
    )
    analysis_included_count = sum(1 for row in rows if row.analysis_status == "included")
    reporting_retained_count = sum(
        1 for row in rows if row.reporting_status == "retained"
    )
    return SupplementaryTaxonTableResult(
        output_path=path,
        row_count=len(rows),
        analysis_included_count=analysis_included_count,
        analysis_excluded_count=len(rows) - analysis_included_count,
        reporting_retained_count=reporting_retained_count,
        reporting_dropped_count=len(rows) - reporting_retained_count,
        metadata_column_count=len(metadata_table.columns) - 1,
        trait_column_count=len(traits_table.columns) - 1,
        columns=columns,
        rows=rows,
    )


def write_supplementary_alignment_diagnostics_table(
    path: Path,
    *,
    alignment_path: Path,
    filtered_alignment_path: Path | None = None,
) -> SupplementaryAlignmentDiagnosticsTableResult:
    """Write one supplementary alignment diagnostics table with optional filtering outcomes."""
    original_summary = summarise_fasta(alignment_path)
    original_quality = build_alignment_quality_report(alignment_path)
    original_low_information = assess_alignment_low_information(alignment_path)
    original_ranking = build_sequence_quality_ranking(alignment_path)
    filtered_summary = (
        None if filtered_alignment_path is None else summarise_fasta(filtered_alignment_path)
    )
    filtered_low_information = (
        None
        if filtered_alignment_path is None
        else assess_alignment_low_information(filtered_alignment_path)
    )
    if filtered_summary is not None:
        compare_alignment_summaries(alignment_path, original_summary, filtered_summary)
    sequence_ids = _alignment_sequence_order(original_summary, filtered_summary)
    original_uncertainty_lookup = _sequence_uncertainty_lookup(original_summary)
    original_quality_lookup = _sequence_quality_lookup(original_ranking)
    rows = [
        _build_alignment_row(
            sequence_id=sequence_id,
            original_summary=original_summary,
            original_quality=original_quality,
            original_low_information=original_low_information,
            original_quality_lookup=original_quality_lookup,
            original_uncertainty_lookup=original_uncertainty_lookup,
            filtered_summary=filtered_summary,
            filtered_low_information=filtered_low_information,
        )
        for sequence_id in sequence_ids
    ]
    columns = _alignment_table_columns()
    _write_alignment_rows(path, columns=columns, rows=rows)
    retained_sequence_count = sum(
        1 for row in rows if row.filtering_status == "retained_after_filtering"
    )
    removed_sequence_count = sum(
        1 for row in rows if row.filtering_status == "removed_during_filtering"
    )
    filtered_only_sequence_count = sum(
        1 for row in rows if row.filtering_status == "only_in_filtered_alignment"
    )
    return SupplementaryAlignmentDiagnosticsTableResult(
        output_path=path,
        row_count=len(rows),
        retained_sequence_count=retained_sequence_count,
        removed_sequence_count=removed_sequence_count,
        filtered_only_sequence_count=filtered_only_sequence_count,
        columns=columns,
        rows=rows,
    )


def write_supplementary_tree_diagnostics_table(
    path: Path,
    *,
    tree_path: Path,
) -> SupplementaryTreeDiagnosticsTableResult:
    """Write one supplementary tree diagnostics table with topology and warning summaries."""
    validation = validate_tree_path(
        tree_path,
        allow_duplicates=True,
        allow_negative_branch_lengths=True,
    )
    inspection = inspect_tree_path(tree_path)
    forensic = _build_tree_forensic_review(
        tree_path=tree_path,
        validation=validation,
        inspection=inspection,
    )
    clades = extract_tree_clades(tree_path)
    support_rows = summarize_tree_support(clades)
    branch_lengths = analyze_branch_length_distribution(tree_path)
    branch_stats = summarize_tree_branch_statistics(branch_lengths)
    rows = [
        _build_tree_row(
            tree_path=tree_path,
            validation=validation,
            inspection=inspection,
            forensic=forensic,
            support_rows=support_rows,
            branch_stats=branch_stats,
        )
    ]
    columns = _tree_table_columns()
    _write_tree_rows(path, columns=columns, rows=rows)
    return SupplementaryTreeDiagnosticsTableResult(
        output_path=path,
        row_count=len(rows),
        columns=columns,
        rows=rows,
    )
