from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile

from bijux_phylogenetics.ancestral import (
    ContinuousAncestralReport,
    DiscreteAncestralReport,
    continuous_ancestral_exclusions,
    discrete_ancestral_exclusions,
    reconstruct_continuous_ancestral_states,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.comparative.pgls import run_pgls
from bijux_phylogenetics.comparative.phylogenetic_logistic import (
    summarize_phylogenetic_logistic,
)
from bijux_phylogenetics.comparative.diversification import (
    CladeDiversificationObservation,
    DiversificationRateReport,
    SamplingFractionIssue,
    SamplingFractionReport,
    compare_diversification_models,
    detect_diversification_outlier_clades,
    detect_incomplete_taxon_sampling_metadata,
    estimate_diversification_rate,
)
from bijux_phylogenetics.comparative.regression_model_selection import (
    ComparativeRegressionModelExclusion,
    ComparativeRegressionModelRow,
    ComparativeRegressionModelSelectionReport,
    compare_comparative_regression_models,
)
from bijux_phylogenetics.core.alignment import (
    AlignmentLowInformationReport,
    AlignmentQualityReport,
    AlignmentSummary,
    SequenceQualityRankingReport,
    SequenceQualityRankingRow,
    SequenceUncertaintyProfile,
)
from bijux_phylogenetics.engines.iqtree_artifacts import (
    IqtreeModelCandidate,
    IqtreeModelSelectionSummary,
    parse_iqtree_model_selection_summary,
    resolve_iqtree_model_sidecar,
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
from bijux_phylogenetics.core.metadata import (
    TaxonTable,
    load_taxon_table,
    write_taxon_rows,
)
from bijux_phylogenetics.core.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.core.taxon_workflows import (
    TaxonWorkflowLossReport,
    TaxonWorkflowLossRow,
    build_taxon_workflow_loss_report,
)
from bijux_phylogenetics.io.fasta.cleaning import compare_alignment_summaries
from bijux_phylogenetics.io.fasta.quality import (
    assess_alignment_low_information,
    build_alignment_quality_report,
    build_sequence_quality_ranking,
)
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.reports.tree_package import (
    TreeBranchStatisticsRow,
    TreeSupportRow,
    summarize_tree_branch_statistics,
    summarize_tree_support,
)
from bijux_phylogenetics.trees import (
    analyze_branch_length_distribution,
    compute_reference_tree_clade_support,
    extract_tree_clades,
    TreeSetCladeSupportReport,
    TreeSetCladeSupportRow,
)
from bijux_phylogenetics.io.newick import dumps_newick
from .columns import (
    alignment_table_columns as _alignment_table_columns,
    ancestral_state_table_columns as _ancestral_state_table_columns,
    batch_summary_table_columns as _batch_summary_table_columns,
    clade_support_table_columns as _clade_support_table_columns,
    comparative_model_table_columns as _comparative_model_table_columns,
    diversification_table_columns as _diversification_table_columns,
    model_selection_table_columns as _model_selection_table_columns,
    taxon_table_columns as _table_columns,
    tree_table_columns as _tree_table_columns,
)
from .models import (
    SupplementaryAlignmentDiagnosticsRow,
    SupplementaryAlignmentDiagnosticsTableResult,
    SupplementaryAncestralStateRow,
    SupplementaryAncestralStateTableResult,
    SupplementaryBatchSummaryRow,
    SupplementaryBatchSummaryTableResult,
    SupplementaryCladeSupportRow,
    SupplementaryCladeSupportTableResult,
    SupplementaryComparativeModelRow,
    SupplementaryComparativeModelTableResult,
    SupplementaryDiversificationRow,
    SupplementaryDiversificationTableResult,
    SupplementaryModelSelectionRow,
    SupplementaryModelSelectionTableResult,
    SupplementaryTaxonTableResult,
    SupplementaryTaxonTableRow,
    SupplementaryTreeDiagnosticsRow,
    SupplementaryTreeDiagnosticsTableResult,
)
from .shared import (
    row_lookup as _row_lookup,
    stringify_list as _stringify_list,
    stringify_mapping as _stringify_mapping,
    table_delimiter as _table_delimiter,
    write_dict_rows as _write_dict_rows,
)


def _serialize_row(
    row: SupplementaryTaxonTableRow,
    *,
    metadata_table: TaxonTable,
    traits_table: TaxonTable,
) -> dict[str, object]:
    from .taxon import _serialize_row as serialize_impl

    return serialize_impl(
        row,
        metadata_table=metadata_table,
        traits_table=traits_table,
    )


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
    from .tree_diagnostics import _serialize_tree_row as serialize_impl

    return serialize_impl(row)


def _serialize_clade_support_row(
    row: SupplementaryCladeSupportRow,
) -> dict[str, object]:
    return {
        "tree_source": row.tree_source,
        "comparison_tree_set_source": ""
        if row.comparison_tree_set_source is None
        else row.comparison_tree_set_source,
        "clade_id": row.clade_id,
        "node_kind": row.node_kind,
        "node_label": "" if row.node_label is None else row.node_label,
        "descendant_taxa": _stringify_list(row.descendant_taxa),
        "support": "" if row.support is None else row.support,
        "support_fraction": "" if row.support_fraction is None else row.support_fraction,
        "support_class": row.support_class,
        "support_method": row.support_method,
        "branch_length": "" if row.branch_length is None else row.branch_length,
        "root_depth": "" if row.root_depth is None else row.root_depth,
        "supporting_tree_count": ""
        if row.supporting_tree_count is None
        else row.supporting_tree_count,
        "clade_frequency": ""
        if row.clade_frequency is None
        else row.clade_frequency,
        "support_percent": ""
        if row.support_percent is None
        else row.support_percent,
        "frequency_method": ""
        if row.frequency_method is None
        else row.frequency_method,
        "frequency_status": ""
        if row.frequency_status is None
        else row.frequency_status,
        "frequency_explanation": ""
        if row.frequency_explanation is None
        else row.frequency_explanation,
    }


def _serialize_comparative_model_row(
    row: SupplementaryComparativeModelRow,
) -> dict[str, object]:
    return {
        "tree_source": row.tree_source,
        "traits_source": row.traits_source,
        "response": row.response,
        "formula": row.formula,
        "model_family": row.model_family,
        "selected_criterion": row.selected_criterion,
        "best_formula": row.best_formula,
        "rank": row.rank,
        "selected": row.selected,
        "analysis_taxon_count": row.analysis_taxon_count,
        "excluded_taxon_count": row.excluded_taxon_count,
        "excluded_taxa": _stringify_list(row.excluded_taxa),
        "encoded_columns": _stringify_list(row.encoded_columns),
        "coefficient_name": row.coefficient_name,
        "estimate": row.estimate,
        "standard_error": row.standard_error,
        "test_statistic": row.test_statistic,
        "p_value": row.p_value,
        "lower_95_confidence_interval": row.lower_95_confidence_interval,
        "upper_95_confidence_interval": row.upper_95_confidence_interval,
        "inference_distribution": row.inference_distribution,
        "phylogenetic_parameter_name": row.phylogenetic_parameter_name,
        "phylogenetic_parameter_value": row.phylogenetic_parameter_value,
        "phylogenetic_parameter_estimated": row.phylogenetic_parameter_estimated,
        "log_likelihood": row.log_likelihood,
        "aic": row.aic,
        "aicc": row.aicc,
        "bic": row.bic,
        "delta_aicc": row.delta_aicc,
        "delta_bic": row.delta_bic,
        "akaike_weight": row.akaike_weight,
        "residual_mean": "" if row.residual_mean is None else row.residual_mean,
        "outlier_taxon_count": row.outlier_taxon_count,
        "outlier_taxa": _stringify_list(row.outlier_taxa),
        "max_leverage": "" if row.max_leverage is None else row.max_leverage,
        "max_abs_standardized_residual": ""
        if row.max_abs_standardized_residual is None
        else row.max_abs_standardized_residual,
        "converged": "" if row.converged is None else row.converged,
        "iteration_count": "" if row.iteration_count is None else row.iteration_count,
        "separation_detected": row.separation_detected,
        "diagnostic_warning_count": row.diagnostic_warning_count,
        "diagnostic_warnings": _stringify_list(row.diagnostic_warnings),
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
    from .taxon import _build_supplementary_row as build_row_impl

    return build_row_impl(
        taxon=taxon,
        analysis_taxa=analysis_taxa,
        crosswalk=crosswalk,
        completeness=completeness,
        exclusion=exclusion,
        workflow_loss=workflow_loss,
        metadata_values=metadata_values,
        trait_values=trait_values,
    )


def _sequence_uncertainty_lookup(
    summary: AlignmentSummary,
) -> dict[str, SequenceUncertaintyProfile]:
    from .alignment import _sequence_uncertainty_lookup as lookup_impl

    return lookup_impl(summary)


def _sequence_quality_lookup(
    ranking: SequenceQualityRankingReport,
) -> dict[str, SequenceQualityRankingRow]:
    from .alignment import _sequence_quality_lookup as lookup_impl

    return lookup_impl(ranking)


def _filtering_status(
    *,
    sequence_id: str,
    original_ids: set[str],
    filtered_ids: set[str] | None,
) -> tuple[str, str | None]:
    from .alignment import _filtering_status as filtering_status_impl

    return filtering_status_impl(
        sequence_id=sequence_id,
        original_ids=original_ids,
        filtered_ids=filtered_ids,
    )


def _alignment_sequence_order(
    original_summary: AlignmentSummary,
    filtered_summary: AlignmentSummary | None,
) -> list[str]:
    from .alignment import _alignment_sequence_order as sequence_order_impl

    return sequence_order_impl(original_summary, filtered_summary)


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
    from .alignment import _build_alignment_row as build_row_impl

    return build_row_impl(
        sequence_id=sequence_id,
        original_summary=original_summary,
        original_quality=original_quality,
        original_low_information=original_low_information,
        original_quality_lookup=original_quality_lookup,
        original_uncertainty_lookup=original_uncertainty_lookup,
        filtered_summary=filtered_summary,
        filtered_low_information=filtered_low_information,
    )


def _write_alignment_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryAlignmentDiagnosticsRow],
) -> Path:
    from .alignment import _write_alignment_rows as write_rows_impl

    return write_rows_impl(path, columns=columns, rows=rows)


def _support_counts(rows: list[TreeSupportRow]) -> dict[str, int]:
    from .tree_diagnostics import _support_counts as support_counts_impl

    return support_counts_impl(rows)


def _topology_shape(inspection: TreeInspectionReport) -> str:
    from .tree_diagnostics import _topology_shape as topology_shape_impl

    return topology_shape_impl(inspection)


def _tree_warning_ledger(
    *,
    validation: TreeValidationReport,
    inspection: TreeInspectionReport,
    forensic: TreeForensicReport,
) -> list[str]:
    from .tree_diagnostics import _tree_warning_ledger as warning_ledger_impl

    return warning_ledger_impl(
        validation=validation,
        inspection=inspection,
        forensic=forensic,
    )


def _build_tree_forensic_review(
    *,
    tree_path: Path,
    validation: TreeValidationReport,
    inspection: TreeInspectionReport,
) -> TreeForensicReport:
    from .tree_diagnostics import _build_tree_forensic_review as build_review_impl

    return build_review_impl(
        tree_path=tree_path,
        validation=validation,
        inspection=inspection,
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
    from .tree_diagnostics import _build_tree_row as build_row_impl

    return build_row_impl(
        tree_path=tree_path,
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        support_rows=support_rows,
        branch_stats=branch_stats,
    )


def _write_tree_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryTreeDiagnosticsRow],
) -> Path:
    from .tree_diagnostics import _write_tree_rows as write_rows_impl

    return write_rows_impl(path, columns=columns, rows=rows)


def _clade_support_row_lookup(
    rows: list[TreeSetCladeSupportRow],
) -> dict[tuple[str, ...], TreeSetCladeSupportRow]:
    from .clade_support import _clade_support_row_lookup as lookup_impl

    return lookup_impl(rows)


def _build_clade_support_row(
    *,
    tree_path: Path,
    comparison_tree_set_path: Path | None,
    support_row: TreeSupportRow,
    frequency_row: TreeSetCladeSupportRow | None,
) -> SupplementaryCladeSupportRow:
    from .clade_support import _build_clade_support_row as build_row_impl

    return build_row_impl(
        tree_path=tree_path,
        comparison_tree_set_path=comparison_tree_set_path,
        support_row=support_row,
        frequency_row=frequency_row,
    )


def _write_clade_support_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryCladeSupportRow],
) -> Path:
    from .clade_support import _write_clade_support_rows as write_rows_impl

    return write_rows_impl(path, columns=columns, rows=rows)


def _write_comparative_model_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryComparativeModelRow],
) -> Path:
    return _write_dict_rows(
        path,
        columns=columns,
        rows=[_serialize_comparative_model_row(row) for row in rows],
    )


def _serialize_ancestral_state_row(
    row: SupplementaryAncestralStateRow,
) -> dict[str, str]:
    return {
        "tree_source": row.tree_source,
        "traits_source": row.traits_source,
        "trait": row.trait,
        "reconstruction_kind": row.reconstruction_kind,
        "model": row.model,
        "estimator": "" if row.estimator is None else row.estimator,
        "state_ordering": "" if row.state_ordering is None else row.state_ordering,
        "root_prior_mode": (
            "" if row.root_prior_mode is None else row.root_prior_mode
        ),
        "fixed_root_state": (
            "" if row.fixed_root_state is None else row.fixed_root_state
        ),
        "alpha": "" if row.alpha is None else str(row.alpha),
        "analysis_taxon_count": str(row.analysis_taxon_count),
        "excluded_taxon_count": str(row.excluded_taxon_count),
        "excluded_taxa": _stringify_list(row.excluded_taxa),
        "warning_count": str(row.warning_count),
        "warnings": _stringify_list(row.warnings),
        "node": row.node,
        "node_name": "" if row.node_name is None else row.node_name,
        "descendant_taxa": _stringify_list(row.descendant_taxa),
        "descendant_taxon_count": str(row.descendant_taxon_count),
        "estimate_value": "" if row.estimate_value is None else str(row.estimate_value),
        "most_likely_state": (
            "" if row.most_likely_state is None else row.most_likely_state
        ),
        "state_set": _stringify_list(row.state_set),
        "state_probabilities": json.dumps(row.state_probabilities, sort_keys=True),
        "standard_error": "" if row.standard_error is None else str(row.standard_error),
        "lower_95_interval": (
            "" if row.lower_95_interval is None else str(row.lower_95_interval)
        ),
        "upper_95_interval": (
            "" if row.upper_95_interval is None else str(row.upper_95_interval)
        ),
        "confidence": str(row.confidence),
        "ambiguous": "" if row.ambiguous is None else str(row.ambiguous).lower(),
        "unstable": str(row.unstable).lower(),
        "interpretation": row.interpretation,
        "downstream_risks": _stringify_list(row.downstream_risks),
    }


def _write_ancestral_state_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryAncestralStateRow],
) -> Path:
    return _write_dict_rows(
        path,
        columns=columns,
        rows=[_serialize_ancestral_state_row(row) for row in rows],
    )


def _split_batch_values(value: str) -> list[str]:
    stripped = value.strip()
    if not stripped:
        return []
    for separator in ("|", "; ", ";"):
        if separator in stripped:
            return [item for item in stripped.split(separator) if item]
    return [stripped]


def _read_bundle_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=_table_delimiter(path)))


def _read_bundle_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _optional_int(value: str | None) -> int | None:
    if value is None:
        return None
    stripped = value.strip()
    if not stripped:
        return None
    return int(stripped)


def _serialize_batch_summary_row(
    row: SupplementaryBatchSummaryRow,
) -> dict[str, str]:
    return {
        "row_scope": row.row_scope,
        "dataset_id": row.dataset_id,
        "dataset_label": row.dataset_label,
        "workflow_status": row.workflow_status,
        "variant_id": "" if row.variant_id is None else row.variant_id,
        "label": "" if row.label is None else row.label,
        "execution_mode": "" if row.execution_mode is None else row.execution_mode,
        "task_status": "" if row.task_status is None else row.task_status,
        "job_status": "" if row.job_status is None else row.job_status,
        "output_freshness_status": (
            "" if row.output_freshness_status is None else row.output_freshness_status
        ),
        "recovery_action": "" if row.recovery_action is None else row.recovery_action,
        "merge_status": "" if row.merge_status is None else row.merge_status,
        "evidence_status": "" if row.evidence_status is None else row.evidence_status,
        "reproducibility_status": (
            "" if row.reproducibility_status is None else row.reproducibility_status
        ),
        "selected_model": "" if row.selected_model is None else row.selected_model,
        "output_root": "" if row.output_root is None else row.output_root,
        "task_log_path": "" if row.task_log_path is None else row.task_log_path,
        "evidence_json_path": (
            "" if row.evidence_json_path is None else row.evidence_json_path
        ),
        "evidence_html_path": (
            "" if row.evidence_html_path is None else row.evidence_html_path
        ),
        "variant_count": "" if row.variant_count is None else str(row.variant_count),
        "successful_variant_count": (
            ""
            if row.successful_variant_count is None
            else str(row.successful_variant_count)
        ),
        "failed_variant_count": (
            "" if row.failed_variant_count is None else str(row.failed_variant_count)
        ),
        "output_file_count": (
            "" if row.output_file_count is None else str(row.output_file_count)
        ),
        "output_byte_count": (
            "" if row.output_byte_count is None else str(row.output_byte_count)
        ),
        "artifact_file_count": (
            "" if row.artifact_file_count is None else str(row.artifact_file_count)
        ),
        "linked_artifact_count": (
            ""
            if row.linked_artifact_count is None
            else str(row.linked_artifact_count)
        ),
        "linked_artifact_bytes": (
            ""
            if row.linked_artifact_bytes is None
            else str(row.linked_artifact_bytes)
        ),
        "issue_count": str(row.issue_count),
        "issues": _stringify_list(row.issues),
        "error_code": "" if row.error_code is None else row.error_code,
        "error_message": "" if row.error_message is None else row.error_message,
        "job_evidence_warning_count": (
            ""
            if row.job_evidence_warning_count is None
            else str(row.job_evidence_warning_count)
        ),
        "warning_count": str(row.warning_count),
        "warnings": _stringify_list(row.warnings),
    }


def _write_batch_summary_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryBatchSummaryRow],
) -> Path:
    return _write_dict_rows(
        path,
        columns=columns,
        rows=[_serialize_batch_summary_row(row) for row in rows],
    )


def _build_taxon_rows(
    *,
    dataset_audit: DatasetAuditReport,
    workflow_loss: TaxonWorkflowLossReport,
    metadata_table: TaxonTable,
    traits_table: TaxonTable,
) -> list[SupplementaryTaxonTableRow]:
    from .taxon import _build_taxon_rows as build_rows_impl

    return build_rows_impl(
        dataset_audit=dataset_audit,
        workflow_loss=workflow_loss,
        metadata_table=metadata_table,
        traits_table=traits_table,
    )


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
    from .taxon import write_supplementary_taxon_table as write_taxon_table_impl

    return write_taxon_table_impl(
        path,
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        alignment_path=alignment_path,
        filtered_alignment_path=filtered_alignment_path,
        inference_tree_path=inference_tree_path,
        reported_taxa_path=reported_taxa_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
    )


def write_supplementary_alignment_diagnostics_table(
    path: Path,
    *,
    alignment_path: Path,
    filtered_alignment_path: Path | None = None,
) -> SupplementaryAlignmentDiagnosticsTableResult:
    from .alignment import (
        write_supplementary_alignment_diagnostics_table as write_alignment_table_impl,
    )

    return write_alignment_table_impl(
        path,
        alignment_path=alignment_path,
        filtered_alignment_path=filtered_alignment_path,
    )


def write_supplementary_tree_diagnostics_table(
    path: Path,
    *,
    tree_path: Path,
) -> SupplementaryTreeDiagnosticsTableResult:
    from .tree_diagnostics import (
        write_supplementary_tree_diagnostics_table as write_tree_table_impl,
    )

    return write_tree_table_impl(path, tree_path=tree_path)


def write_supplementary_clade_support_table(
    path: Path,
    *,
    tree_path: Path,
    comparison_tree_set_path: Path | None = None,
) -> SupplementaryCladeSupportTableResult:
    from .clade_support import (
        write_supplementary_clade_support_table as write_clade_support_table_impl,
    )

    return write_clade_support_table_impl(
        path,
        tree_path=tree_path,
        comparison_tree_set_path=comparison_tree_set_path,
    )


def write_supplementary_model_selection_table(
    path: Path,
    *,
    iqtree_report_path: Path,
    model_sidecar_path: Path | None = None,
) -> SupplementaryModelSelectionTableResult:
    from .model_selection import (
        write_supplementary_model_selection_table as write_model_selection_table_impl,
    )

    return write_model_selection_table_impl(
        path,
        iqtree_report_path=iqtree_report_path,
        model_sidecar_path=model_sidecar_path,
    )


def write_supplementary_diversification_table(
    path: Path,
    *,
    tree_path: Path,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
    sampling_column: str | None = None,
    clade_model: str = "birth-death",
) -> SupplementaryDiversificationTableResult:
    from .diversification import (
        write_supplementary_diversification_table as write_diversification_table_impl,
    )

    return write_diversification_table_impl(
        path,
        tree_path=tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
        clade_model=clade_model,
    )


def _bundle_output_totals(bundle_root: Path) -> tuple[int, int]:
    file_paths = [path for path in bundle_root.rglob("*") if path.is_file()]
    return len(file_paths), sum(path.stat().st_size for path in file_paths)


def _maybe_path(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    return None if not stripped else stripped


def _build_batch_variant_row(
    *,
    dataset_id: str,
    dataset_label: str,
    workflow_status: str,
    task_row: dict[str, str],
    task_record: dict[str, object],
    variant_summary_row: dict[str, str],
    job_status_row: dict[str, str],
    freshness_row: dict[str, str],
    recovery_row: dict[str, str],
    merge_row: dict[str, str],
    evidence_row: dict[str, str],
    reproducibility_row: dict[str, str],
) -> SupplementaryBatchSummaryRow:
    issues = (
        _split_batch_values(merge_row.get("issues", ""))
        + _split_batch_values(reproducibility_row.get("issues", ""))
        + _split_batch_values(freshness_row.get("stale_reason_codes", ""))
        + _split_batch_values(job_status_row.get("output_freshness_reason_codes", ""))
    )
    warnings: list[str] = []
    task_status = _maybe_path(task_row.get("status"))
    job_status = _maybe_path(job_status_row.get("status"))
    freshness_status = _maybe_path(freshness_row.get("freshness_status"))
    recovery_action = _maybe_path(recovery_row.get("recovery_action"))
    merge_status = _maybe_path(merge_row.get("merge_status"))
    evidence_status = _maybe_path(merge_row.get("evidence_status"))
    reproducibility_status = _maybe_path(reproducibility_row.get("status"))
    error_code = _maybe_path(task_row.get("error_code")) or _maybe_path(
        task_record.get("error_code") if isinstance(task_record, dict) else None
    )
    error_message = _maybe_path(
        task_record.get("error_message") if isinstance(task_record, dict) else None
    )
    if task_status is not None and task_status != "succeeded":
        warnings.append(f"task-status:{task_status}")
    if error_code is not None:
        warnings.append(f"task-error:{error_code}")
    if job_status is not None and job_status != "completed":
        warnings.append(f"job-status:{job_status}")
    if freshness_status is not None and freshness_status != "fresh":
        warnings.append(f"freshness-status:{freshness_status}")
    if recovery_action is not None and recovery_action != "no_action":
        warnings.append(f"recovery-action:{recovery_action}")
    if merge_status is not None and merge_status != "merged":
        warnings.append(f"merge-status:{merge_status}")
    if evidence_status is not None and evidence_status != "present":
        warnings.append(f"evidence-status:{evidence_status}")
    if reproducibility_status is not None and reproducibility_status != "passed":
        warnings.append(f"reproducibility-status:{reproducibility_status}")
    job_evidence_warning_count = _optional_int(evidence_row.get("warning_count"))
    if job_evidence_warning_count is not None and job_evidence_warning_count > 0:
        warnings.append(f"job-evidence-warning-count:{job_evidence_warning_count}")
    warnings.extend(issue for issue in issues if issue not in warnings)
    return SupplementaryBatchSummaryRow(
        row_scope="variant",
        dataset_id=dataset_id,
        dataset_label=dataset_label,
        workflow_status=workflow_status,
        variant_id=str(task_row["variant_id"]),
        label=_maybe_path(task_row.get("label")),
        execution_mode=_maybe_path(task_row.get("execution_mode")),
        task_status=task_status,
        job_status=job_status,
        output_freshness_status=freshness_status,
        recovery_action=recovery_action,
        merge_status=merge_status,
        evidence_status=evidence_status,
        reproducibility_status=reproducibility_status,
        selected_model=_maybe_path(variant_summary_row.get("selected_model")),
        output_root=_maybe_path(job_status_row.get("output_root"))
        or _maybe_path(task_record.get("output_root") if isinstance(task_record, dict) else None),
        task_log_path=_maybe_path(task_row.get("log_path"))
        or _maybe_path(job_status_row.get("task_log_path")),
        evidence_json_path=_maybe_path(merge_row.get("evidence_json_path"))
        or _maybe_path(evidence_row.get("evidence_json_path")),
        evidence_html_path=_maybe_path(merge_row.get("evidence_html_path"))
        or _maybe_path(evidence_row.get("evidence_html_path")),
        variant_count=None,
        successful_variant_count=None,
        failed_variant_count=None,
        output_file_count=_optional_int(job_status_row.get("output_file_count")),
        output_byte_count=_optional_int(job_status_row.get("output_byte_count")),
        artifact_file_count=_optional_int(evidence_row.get("artifact_file_count")),
        linked_artifact_count=None,
        linked_artifact_bytes=None,
        issue_count=len(issues),
        issues=issues,
        error_code=error_code,
        error_message=error_message,
        job_evidence_warning_count=job_evidence_warning_count,
        warning_count=len(warnings),
        warnings=warnings,
    )


def _build_batch_dataset_row(
    *,
    bundle_root: Path,
    dataset_id: str,
    dataset_label: str,
    workflow_status: str,
    execution_mode: str,
    workflow_run: dict[str, object],
    workflow_status_summary: dict[str, object],
    failure_recovery_summary: dict[str, object],
    merge_summary: dict[str, object],
    reproducibility_summary: dict[str, object],
    report_manifest: dict[str, object],
) -> SupplementaryBatchSummaryRow:
    total_file_count, total_byte_count = _bundle_output_totals(bundle_root)
    warnings: list[str] = []
    issues: list[str] = []
    active_run_state = str(workflow_status_summary.get("active_run_state", "unknown"))
    overall_recovery_status = str(
        failure_recovery_summary.get("overall_recovery_status", "unknown")
    )
    merge_status = str(merge_summary.get("merge_status", "unknown"))
    reproducibility_status = (
        "passed"
        if bool(reproducibility_summary.get("all_passed", False))
        else "failed"
    )
    failed_variant_count = len(list(workflow_run.get("failed_variants", [])))
    successful_variant_count = len(list(workflow_run.get("successful_variants", [])))
    if workflow_status != "succeeded":
        warnings.append(f"workflow-status:{workflow_status}")
    if active_run_state != "absent":
        warnings.append(f"active-run-state:{active_run_state}")
    if overall_recovery_status != "clean":
        warnings.append(f"recovery-status:{overall_recovery_status}")
    if merge_status != "merge-ready":
        warnings.append(f"merge-status:{merge_status}")
    if reproducibility_status != "passed":
        warnings.append("reproducibility-status:failed")
    failed_job_count = int(workflow_status_summary.get("failed_job_count", 0))
    stale_job_count = int(workflow_status_summary.get("stale_job_count", 0))
    stale_output_job_count = int(workflow_status_summary.get("stale_output_job_count", 0))
    if failed_job_count > 0:
        issues.append(f"failed-job-count:{failed_job_count}")
    if stale_job_count > 0:
        issues.append(f"stale-job-count:{stale_job_count}")
    if stale_output_job_count > 0:
        issues.append(f"stale-output-job-count:{stale_output_job_count}")
    if int(reproducibility_summary.get("failed_variant_count", 0)) > 0:
        issues.append(
            "reproducibility-failed-variant-count:"
            + str(reproducibility_summary["failed_variant_count"])
        )
    warnings.extend(issue for issue in issues if issue not in warnings)
    linked_artifacts = report_manifest.get("linked_artifacts", {})
    linked_artifact_bytes = 0
    if isinstance(linked_artifacts, dict):
        linked_artifact_bytes = sum(
            int(payload.get("byte_count", 0))
            for payload in linked_artifacts.values()
            if isinstance(payload, dict)
        )
    return SupplementaryBatchSummaryRow(
        row_scope="dataset",
        dataset_id=dataset_id,
        dataset_label=dataset_label,
        workflow_status=workflow_status,
        variant_id=None,
        label=dataset_label,
        execution_mode=execution_mode,
        task_status=None,
        job_status=None,
        output_freshness_status=None,
        recovery_action=overall_recovery_status,
        merge_status=merge_status,
        evidence_status=None,
        reproducibility_status=reproducibility_status,
        selected_model=None,
        output_root=".",
        task_log_path=None,
        evidence_json_path=None,
        evidence_html_path=None,
        variant_count=int(workflow_run.get("variant_count", 0)),
        successful_variant_count=successful_variant_count,
        failed_variant_count=failed_variant_count,
        output_file_count=total_file_count,
        output_byte_count=total_byte_count,
        artifact_file_count=None,
        linked_artifact_count=int(report_manifest.get("linked_artifact_count", 0)),
        linked_artifact_bytes=linked_artifact_bytes,
        issue_count=len(issues),
        issues=issues,
        error_code=None,
        error_message=None,
        job_evidence_warning_count=None,
        warning_count=len(warnings),
        warnings=warnings,
    )


def write_supplementary_batch_summary_table(
    path: Path,
    *,
    workflow_bundle_root: Path,
) -> SupplementaryBatchSummaryTableResult:
    """Write one supplementary batch summary table from a written workflow bundle."""
    bundle_root = workflow_bundle_root.resolve()
    workflow_run = _read_bundle_json(
        bundle_root / "rabies-method-sensitivity-panel.run.json"
    )
    workflow_manifest = _read_bundle_json(
        bundle_root / "rabies-method-sensitivity.manifest.json"
    )
    report_manifest = _read_bundle_json(
        bundle_root
        / "report-artifacts"
        / "rabies-method-sensitivity-report.manifest.json"
    )
    workflow_summary_rows = _read_bundle_rows(bundle_root / "workflow-summary.tsv")
    task_rows = _read_bundle_rows(bundle_root / "parallel-execution-summary.tsv")
    variant_summary_rows = _read_bundle_rows(bundle_root / "variant-summary.tsv")
    job_status_rows = _read_bundle_rows(bundle_root / "slurm-job-status.tsv")
    freshness_rows = _read_bundle_rows(bundle_root / "slurm-output-freshness.tsv")
    recovery_rows = _read_bundle_rows(bundle_root / "slurm-failure-recovery-jobs.tsv")
    merge_rows = _read_bundle_rows(bundle_root / "slurm-merge-variants.tsv")
    evidence_rows = _read_bundle_rows(bundle_root / "slurm-job-evidence.tsv")
    reproducibility_rows = _read_bundle_rows(bundle_root / "reproducibility-variants.tsv")
    workflow_status_summary = _read_bundle_json(bundle_root / "slurm-workflow-status.json")
    failure_recovery_summary = _read_bundle_json(
        bundle_root / "slurm-failure-recovery-report.json"
    )
    merge_summary = _read_bundle_json(bundle_root / "slurm-merge-report.json")
    reproducibility_summary = _read_bundle_json(bundle_root / "reproducibility-audit.json")

    workflow_summary = workflow_summary_rows[0]
    dataset_id = str(workflow_summary["dataset_id"])
    dataset_label = str(workflow_manifest.get("label", dataset_id))
    task_records = workflow_run.get("task_records", [])
    task_record_by_variant = {
        str(record["variant_id"]): record
        for record in task_records
        if isinstance(record, dict)
    }
    task_row_by_variant = {str(row["variant_id"]): row for row in task_rows}
    variant_summary_by_variant = {
        str(row["variant_id"]): row for row in variant_summary_rows
    }
    job_status_by_variant = {str(row["variant_id"]): row for row in job_status_rows}
    freshness_by_variant = {str(row["variant_id"]): row for row in freshness_rows}
    recovery_by_variant = {str(row["variant_id"]): row for row in recovery_rows}
    merge_by_variant = {str(row["variant_id"]): row for row in merge_rows}
    evidence_by_variant = {str(row["variant_id"]): row for row in evidence_rows}
    reproducibility_by_variant = {
        str(row["variant_id"]): row for row in reproducibility_rows
    }
    workflow_status = str(workflow_run.get("status", "unknown"))
    execution_mode = str(workflow_run.get("execution_mode", "unknown"))
    variant_ids = [
        str(record["variant_id"])
        for record in task_records
        if isinstance(record, dict) and "variant_id" in record
    ]
    rows = [
        _build_batch_dataset_row(
            bundle_root=bundle_root,
            dataset_id=dataset_id,
            dataset_label=dataset_label,
            workflow_status=workflow_status,
            execution_mode=execution_mode,
            workflow_run=workflow_run,
            workflow_status_summary=workflow_status_summary,
            failure_recovery_summary=failure_recovery_summary,
            merge_summary=merge_summary,
            reproducibility_summary=reproducibility_summary,
            report_manifest=report_manifest,
        )
    ]
    rows.extend(
        _build_batch_variant_row(
            dataset_id=dataset_id,
            dataset_label=dataset_label,
            workflow_status=workflow_status,
            task_row=task_row_by_variant[variant_id],
            task_record=task_record_by_variant[variant_id],
            variant_summary_row=variant_summary_by_variant[variant_id],
            job_status_row=job_status_by_variant[variant_id],
            freshness_row=freshness_by_variant[variant_id],
            recovery_row=recovery_by_variant[variant_id],
            merge_row=merge_by_variant[variant_id],
            evidence_row=evidence_by_variant[variant_id],
            reproducibility_row=reproducibility_by_variant[variant_id],
        )
        for variant_id in variant_ids
    )
    columns = _batch_summary_table_columns()
    _write_batch_summary_rows(path, columns=columns, rows=rows)
    return SupplementaryBatchSummaryTableResult(
        output_path=path,
        row_count=len(rows),
        dataset_row_count=sum(1 for row in rows if row.row_scope == "dataset"),
        variant_row_count=sum(1 for row in rows if row.row_scope == "variant"),
        workflow_status=workflow_status,
        warning_count=sum(row.warning_count for row in rows),
        columns=columns,
        rows=rows,
    )


def _serialize_comparative_exclusion(
    row: ComparativeRegressionModelExclusion,
) -> str:
    if row.missing_columns:
        return f"{row.taxon}:{row.reason}:{','.join(row.missing_columns)}"
    return f"{row.taxon}:{row.reason}"


def _prepare_shared_comparative_inputs(
    *,
    tree_path: Path,
    traits_path: Path,
    taxon_column: str | None,
    analysis_taxa: list[str],
) -> tuple[Path, Path]:
    reduced_tree, _ = prune_tree_to_requested_taxa(tree_path, analysis_taxa)
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    reduced_rows = [rows_by_taxon[taxon] for taxon in analysis_taxa]
    temp_dir = Path(tempfile.mkdtemp(prefix="bijux-comparative-supplement-"))
    reduced_tree_path = temp_dir / "comparative-model-tree.nwk"
    reduced_table_path = temp_dir / "comparative-model-traits.tsv"
    reduced_tree_path.write_text(dumps_newick(reduced_tree) + "\n", encoding="utf-8")
    write_taxon_rows(reduced_table_path, columns=table.columns, rows=reduced_rows)
    return reduced_tree_path, reduced_table_path


def _build_pgls_comparative_rows(
    *,
    report: ComparativeRegressionModelSelectionReport,
    model_row: ComparativeRegressionModelRow,
    reduced_tree_path: Path,
    reduced_table_path: Path,
    taxon_column: str,
    lambda_value: float | str,
) -> list[SupplementaryComparativeModelRow]:
    fitted = run_pgls(
        reduced_tree_path,
        reduced_table_path,
        formula=model_row.formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    outlier_taxa = [row.taxon for row in fitted.diagnostics.outlier_taxa]
    max_leverage = max(
        (row.leverage for row in fitted.diagnostics.leverage_rows),
        default=None,
    )
    max_abs_standardized_residual = max(
        (abs(row.standardized_residual) for row in fitted.diagnostics.leverage_rows),
        default=None,
    )
    excluded_taxa = [
        _serialize_comparative_exclusion(row) for row in report.excluded_taxa
    ]
    return [
        SupplementaryComparativeModelRow(
            tree_source=str(report.tree_path),
            traits_source=str(report.traits_path),
            response=report.response,
            formula=model_row.formula,
            model_family=model_row.model_family,
            selected_criterion=report.selected_criterion,
            best_formula=report.best_formula,
            rank=model_row.rank,
            selected=model_row.selected,
            analysis_taxon_count=len(report.analysis_taxa),
            excluded_taxon_count=len(report.excluded_taxa),
            excluded_taxa=excluded_taxa,
            encoded_columns=list(model_row.encoded_columns),
            coefficient_name=coefficient.name,
            estimate=coefficient.estimate,
            standard_error=coefficient.standard_error,
            test_statistic=coefficient.test_statistic,
            p_value=coefficient.p_value,
            lower_95_confidence_interval=coefficient.lower_95_confidence_interval,
            upper_95_confidence_interval=coefficient.upper_95_confidence_interval,
            inference_distribution=coefficient.inference_distribution,
            phylogenetic_parameter_name=model_row.phylogenetic_parameter_name,
            phylogenetic_parameter_value=fitted.lambda_fit.lambda_value,
            phylogenetic_parameter_estimated=fitted.lambda_fit.mode == "estimated",
            log_likelihood=fitted.log_likelihood,
            aic=model_row.aic,
            aicc=model_row.aicc,
            bic=model_row.bic,
            delta_aicc=model_row.delta_aicc,
            delta_bic=model_row.delta_bic,
            akaike_weight=model_row.akaike_weight,
            residual_mean=fitted.diagnostics.residual_mean,
            outlier_taxon_count=len(outlier_taxa),
            outlier_taxa=outlier_taxa,
            max_leverage=max_leverage,
            max_abs_standardized_residual=max_abs_standardized_residual,
            converged=None,
            iteration_count=None,
            separation_detected=False,
            diagnostic_warning_count=0,
            diagnostic_warnings=[],
        )
        for coefficient in fitted.coefficients
    ]


def _build_logistic_comparative_rows(
    *,
    report: ComparativeRegressionModelSelectionReport,
    model_row: ComparativeRegressionModelRow,
    reduced_tree_path: Path,
    reduced_table_path: Path,
    taxon_column: str,
    lambda_value: float,
) -> list[SupplementaryComparativeModelRow]:
    fitted = summarize_phylogenetic_logistic(
        reduced_tree_path,
        reduced_table_path,
        formula=model_row.formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    residual_mean = (
        None
        if not fitted.fitted_rows
        else sum(row.residual for row in fitted.fitted_rows) / len(fitted.fitted_rows)
    )
    excluded_taxa = [
        _serialize_comparative_exclusion(row) for row in report.excluded_taxa
    ]
    warning_messages = [warning.message for warning in fitted.warnings]
    return [
        SupplementaryComparativeModelRow(
            tree_source=str(report.tree_path),
            traits_source=str(report.traits_path),
            response=report.response,
            formula=model_row.formula,
            model_family=model_row.model_family,
            selected_criterion=report.selected_criterion,
            best_formula=report.best_formula,
            rank=model_row.rank,
            selected=model_row.selected,
            analysis_taxon_count=len(report.analysis_taxa),
            excluded_taxon_count=len(report.excluded_taxa),
            excluded_taxa=excluded_taxa,
            encoded_columns=list(model_row.encoded_columns),
            coefficient_name=coefficient.name,
            estimate=coefficient.estimate,
            standard_error=coefficient.standard_error,
            test_statistic=coefficient.test_statistic,
            p_value=coefficient.p_value,
            lower_95_confidence_interval=coefficient.lower_95_confidence_interval,
            upper_95_confidence_interval=coefficient.upper_95_confidence_interval,
            inference_distribution=coefficient.inference_distribution,
            phylogenetic_parameter_name=model_row.phylogenetic_parameter_name,
            phylogenetic_parameter_value=fitted.lambda_value,
            phylogenetic_parameter_estimated=False,
            log_likelihood=fitted.binomial_log_likelihood,
            aic=model_row.aic,
            aicc=model_row.aicc,
            bic=model_row.bic,
            delta_aicc=model_row.delta_aicc,
            delta_bic=model_row.delta_bic,
            akaike_weight=model_row.akaike_weight,
            residual_mean=residual_mean,
            outlier_taxon_count=0,
            outlier_taxa=[],
            max_leverage=None,
            max_abs_standardized_residual=None,
            converged=fitted.converged,
            iteration_count=fitted.iteration_count,
            separation_detected=fitted.separation_detected,
            diagnostic_warning_count=len(fitted.warnings),
            diagnostic_warnings=warning_messages,
        )
        for coefficient in fitted.coefficients
    ]


def write_supplementary_comparative_model_table(
    path: Path,
    *,
    tree_path: Path,
    traits_path: Path,
    formulas: list[str],
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> SupplementaryComparativeModelTableResult:
    """Write one coefficient-level supplementary comparative-model table."""
    report = compare_comparative_regression_models(
        tree_path,
        traits_path,
        formulas=formulas,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    reduced_tree_path, reduced_table_path = _prepare_shared_comparative_inputs(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=taxon_column,
        analysis_taxa=report.analysis_taxa,
    )
    rows: list[SupplementaryComparativeModelRow] = []
    for model_row in report.rows:
        if model_row.model_family == "pgls":
            rows.extend(
                _build_pgls_comparative_rows(
                    report=report,
                    model_row=model_row,
                    reduced_tree_path=reduced_tree_path,
                    reduced_table_path=reduced_table_path,
                    taxon_column=report.taxon_column,
                    lambda_value=lambda_value,
                )
            )
            continue
        rows.extend(
            _build_logistic_comparative_rows(
                report=report,
                model_row=model_row,
                reduced_tree_path=reduced_tree_path,
                reduced_table_path=reduced_table_path,
                taxon_column=report.taxon_column,
                lambda_value=float(lambda_value),
            )
        )
    columns = _comparative_model_table_columns()
    _write_comparative_model_rows(path, columns=columns, rows=rows)
    return SupplementaryComparativeModelTableResult(
        output_path=path,
        row_count=len(rows),
        model_count=len(report.rows),
        selected_formula=report.best_formula,
        selected_criterion=report.selected_criterion,
        excluded_taxon_count=len(report.excluded_taxa),
        columns=columns,
        rows=rows,
    )


def _build_continuous_ancestral_state_rows(
    report: ContinuousAncestralReport,
) -> list[SupplementaryAncestralStateRow]:
    excluded_taxa = [
        f"{row.taxon}:{row.reason}" for row in continuous_ancestral_exclusions(report)
    ]
    return [
        SupplementaryAncestralStateRow(
            tree_source=str(report.tree_path),
            traits_source=str(report.traits_path),
            trait=report.trait,
            reconstruction_kind="continuous",
            model=report.model,
            estimator=report.estimator,
            state_ordering=None,
            root_prior_mode=None,
            fixed_root_state=None,
            alpha=report.alpha,
            analysis_taxon_count=report.taxon_count,
            excluded_taxon_count=len(excluded_taxa),
            excluded_taxa=excluded_taxa,
            warning_count=len(report.warnings),
            warnings=list(report.warnings),
            node=estimate.node,
            node_name=estimate.node_name,
            descendant_taxa=list(estimate.descendant_taxa),
            descendant_taxon_count=len(estimate.descendant_taxa),
            estimate_value=estimate.estimate,
            most_likely_state=None,
            state_set=[],
            state_probabilities={},
            standard_error=estimate.standard_error,
            lower_95_interval=estimate.lower_95_interval,
            upper_95_interval=estimate.upper_95_interval,
            confidence=estimate.confidence,
            ambiguous=None,
            unstable=estimate.unstable,
            interpretation=estimate.interpretation,
            downstream_risks=list(estimate.downstream_risks),
        )
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def _build_discrete_ancestral_state_rows(
    report: DiscreteAncestralReport,
) -> list[SupplementaryAncestralStateRow]:
    excluded_taxa = [
        f"{row.taxon}:{row.reason}" for row in discrete_ancestral_exclusions(report)
    ]
    return [
        SupplementaryAncestralStateRow(
            tree_source=str(report.tree_path),
            traits_source=str(report.traits_path),
            trait=report.trait,
            reconstruction_kind="discrete",
            model=report.model,
            estimator=None,
            state_ordering=report.state_ordering,
            root_prior_mode=report.root_prior_mode,
            fixed_root_state=report.fixed_root_state,
            alpha=None,
            analysis_taxon_count=report.taxon_count,
            excluded_taxon_count=len(excluded_taxa),
            excluded_taxa=excluded_taxa,
            warning_count=len(report.warnings),
            warnings=list(report.warnings),
            node=estimate.node,
            node_name=estimate.node_name,
            descendant_taxa=list(estimate.descendant_taxa),
            descendant_taxon_count=len(estimate.descendant_taxa),
            estimate_value=None,
            most_likely_state=estimate.most_likely_state,
            state_set=list(estimate.state_set),
            state_probabilities=dict(estimate.state_probabilities),
            standard_error=None,
            lower_95_interval=None,
            upper_95_interval=None,
            confidence=estimate.confidence,
            ambiguous=estimate.ambiguous,
            unstable=estimate.unstable,
            interpretation=estimate.interpretation,
            downstream_risks=list(estimate.downstream_risks),
        )
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def write_supplementary_ancestral_state_table(
    path: Path,
    *,
    tree_path: Path,
    traits_path: Path,
    trait: str,
    reconstruction_kind: str,
    taxon_column: str | None = None,
    model: str | None = None,
    estimator: str | None = None,
    alpha: float = 1.0,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    root_prior_mode: str = "equal",
    fixed_root_state: str | None = None,
) -> SupplementaryAncestralStateTableResult:
    """Write one supplementary internal-node ancestral-state table."""
    if reconstruction_kind == "continuous":
        report = reconstruct_continuous_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model="brownian" if model is None else model,
            estimator=estimator,
            alpha=alpha,
        )
        rows = _build_continuous_ancestral_state_rows(report)
        resolved_model = report.model
        unstable_node_count = len(report.unstable_nodes)
    elif reconstruction_kind == "discrete":
        report = reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model="fitch" if model is None else model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
            root_prior_mode=root_prior_mode,
            fixed_root_state=fixed_root_state,
        )
        rows = _build_discrete_ancestral_state_rows(report)
        resolved_model = report.model
        unstable_node_count = len(report.unstable_nodes)
    else:
        raise ValueError(
            "reconstruction_kind must be 'continuous' or 'discrete'"
        )
    columns = _ancestral_state_table_columns()
    _write_ancestral_state_rows(path, columns=columns, rows=rows)
    return SupplementaryAncestralStateTableResult(
        output_path=path,
        row_count=len(rows),
        reconstruction_kind=reconstruction_kind,
        model=resolved_model,
        analysis_taxon_count=0 if not rows else rows[0].analysis_taxon_count,
        excluded_taxon_count=0 if not rows else rows[0].excluded_taxon_count,
        unstable_node_count=unstable_node_count,
        columns=columns,
        rows=rows,
    )
