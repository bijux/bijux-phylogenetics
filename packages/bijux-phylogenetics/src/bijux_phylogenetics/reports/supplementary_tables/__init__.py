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


def _clade_support_table_columns() -> list[str]:
    return [
        "tree_source",
        "comparison_tree_set_source",
        "clade_id",
        "node_kind",
        "node_label",
        "descendant_taxa",
        "support",
        "support_fraction",
        "support_class",
        "support_method",
        "branch_length",
        "root_depth",
        "supporting_tree_count",
        "clade_frequency",
        "support_percent",
        "frequency_method",
        "frequency_status",
        "frequency_explanation",
    ]


def _model_selection_table_columns() -> list[str]:
    return [
        "iqtree_report_source",
        "model_sidecar_source",
        "rank",
        "model",
        "log_likelihood",
        "parameter_count",
        "aic",
        "aicc",
        "bic",
        "best_aic",
        "best_aicc",
        "best_bic",
        "selected_model",
        "selected_model_name",
        "selected_criterion",
    ]


def _comparative_model_table_columns() -> list[str]:
    return [
        "tree_source",
        "traits_source",
        "response",
        "formula",
        "model_family",
        "selected_criterion",
        "best_formula",
        "rank",
        "selected",
        "analysis_taxon_count",
        "excluded_taxon_count",
        "excluded_taxa",
        "encoded_columns",
        "coefficient_name",
        "estimate",
        "standard_error",
        "test_statistic",
        "p_value",
        "lower_95_confidence_interval",
        "upper_95_confidence_interval",
        "inference_distribution",
        "phylogenetic_parameter_name",
        "phylogenetic_parameter_value",
        "phylogenetic_parameter_estimated",
        "log_likelihood",
        "aic",
        "aicc",
        "bic",
        "delta_aicc",
        "delta_bic",
        "akaike_weight",
        "residual_mean",
        "outlier_taxon_count",
        "outlier_taxa",
        "max_leverage",
        "max_abs_standardized_residual",
        "converged",
        "iteration_count",
        "separation_detected",
        "diagnostic_warning_count",
        "diagnostic_warnings",
    ]


def _ancestral_state_table_columns() -> list[str]:
    return [
        "tree_source",
        "traits_source",
        "trait",
        "reconstruction_kind",
        "model",
        "estimator",
        "state_ordering",
        "root_prior_mode",
        "fixed_root_state",
        "alpha",
        "analysis_taxon_count",
        "excluded_taxon_count",
        "excluded_taxa",
        "warning_count",
        "warnings",
        "node",
        "node_name",
        "descendant_taxa",
        "descendant_taxon_count",
        "estimate_value",
        "most_likely_state",
        "state_set",
        "state_probabilities",
        "standard_error",
        "lower_95_interval",
        "upper_95_interval",
        "confidence",
        "ambiguous",
        "unstable",
        "interpretation",
        "downstream_risks",
    ]


def _diversification_table_columns() -> list[str]:
    return [
        "tree_source",
        "metadata_source",
        "clade_model",
        "better_model",
        "node",
        "node_name",
        "descendant_taxa",
        "tip_count",
        "crown_age",
        "clade_diversification_rate",
        "clade_rate_z_score",
        "clade_classification",
        "global_diversification_rate",
        "yule_log_likelihood",
        "yule_aic",
        "yule_corrected_tip_count",
        "yule_sampling_fraction",
        "yule_net_diversification_rate",
        "yule_relative_extinction",
        "birth_death_log_likelihood",
        "birth_death_aic",
        "birth_death_corrected_tip_count",
        "birth_death_sampling_fraction",
        "birth_death_net_diversification_rate",
        "birth_death_relative_extinction",
        "sampling_metadata_complete",
        "sampling_column",
        "sampling_fraction",
        "sampling_heterogeneous",
        "sampling_missing_taxa",
        "sampling_invalid_rows",
        "warning_count",
        "warnings",
    ]


def _batch_summary_table_columns() -> list[str]:
    return [
        "row_scope",
        "dataset_id",
        "dataset_label",
        "workflow_status",
        "variant_id",
        "label",
        "execution_mode",
        "task_status",
        "job_status",
        "output_freshness_status",
        "recovery_action",
        "merge_status",
        "evidence_status",
        "reproducibility_status",
        "selected_model",
        "output_root",
        "task_log_path",
        "evidence_json_path",
        "evidence_html_path",
        "variant_count",
        "successful_variant_count",
        "failed_variant_count",
        "output_file_count",
        "output_byte_count",
        "artifact_file_count",
        "linked_artifact_count",
        "linked_artifact_bytes",
        "issue_count",
        "issues",
        "error_code",
        "error_message",
        "job_evidence_warning_count",
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


def _serialize_model_selection_row(
    row: SupplementaryModelSelectionRow,
) -> dict[str, object]:
    return {
        "iqtree_report_source": row.iqtree_report_source,
        "model_sidecar_source": ""
        if row.model_sidecar_source is None
        else row.model_sidecar_source,
        "rank": row.rank,
        "model": row.model,
        "log_likelihood": row.log_likelihood,
        "parameter_count": "" if row.parameter_count is None else row.parameter_count,
        "aic": row.aic,
        "aicc": row.aicc,
        "bic": row.bic,
        "best_aic": row.best_aic,
        "best_aicc": row.best_aicc,
        "best_bic": row.best_bic,
        "selected_model": row.selected_model,
        "selected_model_name": row.selected_model_name,
        "selected_criterion": ""
        if row.selected_criterion is None
        else row.selected_criterion,
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


def _clade_support_row_lookup(
    rows: list[TreeSetCladeSupportRow],
) -> dict[tuple[str, ...], TreeSetCladeSupportRow]:
    return {tuple(row.descendant_taxa): row for row in rows}


def _build_clade_support_row(
    *,
    tree_path: Path,
    comparison_tree_set_path: Path | None,
    support_row: TreeSupportRow,
    frequency_row: TreeSetCladeSupportRow | None,
) -> SupplementaryCladeSupportRow:
    return SupplementaryCladeSupportRow(
        tree_source=str(tree_path),
        comparison_tree_set_source=(
            None
            if comparison_tree_set_path is None
            else str(comparison_tree_set_path)
        ),
        clade_id=support_row.node,
        node_kind=support_row.node_kind,
        node_label=support_row.node_label,
        descendant_taxa=list(support_row.descendant_taxa),
        support=support_row.support,
        support_fraction=support_row.support_fraction,
        support_class=support_row.support_class,
        support_method="tree-label",
        branch_length=support_row.branch_length,
        root_depth=support_row.root_depth,
        supporting_tree_count=(
            None if frequency_row is None else frequency_row.supporting_tree_count
        ),
        clade_frequency=None if frequency_row is None else frequency_row.clade_frequency,
        support_percent=None if frequency_row is None else frequency_row.support_percent,
        frequency_method=(
            None if frequency_row is None else "reference-tree-clade-frequency"
        ),
        frequency_status=None if frequency_row is None else frequency_row.support_status,
        frequency_explanation=(
            None if frequency_row is None else frequency_row.explanation
        ),
    )


def _write_clade_support_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryCladeSupportRow],
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
            writer.writerow(_serialize_clade_support_row(row))
    return path


def _write_model_selection_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryModelSelectionRow],
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
            writer.writerow(_serialize_model_selection_row(row))
    return path


def _write_comparative_model_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryComparativeModelRow],
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
            writer.writerow(_serialize_comparative_model_row(row))
    return path


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
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter=_table_delimiter(path),
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(_serialize_ancestral_state_row(row))
    return path


def _serialize_sampling_issue(row: SamplingFractionIssue) -> str:
    raw_value = row.raw_value if row.raw_value else "<missing>"
    return f"{row.taxon}:{row.code}:{raw_value}"


def _serialize_diversification_row(
    row: SupplementaryDiversificationRow,
) -> dict[str, str]:
    return {
        "tree_source": row.tree_source,
        "metadata_source": "" if row.metadata_source is None else row.metadata_source,
        "clade_model": row.clade_model,
        "better_model": row.better_model,
        "node": row.node,
        "node_name": "" if row.node_name is None else row.node_name,
        "descendant_taxa": _stringify_list(row.descendant_taxa),
        "tip_count": str(row.tip_count),
        "crown_age": str(row.crown_age),
        "clade_diversification_rate": str(row.clade_diversification_rate),
        "clade_rate_z_score": str(row.clade_rate_z_score),
        "clade_classification": row.clade_classification,
        "global_diversification_rate": str(row.global_diversification_rate),
        "yule_log_likelihood": str(row.yule_log_likelihood),
        "yule_aic": str(row.yule_aic),
        "yule_corrected_tip_count": str(row.yule_corrected_tip_count),
        "yule_sampling_fraction": str(row.yule_sampling_fraction),
        "yule_net_diversification_rate": str(row.yule_net_diversification_rate),
        "yule_relative_extinction": str(row.yule_relative_extinction),
        "birth_death_log_likelihood": str(row.birth_death_log_likelihood),
        "birth_death_aic": str(row.birth_death_aic),
        "birth_death_corrected_tip_count": str(row.birth_death_corrected_tip_count),
        "birth_death_sampling_fraction": str(row.birth_death_sampling_fraction),
        "birth_death_net_diversification_rate": str(
            row.birth_death_net_diversification_rate
        ),
        "birth_death_relative_extinction": str(row.birth_death_relative_extinction),
        "sampling_metadata_complete": (
            ""
            if row.sampling_metadata_complete is None
            else str(row.sampling_metadata_complete).lower()
        ),
        "sampling_column": "" if row.sampling_column is None else row.sampling_column,
        "sampling_fraction": (
            "" if row.sampling_fraction is None else str(row.sampling_fraction)
        ),
        "sampling_heterogeneous": (
            ""
            if row.sampling_heterogeneous is None
            else str(row.sampling_heterogeneous).lower()
        ),
        "sampling_missing_taxa": _stringify_list(row.sampling_missing_taxa),
        "sampling_invalid_rows": _stringify_list(row.sampling_invalid_rows),
        "warning_count": str(row.warning_count),
        "warnings": _stringify_list(row.warnings),
    }


def _write_diversification_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryDiversificationRow],
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
            writer.writerow(_serialize_diversification_row(row))
    return path


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
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter=_table_delimiter(path),
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(_serialize_batch_summary_row(row))
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


def write_supplementary_clade_support_table(
    path: Path,
    *,
    tree_path: Path,
    comparison_tree_set_path: Path | None = None,
) -> SupplementaryCladeSupportTableResult:
    """Write one supplementary clade-support table from a reference tree and optional tree set."""
    clades = extract_tree_clades(tree_path)
    support_rows = summarize_tree_support(clades)
    frequency_report: TreeSetCladeSupportReport | None = None
    frequency_rows: dict[tuple[str, ...], TreeSetCladeSupportRow] = {}
    if comparison_tree_set_path is not None:
        frequency_report = compute_reference_tree_clade_support(
            tree_path,
            comparison_tree_set_path,
        )
        frequency_rows = _clade_support_row_lookup(frequency_report.rows)
    rows = [
        _build_clade_support_row(
            tree_path=tree_path,
            comparison_tree_set_path=comparison_tree_set_path,
            support_row=support_row,
            frequency_row=frequency_rows.get(tuple(support_row.descendant_taxa)),
        )
        for support_row in support_rows
    ]
    columns = _clade_support_table_columns()
    _write_clade_support_rows(path, columns=columns, rows=rows)
    return SupplementaryCladeSupportTableResult(
        output_path=path,
        row_count=len(rows),
        supported_clade_count=sum(1 for row in rows if row.support is not None),
        frequency_scored_clade_count=sum(
            1 for row in rows if row.clade_frequency is not None
        ),
        frequency_partial_support_count=sum(
            1 for row in rows if row.frequency_status == "partial-support"
        ),
        frequency_absent_clade_count=sum(
            1 for row in rows if row.frequency_status == "absent"
        ),
        frequency_unscored_clade_count=sum(
            1 for row in rows if row.frequency_status == "not-counted"
        ),
        columns=columns,
        rows=rows,
    )


def _resolve_model_sidecar_path(
    iqtree_report_path: Path,
    model_sidecar_path: Path | None,
) -> Path | None:
    if model_sidecar_path is not None:
        return model_sidecar_path
    return resolve_iqtree_model_sidecar(iqtree_report_path.with_suffix(""))


def _build_model_selection_row(
    *,
    iqtree_report_path: Path,
    model_sidecar_path: Path | None,
    candidate: IqtreeModelCandidate,
    summary: IqtreeModelSelectionSummary,
) -> SupplementaryModelSelectionRow:
    return SupplementaryModelSelectionRow(
        iqtree_report_source=str(iqtree_report_path),
        model_sidecar_source=(
            None if model_sidecar_path is None else str(model_sidecar_path)
        ),
        rank=candidate.rank,
        model=candidate.model,
        log_likelihood=candidate.log_likelihood,
        parameter_count=candidate.parameter_count,
        aic=candidate.aic,
        aicc=candidate.aicc,
        bic=candidate.bic,
        best_aic=candidate.model == summary.best_model_aic,
        best_aicc=candidate.model == summary.best_model_aicc,
        best_bic=candidate.model == summary.best_model_bic,
        selected_model=candidate.model == summary.selected_model,
        selected_model_name=summary.selected_model or candidate.model,
        selected_criterion=summary.selected_criterion,
    )


def write_supplementary_model_selection_table(
    path: Path,
    *,
    iqtree_report_path: Path,
    model_sidecar_path: Path | None = None,
) -> SupplementaryModelSelectionTableResult:
    """Write one supplementary model-selection table from parsed IQ-TREE artifacts."""
    resolved_sidecar_path = _resolve_model_sidecar_path(
        iqtree_report_path,
        model_sidecar_path,
    )
    summary = parse_iqtree_model_selection_summary(
        iqtree_report_path=iqtree_report_path,
        model_sidecar_path=resolved_sidecar_path,
    )
    if summary is None or summary.selected_model is None:
        raise ValueError(
            "iqtree model-selection artifacts do not expose a selected model"
        )
    if not summary.candidates:
        raise ValueError(
            "iqtree model-selection artifacts do not expose candidate model rows"
        )
    rows = [
        _build_model_selection_row(
            iqtree_report_path=iqtree_report_path,
            model_sidecar_path=resolved_sidecar_path,
            candidate=candidate,
            summary=summary,
        )
        for candidate in summary.candidates
    ]
    columns = _model_selection_table_columns()
    _write_model_selection_rows(path, columns=columns, rows=rows)
    return SupplementaryModelSelectionTableResult(
        output_path=path,
        row_count=len(rows),
        selected_model=summary.selected_model,
        selected_criterion=summary.selected_criterion,
        candidate_count=summary.candidate_count,
        columns=columns,
        rows=rows,
    )


def _diversification_model_by_name(
    reports: list[DiversificationRateReport],
) -> dict[str, DiversificationRateReport]:
    return {report.model: report for report in reports}


def _build_diversification_row(
    *,
    observation: CladeDiversificationObservation,
    tree_path: Path,
    metadata_path: Path | None,
    clade_model: str,
    better_model: str,
    global_diversification_rate: float,
    yule_report: DiversificationRateReport,
    birth_death_report: DiversificationRateReport,
    sampling_report: SamplingFractionReport | None,
    warnings: list[str],
) -> SupplementaryDiversificationRow:
    return SupplementaryDiversificationRow(
        tree_source=str(tree_path),
        metadata_source=None if metadata_path is None else str(metadata_path),
        clade_model=clade_model,
        better_model=better_model,
        node=observation.node,
        node_name=observation.node_name,
        descendant_taxa=list(observation.descendant_taxa),
        tip_count=observation.tip_count,
        crown_age=observation.crown_age,
        clade_diversification_rate=observation.diversification_rate,
        clade_rate_z_score=observation.z_score,
        clade_classification=observation.classification,
        global_diversification_rate=global_diversification_rate,
        yule_log_likelihood=yule_report.log_likelihood,
        yule_aic=yule_report.aic,
        yule_corrected_tip_count=yule_report.corrected_tip_count,
        yule_sampling_fraction=yule_report.sampling_fraction,
        yule_net_diversification_rate=yule_report.net_diversification_rate,
        yule_relative_extinction=yule_report.relative_extinction,
        birth_death_log_likelihood=birth_death_report.log_likelihood,
        birth_death_aic=birth_death_report.aic,
        birth_death_corrected_tip_count=birth_death_report.corrected_tip_count,
        birth_death_sampling_fraction=birth_death_report.sampling_fraction,
        birth_death_net_diversification_rate=(
            birth_death_report.net_diversification_rate
        ),
        birth_death_relative_extinction=birth_death_report.relative_extinction,
        sampling_metadata_complete=(
            None if sampling_report is None else sampling_report.complete
        ),
        sampling_column=None if sampling_report is None else sampling_report.sampling_column,
        sampling_fraction=(
            None if sampling_report is None else sampling_report.sampling_fraction
        ),
        sampling_heterogeneous=(
            None if sampling_report is None else sampling_report.heterogeneous_values
        ),
        sampling_missing_taxa=(
            [] if sampling_report is None else list(sampling_report.missing_taxa)
        ),
        sampling_invalid_rows=(
            []
            if sampling_report is None
            else [
                _serialize_sampling_issue(issue)
                for issue in sampling_report.invalid_rows
            ]
        ),
        warning_count=len(warnings),
        warnings=list(warnings),
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
    """Write one supplementary diversification table over clade and model evidence."""
    if clade_model not in {"yule", "birth-death"}:
        raise ValueError(
            "clade_model must be 'yule' or 'birth-death'"
        )
    sampling_report = (
        None
        if metadata_path is None
        else detect_incomplete_taxon_sampling_metadata(
            tree_path,
            metadata_path,
            taxon_column=taxon_column,
            sampling_column=sampling_column,
        )
    )
    comparison = compare_diversification_models(
        tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
    )
    diversification_reports = [
        estimate_diversification_rate(
            tree_path,
            metadata_path=metadata_path,
            taxon_column=taxon_column,
            sampling_column=sampling_column,
            model="yule",
        ),
        estimate_diversification_rate(
            tree_path,
            metadata_path=metadata_path,
            taxon_column=taxon_column,
            sampling_column=sampling_column,
            model="birth-death",
        ),
    ]
    reports_by_name = _diversification_model_by_name(diversification_reports)
    clade_report = detect_diversification_outlier_clades(
        tree_path,
        model=clade_model,
    )
    warnings = sorted(
        set(
            clade_report.warnings
            + reports_by_name["yule"].warnings
            + reports_by_name["birth-death"].warnings
            + ([] if sampling_report is None else sampling_report.warnings)
        )
    )
    rows = [
        _build_diversification_row(
            observation=observation,
            tree_path=tree_path,
            metadata_path=metadata_path,
            clade_model=clade_model,
            better_model=comparison.better_model,
            global_diversification_rate=clade_report.global_rate,
            yule_report=reports_by_name["yule"],
            birth_death_report=reports_by_name["birth-death"],
            sampling_report=sampling_report,
            warnings=warnings,
        )
        for observation in clade_report.observations
    ]
    columns = _diversification_table_columns()
    _write_diversification_rows(path, columns=columns, rows=rows)
    return SupplementaryDiversificationTableResult(
        output_path=path,
        row_count=len(rows),
        better_model=comparison.better_model,
        clade_model=clade_model,
        high_clade_count=len(clade_report.high_diversification_clades),
        low_clade_count=len(clade_report.low_diversification_clades),
        warning_count=len(warnings),
        sampling_metadata_complete=(
            None if sampling_report is None else sampling_report.complete
        ),
        columns=columns,
        rows=rows,
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
