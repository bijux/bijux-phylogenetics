from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


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


@dataclass(frozen=True, slots=True)
class SupplementaryCladeSupportRow:
    """One reviewer-facing clade-support row for one clade in one reference tree."""

    tree_source: str
    comparison_tree_set_source: str | None
    clade_id: str
    node_kind: str
    node_label: str | None
    descendant_taxa: list[str]
    support: float | None
    support_fraction: float | None
    support_class: str
    support_method: str
    branch_length: float | None
    root_depth: float | None
    supporting_tree_count: int | None
    clade_frequency: float | None
    support_percent: float | None
    frequency_method: str | None
    frequency_status: str | None
    frequency_explanation: str | None


@dataclass(slots=True)
class SupplementaryCladeSupportTableResult:
    output_path: Path
    row_count: int
    supported_clade_count: int
    frequency_scored_clade_count: int
    frequency_partial_support_count: int
    frequency_absent_clade_count: int
    frequency_unscored_clade_count: int
    columns: list[str]
    rows: list[SupplementaryCladeSupportRow]


@dataclass(frozen=True, slots=True)
class SupplementaryModelSelectionRow:
    """One reviewer-facing row for one candidate substitution model."""

    iqtree_report_source: str
    model_sidecar_source: str | None
    rank: int
    model: str
    log_likelihood: float
    parameter_count: int | None
    aic: float
    aicc: float
    bic: float
    best_aic: bool
    best_aicc: bool
    best_bic: bool
    selected_model: bool
    selected_model_name: str
    selected_criterion: str | None


@dataclass(slots=True)
class SupplementaryModelSelectionTableResult:
    output_path: Path
    row_count: int
    selected_model: str
    selected_criterion: str | None
    candidate_count: int
    columns: list[str]
    rows: list[SupplementaryModelSelectionRow]


@dataclass(frozen=True, slots=True)
class SupplementaryComparativeModelRow:
    """One reviewer-facing coefficient row for one comparative candidate model."""

    tree_source: str
    traits_source: str
    response: str
    formula: str
    model_family: str
    selected_criterion: str
    best_formula: str
    rank: int
    selected: bool
    analysis_taxon_count: int
    excluded_taxon_count: int
    excluded_taxa: list[str]
    encoded_columns: list[str]
    coefficient_name: str
    estimate: float
    standard_error: float
    test_statistic: float
    p_value: float
    lower_95_confidence_interval: float
    upper_95_confidence_interval: float
    inference_distribution: str
    phylogenetic_parameter_name: str
    phylogenetic_parameter_value: float
    phylogenetic_parameter_estimated: bool
    log_likelihood: float
    aic: float
    aicc: float
    bic: float
    delta_aicc: float
    delta_bic: float
    akaike_weight: float
    residual_mean: float | None
    outlier_taxon_count: int
    outlier_taxa: list[str]
    max_leverage: float | None
    max_abs_standardized_residual: float | None
    converged: bool | None
    iteration_count: int | None
    separation_detected: bool
    diagnostic_warning_count: int
    diagnostic_warnings: list[str]


@dataclass(slots=True)
class SupplementaryComparativeModelTableResult:
    output_path: Path
    row_count: int
    model_count: int
    selected_formula: str
    selected_criterion: str
    excluded_taxon_count: int
    columns: list[str]
    rows: list[SupplementaryComparativeModelRow]


@dataclass(frozen=True, slots=True)
class SupplementaryAncestralStateRow:
    """One reviewer-facing internal-node row for one ancestral reconstruction."""

    tree_source: str
    traits_source: str
    trait: str
    reconstruction_kind: str
    model: str
    estimator: str | None
    state_ordering: str | None
    root_prior_mode: str | None
    fixed_root_state: str | None
    alpha: float | None
    analysis_taxon_count: int
    excluded_taxon_count: int
    excluded_taxa: list[str]
    warning_count: int
    warnings: list[str]
    node: str
    node_name: str | None
    descendant_taxa: list[str]
    descendant_taxon_count: int
    estimate_value: float | None
    most_likely_state: str | None
    state_set: list[str]
    state_probabilities: dict[str, float]
    standard_error: float | None
    lower_95_interval: float | None
    upper_95_interval: float | None
    confidence: float
    ambiguous: bool | None
    unstable: bool
    interpretation: str
    downstream_risks: list[str]


@dataclass(slots=True)
class SupplementaryAncestralStateTableResult:
    output_path: Path
    row_count: int
    reconstruction_kind: str
    model: str
    analysis_taxon_count: int
    excluded_taxon_count: int
    unstable_node_count: int
    columns: list[str]
    rows: list[SupplementaryAncestralStateRow]


@dataclass(frozen=True, slots=True)
class SupplementaryDiversificationRow:
    """One reviewer-facing clade row for diversification review."""

    tree_source: str
    metadata_source: str | None
    clade_model: str
    better_model: str
    node: str
    node_name: str | None
    descendant_taxa: list[str]
    tip_count: int
    crown_age: float
    clade_diversification_rate: float
    clade_rate_z_score: float
    clade_classification: str
    global_diversification_rate: float
    yule_log_likelihood: float
    yule_aic: float
    yule_corrected_tip_count: float
    yule_sampling_fraction: float
    yule_net_diversification_rate: float
    yule_relative_extinction: float
    birth_death_log_likelihood: float
    birth_death_aic: float
    birth_death_corrected_tip_count: float
    birth_death_sampling_fraction: float
    birth_death_net_diversification_rate: float
    birth_death_relative_extinction: float
    sampling_metadata_complete: bool | None
    sampling_column: str | None
    sampling_fraction: float | None
    sampling_heterogeneous: bool | None
    sampling_missing_taxa: list[str]
    sampling_invalid_rows: list[str]
    warning_count: int
    warnings: list[str]


@dataclass(slots=True)
class SupplementaryDiversificationTableResult:
    output_path: Path
    row_count: int
    better_model: str
    clade_model: str
    high_clade_count: int
    low_clade_count: int
    warning_count: int
    sampling_metadata_complete: bool | None
    columns: list[str]
    rows: list[SupplementaryDiversificationRow]


@dataclass(frozen=True, slots=True)
class SupplementaryBatchSummaryRow:
    """One reviewer-facing dataset or variant row over a written batch bundle."""

    row_scope: str
    dataset_id: str
    dataset_label: str
    workflow_status: str
    variant_id: str | None
    label: str | None
    execution_mode: str | None
    task_status: str | None
    job_status: str | None
    output_freshness_status: str | None
    recovery_action: str | None
    merge_status: str | None
    evidence_status: str | None
    reproducibility_status: str | None
    selected_model: str | None
    output_root: str | None
    task_log_path: str | None
    evidence_json_path: str | None
    evidence_html_path: str | None
    variant_count: int | None
    successful_variant_count: int | None
    failed_variant_count: int | None
    output_file_count: int | None
    output_byte_count: int | None
    artifact_file_count: int | None
    linked_artifact_count: int | None
    linked_artifact_bytes: int | None
    issue_count: int
    issues: list[str]
    error_code: str | None
    error_message: str | None
    job_evidence_warning_count: int | None
    warning_count: int
    warnings: list[str]


@dataclass(slots=True)
class SupplementaryBatchSummaryTableResult:
    output_path: Path
    row_count: int
    dataset_row_count: int
    variant_row_count: int
    workflow_status: str
    warning_count: int
    columns: list[str]
    rows: list[SupplementaryBatchSummaryRow]
