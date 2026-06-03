from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.bayesian.posterior_sets.burnin import (
    BurninSensitivityCladeShift,
    BurninSensitivityParameterShift,
)


@dataclass(slots=True)
class CalibrationValidationIssue:
    calibration_id: str
    code: str
    message: str


@dataclass(slots=True)
class ValidatedCalibration:
    calibration_id: str
    target_kind: str
    target_label: str
    taxa: list[str]
    minimum_age: float | None
    maximum_age: float | None
    distribution: str
    valid: bool


@dataclass(slots=True)
class FossilCalibrationValidationReport:
    tree_path: Path
    calibration_path: Path
    tree_taxa: list[str]
    calibration_count: int
    valid_calibration_count: int
    invalid_calibration_count: int
    calibrations: list[ValidatedCalibration]
    issues: list[CalibrationValidationIssue]


@dataclass(slots=True)
class ImpossibleCalibrationConstraintReport:
    tree_path: Path
    calibration_path: Path
    impossible_calibration_ids: list[str]
    issues: list[CalibrationValidationIssue]


@dataclass(slots=True)
class ValidatedTipDate:
    taxon: str
    date: float | None
    valid: bool


@dataclass(slots=True)
class TipDatingValidationIssue:
    taxon: str
    code: str
    message: str


@dataclass(slots=True)
class TipDatingValidationReport:
    tree_path: Path
    tip_dates_path: Path
    alignment_path: Path | None
    taxon_column: str
    date_column: str
    valid_tip_count: int
    invalid_tip_count: int
    missing_tree_taxa: list[str]
    extra_tip_taxa: list[str]
    extra_alignment_taxa: list[str]
    tip_dates: list[ValidatedTipDate]
    issues: list[TipDatingValidationIssue]


@dataclass(slots=True)
class CalibrationDominanceObservation:
    calibration_id: str
    target_label: str
    bounded_span_fraction: float | None
    dominates_root_age: bool
    warning: str | None


@dataclass(slots=True)
class CalibrationDominanceReport:
    tree_path: Path
    calibration_path: Path
    root_age: float
    valid_calibration_count: int
    dominant_calibration_ids: list[str]
    observations: list[CalibrationDominanceObservation]
    warnings: list[str]


@dataclass(slots=True)
class TimeTreeReadinessReport:
    tree_path: Path
    calibration_path: Path | None
    tip_dates_path: Path | None
    decision: str
    rooted: bool
    ultrametric: bool
    branch_length_status: str
    blockers: list[str]
    warnings: list[str]
    calibration_report: FossilCalibrationValidationReport | None
    tip_date_report: TipDatingValidationReport | None
    calibration_dominance: CalibrationDominanceReport | None


@dataclass(slots=True)
class BeastPreparationReport:
    alignment_path: Path
    output_path: Path
    tree_path: Path | None
    calibration_path: Path | None
    tip_dates_path: Path | None
    taxon_count: int
    character_count: int
    inferred_alphabet: str
    beast_data_type: str
    substitution_model: str
    clock_model: str
    tree_prior: str
    starting_tree_source: str
    chain_length: int
    log_every: int
    calibration_count: int
    tip_date_count: int
    warning_count: int
    warnings: list[str]
    log_path: Path
    tree_log_path: Path
    calibrations: list[BeastCalibration]


@dataclass(slots=True)
class BeastAnalysisXmlIssue:
    code: str
    message: str


@dataclass(slots=True)
class BeastAnalysisXmlLogger:
    logger_kind: str
    file_name: str | None
    log_every: int | None


@dataclass(slots=True)
class BeastAnalysisXmlReport:
    path: Path
    beast_version: str | None
    beast_namespace: str | None
    taxon_count: int
    character_count: int
    beast_data_type: str | None
    substitution_model: str | None
    clock_model: str | None
    tree_prior: str | None
    starting_tree_source: str | None
    chain_length: int | None
    state_node_count: int
    logger_count: int
    posterior_log_path: Path | None
    posterior_tree_path: Path | None
    calibration_count: int
    calibration_ids: list[str]
    tip_date_count: int
    tip_date_units: str | None
    tip_date_direction: str | None
    issues: list[BeastAnalysisXmlIssue]
    valid: bool


@dataclass(slots=True)
class BeastCalibration:
    calibration_id: str
    beast_distribution: str
    target_label: str
    lower_bound: float | None
    upper_bound: float | None
    translated: bool
    translation_note: str | None


@dataclass(slots=True)
class BeastLogRow:
    state: int
    values: dict[str, float]


@dataclass(slots=True)
class BeastLogReport:
    path: Path
    row_count: int
    columns: list[str]
    rows: list[BeastLogRow]


@dataclass(slots=True)
class BeastLogParameterSummary:
    parameter: str
    parameter_category: str
    sample_count: int
    effective_sample_size: float
    mean: float
    median: float
    standard_deviation: float
    minimum: float
    maximum: float
    hpd_95_lower: float
    hpd_95_upper: float
    first_half_mean: float
    second_half_mean: float
    standardized_mean_shift: float


@dataclass(slots=True)
class BeastLogSummaryReport:
    path: Path
    burnin_fraction: float
    burnin_row_count: int
    kept_row_count: int
    first_kept_state: int
    last_kept_state: int
    posterior_parameters: list[str]
    likelihood_parameters: list[str]
    prior_parameters: list[str]
    clock_parameters: list[str]
    tree_parameters: list[str]
    other_parameters: list[str]
    parameter_summaries: list[BeastLogParameterSummary]


@dataclass(slots=True)
class BeastPosteriorDecompositionRow:
    state: int
    log_posterior: float
    log_likelihood: float
    log_prior: float
    decomposition_delta: float
    decomposition_valid: bool


@dataclass(slots=True)
class BeastPosteriorDecompositionReport:
    path: Path
    burnin_fraction: float
    burnin_row_count: int
    kept_row_count: int
    first_kept_state: int
    last_kept_state: int
    posterior_term_source: str
    likelihood_term_source: str
    prior_term_source: str
    identity_tolerance: float
    verified: bool
    maximum_absolute_delta: float
    rows: list[BeastPosteriorDecompositionRow]


@dataclass(slots=True)
class BeastPosteriorTreeSample:
    tree_name: str
    state: int | None
    rooted: bool
    tip_names: list[str]
    newick: str
    annotation_key_count: int
    annotation_record_count: int
    annotation_keys: list[str]
    annotation_values: dict[str, str]


@dataclass(slots=True)
class BeastPosteriorClade:
    clade: str
    tree_count: int
    frequency: float


@dataclass(slots=True)
class BeastPosteriorTreeSetReport:
    path: Path
    burnin_fraction: float
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    rooted_tree_count: int
    sampled_states: list[int]
    tip_names: list[str]
    clades: list[BeastPosteriorClade]
    trees: list[BeastPosteriorTreeSample]


@dataclass(slots=True)
class BeastPosteriorConsensusReport:
    source_path: Path
    retained_tree_set_path: Path
    burnin_fraction: float
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    rooted_topology_count: int
    shared_taxa: list[str]
    consensus_newick: str
    clade_frequency_count: int
    annotated_node_count: int
    minimum_posterior_probability: float | None
    maximum_posterior_probability: float | None


@dataclass(slots=True)
class BeastPosteriorTopologyDiversityReport:
    source_path: Path
    retained_tree_set_path: Path
    burnin_fraction: float
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    rooted_topology_count: int
    dominant_topology_frequency: float
    effective_topology_count: float
    pair_count: int
    mean_robinson_foulds_distance: float
    mean_normalized_robinson_foulds_distance: float
    maximum_robinson_foulds_distance: int
    maximum_normalized_robinson_foulds_distance: float
    unstable_clade_count: int


@dataclass(slots=True)
class BeastLogValidationIssue:
    code: str
    message: str
    row: int | None = None
    column: str | None = None


@dataclass(slots=True)
class BeastPosteriorLogValidationReport:
    path: Path
    row_count: int
    state_count: int
    required_columns: list[str]
    observed_columns: list[str]
    missing_columns: list[str]
    issues: list[BeastLogValidationIssue]
    valid: bool


@dataclass(slots=True)
class BeastBurninSensitivitySlice:
    burnin_fraction: float
    burnin_tree_count: int
    kept_tree_count: int
    rooted_topology_count: int
    selected_tree_index: int
    clade_credibility_score: float
    consensus_newick: str
    clade_frequency_count: int
    kept_row_count: int | None
    first_kept_state: int | None
    last_kept_state: int | None
    posterior_mean: float | None
    likelihood_mean: float | None
    tree_height_mean: float | None


@dataclass(slots=True)
class BeastBurninSensitivityReport:
    posterior_tree_path: Path
    log_path: Path | None
    slices: list[BeastBurninSensitivitySlice]
    changed_mcc_count: int
    changed_consensus_count: int
    parameter_shifts: list[BurninSensitivityParameterShift]
    clade_shifts: list[BurninSensitivityCladeShift]
    unstable_parameter_count: int
    unstable_clade_count: int
    warnings: list[str]


@dataclass(slots=True)
class BeastChainMixingIssue:
    path: Path | None
    parameter: str
    code: str
    message: str
    observed_value: float
    threshold: float


@dataclass(slots=True)
class BeastChainMixingReport:
    log_paths: list[Path]
    chain_count: int
    converged: bool
    issues: list[BeastChainMixingIssue]
    chain_summaries: list[BeastConvergenceReport]


@dataclass(slots=True)
class BeastConvergenceReport:
    path: Path
    burnin_fraction: float
    burnin_row_count: int
    sample_count: int
    converged: bool
    ess_threshold: float
    mean_shift_threshold: float
    warnings: list[dict[str, object]]
    parameter_summaries: list[dict[str, object]]
