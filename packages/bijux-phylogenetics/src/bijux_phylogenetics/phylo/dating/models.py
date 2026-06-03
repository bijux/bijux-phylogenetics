from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class LeastSquaresDatingNodeRow:
    """One dated node recovered from a rooted substitution tree."""

    node_id: str
    node_kind: str
    node_label: str | None
    descendant_taxa: list[str]
    estimated_date: float
    fixed_tip_date: bool
    time_height: float


@dataclass(frozen=True, slots=True)
class LeastSquaresDatingBranchRow:
    """One branch residual row under a least-squares dated-tree fit."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    parent_date: float
    child_date: float
    fitted_time_duration: float
    observed_branch_length: float
    fitted_branch_length: float
    residual: float


@dataclass(slots=True)
class LeastSquaresDatingReport:
    """Closed-form least-squares dating fit on one rooted tree and one tip-date table."""

    tree_newick: str
    dated_tree_newick: str
    taxa: list[str]
    tip_count: int
    internal_node_count: int
    branch_count: int
    parameter_count: int
    tree_path: str | None
    metadata_path: str | None
    taxon_column: str
    date_column: str
    minimum_tip_date: float
    maximum_tip_date: float
    root_date: float
    estimated_clock_rate: float
    residual_sum_squares: float
    condition_number: float
    exact_fit: bool
    optimizer_name: str
    converged: bool
    node_rows: list[LeastSquaresDatingNodeRow]
    branch_rows: list[LeastSquaresDatingBranchRow]


@dataclass(frozen=True, slots=True)
class PenalizedLikelihoodDatingNodeRow:
    """One dated node recovered from one penalized likelihood dating fit."""

    node_id: str
    node_kind: str
    node_label: str | None
    descendant_taxa: list[str]
    estimated_date: float
    fixed_tip_date: bool
    time_height: float
    estimated_log_rate: float
    estimated_rate: float


@dataclass(frozen=True, slots=True)
class PenalizedLikelihoodDatingBranchRow:
    """One branch fit row under one penalized likelihood dating run."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    parent_date: float
    child_date: float
    fitted_time_duration: float
    observed_branch_length: float
    observed_log_rate: float
    fitted_log_rate: float
    estimated_branch_rate: float
    fitted_branch_length: float
    data_score_contribution: float
    smoothing_penalty_contribution: float


@dataclass(slots=True)
class PenalizedLikelihoodDatingReport:
    """Penalized likelihood dating fit on one rooted tree and one tip-date table."""

    tree_newick: str
    dated_tree_newick: str
    taxa: list[str]
    tip_count: int
    internal_node_count: int
    branch_count: int
    parameter_count: int
    tree_path: str | None
    metadata_path: str | None
    taxon_column: str
    date_column: str
    minimum_tip_date: float
    maximum_tip_date: float
    root_date: float
    smoothing_parameter: float
    data_score: float
    penalty_score: float
    total_score: float
    condition_number: float
    optimizer_name: str
    optimization_pass_count: int
    function_evaluation_count: int
    converged: bool
    node_rows: list[PenalizedLikelihoodDatingNodeRow]
    branch_rows: list[PenalizedLikelihoodDatingBranchRow]


@dataclass(frozen=True, slots=True)
class DatingCalibrationAnchor:
    """One fixed internal-node calibration resolved onto one rooted tree."""

    calibration_id: str
    target_kind: str
    target_label: str
    descendant_taxa: list[str]
    node_id: str
    node_kind: str
    fixed_date: float


@dataclass(frozen=True, slots=True)
class DatingCalibrationConstraintRow:
    """One calibration constraint resolved onto one rooted node."""

    calibration_id: str
    target_kind: str
    target_label: str
    descendant_taxa: list[str]
    node_id: str
    node_kind: str
    minimum_bound: float | None
    maximum_bound: float | None
    fixed_date: float | None
    contradictory: bool
    issue_codes: list[str]


@dataclass(frozen=True, slots=True)
class DatingCalibrationNodeWindowRow:
    """One calibrated node with aggregated and propagated feasible date bounds."""

    node_id: str
    node_kind: str
    node_label: str | None
    descendant_taxa: list[str]
    calibration_ids: list[str]
    minimum_bound: float | None
    maximum_bound: float | None
    effective_lower_bound: float | None
    effective_upper_bound: float | None
    contradictory: bool
    issue_codes: list[str]


@dataclass(frozen=True, slots=True)
class DatingCalibrationConstraintIssue:
    """One contradiction or validation issue in a dating calibration set."""

    scope_kind: str
    scope_id: str
    code: str
    message: str
    related_node_ids: list[str]
    related_calibration_ids: list[str]


@dataclass(slots=True)
class DatingCalibrationConstraintReport:
    """Dating calibration constraints resolved onto one rooted tree."""

    tree_newick: str
    taxa: list[str]
    tip_count: int
    internal_node_count: int
    tree_path: str
    calibration_path: str
    calibration_count: int
    valid_calibration_count: int
    invalid_calibration_count: int
    resolved_calibration_count: int
    contradictory_calibration_count: int
    contradictory_node_count: int
    feasible: bool
    constraint_rows: list[DatingCalibrationConstraintRow]
    node_window_rows: list[DatingCalibrationNodeWindowRow]
    issue_rows: list[DatingCalibrationConstraintIssue]


@dataclass(frozen=True, slots=True)
class PenalizedLikelihoodCrossValidationCandidateRow:
    """One smoothing-parameter candidate summarized over held-out calibrations."""

    smoothing_parameter: float
    fold_count: int
    mean_absolute_error: float
    mean_squared_error: float
    root_mean_squared_error: float
    max_absolute_error: float
    selected: bool


@dataclass(frozen=True, slots=True)
class PenalizedLikelihoodCrossValidationPredictionRow:
    """One held-out calibration prediction under one smoothing parameter."""

    smoothing_parameter: float
    held_out_calibration_id: str
    held_out_target_label: str
    held_out_descendant_taxa: list[str]
    held_out_node_id: str
    training_calibration_count: int
    observed_date: float
    predicted_date: float
    absolute_error: float
    squared_error: float
    optimization_pass_count: int
    function_evaluation_count: int
    converged: bool


@dataclass(slots=True)
class PenalizedLikelihoodCrossValidationReport:
    """Held-out calibration cross-validation over penalized dating smoothing values."""

    tree_newick: str
    taxa: list[str]
    tip_count: int
    internal_node_count: int
    branch_count: int
    tree_path: str | None
    metadata_path: str | None
    calibration_path: str | None
    taxon_column: str
    date_column: str
    usable_calibration_count: int
    candidate_count: int
    selected_smoothing_parameter: float
    selected_mean_absolute_error: float
    selected_mean_squared_error: float
    selected_root_mean_squared_error: float
    calibration_rows: list[DatingCalibrationAnchor]
    candidate_rows: list[PenalizedLikelihoodCrossValidationCandidateRow]
    prediction_rows: list[PenalizedLikelihoodCrossValidationPredictionRow]
    selected_fit: PenalizedLikelihoodDatingReport


@dataclass(frozen=True, slots=True)
class RelaxedRateBranchSummaryRow:
    """One branch-level rate summary recovered from one dated tree pair."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    substitution_branch_length: float
    dated_time_duration: float
    branch_rate: float
    rate_z_score: float
    outlier: bool


@dataclass(frozen=True, slots=True)
class RelaxedRateBranchOutlier:
    """One ranked branch-rate outlier from one dated tree pair."""

    rank: int
    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    substitution_branch_length: float
    dated_time_duration: float
    branch_rate: float
    rate_z_score: float


@dataclass(slots=True)
class RelaxedRateBranchSummaryReport:
    """Branch-rate summary computed from one substitution tree and one dated tree."""

    substitution_tree_newick: str
    dated_tree_newick: str
    taxa: list[str]
    tip_count: int
    internal_node_count: int
    branch_count: int
    substitution_tree_path: str | None
    dated_tree_path: str | None
    outlier_threshold: float
    mean_branch_rate: float
    standard_deviation_branch_rate: float
    minimum_branch_rate: float
    maximum_branch_rate: float
    outlier_count: int
    branch_rows: list[RelaxedRateBranchSummaryRow]
    outlier_rows: list[RelaxedRateBranchOutlier]
