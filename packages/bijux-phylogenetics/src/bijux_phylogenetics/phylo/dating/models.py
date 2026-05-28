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
