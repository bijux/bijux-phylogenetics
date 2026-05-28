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
