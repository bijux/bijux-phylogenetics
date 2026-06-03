from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class ContinuousAncestralEstimate:
    """One continuous ancestral-state estimate for a tree node."""

    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    estimate: float
    standard_error: float
    lower_95_interval: float
    upper_95_interval: float
    uncertainty_width: float
    confidence: float
    interpretation: str
    unstable: bool
    downstream_risks: list[str]


@dataclass(slots=True)
class ContinuousAncestralBrownianFitDiagnostics:
    """Explicit Brownian covariance and solver diagnostics for one ancestral fit."""

    covariance_model: str
    tree_is_ultrametric: bool
    minimum_root_to_tip_depth: float
    maximum_root_to_tip_depth: float
    minimum_branch_length: float
    maximum_branch_length: float
    covariance_matrix_dimension: int
    covariance_matrix_rank: int
    covariance_singular: bool
    covariance_near_singular: bool
    covariance_positive_definite: bool
    covariance_condition_number: float
    covariance_log_determinant: float | None
    solver_name: str
    solver_regularized: bool
    solver_regularization_epsilon: float | None
    log_likelihood: float
    residual_sigma_squared: float


@dataclass(slots=True)
class ContinuousAncestralOptimizerDiagnostics:
    """Explicit optimizer or closed-form solver diagnostics for one continuous fit."""

    optimizer_name: str
    converged: bool
    iteration_count: int
    function_evaluation_count: int
    convergence_status: str
    message: str | None


@dataclass(slots=True)
class ContinuousAncestralReport:
    """Continuous ancestral-state reconstruction report."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    estimator: str
    alpha: float
    taxon_count: int
    analysis_tree_newick: str
    missing_from_traits_taxa: list[str]
    dropped_missing_taxa: list[str]
    dropped_non_numeric_taxa: list[str]
    warnings: list[str]
    unstable_nodes: list[str]
    weak_support_nodes: list[str]
    brownian_fit_diagnostics: ContinuousAncestralBrownianFitDiagnostics | None
    optimizer_diagnostics: ContinuousAncestralOptimizerDiagnostics | None
    estimates: list[ContinuousAncestralEstimate]


@dataclass(slots=True)
class ContinuousAncestralSummary:
    """Reviewer-facing summary for one continuous ancestral reconstruction."""

    trait: str
    taxon_column: str
    model: str
    estimator: str
    alpha: float
    analyzed_taxon_count: int
    excluded_taxon_count: int
    missing_tip_taxon_count: int
    non_numeric_tip_taxon_count: int
    internal_node_count: int
    unstable_node_count: int
    weak_support_node_count: int
    root_node: str
    root_estimate: float
    root_standard_error: float
    root_lower_95_interval: float
    root_upper_95_interval: float
    tree_is_ultrametric: bool | None
    covariance_near_singular: bool | None
    covariance_condition_number: float | None
    log_likelihood: float | None
    residual_sigma_squared: float | None
    optimizer_name: str | None
    optimizer_converged: bool | None
    optimizer_iteration_count: int | None
    optimizer_function_evaluation_count: int | None
    optimizer_convergence_status: str | None
    warning_count: int


@dataclass(slots=True)
class ContinuousAncestralExclusion:
    """One excluded tip from a continuous ancestral reconstruction."""

    taxon: str
    reason: str
