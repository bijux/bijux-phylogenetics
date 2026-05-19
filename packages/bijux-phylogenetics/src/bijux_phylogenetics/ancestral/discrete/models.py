from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class DiscreteAncestralEstimate:
    """One discrete ancestral-state estimate for a tree node."""

    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    state_set: list[str]
    most_likely_state: str
    state_probabilities: dict[str, float]
    ambiguous: bool
    confidence: float
    interpretation: str
    unstable: bool
    downstream_risks: list[str]


@dataclass(slots=True)
class DiscreteTransitionRateRow:
    """One directed transition rate from a fitted discrete likelihood model."""

    source_state: str
    target_state: str
    transition_allowed: bool
    step_distance: int
    rate: float


@dataclass(slots=True)
class DiscreteOptimizerDiagnostics:
    """Optimizer state for one discrete likelihood fit."""

    optimizer_name: str
    parameter_count: int
    initial_candidate_count: int
    best_initial_scale: float
    converged: bool
    iteration_count: int
    function_evaluation_count: int
    simplex_shrink_count: int
    hit_lower_parameter_bound: bool
    hit_upper_parameter_bound: bool


@dataclass(slots=True)
class DiscreteModelBaselineComparison:
    """Likelihood-model comparison against the equal-rates baseline."""

    baseline_model: str
    baseline_log_likelihood: float
    baseline_parameter_count: int
    baseline_aic: float
    delta_log_likelihood: float
    delta_aic: float
    preferred_model_by_aic: str


@dataclass(slots=True)
class DiscreteRerootingMethodCompatibility:
    """Compatibility of the owned marginal surface with `phytools::rerootingMethod`."""

    comparable: bool
    reference_model: str | None
    reference_root_prior_mode: str | None
    notes: list[str]


@dataclass(slots=True)
class DiscreteLikelihoodFitResult:
    """Internal likelihood fit details for one discrete ancestral reconstruction."""

    estimates: list[DiscreteAncestralEstimate]
    ordered_states: list[str]
    state_order: list[str]
    rerooting_method_compatibility: DiscreteRerootingMethodCompatibility
    log_likelihood: float
    parameter_count: int
    aic: float
    transition_rate_rows: list[DiscreteTransitionRateRow]
    allowed_transition_pairs: list[tuple[str, str]]
    optimizer_diagnostics: DiscreteOptimizerDiagnostics
    overparameterized: bool
    baseline_comparison: DiscreteModelBaselineComparison | None


@dataclass(slots=True)
class DiscreteAncestralReport:
    """Discrete ancestral-state reconstruction report."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    state_ordering: str
    root_prior_mode: str | None
    fixed_root_state: str | None
    ordered_states: list[str]
    taxon_count: int
    observed_states: list[str]
    state_counts: dict[str, int]
    sparse_states: list[str]
    analysis_tree_newick: str
    dropped_missing_taxa: list[str]
    minimal_change_count: int | None
    parsimonious_root_state_count: int | None
    warnings: list[str]
    unstable_nodes: list[str]
    weak_support_nodes: list[str]
    estimates: list[DiscreteAncestralEstimate]
    rerooting_method_compatibility: DiscreteRerootingMethodCompatibility
    log_likelihood: float | None
    parameter_count: int | None
    aic: float | None
    transition_rate_rows: list[DiscreteTransitionRateRow]
    allowed_transition_pairs: list[tuple[str, str]]
    optimizer_diagnostics: DiscreteOptimizerDiagnostics | None
    overparameterized: bool
    baseline_comparison: DiscreteModelBaselineComparison | None


@dataclass(slots=True)
class DiscreteAncestralSummary:
    """Reviewer-facing summary for one discrete ancestral reconstruction."""

    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    root_prior_mode: str | None
    fixed_root_state: str | None
    analyzed_taxon_count: int
    excluded_taxon_count: int
    internal_node_count: int
    ambiguous_internal_node_count: int
    unstable_node_count: int
    weak_support_node_count: int
    observed_state_count: int
    sparse_state_count: int
    minimal_change_count: int | None
    parsimonious_root_state_count: int | None
    root_node: str
    root_most_likely_state: str
    root_confidence: float
    phytools_rerooting_method_comparable: bool
    log_likelihood: float | None
    parameter_count: int | None
    aic: float | None
    optimizer_converged: bool | None
    optimizer_iteration_count: int | None
    optimizer_function_evaluation_count: int | None
    overparameterized: bool
    baseline_model: str | None
    baseline_delta_aic: float | None
    preferred_model_by_aic: str | None
    warning_count: int


@dataclass(slots=True)
class DiscreteAncestralExclusion:
    """One excluded tip from a discrete ancestral reconstruction."""

    taxon: str
    reason: str
