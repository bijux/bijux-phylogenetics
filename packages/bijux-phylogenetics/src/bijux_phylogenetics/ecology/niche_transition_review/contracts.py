from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class NicheStateNodeRow:
    """One internal-node ecological niche reconstruction row."""

    node: str
    node_name: str | None
    descendant_taxa: list[str]
    most_likely_niche: str
    niche_probabilities: dict[str, float]
    confidence: float
    ambiguous: bool
    is_root: bool


@dataclass(frozen=True, slots=True)
class NicheTransitionRateRow:
    """One fitted ecological niche transition-rate row."""

    source_niche: str
    target_niche: str
    transition_allowed: bool
    step_distance: int
    rate: float


@dataclass(frozen=True, slots=True)
class NicheTransitionBranchRow:
    """One branchwise ecological niche transition review row."""

    branch_id: str
    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    branch_length: float | None
    parent_most_likely_niche: str
    child_most_likely_niche: str
    parent_niche_set: list[str]
    child_niche_set: list[str]
    overlapping_niches: list[str]
    changed: bool
    transition: str
    certainty_class: str
    support: float
    strongly_supported: bool
    parent_confidence: float
    child_confidence: float


@dataclass(frozen=True, slots=True)
class NicheTransitionCountRow:
    """One aggregated ecological niche transition count row."""

    transition: str
    source_niche: str
    target_niche: str
    certain_transition_count: int
    uncertain_transition_count: int
    total_transition_count: int
    strongly_supported_transition_count: int


@dataclass(frozen=True, slots=True)
class NicheTransitionCladeRow:
    """One internal-clade ecological niche shift burden row."""

    node: str
    node_name: str | None
    descendant_taxa: list[str]
    descendant_taxon_count: int
    descendant_internal_node_count: int
    changed_branch_count: int
    certain_transition_count: int
    uncertain_transition_count: int
    strongly_supported_transition_count: int
    transition_diversity: int
    dominant_transition: str
    dominant_transition_count: int
    shift_burden_score: float
    contains_repeated_shifts: bool
    rank: int


@dataclass(frozen=True, slots=True)
class NicheTransitionExclusionRow:
    """One excluded metadata row from ecological niche transition analysis."""

    taxon: str
    raw_niche: str
    normalized_niche: str | None
    reason: str
    note: str


@dataclass(frozen=True, slots=True)
class NicheTransitionSummary:
    """Reviewer-facing summary for ecological niche transition analysis."""

    trait: str
    taxon_column: str
    model: str
    internal_model: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    observed_niche_count: int
    internal_node_count: int
    ambiguous_internal_node_count: int
    log_likelihood: float
    parameter_count: int
    aic: float
    transition_rate_row_count: int
    changed_branch_count: int
    certain_transition_count: int
    uncertain_transition_count: int
    strongly_supported_transition_count: int
    clade_shift_row_count: int
    repeated_shift_clade_count: int
    root_niche: str
    root_confidence: float
    warning_count: int


@dataclass(slots=True)
class NicheTransitionReport:
    """Owned ecological niche transition review surface for one tree."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    internal_model: str
    summary: NicheTransitionSummary
    node_rows: list[NicheStateNodeRow]
    rate_rows: list[NicheTransitionRateRow]
    branch_rows: list[NicheTransitionBranchRow]
    count_rows: list[NicheTransitionCountRow]
    clade_rows: list[NicheTransitionCladeRow]
    exclusion_rows: list[NicheTransitionExclusionRow]
    warnings: list[str]
