from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class AncestralTransitionBranchRow:
    """One branchwise ancestral transition review row."""

    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    branch_length: float | None
    parent_most_likely_state: str
    child_most_likely_state: str
    parent_state_set: list[str]
    child_state_set: list[str]
    overlapping_states: list[str]
    changed: bool
    transition: str
    certainty_class: str


@dataclass(frozen=True, slots=True)
class AncestralTransitionCountRow:
    """One transition-pair count summary for a discrete ancestral reconstruction."""

    transition: str
    source_state: str
    target_state: str
    certain_change_count: int
    uncertain_change_count: int
    total_change_count: int


@dataclass(frozen=True, slots=True)
class AncestralTransitionExclusion:
    """One excluded tip from ancestral transition counting."""

    taxon: str
    reason: str


@dataclass(frozen=True, slots=True)
class AncestralTransitionSummary:
    """Reviewer-facing summary for one ancestral transition count report."""

    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    total_branch_count: int
    changed_branch_count: int
    certain_change_count: int
    uncertain_change_count: int
    unique_transition_count: int
    warning_count: int


@dataclass(slots=True)
class AncestralTransitionReport:
    """Discrete ancestral transition counts for one analyzed rooted tree."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    ordered_states: list[str]
    taxon_count: int
    branch_rows: list[AncestralTransitionBranchRow]
    transition_rows: list[AncestralTransitionCountRow]
    exclusions: list[AncestralTransitionExclusion]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class AncestralTransitionTreeRow:
    """One retained tree summary from ancestral transition counting over a tree set."""

    source_tree_index: int
    post_burnin_index: int
    rooted_topology_id: str
    unrooted_topology_id: str
    total_branch_count: int
    changed_branch_count: int
    certain_change_count: int
    uncertain_change_count: int


@dataclass(frozen=True, slots=True)
class AncestralTransitionTreeSetBranchRow:
    """One branchwise ancestral transition row from one retained tree."""

    source_tree_index: int
    post_burnin_index: int
    rooted_topology_id: str
    unrooted_topology_id: str
    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    branch_length: float | None
    parent_most_likely_state: str
    child_most_likely_state: str
    parent_state_set: list[str]
    child_state_set: list[str]
    overlapping_states: list[str]
    changed: bool
    transition: str
    certainty_class: str


@dataclass(frozen=True, slots=True)
class AncestralTransitionTreeSetCountRow:
    """One transition-pair summary across retained trees."""

    transition: str
    source_state: str
    target_state: str
    tree_presence_count: int
    tree_presence_fraction: float
    mean_certain_change_count: float
    mean_uncertain_change_count: float
    mean_total_change_count: float
    minimum_total_change_count: int
    maximum_total_change_count: int
    lower_95_total_change_count: float
    upper_95_total_change_count: float
    stability_class: str


@dataclass(frozen=True, slots=True)
class AncestralTransitionTreeSetSummary:
    """Reviewer-facing summary for ancestral transition counts across a tree set."""

    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    shared_tree_taxon_count: int
    analysis_taxon_count: int
    rooted_topology_count: int
    unrooted_topology_count: int
    transition_pair_count: int
    topology_sensitive_transition_pair_count: int
    uncertainty_sensitive_transition_pair_count: int
    warning_count: int


@dataclass(slots=True)
class AncestralTransitionTreeSetReport:
    """Discrete ancestral transition counts aggregated across a posterior/bootstrap tree set."""

    tree_set_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    ordered_states: list[str]
    burnin_fraction: float
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    shared_tree_taxa: list[str]
    analysis_taxa: list[str]
    rooted_topology_count: int
    unrooted_topology_count: int
    tree_rows: list[AncestralTransitionTreeRow]
    branch_rows: list[AncestralTransitionTreeSetBranchRow]
    transition_rows: list[AncestralTransitionTreeSetCountRow]
    exclusions: list[AncestralTransitionExclusion]
    warnings: list[str]
