from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class AncestralTreeSetTreeRow:
    """One retained tree from an ancestral reconstruction tree-set review."""

    source_tree_index: int
    post_burnin_index: int
    rooted_topology_id: str
    unrooted_topology_id: str
    internal_clade_count: int


@dataclass(frozen=True, slots=True)
class AncestralTreeSetExclusion:
    """One taxon excluded before ancestral tree-set reconstruction."""

    taxon: str
    reason: str


@dataclass(frozen=True, slots=True)
class ContinuousAncestralTreeSetNodeRow:
    """One internal-node continuous ancestral reconstruction from one retained tree."""

    source_tree_index: int
    post_burnin_index: int
    rooted_topology_id: str
    unrooted_topology_id: str
    clade_id: str
    clade_taxa: list[str]
    estimate: float
    standard_error: float
    lower_95_interval: float
    upper_95_interval: float
    confidence: float
    unstable: bool


@dataclass(frozen=True, slots=True)
class ContinuousAncestralTreeSetCladeSummaryRow:
    """One comparable clade summary across retained continuous reconstructions."""

    clade_id: str
    clade_taxa: list[str]
    tree_presence_count: int
    tree_presence_fraction: float
    mean_estimate: float
    median_estimate: float
    standard_deviation: float
    minimum_estimate: float
    maximum_estimate: float
    lower_95_empirical_estimate: float
    upper_95_empirical_estimate: float
    empirical_interval_width: float
    mean_standard_error: float
    unstable_tree_count: int
    unstable_tree_fraction: float
    instability_score: float
    stability_class: str


@dataclass(frozen=True, slots=True)
class ContinuousAncestralTreeSetSummary:
    """Reviewer-facing summary for one continuous ancestral tree-set analysis."""

    trait: str
    taxon_column: str
    model: str
    alpha: float
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    shared_tree_taxon_count: int
    analysis_taxon_count: int
    rooted_topology_count: int
    unrooted_topology_count: int
    clade_summary_count: int
    unstable_clade_count: int
    top_unstable_clade: str | None
    warning_count: int


@dataclass(slots=True)
class ContinuousAncestralTreeSetReport:
    """Continuous ancestral reconstruction stability across a posterior or bootstrap tree set."""

    tree_set_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    alpha: float
    burnin_fraction: float
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    shared_tree_taxa: list[str]
    analysis_taxa: list[str]
    rooted_topology_count: int
    unrooted_topology_count: int
    tree_rows: list[AncestralTreeSetTreeRow]
    node_rows: list[ContinuousAncestralTreeSetNodeRow]
    clade_summaries: list[ContinuousAncestralTreeSetCladeSummaryRow]
    exclusions: list[AncestralTreeSetExclusion]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class DiscreteAncestralTreeSetNodeRow:
    """One internal-node discrete ancestral reconstruction from one retained tree."""

    source_tree_index: int
    post_burnin_index: int
    rooted_topology_id: str
    unrooted_topology_id: str
    clade_id: str
    clade_taxa: list[str]
    most_likely_state: str
    state_set: list[str]
    confidence: float
    ambiguous: bool
    unstable: bool


@dataclass(frozen=True, slots=True)
class DiscreteAncestralTreeSetCladeSummaryRow:
    """One comparable clade summary across retained discrete reconstructions."""

    clade_id: str
    clade_taxa: list[str]
    tree_presence_count: int
    tree_presence_fraction: float
    unique_state_count: int
    dominant_state: str
    dominant_state_tree_count: int
    dominant_state_fraction: float
    ambiguous_tree_count: int
    ambiguous_tree_fraction: float
    unstable_tree_count: int
    unstable_tree_fraction: float
    state_distribution: dict[str, int]
    instability_score: float
    stability_class: str


@dataclass(frozen=True, slots=True)
class DiscreteAncestralTreeSetSummary:
    """Reviewer-facing summary for one discrete ancestral tree-set analysis."""

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
    clade_summary_count: int
    unstable_clade_count: int
    top_unstable_clade: str | None
    warning_count: int


@dataclass(slots=True)
class DiscreteAncestralTreeSetReport:
    """Discrete ancestral reconstruction stability across a posterior or bootstrap tree set."""

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
    tree_rows: list[AncestralTreeSetTreeRow]
    node_rows: list[DiscreteAncestralTreeSetNodeRow]
    clade_summaries: list[DiscreteAncestralTreeSetCladeSummaryRow]
    exclusions: list[AncestralTreeSetExclusion]
    warnings: list[str]
