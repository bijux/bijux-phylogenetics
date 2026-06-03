from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.biogeography.state_models import GeographicExcludedTaxonRow


@dataclass(frozen=True, slots=True)
class GeographicMigrationEventRow:
    """One inferred geographic movement event on one analyzed rooted tree."""

    branch_id: str
    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    branch_length: float
    parent_depth: float
    child_depth: float
    midpoint_depth: float
    source_region: str
    target_region: str
    support: float
    strongly_supported: bool
    confidence_class: str


@dataclass(frozen=True, slots=True)
class GeographicMigrationEventSummary:
    """Reviewer-facing summary for geographic movement events on one tree."""

    trait: str
    taxon_column: str
    model: str
    internal_model: str
    likelihood_method: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    tree_depth: float
    event_count: int
    strongly_supported_event_count: int
    mean_event_support: float
    earliest_midpoint_depth: float | None
    latest_midpoint_depth: float | None
    warning_count: int


@dataclass(slots=True)
class GeographicMigrationEventReport:
    """Owned geographic movement-event review surface for one rooted tree."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    internal_model: str
    likelihood_method: str
    summary: GeographicMigrationEventSummary
    event_rows: list[GeographicMigrationEventRow]
    exclusion_rows: list[GeographicExcludedTaxonRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class GeographicMigrationTreeRow:
    """One retained tree summary from geographic movement-event tree-set review."""

    source_tree_index: int
    post_burnin_index: int
    rooted_topology_id: str
    unrooted_topology_id: str
    event_count: int
    strongly_supported_event_count: int


@dataclass(frozen=True, slots=True)
class GeographicMigrationTreeSetEventRow:
    """One inferred geographic movement event from one retained tree."""

    source_tree_index: int
    post_burnin_index: int
    rooted_topology_id: str
    unrooted_topology_id: str
    branch_id: str
    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    branch_length: float
    parent_depth: float
    child_depth: float
    midpoint_depth: float
    source_region: str
    target_region: str
    support: float
    strongly_supported: bool
    confidence_class: str


@dataclass(frozen=True, slots=True)
class GeographicMigrationTreeSetEventSummaryRow:
    """One comparable event summary across retained trees."""

    branch_id: str
    child_descendant_taxa: list[str]
    source_region: str
    target_region: str
    tree_presence_count: int
    tree_presence_fraction: float
    strongly_supported_tree_count: int
    strongly_supported_tree_fraction: float
    mean_support: float
    lower_95_midpoint_depth: float
    upper_95_midpoint_depth: float
    minimum_parent_depth: float
    maximum_child_depth: float
    stability_class: str


@dataclass(frozen=True, slots=True)
class GeographicMigrationTreeSetSummary:
    """Reviewer-facing summary for geographic movement events across a tree set."""

    trait: str
    taxon_column: str
    model: str
    internal_model: str
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    shared_tree_taxon_count: int
    analysis_taxon_count: int
    rooted_topology_count: int
    unrooted_topology_count: int
    event_row_count: int
    event_summary_count: int
    topology_sensitive_event_count: int
    low_support_event_count: int
    excluded_taxon_count: int
    warning_count: int


@dataclass(slots=True)
class GeographicMigrationTreeSetReport:
    """Geographic movement-event stability across a posterior or bootstrap tree set."""

    tree_set_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    internal_model: str
    burnin_fraction: float
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    shared_tree_taxa: list[str]
    analysis_taxa: list[str]
    rooted_topology_count: int
    unrooted_topology_count: int
    summary: GeographicMigrationTreeSetSummary
    tree_rows: list[GeographicMigrationTreeRow]
    event_rows: list[GeographicMigrationTreeSetEventRow]
    event_summaries: list[GeographicMigrationTreeSetEventSummaryRow]
    exclusion_rows: list[GeographicExcludedTaxonRow]
    warnings: list[str]
