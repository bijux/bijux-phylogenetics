from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class TreeTransformationSummary:
    """Before/after summary for a tree transformation."""

    transformation: str
    original_tip_count: int
    transformed_tip_count: int
    retained_taxa: list[str]
    removed_taxa: list[str]
    added_taxa: list[str]
    original_internal_node_count: int
    transformed_internal_node_count: int
    nodes_changed: list[str]
    original_total_branch_length: float
    transformed_total_branch_length: float
    branch_length_delta: float
    branch_lengths_affected: list[str]


@dataclass(slots=True)
class CladeExtractionReport:
    """Explicit record of subtree extraction for a named clade."""

    tree_path: Path
    clade_name: str
    tip_count: int
    taxa: list[str]
    retained_all_requested_descendants: bool
    missing_requested_descendants: list[str]
    unexpected_retained_taxa: list[str]
    summary: TreeTransformationSummary


@dataclass(slots=True)
class SubtreeExtractionReport:
    """Explicit record of subtree extraction by node identity or taxa."""

    tree_path: Path
    selector_kind: str
    requested_node_id: int | None
    matched_node_id: int
    requested_taxa: list[str]
    matched_node_name: str | None
    tip_count: int
    taxa: list[str]
    retained_all_requested_descendants: bool
    missing_requested_descendants: list[str]
    unexpected_retained_taxa: list[str]
    summary: TreeTransformationSummary


@dataclass(slots=True)
class TreeMrcaReport:
    """Explicit record of a tree MRCA resolution request."""

    tree_path: Path
    requested_taxa: list[str]
    unique_requested_taxa: list[str]
    duplicate_requested_taxa: list[str]
    matched_node_id: int
    matched_node_name: str | None
    matched_taxa: list[str]
    matched_extra_taxa: list[str]
    matched_tip_count: int
    is_root: bool
    rooted: bool | None


@dataclass(slots=True)
class TreeMonophylyReport:
    """Explicit record of one tree monophyly assessment."""

    tree_path: Path
    requested_taxa: list[str]
    unique_requested_taxa: list[str]
    duplicate_requested_taxa: list[str]
    missing_requested_taxa: list[str]
    present_requested_taxa: list[str]
    reroot: bool
    rooted: bool | None
    monophyletic: bool
    complementary_clade_used: bool
    matched_node_id: int | None
    matched_node_name: str | None
    matched_taxa: list[str]
    matched_extra_taxa: list[str]
    matched_tip_count: int
    is_root: bool | None


@dataclass(slots=True)
class BranchCollapseReport:
    """Explicit record of internal branches collapsed by a length threshold."""

    tree_path: Path
    threshold: float
    collapsed_clades: list[str]
    topology_preserved: bool
    summary: TreeTransformationSummary


@dataclass(slots=True)
class TreeOrderingReport:
    """Explicit record of a deterministic child-ordering transform."""

    tree_path: Path
    strategy: str
    tip_order: list[str]
    rooted_topology_preserved: bool
    unrooted_topology_preserved: bool
    summary: TreeTransformationSummary


@dataclass(slots=True)
class TreeRootingReport:
    """Explicit record of a tree rooting transform."""

    tree_path: Path
    strategy: str
    requested_taxa: list[str]
    matched_taxa: list[str]
    absent_taxa: list[str]
    ingroup_taxa: list[str]
    outgroup_monophyletic: bool | None
    outgroup_mrca_taxa: list[str]
    outgroup_mrca_extra_taxa: list[str]
    rooted_outgroup_taxa: list[str]
    rooted_ingroup_taxa: list[str]
    tip_order: list[str]
    warnings: list[str]
    midpoint_anchor_taxa: list[str]
    midpoint_path_length: float | None
    midpoint_distance_from_anchor: float | None
    midpoint_position_kind: str | None
    midpoint_anchor_side_taxa: list[str]
    midpoint_opposite_side_taxa: list[str]
    midpoint_suitable: bool | None
    summary: TreeTransformationSummary


@dataclass(slots=True)
class _TopologyComparison:
    topology_equal: bool
    same_unrooted_topology: bool


@dataclass(slots=True)
class RandomBifurcatingTreeReport:
    """Explicit record of one seeded random bifurcating tree generation."""

    algorithm: str
    seed: int
    branch_length_policy: str
    requested_taxa: list[str]
    tip_order: list[str]
    tip_count: int
    internal_node_count: int
    rooted: bool | None
    strictly_bifurcating: bool
    all_requested_taxa_present_once: bool
    missing_requested_taxa: list[str]
    duplicate_generated_taxa: list[str]
    unexpected_generated_taxa: list[str]
    validation_errors: list[str]
    tree_newick: str


@dataclass(frozen=True, slots=True)
class StepwiseAdditionCandidateScore:
    """One evaluated insertion edge for one greedy stepwise-addition step."""

    branch_id: str
    descendant_taxa: list[str]
    score: float
    candidate_tree_newick: str


@dataclass(frozen=True, slots=True)
class StepwiseAdditionTraceRow:
    """Trace of one taxon insertion during greedy stepwise addition."""

    step_index: int
    taxon: str
    inserted_taxa: list[str]
    tested_edge_rows: list[StepwiseAdditionCandidateScore]
    best_edge_id: str
    best_edge_descendant_taxa: list[str]
    best_score: float
    selected_tree_newick: str


@dataclass(slots=True)
class StepwiseAdditionTreeReport:
    """Explicit record of one greedy stepwise-addition tree build."""

    algorithm: str
    objective_name: str
    objective_direction: str
    insertion_order: list[str]
    starting_taxa: list[str]
    tip_order: list[str]
    tip_count: int
    internal_node_count: int
    rooted: bool | None
    strictly_bifurcating: bool
    all_requested_taxa_present_once: bool
    missing_requested_taxa: list[str]
    duplicate_generated_taxa: list[str]
    unexpected_generated_taxa: list[str]
    validation_errors: list[str]
    final_score: float
    trace_rows: list[StepwiseAdditionTraceRow]
    tree_newick: str
