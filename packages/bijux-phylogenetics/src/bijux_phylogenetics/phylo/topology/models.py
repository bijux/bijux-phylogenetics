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
class RootedNniNeighborRow:
    """One generated rooted NNI neighbor with its defining move metadata."""

    neighbor_index: int
    parent_node_id: str
    child_node_id: str
    sibling_node_id: str
    exchanged_child_node_id: str
    pivot_branch_id: str
    sibling_clade_id: str
    exchanged_clade_id: str
    neighbor_tree_newick: str
    neighbor_topology_fingerprint: str
    tip_order: list[str]
    validation_errors: list[str]


@dataclass(slots=True)
class RootedNniNeighborhoodReport:
    """Explicit record of rooted NNI neighbors for one rooted binary tree."""

    algorithm: str
    input_tree_path: Path | None
    input_tree_newick: str
    tip_count: int
    internal_node_count: int
    rooted: bool | None
    strictly_bifurcating: bool
    expected_neighbor_count: int
    generated_neighbor_count: int
    unique_neighbor_topology_count: int
    duplicate_neighbor_topologies: list[str]
    missing_tip_taxa: list[str]
    unexpected_tip_taxa: list[str]
    input_validation_errors: list[str]
    neighbor_rows: list[RootedNniNeighborRow]


@dataclass(slots=True)
class RootedNniMoveApplicationReport:
    """Explicit record of one rooted NNI move, its reverse, and preserved payloads."""

    algorithm: str
    input_tree_path: Path | None
    input_tree_newick: str
    input_topology_fingerprint: str
    selected_move_index: int
    available_move_count: int
    selected_parent_node_id: str
    selected_child_node_id: str
    selected_sibling_node_id: str
    selected_exchanged_child_node_id: str
    selected_pivot_branch_id: str
    selected_sibling_clade_id: str
    selected_exchanged_clade_id: str
    moved_tree_newick: str
    moved_topology_fingerprint: str
    moved_topology_changed: bool
    reverse_parent_node_id: str
    reverse_child_node_id: str
    reverse_sibling_node_id: str
    reverse_exchanged_child_node_id: str
    reverse_pivot_branch_id: str
    reverse_sibling_clade_id: str
    reverse_exchanged_clade_id: str
    reversed_tree_newick: str
    reversed_topology_fingerprint: str
    reverse_restores_original_topology: bool
    tip_count: int
    internal_node_count: int
    rooted: bool | None
    strictly_bifurcating: bool
    missing_tip_taxa: list[str]
    unexpected_tip_taxa: list[str]
    moved_validation_errors: list[str]
    reversed_validation_errors: list[str]
    node_names_preserved: bool
    node_metadata_preserved: bool
    edge_metadata_preserved: bool
    branch_lengths_preserved: bool
    total_branch_length_preserved: bool


@dataclass(frozen=True, slots=True)
class RootedSprNeighborRow:
    """One rooted SPR neighbor with a representative move and multiplicity."""

    neighbor_index: int
    representative_pruned_node_id: str
    representative_pruned_clade_id: str
    representative_pruned_descendant_taxa: list[str]
    representative_regraft_target_branch_id: str
    representative_regraft_target_descendant_taxa: list[str] | None
    supporting_move_count: int
    neighbor_tree_newick: str
    neighbor_topology_fingerprint: str
    tip_order: list[str]
    validation_errors: list[str]


@dataclass(frozen=True, slots=True)
class RootedSprEnumerationBudget:
    """Explicit rooted SPR enumeration limits over prune nodes and regraft targets."""

    max_pruned_clade_count: int | None = None
    max_regraft_target_count_per_pruned_clade: int | None = None


@dataclass(slots=True)
class RootedSprNeighborhoodReport:
    """Explicit record of rooted SPR neighbors for one binary-root tree."""

    algorithm: str
    input_tree_path: Path | None
    input_tree_newick: str
    tip_count: int
    internal_node_count: int
    rooted: bool | None
    strictly_bifurcating: bool
    max_pruned_clade_count: int | None
    max_regraft_target_count_per_pruned_clade: int | None
    skipped_pruned_clade_count: int
    skipped_regraft_target_count: int
    generated_move_candidate_count: int
    identity_move_candidate_count: int
    self_regraft_candidate_count: int
    generated_neighbor_count: int
    unique_neighbor_topology_count: int
    duplicate_move_neighbor_topologies: list[str]
    missing_tip_taxa: list[str]
    unexpected_tip_taxa: list[str]
    input_validation_errors: list[str]
    neighbor_rows: list[RootedSprNeighborRow]


@dataclass(slots=True)
class RootedSprMoveApplicationReport:
    """Explicit record of one rooted SPR move application over a binary rooted tree."""

    algorithm: str
    input_tree_path: Path | None
    input_tree_newick: str
    input_topology_fingerprint: str
    selected_move_index: int
    available_move_count: int
    max_pruned_clade_count: int | None
    max_regraft_target_count_per_pruned_clade: int | None
    selected_pruned_node_id: str
    selected_pruned_clade_id: str
    selected_pruned_descendant_taxa: list[str]
    selected_regraft_target_branch_id: str
    selected_regraft_target_descendant_taxa: list[str] | None
    moved_tree_newick: str
    moved_topology_fingerprint: str
    moved_topology_changed: bool
    tip_count: int
    internal_node_count: int
    rooted: bool | None
    strictly_bifurcating: bool
    missing_tip_taxa: list[str]
    unexpected_tip_taxa: list[str]
    moved_validation_errors: list[str]
    affected_clade_ids: list[str]
    pruned_edge_id: str
    regraft_edge_id: str


@dataclass(frozen=True, slots=True)
class RootedTbrNeighborRow:
    """One rooted TBR neighbor with a representative cut and reconnection."""

    neighbor_index: int
    representative_cut_edge_id: str
    representative_cut_descendant_taxa: list[str]
    representative_left_attachment_branch_id: str
    representative_left_attachment_descendant_taxa: list[str]
    representative_right_attachment_branch_id: str
    representative_right_attachment_descendant_taxa: list[str]
    supporting_reconnection_count: int
    neighbor_tree_newick: str
    neighbor_topology_fingerprint: str
    tip_order: list[str]
    validation_errors: list[str]


@dataclass(slots=True)
class RootedTbrNeighborhoodReport:
    """Explicit record of rooted TBR neighbors for one binary-root tree."""

    algorithm: str
    input_tree_path: Path | None
    input_tree_newick: str
    tip_count: int
    internal_node_count: int
    rooted: bool | None
    strictly_bifurcating: bool
    generated_cut_edge_count: int
    generated_reconnection_count: int
    identity_reconnection_count: int
    generated_neighbor_count: int
    unique_neighbor_topology_count: int
    duplicate_reconnection_neighbor_topologies: list[str]
    missing_tip_taxa: list[str]
    unexpected_tip_taxa: list[str]
    input_validation_errors: list[str]
    neighbor_rows: list[RootedTbrNeighborRow]


@dataclass(slots=True)
class RootedTbrMoveApplicationReport:
    """Explicit record of one rooted TBR move application over a binary rooted tree."""

    algorithm: str
    input_tree_path: Path | None
    input_tree_newick: str
    input_topology_fingerprint: str
    selected_move_index: int
    available_move_count: int
    selected_cut_parent_node_id: str
    selected_cut_child_node_id: str
    selected_cut_edge_id: str
    selected_cut_descendant_taxa: list[str]
    left_component_tip_count: int
    right_component_tip_count: int
    selected_left_attachment_branch_id: str
    selected_left_attachment_descendant_taxa: list[str]
    selected_right_attachment_branch_id: str
    selected_right_attachment_descendant_taxa: list[str]
    moved_tree_newick: str
    moved_topology_fingerprint: str
    moved_topology_changed: bool
    reverse_move_available: bool
    reverse_available_move_count: int
    reverse_cut_edge_id: str | None
    reverse_left_attachment_branch_id: str | None
    reverse_right_attachment_branch_id: str | None
    tip_count: int
    internal_node_count: int
    rooted: bool | None
    strictly_bifurcating: bool
    missing_tip_taxa: list[str]
    unexpected_tip_taxa: list[str]
    moved_validation_errors: list[str]


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
