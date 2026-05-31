from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import json
from pathlib import Path

from bijux_phylogenetics.io.newick import load_newick, loads_newick, write_newick
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_clade_id,
    rooted_topology_fingerprint,
)
from bijux_phylogenetics.phylo.topology.models import (
    RootedTbrMoveApplicationReport,
    RootedTbrNeighborhoodReport,
    RootedTbrNeighborRow,
)

from .affected_subtrees import summarize_affected_subtrees
from .neighborhood_summary import (
    summarize_topology_neighborhood,
    write_topology_neighborhood_summary_table,
)
from .tree import PhyloTree, TreeNode, descendant_taxa

ROOTED_TBR_INTERFACE_BRANCH_ID = "interface"


@dataclass(frozen=True, slots=True)
class _RootedTbrCutEdgeCandidate:
    parent_node_id: str
    child_node_id: str
    cut_edge_id: str
    cut_descendant_taxa: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _RootedTbrAttachmentBranchCandidate:
    left_node_id: str
    right_node_id: str
    branch_id: str
    descendant_taxa: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RootedTbrMoveCandidate:
    """One deterministic rooted TBR reconnection over a binary rooted tree."""

    cut_parent_node_id: str
    cut_child_node_id: str
    cut_edge_id: str
    cut_descendant_taxa: tuple[str, ...]
    left_attachment_left_node_id: str
    left_attachment_right_node_id: str
    left_attachment_branch_id: str
    left_attachment_descendant_taxa: tuple[str, ...]
    right_attachment_left_node_id: str
    right_attachment_right_node_id: str
    right_attachment_branch_id: str
    right_attachment_descendant_taxa: tuple[str, ...]


@dataclass(slots=True)
class _RootedTbrGraph:
    node_labels: dict[str, str | None]
    adjacency: dict[str, dict[str, float | None]]
    interface_node_ids: set[str]

    def copy(self) -> _RootedTbrGraph:
        return _RootedTbrGraph(
            node_labels=dict(self.node_labels),
            adjacency={
                node_id: dict(neighbors)
                for node_id, neighbors in self.adjacency.items()
            },
            interface_node_ids=set(self.interface_node_ids),
        )

    def add_node(self, node_id: str, label: str | None) -> None:
        self.node_labels[node_id] = label
        self.adjacency.setdefault(node_id, {})

    def add_edge(
        self, left_node_id: str, right_node_id: str, length: float | None
    ) -> None:
        self.adjacency.setdefault(left_node_id, {})[right_node_id] = length
        self.adjacency.setdefault(right_node_id, {})[left_node_id] = length

    def remove_edge(self, left_node_id: str, right_node_id: str) -> float | None:
        length = self.adjacency[left_node_id].pop(right_node_id)
        self.adjacency[right_node_id].pop(left_node_id)
        return length

    def remove_node(self, node_id: str) -> None:
        for neighbor_id in list(self.adjacency.get(node_id, {})):
            self.adjacency[neighbor_id].pop(node_id, None)
        self.adjacency.pop(node_id, None)
        self.node_labels.pop(node_id, None)
        self.interface_node_ids.discard(node_id)


def rooted_tbr_node_sort_key(node: TreeNode) -> tuple[int, tuple[str, ...]]:
    """Sort rooted TBR cut edges deterministically by descendant taxa."""
    descendants = tuple(descendant_taxa(node))
    return (len(descendants), descendants)


def rooted_tbr_clade_id(node: TreeNode) -> str:
    """Render one rooted TBR edge endpoint as a stable descendant-clade identifier."""
    return canonical_clade_id(frozenset(descendant_taxa(node)))


def require_rooted_tbr_node_id(node: TreeNode) -> str:
    """Require one refreshed stable node identifier on a rooted TBR tree."""
    if node.node_id is None:
        raise AssertionError(
            "rooted TBR enumeration requires refreshed node identities"
        )
    return node.node_id


def validate_rooted_tbr_tree(tree: PhyloTree) -> None:
    """Require one structurally valid strictly bifurcating binary-root tree."""
    validation_errors = tree.validation_errors()
    if validation_errors:
        raise ValueError(
            "rooted TBR enumeration requires a structurally valid tree: "
            + "; ".join(validation_errors)
        )
    if len(tree.root.children) != 2:
        raise ValueError("rooted TBR enumeration requires a binary root")
    invalid_internal_nodes = [
        node.node_id
        for node in tree.iter_internal_nodes(order="preorder")
        if len(node.children) != 2
    ]
    if invalid_internal_nodes:
        raise ValueError("rooted TBR enumeration requires a strictly bifurcating tree")


def enumerate_rooted_tbr_neighbors(
    tree: PhyloTree | Path,
) -> RootedTbrNeighborhoodReport:
    """Enumerate rooted TBR neighbors exactly once by unique non-identity topology."""
    resolved_tree, input_tree_path = _resolve_rooted_tbr_tree(tree)
    validate_rooted_tbr_tree(resolved_tree)
    input_tip_taxa = set(resolved_tree.tip_names)
    input_topology_fingerprint = rooted_topology_fingerprint(resolved_tree)
    generated_cut_edge_count = 0
    generated_reconnection_count = 0
    identity_reconnection_count = 0
    duplicate_reconnection_counts: dict[str, int] = {}
    neighbor_row_by_fingerprint: dict[str, RootedTbrNeighborRow] = {}
    for cut_candidate in iter_rooted_tbr_cut_edge_candidates(resolved_tree):
        generated_cut_edge_count += 1
        left_component, right_component = bisect_rooted_tbr_edge(
            resolved_tree,
            cut_candidate.parent_node_id,
            cut_candidate.child_node_id,
        )
        for left_attachment in iter_rooted_tbr_attachment_branch_candidates(
            left_component
        ):
            for right_attachment in iter_rooted_tbr_attachment_branch_candidates(
                right_component
            ):
                generated_reconnection_count += 1
                neighbor_tree = reconnect_rooted_tbr_components(
                    left_component,
                    right_component,
                    left_attachment=left_attachment,
                    right_attachment=right_attachment,
                    reconnection_index=generated_reconnection_count,
                )
                topology_fingerprint = rooted_topology_fingerprint(neighbor_tree)
                if topology_fingerprint == input_topology_fingerprint:
                    identity_reconnection_count += 1
                    continue
                duplicate_reconnection_counts[topology_fingerprint] = (
                    duplicate_reconnection_counts.get(topology_fingerprint, 0) + 1
                )
                if topology_fingerprint in neighbor_row_by_fingerprint:
                    representative_row = neighbor_row_by_fingerprint[
                        topology_fingerprint
                    ]
                    neighbor_row_by_fingerprint[topology_fingerprint] = (
                        RootedTbrNeighborRow(
                            neighbor_index=representative_row.neighbor_index,
                            representative_cut_edge_id=representative_row.representative_cut_edge_id,
                            representative_cut_descendant_taxa=representative_row.representative_cut_descendant_taxa,
                            representative_left_attachment_branch_id=representative_row.representative_left_attachment_branch_id,
                            representative_left_attachment_descendant_taxa=representative_row.representative_left_attachment_descendant_taxa,
                            representative_right_attachment_branch_id=representative_row.representative_right_attachment_branch_id,
                            representative_right_attachment_descendant_taxa=representative_row.representative_right_attachment_descendant_taxa,
                            supporting_reconnection_count=(
                                representative_row.supporting_reconnection_count + 1
                            ),
                            neighbor_tree_newick=representative_row.neighbor_tree_newick,
                            neighbor_topology_fingerprint=representative_row.neighbor_topology_fingerprint,
                            tip_order=representative_row.tip_order,
                            validation_errors=representative_row.validation_errors,
                        )
                    )
                    continue
                neighbor_row_by_fingerprint[topology_fingerprint] = (
                    RootedTbrNeighborRow(
                        neighbor_index=len(neighbor_row_by_fingerprint) + 1,
                        representative_cut_edge_id=cut_candidate.cut_edge_id,
                        representative_cut_descendant_taxa=list(
                            cut_candidate.cut_descendant_taxa
                        ),
                        representative_left_attachment_branch_id=left_attachment.branch_id,
                        representative_left_attachment_descendant_taxa=list(
                            left_attachment.descendant_taxa
                        ),
                        representative_right_attachment_branch_id=right_attachment.branch_id,
                        representative_right_attachment_descendant_taxa=list(
                            right_attachment.descendant_taxa
                        ),
                        supporting_reconnection_count=1,
                        neighbor_tree_newick=neighbor_tree.to_newick(),
                        neighbor_topology_fingerprint=topology_fingerprint,
                        tip_order=neighbor_tree.tip_names,
                        validation_errors=neighbor_tree.validation_errors(),
                    )
                )
    neighbor_rows = list(neighbor_row_by_fingerprint.values())
    report = RootedTbrNeighborhoodReport(
        algorithm="rooted-tbr-neighbor-enumeration",
        input_tree_path=input_tree_path,
        input_tree_newick=resolved_tree.to_newick(),
        tip_count=resolved_tree.tip_count,
        internal_node_count=resolved_tree.internal_node_count,
        rooted=resolved_tree.rooted,
        strictly_bifurcating=True,
        generated_cut_edge_count=generated_cut_edge_count,
        generated_reconnection_count=generated_reconnection_count,
        identity_reconnection_count=identity_reconnection_count,
        generated_neighbor_count=len(neighbor_rows),
        unique_neighbor_topology_count=len(neighbor_row_by_fingerprint),
        duplicate_reconnection_neighbor_topologies=sorted(
            fingerprint
            for fingerprint, count in duplicate_reconnection_counts.items()
            if count > 1
        ),
        missing_tip_taxa=sorted(
            taxon
            for taxon in input_tip_taxa
            if any(taxon not in row.tip_order for row in neighbor_rows)
        ),
        unexpected_tip_taxa=sorted(
            {
                taxon
                for row in neighbor_rows
                for taxon in row.tip_order
                if taxon not in input_tip_taxa
            }
        ),
        input_validation_errors=[],
        neighbor_rows=neighbor_rows,
    )
    _validate_rooted_tbr_neighbor_report(report)
    return report


def iter_rooted_tbr_move_candidates(tree: PhyloTree):
    """Yield deterministic legal rooted TBR reconnection moves for one rooted tree."""
    input_topology_fingerprint = rooted_topology_fingerprint(tree)
    for cut_candidate in iter_rooted_tbr_cut_edge_candidates(tree):
        left_component, right_component = bisect_rooted_tbr_edge(
            tree,
            cut_candidate.parent_node_id,
            cut_candidate.child_node_id,
        )
        for left_attachment in iter_rooted_tbr_attachment_branch_candidates(
            left_component
        ):
            for right_attachment in iter_rooted_tbr_attachment_branch_candidates(
                right_component
            ):
                candidate = RootedTbrMoveCandidate(
                    cut_parent_node_id=cut_candidate.parent_node_id,
                    cut_child_node_id=cut_candidate.child_node_id,
                    cut_edge_id=cut_candidate.cut_edge_id,
                    cut_descendant_taxa=cut_candidate.cut_descendant_taxa,
                    left_attachment_left_node_id=left_attachment.left_node_id,
                    left_attachment_right_node_id=left_attachment.right_node_id,
                    left_attachment_branch_id=left_attachment.branch_id,
                    left_attachment_descendant_taxa=left_attachment.descendant_taxa,
                    right_attachment_left_node_id=right_attachment.left_node_id,
                    right_attachment_right_node_id=right_attachment.right_node_id,
                    right_attachment_branch_id=right_attachment.branch_id,
                    right_attachment_descendant_taxa=right_attachment.descendant_taxa,
                )
                moved_tree = apply_rooted_tbr_move(tree, candidate)
                if (
                    rooted_topology_fingerprint(moved_tree)
                    == input_topology_fingerprint
                ):
                    continue
                yield candidate


def apply_rooted_tbr_move(
    tree: PhyloTree,
    candidate: RootedTbrMoveCandidate,
) -> PhyloTree:
    """Return one copied rooted tree with the selected TBR reconnection applied."""
    left_component, right_component = bisect_rooted_tbr_edge(
        tree,
        candidate.cut_parent_node_id,
        candidate.cut_child_node_id,
    )
    return reconnect_rooted_tbr_components(
        left_component,
        right_component,
        left_attachment=_RootedTbrAttachmentBranchCandidate(
            left_node_id=candidate.left_attachment_left_node_id,
            right_node_id=candidate.left_attachment_right_node_id,
            branch_id=candidate.left_attachment_branch_id,
            descendant_taxa=candidate.left_attachment_descendant_taxa,
        ),
        right_attachment=_RootedTbrAttachmentBranchCandidate(
            left_node_id=candidate.right_attachment_left_node_id,
            right_node_id=candidate.right_attachment_right_node_id,
            branch_id=candidate.right_attachment_branch_id,
            descendant_taxa=candidate.right_attachment_descendant_taxa,
        ),
        reconnection_index=1,
    )


def resolve_rooted_tbr_move_candidate(
    tree: PhyloTree | Path,
    move_index: int,
) -> tuple[RootedTbrMoveCandidate, int]:
    """Resolve one 1-indexed rooted TBR move candidate from a validated tree."""
    resolved_tree, _input_tree_path = _resolve_rooted_tbr_tree(tree)
    validate_rooted_tbr_tree(resolved_tree)
    candidates = list(iter_rooted_tbr_move_candidates(resolved_tree))
    if move_index < 1 or move_index > len(candidates):
        raise ValueError(
            f"rooted TBR move index must be between 1 and {len(candidates)}"
        )
    return candidates[move_index - 1], len(candidates)


def summarize_rooted_tbr_move_application(
    tree: PhyloTree | Path,
    move_index: int,
) -> RootedTbrMoveApplicationReport:
    """Apply one rooted TBR move and report the resulting valid transformed tree."""
    resolved_tree, input_tree_path = _resolve_rooted_tbr_tree(tree)
    validate_rooted_tbr_tree(resolved_tree)
    selected_candidate, available_move_count = resolve_rooted_tbr_move_candidate(
        resolved_tree,
        move_index,
    )
    moved_tree = apply_rooted_tbr_move(resolved_tree, selected_candidate)
    input_tip_taxa = set(resolved_tree.tip_names)
    moved_tip_taxa = set(moved_tree.tip_names)
    input_topology_fingerprint = rooted_topology_fingerprint(resolved_tree)
    moved_topology_fingerprint = rooted_topology_fingerprint(moved_tree)
    reverse_candidates = [
        candidate
        for candidate in iter_rooted_tbr_move_candidates(moved_tree)
        if rooted_topology_fingerprint(apply_rooted_tbr_move(moved_tree, candidate))
        == input_topology_fingerprint
    ]
    reverse_candidate = reverse_candidates[0] if reverse_candidates else None
    report = RootedTbrMoveApplicationReport(
        algorithm="rooted-tbr-move-application",
        input_tree_path=input_tree_path,
        input_tree_newick=resolved_tree.to_newick(),
        input_topology_fingerprint=input_topology_fingerprint,
        selected_move_index=move_index,
        available_move_count=available_move_count,
        selected_cut_parent_node_id=selected_candidate.cut_parent_node_id,
        selected_cut_child_node_id=selected_candidate.cut_child_node_id,
        selected_cut_edge_id=selected_candidate.cut_edge_id,
        selected_cut_descendant_taxa=list(selected_candidate.cut_descendant_taxa),
        left_component_tip_count=len(selected_candidate.cut_descendant_taxa),
        right_component_tip_count=(
            resolved_tree.tip_count - len(selected_candidate.cut_descendant_taxa)
        ),
        selected_left_attachment_branch_id=selected_candidate.left_attachment_branch_id,
        selected_left_attachment_descendant_taxa=list(
            selected_candidate.left_attachment_descendant_taxa
        ),
        selected_right_attachment_branch_id=(
            selected_candidate.right_attachment_branch_id
        ),
        selected_right_attachment_descendant_taxa=list(
            selected_candidate.right_attachment_descendant_taxa
        ),
        moved_tree_newick=moved_tree.to_newick(),
        moved_topology_fingerprint=moved_topology_fingerprint,
        moved_topology_changed=moved_topology_fingerprint != input_topology_fingerprint,
        reverse_move_available=bool(reverse_candidates),
        reverse_available_move_count=len(reverse_candidates),
        reverse_cut_edge_id=(
            None if reverse_candidate is None else reverse_candidate.cut_edge_id
        ),
        reverse_left_attachment_branch_id=(
            None
            if reverse_candidate is None
            else reverse_candidate.left_attachment_branch_id
        ),
        reverse_right_attachment_branch_id=(
            None
            if reverse_candidate is None
            else reverse_candidate.right_attachment_branch_id
        ),
        tip_count=resolved_tree.tip_count,
        internal_node_count=resolved_tree.internal_node_count,
        rooted=resolved_tree.rooted,
        strictly_bifurcating=True,
        missing_tip_taxa=sorted(input_tip_taxa - moved_tip_taxa),
        unexpected_tip_taxa=sorted(moved_tip_taxa - input_tip_taxa),
        moved_validation_errors=moved_tree.validation_errors(),
        affected_subtree_report=summarize_affected_subtrees(
            resolved_tree,
            moved_tree,
        ),
    )
    _validate_rooted_tbr_move_application_report(report)
    return report


def iter_rooted_tbr_cut_edge_candidates(tree: PhyloTree):
    """Yield deterministic internal cut edges for rooted TBR enumeration."""
    sorted_internal_nodes = sorted(
        (
            node
            for node in tree.iter_internal_nodes(order="preorder")
            if node is not tree.root
        ),
        key=rooted_tbr_node_sort_key,
    )
    for node in sorted_internal_nodes:
        parent = node.parent
        if parent is None:
            raise AssertionError("rooted TBR cut nodes must not be the root")
        yield _RootedTbrCutEdgeCandidate(
            parent_node_id=require_rooted_tbr_node_id(parent),
            child_node_id=require_rooted_tbr_node_id(node),
            cut_edge_id=rooted_tbr_clade_id(node),
            cut_descendant_taxa=tuple(sorted(node.descendant_taxa)),
        )


def bisect_rooted_tbr_edge(
    tree: PhyloTree,
    parent_node_id: str,
    child_node_id: str,
) -> tuple[_RootedTbrGraph, _RootedTbrGraph]:
    """Cut one internal edge and return the two augmented components it induces."""
    graph = _build_rooted_tbr_graph(tree)
    cut_edge_length = graph.remove_edge(parent_node_id, child_node_id)
    left_component_node_ids = _graph_component_node_ids(graph, child_node_id)
    right_component_node_ids = _graph_component_node_ids(graph, parent_node_id)
    left_component = _graph_induced_subgraph(graph, left_component_node_ids)
    right_component = _graph_induced_subgraph(graph, right_component_node_ids)
    left_interface_node_id = f"{child_node_id}::tbr-interface"
    right_interface_node_id = f"{parent_node_id}::tbr-interface"
    left_component.add_node(left_interface_node_id, None)
    right_component.add_node(right_interface_node_id, None)
    left_component.interface_node_ids.add(left_interface_node_id)
    right_component.interface_node_ids.add(right_interface_node_id)
    left_component.add_edge(left_interface_node_id, child_node_id, cut_edge_length)
    right_component.add_edge(right_interface_node_id, parent_node_id, cut_edge_length)
    return left_component, right_component


def iter_rooted_tbr_attachment_branch_candidates(component: _RootedTbrGraph):
    """Yield all legal TBR attachment branches for one cut component."""
    component_tip_taxa = tuple(sorted(_graph_tip_taxa(component)))
    candidates = [
        _RootedTbrAttachmentBranchCandidate(
            left_node_id=left_node_id,
            right_node_id=right_node_id,
            branch_id=branch_id,
            descendant_taxa=descendant_taxa,
        )
        for left_node_id, right_node_id, branch_id, descendant_taxa in _graph_edges_with_labels(
            component,
            component_tip_taxa,
        )
    ]
    return iter(
        sorted(
            candidates,
            key=lambda candidate: (
                candidate.branch_id != ROOTED_TBR_INTERFACE_BRANCH_ID,
                len(candidate.descendant_taxa),
                candidate.descendant_taxa,
                tuple(sorted((candidate.left_node_id, candidate.right_node_id))),
            ),
        )
    )


def reconnect_rooted_tbr_components(
    left_component: _RootedTbrGraph,
    right_component: _RootedTbrGraph,
    *,
    left_attachment: _RootedTbrAttachmentBranchCandidate,
    right_attachment: _RootedTbrAttachmentBranchCandidate,
    reconnection_index: int,
) -> PhyloTree:
    """Reconnect two cut components across selected attachment edges."""
    merged_graph = _merge_rooted_tbr_graphs(left_component, right_component)
    left_bridge_node_id = _subdivide_rooted_tbr_edge(
        merged_graph,
        left_attachment.left_node_id,
        left_attachment.right_node_id,
        subdivision_node_id=f"left-bridge-{reconnection_index}",
    )
    right_bridge_node_id = _subdivide_rooted_tbr_edge(
        merged_graph,
        right_attachment.left_node_id,
        right_attachment.right_node_id,
        subdivision_node_id=f"right-bridge-{reconnection_index}",
    )
    merged_graph.add_node(left_bridge_node_id, None)
    merged_graph.add_node(right_bridge_node_id, None)
    merged_graph.add_edge(left_bridge_node_id, right_bridge_node_id, None)
    for interface_node_id in list(merged_graph.interface_node_ids):
        merged_graph.remove_node(interface_node_id)
    _prune_rooted_tbr_graph_dead_ends(
        merged_graph,
        preserved_node_ids={left_bridge_node_id, right_bridge_node_id},
    )
    return _root_rooted_tbr_graph_on_edge(
        merged_graph,
        left_bridge_node_id,
        right_bridge_node_id,
    )


def _resolve_rooted_tbr_tree(
    tree: PhyloTree | Path,
) -> tuple[PhyloTree, Path | None]:
    resolved_tree_path = tree if isinstance(tree, Path) else None
    resolved_tree = load_newick(tree) if isinstance(tree, Path) else tree.copy()
    return resolved_tree.refresh(), resolved_tree_path


def _validate_rooted_tbr_neighbor_report(
    report: RootedTbrNeighborhoodReport,
) -> None:
    if report.generated_neighbor_count != report.unique_neighbor_topology_count:
        raise ValueError(
            "rooted TBR enumeration did not collapse neighbors to unique topologies"
        )
    if report.missing_tip_taxa or report.unexpected_tip_taxa:
        raise ValueError("rooted TBR enumeration changed the input taxon set")
    invalid_neighbor_rows = [
        row.neighbor_index for row in report.neighbor_rows if row.validation_errors
    ]
    if invalid_neighbor_rows:
        raise ValueError("rooted TBR enumeration generated invalid neighbor trees")


def _validate_rooted_tbr_move_application_report(
    report: RootedTbrMoveApplicationReport,
) -> None:
    if not report.moved_topology_changed:
        raise ValueError("rooted TBR move application did not change the topology")
    if report.missing_tip_taxa or report.unexpected_tip_taxa:
        raise ValueError("rooted TBR move application changed the input taxon set")
    if report.moved_validation_errors:
        raise ValueError("rooted TBR move application generated an invalid tree")


def write_rooted_tbr_neighbor_table(
    path: Path,
    report: RootedTbrNeighborhoodReport,
) -> Path:
    """Write one row per unique rooted TBR neighbor topology."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "neighbor_index",
                "representative_cut_edge_id",
                "representative_cut_descendant_taxa",
                "representative_left_attachment_branch_id",
                "representative_left_attachment_descendant_taxa",
                "representative_right_attachment_branch_id",
                "representative_right_attachment_descendant_taxa",
                "supporting_reconnection_count",
                "neighbor_topology_fingerprint",
                "tip_order",
                "validation_errors",
                "neighbor_tree_newick",
            ]
        )
    ]
    for row in report.neighbor_rows:
        lines.append(
            "\t".join(
                [
                    str(row.neighbor_index),
                    row.representative_cut_edge_id,
                    ",".join(row.representative_cut_descendant_taxa),
                    row.representative_left_attachment_branch_id,
                    ",".join(row.representative_left_attachment_descendant_taxa),
                    row.representative_right_attachment_branch_id,
                    ",".join(row.representative_right_attachment_descendant_taxa),
                    str(row.supporting_reconnection_count),
                    row.neighbor_topology_fingerprint,
                    ",".join(row.tip_order),
                    ",".join(row.validation_errors),
                    row.neighbor_tree_newick,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_rooted_tbr_run_json(
    path: Path,
    report: RootedTbrNeighborhoodReport,
) -> Path:
    """Write one machine-readable rooted TBR neighborhood payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "algorithm": report.algorithm,
        "input_tree_path": (
            None if report.input_tree_path is None else str(report.input_tree_path)
        ),
        "input_tree_newick": report.input_tree_newick,
        "tip_count": report.tip_count,
        "internal_node_count": report.internal_node_count,
        "rooted": report.rooted,
        "strictly_bifurcating": report.strictly_bifurcating,
        "generated_cut_edge_count": report.generated_cut_edge_count,
        "generated_reconnection_count": report.generated_reconnection_count,
        "identity_reconnection_count": report.identity_reconnection_count,
        "generated_neighbor_count": report.generated_neighbor_count,
        "unique_neighbor_topology_count": report.unique_neighbor_topology_count,
        "duplicate_reconnection_neighbor_topologies": (
            report.duplicate_reconnection_neighbor_topologies
        ),
        "missing_tip_taxa": report.missing_tip_taxa,
        "unexpected_tip_taxa": report.unexpected_tip_taxa,
        "input_validation_errors": report.input_validation_errors,
        "neighbor_rows": [
            {
                "neighbor_index": row.neighbor_index,
                "representative_cut_edge_id": row.representative_cut_edge_id,
                "representative_cut_descendant_taxa": (
                    row.representative_cut_descendant_taxa
                ),
                "representative_left_attachment_branch_id": (
                    row.representative_left_attachment_branch_id
                ),
                "representative_left_attachment_descendant_taxa": (
                    row.representative_left_attachment_descendant_taxa
                ),
                "representative_right_attachment_branch_id": (
                    row.representative_right_attachment_branch_id
                ),
                "representative_right_attachment_descendant_taxa": (
                    row.representative_right_attachment_descendant_taxa
                ),
                "supporting_reconnection_count": row.supporting_reconnection_count,
                "neighbor_tree_newick": row.neighbor_tree_newick,
                "neighbor_topology_fingerprint": row.neighbor_topology_fingerprint,
                "tip_order": row.tip_order,
                "validation_errors": row.validation_errors,
            }
            for row in report.neighbor_rows
        ],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def write_rooted_tbr_artifacts(
    out_dir: Path,
    report: RootedTbrNeighborhoodReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one rooted TBR enumeration run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    input_tree_path = write_newick(
        out_dir / "input_tree.nwk",
        loads_newick(report.input_tree_newick),
    )
    neighbors_path = write_rooted_tbr_neighbor_table(out_dir / "neighbors.tsv", report)
    summary_path = write_topology_neighborhood_summary_table(
        out_dir / "summary.tsv",
        [summarize_topology_neighborhood(report)],
    )
    run_json_path = write_rooted_tbr_run_json(out_dir / "run.json", report)
    return {
        "input_tree_path": input_tree_path,
        "neighbors_path": neighbors_path,
        "summary_path": summary_path,
        "run_json_path": run_json_path,
    }


def write_rooted_tbr_move_run_json(
    path: Path,
    report: RootedTbrMoveApplicationReport,
) -> Path:
    """Write one machine-readable rooted TBR move-application payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "algorithm": report.algorithm,
        "input_tree_path": (
            None if report.input_tree_path is None else str(report.input_tree_path)
        ),
        "input_tree_newick": report.input_tree_newick,
        "input_topology_fingerprint": report.input_topology_fingerprint,
        "selected_move_index": report.selected_move_index,
        "available_move_count": report.available_move_count,
        "selected_cut_parent_node_id": report.selected_cut_parent_node_id,
        "selected_cut_child_node_id": report.selected_cut_child_node_id,
        "selected_cut_edge_id": report.selected_cut_edge_id,
        "selected_cut_descendant_taxa": report.selected_cut_descendant_taxa,
        "left_component_tip_count": report.left_component_tip_count,
        "right_component_tip_count": report.right_component_tip_count,
        "selected_left_attachment_branch_id": (
            report.selected_left_attachment_branch_id
        ),
        "selected_left_attachment_descendant_taxa": (
            report.selected_left_attachment_descendant_taxa
        ),
        "selected_right_attachment_branch_id": (
            report.selected_right_attachment_branch_id
        ),
        "selected_right_attachment_descendant_taxa": (
            report.selected_right_attachment_descendant_taxa
        ),
        "moved_tree_newick": report.moved_tree_newick,
        "moved_topology_fingerprint": report.moved_topology_fingerprint,
        "moved_topology_changed": report.moved_topology_changed,
        "reverse_move_available": report.reverse_move_available,
        "reverse_available_move_count": report.reverse_available_move_count,
        "reverse_cut_edge_id": report.reverse_cut_edge_id,
        "reverse_left_attachment_branch_id": report.reverse_left_attachment_branch_id,
        "reverse_right_attachment_branch_id": (
            report.reverse_right_attachment_branch_id
        ),
        "tip_count": report.tip_count,
        "internal_node_count": report.internal_node_count,
        "rooted": report.rooted,
        "strictly_bifurcating": report.strictly_bifurcating,
        "missing_tip_taxa": report.missing_tip_taxa,
        "unexpected_tip_taxa": report.unexpected_tip_taxa,
        "moved_validation_errors": report.moved_validation_errors,
        "affected_subtrees": {
            "original_branch_clade_ids": (
                report.affected_subtree_report.original_branch_clade_ids
            ),
            "moved_branch_clade_ids": (
                report.affected_subtree_report.moved_branch_clade_ids
            ),
            "retired_branch_clade_ids": (
                report.affected_subtree_report.retired_branch_clade_ids
            ),
            "introduced_branch_clade_ids": (
                report.affected_subtree_report.introduced_branch_clade_ids
            ),
            "affected_branch_clade_ids": (
                report.affected_subtree_report.affected_branch_clade_ids
            ),
            "unaffected_branch_clade_ids": (
                report.affected_subtree_report.unaffected_branch_clade_ids
            ),
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_rooted_tbr_move_artifacts(
    out_dir: Path,
    report: RootedTbrMoveApplicationReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one rooted TBR move application."""
    out_dir.mkdir(parents=True, exist_ok=True)
    input_tree_path = write_newick(
        out_dir / "input_tree.nwk",
        loads_newick(report.input_tree_newick),
    )
    moved_tree_path = write_newick(
        out_dir / "moved_tree.nwk",
        loads_newick(report.moved_tree_newick),
    )
    run_json_path = write_rooted_tbr_move_run_json(out_dir / "run.json", report)
    return {
        "input_tree_path": input_tree_path,
        "moved_tree_path": moved_tree_path,
        "run_json_path": run_json_path,
    }


def _build_rooted_tbr_graph(tree: PhyloTree) -> _RootedTbrGraph:
    working_tree = tree.copy().refresh()
    node_labels: dict[str, str | None] = {}
    adjacency: dict[str, dict[str, float | None]] = {}
    for node in working_tree.iter_nodes(order="preorder"):
        node_id = require_rooted_tbr_node_id(node)
        node_labels[node_id] = node.name
        adjacency.setdefault(node_id, {})
    for parent, child in working_tree.iter_edges():
        parent_node_id = require_rooted_tbr_node_id(parent)
        child_node_id = require_rooted_tbr_node_id(child)
        adjacency[parent_node_id][child_node_id] = child.branch_length
        adjacency[child_node_id][parent_node_id] = child.branch_length
    return _RootedTbrGraph(
        node_labels=node_labels,
        adjacency=adjacency,
        interface_node_ids=set(),
    )


def _graph_component_node_ids(
    graph: _RootedTbrGraph,
    start_node_id: str,
) -> set[str]:
    visited_node_ids: set[str] = set()
    pending_node_ids = deque([start_node_id])
    while pending_node_ids:
        node_id = pending_node_ids.popleft()
        if node_id in visited_node_ids:
            continue
        visited_node_ids.add(node_id)
        pending_node_ids.extend(
            neighbor_id
            for neighbor_id in graph.adjacency[node_id]
            if neighbor_id not in visited_node_ids
        )
    return visited_node_ids


def _graph_induced_subgraph(
    graph: _RootedTbrGraph,
    included_node_ids: set[str],
) -> _RootedTbrGraph:
    subgraph = _RootedTbrGraph(node_labels={}, adjacency={}, interface_node_ids=set())
    for node_id in included_node_ids:
        subgraph.add_node(node_id, graph.node_labels[node_id])
    for node_id in included_node_ids:
        for neighbor_id, length in graph.adjacency[node_id].items():
            if neighbor_id not in included_node_ids or neighbor_id < node_id:
                continue
            subgraph.add_edge(node_id, neighbor_id, length)
    return subgraph


def _graph_tip_taxa(graph: _RootedTbrGraph) -> list[str]:
    return sorted(
        label
        for node_id, label in graph.node_labels.items()
        if label is not None and node_id not in graph.interface_node_ids
    )


def _graph_edges_with_labels(
    graph: _RootedTbrGraph,
    component_tip_taxa: tuple[str, ...],
):
    seen_edges: set[tuple[str, str]] = set()
    for left_node_id, neighbors in graph.adjacency.items():
        for right_node_id in neighbors:
            edge_key = tuple(sorted((left_node_id, right_node_id)))
            if edge_key in seen_edges:
                continue
            seen_edges.add(edge_key)
            if (
                left_node_id in graph.interface_node_ids
                or right_node_id in graph.interface_node_ids
            ):
                yield (
                    left_node_id,
                    right_node_id,
                    ROOTED_TBR_INTERFACE_BRANCH_ID,
                    component_tip_taxa,
                )
                continue
            left_taxa = _graph_side_tip_taxa(
                graph, left_node_id, blocked_node_id=right_node_id
            )
            right_taxa = _graph_side_tip_taxa(
                graph, right_node_id, blocked_node_id=left_node_id
            )
            non_empty_taxa_sides = [
                tuple(taxa) for taxa in (left_taxa, right_taxa) if taxa
            ]
            selected_taxa = min(
                non_empty_taxa_sides if non_empty_taxa_sides else [()],
                key=lambda taxa: (len(taxa), taxa),
            )
            yield (
                left_node_id,
                right_node_id,
                canonical_clade_id(frozenset(selected_taxa)),
                selected_taxa,
            )


def _graph_side_tip_taxa(
    graph: _RootedTbrGraph,
    start_node_id: str,
    *,
    blocked_node_id: str,
) -> list[str]:
    visited_node_ids: set[str] = set()
    pending_node_ids = deque([start_node_id])
    collected_taxa: list[str] = []
    while pending_node_ids:
        node_id = pending_node_ids.popleft()
        if node_id in visited_node_ids:
            continue
        visited_node_ids.add(node_id)
        node_label = graph.node_labels[node_id]
        if node_label is not None and node_id not in graph.interface_node_ids:
            collected_taxa.append(node_label)
        pending_node_ids.extend(
            neighbor_id
            for neighbor_id in graph.adjacency[node_id]
            if neighbor_id != blocked_node_id and neighbor_id not in visited_node_ids
        )
    return sorted(collected_taxa)


def _merge_rooted_tbr_graphs(
    left_graph: _RootedTbrGraph,
    right_graph: _RootedTbrGraph,
) -> _RootedTbrGraph:
    merged_graph = left_graph.copy()
    for node_id, label in right_graph.node_labels.items():
        merged_graph.add_node(node_id, label)
    for node_id, neighbors in right_graph.adjacency.items():
        for neighbor_id, length in neighbors.items():
            if neighbor_id < node_id:
                continue
            merged_graph.add_edge(node_id, neighbor_id, length)
    merged_graph.interface_node_ids.update(right_graph.interface_node_ids)
    return merged_graph


def _subdivide_rooted_tbr_edge(
    graph: _RootedTbrGraph,
    left_node_id: str,
    right_node_id: str,
    *,
    subdivision_node_id: str,
) -> str:
    edge_length = graph.remove_edge(left_node_id, right_node_id)
    graph.add_node(subdivision_node_id, None)
    if edge_length is None:
        left_length = None
        right_length = None
    else:
        left_length = edge_length / 2.0
        right_length = edge_length / 2.0
    graph.add_edge(left_node_id, subdivision_node_id, left_length)
    graph.add_edge(subdivision_node_id, right_node_id, right_length)
    return subdivision_node_id


def _root_rooted_tbr_graph_on_edge(
    graph: _RootedTbrGraph,
    left_bridge_node_id: str,
    right_bridge_node_id: str,
) -> PhyloTree:
    left_subtree = _clone_rooted_tbr_graph_component(
        graph,
        left_bridge_node_id,
        from_neighbor_id=right_bridge_node_id,
        incoming_length=None,
    )
    right_subtree = _clone_rooted_tbr_graph_component(
        graph,
        right_bridge_node_id,
        from_neighbor_id=left_bridge_node_id,
        incoming_length=None,
    )
    return PhyloTree(
        root=TreeNode(children=[left_subtree, right_subtree]),
        rooted=True,
    ).refresh()


def _prune_rooted_tbr_graph_dead_ends(
    graph: _RootedTbrGraph,
    *,
    preserved_node_ids: set[str],
) -> None:
    changed = True
    while changed:
        changed = False
        for node_id in list(graph.node_labels):
            if node_id in preserved_node_ids:
                continue
            if graph.node_labels[node_id] is not None:
                continue
            if len(graph.adjacency[node_id]) != 1:
                continue
            graph.remove_node(node_id)
            changed = True


def _clone_rooted_tbr_graph_component(
    graph: _RootedTbrGraph,
    node_id: str,
    *,
    from_neighbor_id: str,
    incoming_length: float | None,
) -> TreeNode:
    neighbor_ids = [
        neighbor_id
        for neighbor_id in graph.adjacency[node_id]
        if neighbor_id != from_neighbor_id
    ]
    if neighbor_ids and len(neighbor_ids) == 1 and graph.node_labels[node_id] is None:
        next_neighbor_id = neighbor_ids[0]
        return _clone_rooted_tbr_graph_component(
            graph,
            next_neighbor_id,
            from_neighbor_id=node_id,
            incoming_length=_combine_optional_lengths(
                incoming_length,
                graph.adjacency[node_id][next_neighbor_id],
            ),
        )
    clone = TreeNode(name=graph.node_labels[node_id], branch_length=incoming_length)
    sorted_neighbor_ids = sorted(
        neighbor_ids,
        key=lambda neighbor_id: (
            len(
                _graph_side_tip_taxa(
                    graph,
                    neighbor_id,
                    blocked_node_id=node_id,
                )
            ),
            tuple(
                _graph_side_tip_taxa(
                    graph,
                    neighbor_id,
                    blocked_node_id=node_id,
                )
            ),
        ),
    )
    for neighbor_id in sorted_neighbor_ids:
        clone.append_child(
            _clone_rooted_tbr_graph_component(
                graph,
                neighbor_id,
                from_neighbor_id=node_id,
                incoming_length=graph.adjacency[node_id][neighbor_id],
            )
        )
    return clone


def _combine_optional_lengths(
    left_length: float | None,
    right_length: float | None,
) -> float | None:
    if left_length is None and right_length is None:
        return None
    return (left_length or 0.0) + (right_length or 0.0)
