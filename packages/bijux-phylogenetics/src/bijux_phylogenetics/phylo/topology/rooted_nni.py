from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from bijux_phylogenetics.io.newick import load_newick, loads_newick, write_newick
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_clade_id,
    rooted_topology_fingerprint,
)
from bijux_phylogenetics.phylo.topology.models import (
    RootedNniNeighborRow,
    RootedNniNeighborhoodReport,
)

from .tree import PhyloTree, TreeNode, descendant_taxa


@dataclass(frozen=True, slots=True)
class RootedNniMoveCandidate:
    """One deterministic rooted NNI move over a binary tree."""

    parent_node_id: str
    child_node_id: str
    sibling_node_id: str
    exchanged_child_node_id: str
    pivot_branch_id: str
    sibling_clade_id: str
    exchanged_clade_id: str


def rooted_nni_node_sort_key(node: TreeNode) -> tuple[int, tuple[str, ...]]:
    """Sort candidate NNI branches deterministically by descendant taxa."""
    descendants = tuple(descendant_taxa(node))
    return (len(descendants), descendants)


def rooted_nni_clade_id(node: TreeNode) -> str:
    """Render one branch endpoint as a stable descendant-clade identifier."""
    return canonical_clade_id(frozenset(descendant_taxa(node)))


def require_rooted_nni_node_id(node: TreeNode) -> str:
    """Require one refreshed stable node identifier on a rooted NNI tree."""
    if node.node_id is None:
        raise AssertionError("rooted NNI search requires refreshed node identities")
    return node.node_id


def iter_rooted_nni_move_candidates(tree: PhyloTree):
    """Yield deterministic rooted NNI candidates for one rooted binary tree."""
    for parent in tree.iter_internal_nodes(order="preorder"):
        if len(parent.children) != 2:
            continue
        sorted_parent_children = sorted(parent.children, key=rooted_nni_node_sort_key)
        for child in sorted_parent_children:
            if child.is_leaf() or len(child.children) != 2:
                continue
            sibling = next(
                candidate
                for candidate in sorted_parent_children
                if candidate is not child
            )
            sorted_child_children = sorted(
                child.children,
                key=rooted_nni_node_sort_key,
            )
            for exchanged_child in sorted_child_children:
                yield RootedNniMoveCandidate(
                    parent_node_id=require_rooted_nni_node_id(parent),
                    child_node_id=require_rooted_nni_node_id(child),
                    sibling_node_id=require_rooted_nni_node_id(sibling),
                    exchanged_child_node_id=require_rooted_nni_node_id(
                        exchanged_child
                    ),
                    pivot_branch_id=rooted_nni_clade_id(child),
                    sibling_clade_id=rooted_nni_clade_id(sibling),
                    exchanged_clade_id=rooted_nni_clade_id(exchanged_child),
                )


def apply_rooted_nni_move(
    tree: PhyloTree,
    candidate: RootedNniMoveCandidate,
) -> PhyloTree:
    """Return one copied rooted tree with the selected NNI move applied."""
    swapped_tree = tree.copy().refresh()
    parent = swapped_tree.node_by_id(candidate.parent_node_id)
    child = swapped_tree.node_by_id(candidate.child_node_id)
    sibling = swapped_tree.node_by_id(candidate.sibling_node_id)
    exchanged_child = swapped_tree.node_by_id(candidate.exchanged_child_node_id)
    remaining_child = next(
        branch for branch in child.children if branch is not exchanged_child
    )
    child.replace_children([remaining_child, sibling])
    parent.replace_children([child, exchanged_child])
    return swapped_tree.refresh()


def validate_rooted_nni_tree(tree: PhyloTree) -> None:
    """Require one structurally valid strictly bifurcating binary-root tree."""
    validation_errors = tree.validation_errors()
    if validation_errors:
        raise ValueError(
            "rooted NNI enumeration requires a structurally valid tree: "
            + "; ".join(validation_errors)
        )
    if len(tree.root.children) != 2:
        raise ValueError("rooted NNI enumeration requires a binary root")
    invalid_internal_nodes = [
        node.node_id
        for node in tree.iter_internal_nodes(order="preorder")
        if len(node.children) != 2
    ]
    if invalid_internal_nodes:
        raise ValueError("rooted NNI enumeration requires a strictly bifurcating tree")


def expected_rooted_nni_neighbor_count(tree: PhyloTree) -> int:
    """Return the exact rooted NNI neighbor count for one strict binary tree."""
    validate_rooted_nni_tree(tree)
    return max(0, 2 * (tree.internal_node_count - 1))


def enumerate_rooted_nni_neighbors(
    tree: PhyloTree | Path,
) -> RootedNniNeighborhoodReport:
    """Enumerate all legal rooted NNI neighbors exactly once for one tree."""
    resolved_tree, input_tree_path = _resolve_rooted_nni_tree(tree)
    validate_rooted_nni_tree(resolved_tree)
    input_tip_taxa = set(resolved_tree.tip_names)
    neighbor_rows: list[RootedNniNeighborRow] = []
    topology_counts: dict[str, int] = {}
    for neighbor_index, candidate in enumerate(
        iter_rooted_nni_move_candidates(resolved_tree),
        start=1,
    ):
        neighbor_tree = apply_rooted_nni_move(resolved_tree, candidate)
        topology_fingerprint = rooted_topology_fingerprint(neighbor_tree)
        topology_counts[topology_fingerprint] = (
            topology_counts.get(topology_fingerprint, 0) + 1
        )
        neighbor_rows.append(
            RootedNniNeighborRow(
                neighbor_index=neighbor_index,
                parent_node_id=candidate.parent_node_id,
                child_node_id=candidate.child_node_id,
                sibling_node_id=candidate.sibling_node_id,
                exchanged_child_node_id=candidate.exchanged_child_node_id,
                pivot_branch_id=candidate.pivot_branch_id,
                sibling_clade_id=candidate.sibling_clade_id,
                exchanged_clade_id=candidate.exchanged_clade_id,
                neighbor_tree_newick=neighbor_tree.to_newick(),
                neighbor_topology_fingerprint=topology_fingerprint,
                tip_order=neighbor_tree.tip_names,
                validation_errors=neighbor_tree.validation_errors(),
            )
        )
    duplicate_neighbor_topologies = sorted(
        fingerprint
        for fingerprint, count in topology_counts.items()
        if count > 1
    )
    report = RootedNniNeighborhoodReport(
        algorithm="rooted-nni-neighbor-enumeration",
        input_tree_path=input_tree_path,
        input_tree_newick=resolved_tree.to_newick(),
        tip_count=resolved_tree.tip_count,
        internal_node_count=resolved_tree.internal_node_count,
        rooted=resolved_tree.rooted,
        strictly_bifurcating=True,
        expected_neighbor_count=expected_rooted_nni_neighbor_count(resolved_tree),
        generated_neighbor_count=len(neighbor_rows),
        unique_neighbor_topology_count=len(topology_counts),
        duplicate_neighbor_topologies=duplicate_neighbor_topologies,
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
    _validate_rooted_nni_neighbor_report(report)
    return report


def _resolve_rooted_nni_tree(
    tree: PhyloTree | Path,
) -> tuple[PhyloTree, Path | None]:
    resolved_tree_path = tree if isinstance(tree, Path) else None
    resolved_tree = load_newick(tree) if isinstance(tree, Path) else tree.copy()
    return resolved_tree.refresh(), resolved_tree_path


def _validate_rooted_nni_neighbor_report(
    report: RootedNniNeighborhoodReport,
) -> None:
    if report.generated_neighbor_count != report.expected_neighbor_count:
        raise ValueError(
            "rooted NNI enumeration did not generate the expected number of neighbors"
        )
    if report.duplicate_neighbor_topologies:
        raise ValueError("rooted NNI enumeration generated duplicate neighbor topologies")
    if report.missing_tip_taxa or report.unexpected_tip_taxa:
        raise ValueError("rooted NNI enumeration changed the input taxon set")
    invalid_neighbor_rows = [
        row.neighbor_index for row in report.neighbor_rows if row.validation_errors
    ]
    if invalid_neighbor_rows:
        raise ValueError("rooted NNI enumeration generated invalid neighbor trees")


def write_rooted_nni_neighbor_table(
    path: Path,
    report: RootedNniNeighborhoodReport,
) -> Path:
    """Write one row per rooted NNI neighbor with its defining move metadata."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "neighbor_index",
                "parent_node_id",
                "child_node_id",
                "sibling_node_id",
                "exchanged_child_node_id",
                "pivot_branch_id",
                "sibling_clade_id",
                "exchanged_clade_id",
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
                    row.parent_node_id,
                    row.child_node_id,
                    row.sibling_node_id,
                    row.exchanged_child_node_id,
                    row.pivot_branch_id,
                    row.sibling_clade_id,
                    row.exchanged_clade_id,
                    row.neighbor_topology_fingerprint,
                    ",".join(row.tip_order),
                    ",".join(row.validation_errors),
                    row.neighbor_tree_newick,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_rooted_nni_run_json(
    path: Path,
    report: RootedNniNeighborhoodReport,
) -> Path:
    """Write one machine-readable rooted NNI neighborhood payload."""
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
        "expected_neighbor_count": report.expected_neighbor_count,
        "generated_neighbor_count": report.generated_neighbor_count,
        "unique_neighbor_topology_count": report.unique_neighbor_topology_count,
        "duplicate_neighbor_topologies": report.duplicate_neighbor_topologies,
        "missing_tip_taxa": report.missing_tip_taxa,
        "unexpected_tip_taxa": report.unexpected_tip_taxa,
        "input_validation_errors": report.input_validation_errors,
        "neighbor_rows": [
            {
                "neighbor_index": row.neighbor_index,
                "parent_node_id": row.parent_node_id,
                "child_node_id": row.child_node_id,
                "sibling_node_id": row.sibling_node_id,
                "exchanged_child_node_id": row.exchanged_child_node_id,
                "pivot_branch_id": row.pivot_branch_id,
                "sibling_clade_id": row.sibling_clade_id,
                "exchanged_clade_id": row.exchanged_clade_id,
                "neighbor_tree_newick": row.neighbor_tree_newick,
                "neighbor_topology_fingerprint": row.neighbor_topology_fingerprint,
                "tip_order": row.tip_order,
                "validation_errors": row.validation_errors,
            }
            for row in report.neighbor_rows
        ],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def write_rooted_nni_artifacts(
    out_dir: Path,
    report: RootedNniNeighborhoodReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one rooted NNI enumeration run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    input_tree_path = write_newick(
        out_dir / "input_tree.nwk",
        loads_newick(report.input_tree_newick),
    )
    neighbors_path = write_rooted_nni_neighbor_table(out_dir / "neighbors.tsv", report)
    run_json_path = write_rooted_nni_run_json(out_dir / "run.json", report)
    return {
        "input_tree_path": input_tree_path,
        "neighbors_path": neighbors_path,
        "run_json_path": run_json_path,
    }
