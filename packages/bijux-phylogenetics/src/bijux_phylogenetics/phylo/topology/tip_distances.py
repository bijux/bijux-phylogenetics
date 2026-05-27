from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError

MissingBranchLengthPolicy = str


@dataclass(slots=True)
class TipDistanceMatrixRow:
    """One long-form pairwise tip-distance entry."""

    left_identifier: str
    right_identifier: str
    distance: float


@dataclass(slots=True)
class TipDistanceMatrixReport:
    """Deterministic pairwise tip-distance matrix for one tree."""

    tree_path: Path | None
    identifiers: list[str]
    rooted: bool | None
    matrix: list[list[float]]
    pairs: list[TipDistanceMatrixRow]
    missing_branch_length_policy: MissingBranchLengthPolicy
    complete_branch_lengths: bool
    branch_length_count: int
    expected_branch_length_count: int
    pair_count: int
    diagonal_zero: bool
    symmetric: bool


def _root_to_tip_paths(tree: PhyloTree) -> dict[str, list[TreeNode]]:
    paths: dict[str, list[TreeNode]] = {}

    def visit(node: TreeNode, path: list[TreeNode]) -> None:
        next_path = [*path, node]
        if node.is_leaf():
            if node.name is not None:
                paths[node.name] = next_path
            return
        for child in node.children:
            visit(child, next_path)

    visit(tree.root, [])
    return paths


def _common_prefix_length(left: list[TreeNode], right: list[TreeNode]) -> int:
    length = 0
    for left_node, right_node in zip(left, right, strict=False):
        if left_node is not right_node:
            break
        length += 1
    return length


def _effective_branch_length(
    node: TreeNode,
    *,
    policy: MissingBranchLengthPolicy,
) -> float:
    if node.branch_length is not None:
        return node.branch_length
    if policy == "unit-length":
        return 1.0
    raise InvalidBranchLengthError(
        "tree requires complete branch lengths for tip-distance calculations",
        code="tip_distance_missing_branch_lengths",
        details={"missing_branch_length_policy": policy},
    )


def _interpreted_rooted_state(tree: PhyloTree) -> bool:
    if tree.rooted is True:
        return True
    return len(tree.root.children) == 2


def _root_distance_lookup(
    tree: PhyloTree,
    *,
    policy: MissingBranchLengthPolicy,
) -> dict[str, float]:
    distances: dict[str, float] = {tree.root.node_id or "": 0.0}

    def visit(node: TreeNode) -> None:
        base_distance = distances[node.node_id or ""]
        for child in node.children:
            distances[child.node_id or ""] = base_distance + _effective_branch_length(
                child,
                policy=policy,
            )
            visit(child)

    visit(tree.root)
    return distances


def summarize_tree_tip_distances(
    tree: PhyloTree,
    *,
    tree_path: Path | None = None,
    missing_branch_length_policy: MissingBranchLengthPolicy = "error",
) -> TipDistanceMatrixReport:
    """Compute one deterministic pairwise tip-distance matrix for one in-memory tree."""
    if missing_branch_length_policy not in {"error", "unit-length"}:
        raise ValueError(
            "missing_branch_length_policy must be 'error' or 'unit-length'"
        )
    identifiers = list(tree.tip_names)
    expected_branch_length_count = sum(
        1 for node in tree.iter_nodes() if node is not tree.root
    )
    branch_length_count = sum(
        1
        for node in tree.iter_nodes()
        if node is not tree.root and node.branch_length is not None
    )
    complete_branch_lengths = branch_length_count == expected_branch_length_count
    if not complete_branch_lengths and missing_branch_length_policy == "error":
        raise InvalidBranchLengthError(
            "tree requires complete branch lengths for tip-distance calculations",
            code="tip_distance_missing_branch_lengths",
            details={
                "missing_branch_length_policy": missing_branch_length_policy,
                "branch_length_count": branch_length_count,
                "expected_branch_length_count": expected_branch_length_count,
            },
        )

    paths = _root_to_tip_paths(tree)
    root_distances = _root_distance_lookup(
        tree,
        policy=missing_branch_length_policy,
    )

    matrix: list[list[float]] = []
    pairs: list[TipDistanceMatrixRow] = []
    for left_identifier in identifiers:
        left_path = paths[left_identifier]
        left_root_distance = root_distances[left_path[-1].node_id or ""]
        row: list[float] = []
        for right_identifier in identifiers:
            if left_identifier == right_identifier:
                distance = 0.0
            else:
                right_path = paths[right_identifier]
                right_root_distance = root_distances[right_path[-1].node_id or ""]
                prefix_length = _common_prefix_length(left_path, right_path)
                mrca_distance = root_distances[
                    left_path[prefix_length - 1].node_id or ""
                ]
                distance = (
                    left_root_distance + right_root_distance - (2.0 * mrca_distance)
                )
            distance = round(distance, 15)
            row.append(distance)
            pairs.append(
                TipDistanceMatrixRow(
                    left_identifier=left_identifier,
                    right_identifier=right_identifier,
                    distance=distance,
                )
            )
        matrix.append(row)

    diagonal_zero = all(
        matrix[index][index] == 0.0 for index in range(len(identifiers))
    )
    symmetric = all(
        matrix[left_index][right_index] == matrix[right_index][left_index]
        for left_index in range(len(identifiers))
        for right_index in range(len(identifiers))
    )
    return TipDistanceMatrixReport(
        tree_path=tree_path,
        identifiers=identifiers,
        rooted=_interpreted_rooted_state(tree),
        matrix=matrix,
        pairs=pairs,
        missing_branch_length_policy=missing_branch_length_policy,
        complete_branch_lengths=complete_branch_lengths,
        branch_length_count=branch_length_count,
        expected_branch_length_count=expected_branch_length_count,
        pair_count=len(pairs),
        diagonal_zero=diagonal_zero,
        symmetric=symmetric,
    )


def compute_tree_tip_distance_matrix(
    path: Path,
    *,
    missing_branch_length_policy: MissingBranchLengthPolicy = "error",
) -> TipDistanceMatrixReport:
    """Compute a deterministic pairwise tip-distance matrix for one tree path."""
    tree = load_tree(path)
    return summarize_tree_tip_distances(
        tree,
        tree_path=path,
        missing_branch_length_policy=missing_branch_length_policy,
    )


def write_tree_tip_distance_matrix(
    path: Path,
    report: TipDistanceMatrixReport,
) -> Path:
    """Write one deterministic wide pairwise tip-distance matrix."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["taxon\t" + "\t".join(report.identifiers)]
    for left_index, left_identifier in enumerate(report.identifiers):
        values = [
            format(report.matrix[left_index][right_index], ".15g")
            for right_index in range(len(report.identifiers))
        ]
        lines.append(f"{left_identifier}\t" + "\t".join(values))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_tree_tip_distance_long_table(
    path: Path,
    report: TipDistanceMatrixReport,
) -> Path:
    """Write one deterministic long-form tip-distance table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["left_identifier\tright_identifier\tdistance"]
    for row in report.pairs:
        lines.append(
            "\t".join(
                [
                    row.left_identifier,
                    row.right_identifier,
                    format(row.distance, ".15g"),
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
