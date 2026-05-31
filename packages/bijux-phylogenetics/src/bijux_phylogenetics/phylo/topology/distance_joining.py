from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import InvalidDistanceMatrixError

TIE_TOLERANCE = 1e-12
ZERO_TOLERANCE = 1e-15

ClusterKey = tuple[str, ...]


@dataclass(slots=True)
class ActiveDistanceJoinCluster:
    """One active cluster tracked during an NJ-family distance-joining build."""

    key: ClusterKey
    node: TreeNode


def cluster_pair_key(
    left: ClusterKey,
    right: ClusterKey,
) -> tuple[ClusterKey, ClusterKey]:
    """Return one stable unordered key for a cluster-pair lookup."""
    return (left, right) if left < right else (right, left)


def normalized_branch_length(value: float) -> float:
    """Normalize tiny floating-point residue to exact zero branch length."""
    if abs(value) <= ZERO_TOLERANCE:
        return 0.0
    return value


def require_distance(
    distance_lookup: dict[tuple[str, str], float],
    left_identifier: str,
    right_identifier: str,
) -> float:
    """Require one finite directional distance from a full lookup."""
    try:
        value = distance_lookup[(left_identifier, right_identifier)]
    except KeyError as error:
        raise InvalidDistanceMatrixError(
            f"distance matrix is missing pair {left_identifier}/{right_identifier}"
        ) from error
    if not math.isfinite(value):
        raise InvalidDistanceMatrixError(
            f"distance matrix contains non-finite distance for {left_identifier}/{right_identifier}"
        )
    return value


def validate_distance_lookup(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> None:
    """Validate one full symmetric distance lookup before NJ-family joining."""
    if len(identifiers) < 2:
        raise InvalidDistanceMatrixError(
            "distance matrix must contain at least two taxa"
        )
    if len(set(identifiers)) != len(identifiers):
        raise InvalidDistanceMatrixError(
            "distance matrix contains duplicated taxon labels"
        )
    for left_identifier in identifiers:
        diagonal = require_distance(distance_lookup, left_identifier, left_identifier)
        if not math.isclose(diagonal, 0.0, abs_tol=TIE_TOLERANCE, rel_tol=0.0):
            raise InvalidDistanceMatrixError(
                "distance matrix has nonzero diagonal entries"
            )
        if diagonal < 0.0:
            raise InvalidDistanceMatrixError(
                "distance matrix contains negative distances"
            )
        for right_identifier in identifiers:
            if left_identifier >= right_identifier:
                continue
            left_to_right = require_distance(
                distance_lookup,
                left_identifier,
                right_identifier,
            )
            right_to_left = require_distance(
                distance_lookup,
                right_identifier,
                left_identifier,
            )
            if left_to_right < 0.0 or right_to_left < 0.0:
                raise InvalidDistanceMatrixError(
                    "distance matrix contains negative distances"
                )
            if not math.isclose(
                left_to_right,
                right_to_left,
                rel_tol=TIE_TOLERANCE,
                abs_tol=TIE_TOLERANCE,
            ):
                raise InvalidDistanceMatrixError(
                    "distance matrix contains asymmetric directional entries"
                )


def distance_sums(
    active_distances: dict[tuple[ClusterKey, ClusterKey], float],
    active_keys: list[ClusterKey],
) -> dict[ClusterKey, float]:
    """Return the row sums for one NJ-family active distance matrix."""
    return {
        key: sum(
            active_distances[cluster_pair_key(key, other_key)]
            for other_key in active_keys
            if other_key != key
        )
        for key in active_keys
    }


def choose_join_pair(
    active_distances: dict[tuple[ClusterKey, ClusterKey], float],
    active_keys: list[ClusterKey],
    *,
    method_name: str,
) -> tuple[ClusterKey, ClusterKey]:
    """Choose one deterministic NJ-family join pair from the active Q matrix."""
    row_sums = distance_sums(active_distances, active_keys)
    active_count = len(active_keys)
    best_pair: tuple[ClusterKey, ClusterKey] | None = None
    best_q_value: float | None = None
    for left_index, left_key in enumerate(active_keys):
        for right_key in active_keys[left_index + 1 :]:
            pair_key = cluster_pair_key(left_key, right_key)
            q_value = (
                (active_count - 2) * active_distances[pair_key]
                - row_sums[left_key]
                - row_sums[right_key]
            )
            if best_q_value is None or q_value < best_q_value - TIE_TOLERANCE:
                best_pair = pair_key
                best_q_value = q_value
                continue
            if (
                best_q_value is not None
                and math.isclose(
                    q_value,
                    best_q_value,
                    rel_tol=TIE_TOLERANCE,
                    abs_tol=TIE_TOLERANCE,
                )
                and pair_key < best_pair
            ):
                best_pair = pair_key
    if best_pair is None:
        raise InvalidDistanceMatrixError(
            f"distance matrix did not yield a valid {method_name.upper()} join pair"
        )
    return best_pair


def build_three_taxon_join_tree(
    active_clusters: dict[ClusterKey, ActiveDistanceJoinCluster],
    active_distances: dict[tuple[ClusterKey, ClusterKey], float],
    *,
    join_counter: int,
) -> PhyloTree:
    """Build the final three-taxon NJ-family tree from one reduced matrix."""
    left_key, middle_key, right_key = sorted(active_clusters)
    left_middle = active_distances[cluster_pair_key(left_key, middle_key)]
    left_right = active_distances[cluster_pair_key(left_key, right_key)]
    middle_right = active_distances[cluster_pair_key(middle_key, right_key)]
    left_length = (left_middle + left_right - middle_right) / 2.0
    middle_length = (left_middle + middle_right - left_right) / 2.0
    right_length = (left_right + middle_right - left_middle) / 2.0
    active_clusters[left_key].node.branch_length = normalized_branch_length(left_length)
    active_clusters[middle_key].node.branch_length = normalized_branch_length(
        middle_length
    )
    active_clusters[right_key].node.branch_length = normalized_branch_length(
        right_length
    )
    root = TreeNode(
        name=f"Inner{join_counter}",
        children=[
            active_clusters[left_key].node,
            active_clusters[middle_key].node,
            active_clusters[right_key].node,
        ],
    )
    return PhyloTree(root=root, source_format="newick", rooted=False)


def build_two_taxon_join_tree(
    active_clusters: dict[ClusterKey, ActiveDistanceJoinCluster],
    active_distances: dict[tuple[ClusterKey, ClusterKey], float],
    *,
    join_counter: int,
) -> PhyloTree:
    """Build the final two-taxon NJ-family tree from one reduced matrix."""
    left_key, right_key = sorted(active_clusters)
    pair_distance = active_distances[cluster_pair_key(left_key, right_key)]
    left_length = normalized_branch_length(pair_distance / 2.0)
    right_length = normalized_branch_length(pair_distance / 2.0)
    active_clusters[left_key].node.branch_length = left_length
    active_clusters[right_key].node.branch_length = right_length
    root = TreeNode(
        name=f"Inner{join_counter}",
        children=[active_clusters[left_key].node, active_clusters[right_key].node],
    )
    return PhyloTree(root=root, source_format="newick", rooted=False)
