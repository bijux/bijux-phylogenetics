from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.errors import InvalidDistanceMatrixError

_TIE_TOLERANCE = 1e-12
_ZERO_TOLERANCE = 1e-15


ClusterKey = tuple[str, ...]


@dataclass(slots=True)
class _ActiveCluster:
    key: ClusterKey
    node: TreeNode


def _cluster_pair_key(
    left: ClusterKey, right: ClusterKey
) -> tuple[ClusterKey, ClusterKey]:
    return (left, right) if left < right else (right, left)


def _normalized_branch_length(value: float) -> float:
    if abs(value) <= _ZERO_TOLERANCE:
        return 0.0
    return value


def _require_distance(
    distance_lookup: dict[tuple[str, str], float],
    left_identifier: str,
    right_identifier: str,
) -> float:
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


def _validate_distance_lookup(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> None:
    if len(identifiers) < 2:
        raise InvalidDistanceMatrixError(
            "distance matrix must contain at least two taxa"
        )
    if len(set(identifiers)) != len(identifiers):
        raise InvalidDistanceMatrixError(
            "distance matrix contains duplicated taxon labels"
        )
    for left_identifier in identifiers:
        diagonal = _require_distance(distance_lookup, left_identifier, left_identifier)
        if not math.isclose(diagonal, 0.0, abs_tol=_TIE_TOLERANCE, rel_tol=0.0):
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
            left_to_right = _require_distance(
                distance_lookup,
                left_identifier,
                right_identifier,
            )
            right_to_left = _require_distance(
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
                rel_tol=_TIE_TOLERANCE,
                abs_tol=_TIE_TOLERANCE,
            ):
                raise InvalidDistanceMatrixError(
                    "distance matrix contains asymmetric directional entries"
                )


def _distance_sums(
    active_distances: dict[tuple[ClusterKey, ClusterKey], float],
    active_keys: list[ClusterKey],
) -> dict[ClusterKey, float]:
    return {
        key: sum(
            active_distances[_cluster_pair_key(key, other_key)]
            for other_key in active_keys
            if other_key != key
        )
        for key in active_keys
    }


def _choose_join_pair(
    active_distances: dict[tuple[ClusterKey, ClusterKey], float],
    active_keys: list[ClusterKey],
) -> tuple[ClusterKey, ClusterKey]:
    distance_sums = _distance_sums(active_distances, active_keys)
    active_count = len(active_keys)
    best_pair: tuple[ClusterKey, ClusterKey] | None = None
    best_q_value: float | None = None
    for left_index, left_key in enumerate(active_keys):
        for right_key in active_keys[left_index + 1 :]:
            pair_key = _cluster_pair_key(left_key, right_key)
            q_value = (
                (active_count - 2) * active_distances[pair_key]
                - distance_sums[left_key]
                - distance_sums[right_key]
            )
            if best_q_value is None or q_value < best_q_value - _TIE_TOLERANCE:
                best_pair = pair_key
                best_q_value = q_value
                continue
            if (
                best_q_value is not None
                and math.isclose(
                    q_value,
                    best_q_value,
                    rel_tol=_TIE_TOLERANCE,
                    abs_tol=_TIE_TOLERANCE,
                )
                and pair_key < best_pair
            ):
                best_pair = pair_key
    if best_pair is None:
        raise InvalidDistanceMatrixError(
            "distance matrix did not yield a valid NJ join pair"
        )
    return best_pair


def _join_clusters(
    active_clusters: dict[ClusterKey, _ActiveCluster],
    active_distances: dict[tuple[ClusterKey, ClusterKey], float],
    *,
    join_counter: int,
) -> int:
    active_keys = sorted(active_clusters)
    left_key, right_key = _choose_join_pair(active_distances, active_keys)
    distance_sums = _distance_sums(active_distances, active_keys)
    active_count = len(active_keys)
    pair_distance = active_distances[_cluster_pair_key(left_key, right_key)]
    left_length = 0.5 * pair_distance + (
        (distance_sums[left_key] - distance_sums[right_key]) / (2 * (active_count - 2))
    )
    right_length = pair_distance - left_length
    active_clusters[left_key].node.branch_length = _normalized_branch_length(
        left_length
    )
    active_clusters[right_key].node.branch_length = _normalized_branch_length(
        right_length
    )
    merged_key = tuple(sorted(left_key + right_key))
    merged_node = TreeNode(
        name=f"Inner{join_counter}",
        children=[
            active_clusters[left_key].node,
            active_clusters[right_key].node,
        ],
    )
    updated_distances = {
        pair_key: value
        for pair_key, value in active_distances.items()
        if left_key not in pair_key and right_key not in pair_key
    }
    for other_key in active_keys:
        if other_key in {left_key, right_key}:
            continue
        updated_distances[_cluster_pair_key(merged_key, other_key)] = (
            active_distances[_cluster_pair_key(left_key, other_key)]
            + active_distances[_cluster_pair_key(right_key, other_key)]
            - pair_distance
        ) / 2.0
    del active_clusters[left_key]
    del active_clusters[right_key]
    active_clusters[merged_key] = _ActiveCluster(key=merged_key, node=merged_node)
    active_distances.clear()
    active_distances.update(updated_distances)
    return join_counter + 1


def _build_three_taxon_tree(
    active_clusters: dict[ClusterKey, _ActiveCluster],
    active_distances: dict[tuple[ClusterKey, ClusterKey], float],
    *,
    join_counter: int,
) -> PhyloTree:
    left_key, middle_key, right_key = sorted(active_clusters)
    left_middle = active_distances[_cluster_pair_key(left_key, middle_key)]
    left_right = active_distances[_cluster_pair_key(left_key, right_key)]
    middle_right = active_distances[_cluster_pair_key(middle_key, right_key)]
    left_length = (left_middle + left_right - middle_right) / 2.0
    middle_length = (left_middle + middle_right - left_right) / 2.0
    right_length = (left_right + middle_right - left_middle) / 2.0
    active_clusters[left_key].node.branch_length = _normalized_branch_length(
        left_length
    )
    active_clusters[middle_key].node.branch_length = _normalized_branch_length(
        middle_length
    )
    active_clusters[right_key].node.branch_length = _normalized_branch_length(
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


def _build_two_taxon_tree(
    active_clusters: dict[ClusterKey, _ActiveCluster],
    active_distances: dict[tuple[ClusterKey, ClusterKey], float],
    *,
    join_counter: int,
) -> PhyloTree:
    left_key, right_key = sorted(active_clusters)
    pair_distance = active_distances[_cluster_pair_key(left_key, right_key)]
    left_length = _normalized_branch_length(pair_distance / 2.0)
    right_length = _normalized_branch_length(pair_distance / 2.0)
    active_clusters[left_key].node.branch_length = left_length
    active_clusters[right_key].node.branch_length = right_length
    root = TreeNode(
        name=f"Inner{join_counter}",
        children=[active_clusters[left_key].node, active_clusters[right_key].node],
    )
    return PhyloTree(root=root, source_format="newick", rooted=False)


def build_neighbor_joining_tree(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> PhyloTree:
    """Build one deterministic Neighbor-Joining tree from a full distance lookup."""
    _validate_distance_lookup(identifiers, distance_lookup)
    active_clusters = {
        (identifier,): _ActiveCluster(key=(identifier,), node=TreeNode(name=identifier))
        for identifier in identifiers
    }
    active_distances: dict[tuple[ClusterKey, ClusterKey], float] = {}
    for left_index, left_identifier in enumerate(identifiers):
        for right_identifier in identifiers[left_index + 1 :]:
            active_distances[
                _cluster_pair_key((left_identifier,), (right_identifier,))
            ] = distance_lookup[(left_identifier, right_identifier)]
    join_counter = 1
    while len(active_clusters) > 3:
        join_counter = _join_clusters(
            active_clusters,
            active_distances,
            join_counter=join_counter,
        )
    if len(active_clusters) == 3:
        return _build_three_taxon_tree(
            active_clusters,
            active_distances,
            join_counter=join_counter,
        )
    return _build_two_taxon_tree(
        active_clusters,
        active_distances,
        join_counter=join_counter,
    )
