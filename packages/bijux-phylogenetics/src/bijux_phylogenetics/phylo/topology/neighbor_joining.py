from __future__ import annotations

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .distance_joining import (
    ActiveDistanceJoinCluster as _ActiveCluster,
)
from .distance_joining import (
    ClusterKey,
    build_three_taxon_join_tree,
    build_two_taxon_join_tree,
    choose_join_pair,
)
from .distance_joining import (
    cluster_pair_key as _cluster_pair_key,
)
from .distance_joining import (
    distance_sums as _distance_sums,
)
from .distance_joining import (
    normalized_branch_length as _normalized_branch_length,
)
from .distance_joining import (
    validate_distance_lookup as _validate_distance_lookup,
)


def _join_clusters(
    active_clusters: dict[ClusterKey, _ActiveCluster],
    active_distances: dict[tuple[ClusterKey, ClusterKey], float],
    *,
    join_counter: int,
) -> int:
    active_keys = sorted(active_clusters)
    left_key, right_key = choose_join_pair(
        active_distances,
        active_keys,
        method_name="nj",
    )
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
        return build_three_taxon_join_tree(
            active_clusters,
            active_distances,
            join_counter=join_counter,
        )
    return build_two_taxon_join_tree(
        active_clusters,
        active_distances,
        join_counter=join_counter,
    )
